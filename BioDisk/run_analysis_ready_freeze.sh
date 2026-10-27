#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
READY_ROOT="${PROJECT_ROOT}/02_analysis_ready"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${READY_ROOT}/freeze_${STAMP}"
PY_SCRIPT="${RUN_DIR}/build_analysis_ready_freeze.py"
LOG_FILE="${RUN_DIR}/analysis_ready_freeze.log"

pointer_value() {
  local pointer_file="$1"
  local pointer_key="$2"
  if [[ ! -s "${pointer_file}" ]]; then
    echo "ERROR=pointer_missing:${pointer_file}" >&2
    return 1
  fi
  sed -n "s/^${pointer_key}=//p" "${pointer_file}" | tail -n 1
}

GEO_FREEZE_RUN="$(pointer_value "${AUDIT_ROOT}/latest_geo_freeze.txt" "LATEST_GEO_FREEZE")"
CANDIDATE_AUDIT_RUN="$(pointer_value "${AUDIT_ROOT}/latest_candidate_gene_audit.txt" "LATEST_CANDIDATE_GENE_AUDIT")"
PROMOTER_MAP_RUN="$(pointer_value "${AUDIT_ROOT}/latest_acox2_promoter_map.txt" "LATEST_ACOX2_PROMOTER_MAP")"

mkdir -p "${RUN_DIR}"

cat > "${PY_SCRIPT}" <<'PY'
#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import math
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO

from openpyxl import load_workbook


PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
RUN_DIR = Path(os.environ["RUN_DIR"])
GEO_FREEZE_RUN = Path(os.environ["GEO_FREEZE_RUN"])
CANDIDATE_AUDIT_RUN = Path(os.environ["CANDIDATE_AUDIT_RUN"])
PROMOTER_MAP_RUN = Path(os.environ["PROMOTER_MAP_RUN"])
RAW_DIR = PROJECT_ROOT / "00_raw_data"

RNA_FILE = RAW_DIR / "TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz"
METH_FILE = RAW_DIR / "TCGA_BRCA/TCGA.BRCA.HumanMethylation450.tsv.gz"
CDR_FILE = RAW_DIR / "TCGA_BRCA/TCGA-CDR-SupplementalTableS1.xlsx"
GSE70947_FILE = RAW_DIR / "GEO/GSE70947/GSE70947_series_matrix.txt.gz"
GSE65194_FILE = RAW_DIR / "GEO/GSE65194/GSE65194_series_matrix.txt.gz"

CANDIDATE_MANIFEST = CANDIDATE_AUDIT_RUN / "candidate_gene_manifest.tsv"
GEO_AVAILABILITY = CANDIDATE_AUDIT_RUN / "geo_candidate_gene_availability.tsv"
GEO_SAMPLE_SHEET = GEO_FREEZE_RUN / "geo_frozen_sample_sheet.tsv"
GEO_PAIR_SHEET = GEO_FREEZE_RUN / "gse70947_frozen_pairs.tsv"
PROMOTER_PROBES_FILE = PROMOTER_MAP_RUN / "acox2_450k_promoter_TSS1500_TSS200.tsv"
REGULATORY_PROBES_FILE = PROMOTER_MAP_RUN / "acox2_450k_regulatory_5prime.tsv"

EXPECTED = {
    "RNA_SAMPLES": int(os.environ.get("EXP_RNA_SAMPLES", "1218")),
    "RNA_GENES": int(os.environ.get("EXP_RNA_GENES", "8")),
    "RNA_PRIMARY": int(os.environ.get("EXP_RNA_PRIMARY", "1097")),
    "RNA_NORMAL": int(os.environ.get("EXP_RNA_NORMAL", "114")),
    "METH_SAMPLES": int(os.environ.get("EXP_METH_SAMPLES", "888")),
    "PROMOTER_PROBES": int(os.environ.get("EXP_PROMOTER_PROBES", "3")),
    "REGULATORY_PROBES": int(os.environ.get("EXP_REGULATORY_PROBES", "5")),
    "RNA_METH_EXACT": int(os.environ.get("EXP_RNA_METH_EXACT", "873")),
    "JOINT_PRIMARY": int(os.environ.get("EXP_JOINT_PRIMARY", "783")),
    "JOINT_NORMAL": int(os.environ.get("EXP_JOINT_NORMAL", "85")),
    "CLINICAL_RNA_EXACT": int(os.environ.get("EXP_CLINICAL_RNA_EXACT", "1218")),
    "CDR_BRCA_ROWS": int(os.environ.get("EXP_CDR_BRCA_ROWS", "1097")),
    "GEO_INCLUDED_TISSUES": int(os.environ.get("EXP_GEO_INCLUDED_TISSUES", "460")),
    "GSE70947_INCLUDED": int(os.environ.get("EXP_GSE70947_INCLUDED", "296")),
    "GSE65194_INCLUDED": int(os.environ.get("EXP_GSE65194_INCLUDED", "164")),
    "GSE70947_VALID_PAIRS": int(os.environ.get("EXP_GSE70947_VALID_PAIRS", "148")),
}

MISSING_TOKENS = {"", "na", "nan", "null", "--", "[not available]", "[not applicable]"}
CLINICAL_FIELDS = {
    "sample_id": ["sampleID", "sample_id", "SampleID"],
    "sample_type": ["sample_type", "sample_type.samples"],
    "PAM50": ["PAM50_mRNA_nature2012", "PAM50"],
    "AJCC_stage": ["AJCC_Stage_nature2012", "AJCC_stage", "pathologic_stage"],
    "ER_status": ["ER_Status_nature2012", "ER_status"],
    "PR_status": ["PR_Status_nature2012", "PR_status"],
    "HER2_status": ["HER2_Final_Status_nature2012", "HER2_status"],
}


def strip_quotes(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        text = text[1:-1]
    return text


def clean_value(value: object) -> str:
    text = strip_quotes(value)
    return "NA" if text.lower() in MISSING_TOKENS else text


def open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace", newline="")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {strip_quotes(key): strip_quotes(value) for key, value in row.items() if key is not None}
            for row in reader
        ]


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str], *, gz: bool = False) -> None:
    opener = gzip.open if gz else Path.open
    kwargs = {"mode": "wt", "encoding": "utf-8", "newline": ""} if gz else {
        "mode": "w",
        "encoding": "utf-8",
        "newline": "",
    }
    with opener(path, **kwargs) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_code(sample_id: str) -> str:
    match = re.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2})", sample_id.upper())
    return match.group(1) if match else "NA"


def patient_id(sample_id: str) -> str:
    fields = sample_id.upper().split("-")
    return "-".join(fields[:3]) if len(fields) >= 3 and fields[0] == "TCGA" else "NA"


def barcode_group(code: str) -> str:
    if code == "01":
        return "Primary_Tumor"
    if code == "06":
        return "Metastatic"
    if code == "11":
        return "Solid_Tissue_Normal"
    return "Other"


def read_wide_targets(path: Path, targets: set[str], *, upper_ids: bool) -> tuple[list[str], dict[str, list[str]], int]:
    rows: dict[str, list[str]] = {}
    feature_count = 0
    with open_text(path) as handle:
        header = [strip_quotes(item) for item in handle.readline().rstrip("\r\n").split("\t")]
        if len(header) < 2:
            raise RuntimeError(f"Invalid wide matrix header: {path}")
        samples = header[1:]
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if not fields:
                continue
            feature_count += 1
            raw_id = strip_quotes(fields[0])
            feature_id = raw_id.upper() if upper_ids else raw_id
            if feature_id in targets:
                if feature_id in rows:
                    raise RuntimeError(f"Duplicate target feature {feature_id} in {path}")
                values = [strip_quotes(value) for value in fields[1:]]
                if len(values) != len(samples):
                    raise RuntimeError(f"Value count mismatch for {feature_id} in {path}")
                rows[feature_id] = values
    return samples, rows, feature_count


def write_matrix(path: Path, samples: list[str], rows: dict[str, list[str]], order: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_id", *samples])
        for feature in order:
            writer.writerow([feature, *rows[feature]])


def discover_clinical_file() -> tuple[Path, dict[str, str]]:
    excluded = {RNA_FILE.resolve(), METH_FILE.resolve()}
    candidates: list[tuple[int, Path, dict[str, str]]] = []
    for path in sorted((RAW_DIR / "TCGA_BRCA").rglob("*")):
        if not path.is_file() or path.resolve() in excluded or path.suffix.lower() in {".xlsx", ".xls"}:
            continue
        if path.stat().st_size == 0:
            continue
        try:
            with open_text(path) as handle:
                header = [strip_quotes(item) for item in handle.readline().rstrip("\r\n").split("\t")]
        except (OSError, UnicodeError):
            continue
        resolved: dict[str, str] = {}
        for canonical, aliases in CLINICAL_FIELDS.items():
            hit = next((alias for alias in aliases if alias in header), "")
            if hit:
                resolved[canonical] = hit
        if "sample_id" in resolved:
            candidates.append((len(resolved), path, resolved))
    if not candidates:
        raise RuntimeError("No TCGA clinical matrix containing sampleID was found")
    candidates.sort(key=lambda item: (item[0], item[1].stat().st_size), reverse=True)
    score, path, resolved = candidates[0]
    if score < len(CLINICAL_FIELDS):
        missing = sorted(set(CLINICAL_FIELDS) - set(resolved))
        raise RuntimeError(f"Clinical fields unresolved in {path}: {missing}")
    return path, resolved


def read_clinical(path: Path, resolved: dict[str, str]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"Clinical matrix has no header: {path}")
        for raw in reader:
            row = {strip_quotes(key): strip_quotes(value) for key, value in raw.items() if key is not None}
            sid = clean_value(row.get(resolved["sample_id"], "NA"))
            if sid == "NA":
                continue
            result[sid] = {
                canonical: clean_value(row.get(source, "NA"))
                for canonical, source in resolved.items()
            }
    return result


def read_cdr(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        for sheet in workbook.worksheets:
            preview = list(sheet.iter_rows(min_row=1, max_row=12, values_only=True))
            for offset, values in enumerate(preview, start=1):
                header = [strip_quotes(value) for value in values]
                if "bcr_patient_barcode" not in header or "type" not in header:
                    continue
                records: list[dict[str, str]] = []
                for data in sheet.iter_rows(min_row=offset + 1, values_only=True):
                    row = {
                        header[index]: clean_value(value)
                        for index, value in enumerate(data)
                        if index < len(header) and header[index]
                    }
                    if row.get("type", "").upper() == "BRCA":
                        records.append(row)
                return header, records
    finally:
        workbook.close()
    raise RuntimeError("TCGA-CDR BRCA table was not found")


def extract_geo_matrix(path: Path, targets: set[str]) -> tuple[list[str], dict[str, list[str]], int]:
    samples: list[str] = []
    rows: dict[str, list[str]] = {}
    feature_count = 0
    in_table = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            fields = [strip_quotes(value) for value in line.split("\t")]
            if not samples:
                samples = fields[1:]
                continue
            feature_count += 1
            probe = fields[0]
            if probe in targets:
                values = fields[1:]
                if len(values) != len(samples):
                    raise RuntimeError(f"GEO value count mismatch for {probe} in {path}")
                rows[probe] = values
    return samples, rows, feature_count


def valid_pair_flag(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


required_files = [
    RNA_FILE,
    METH_FILE,
    CDR_FILE,
    GSE70947_FILE,
    GSE65194_FILE,
    CANDIDATE_MANIFEST,
    GEO_AVAILABILITY,
    GEO_SAMPLE_SHEET,
    GEO_PAIR_SHEET,
    PROMOTER_PROBES_FILE,
    REGULATORY_PROBES_FILE,
]
for required in required_files:
    if not required.is_file() or required.stat().st_size == 0:
        raise RuntimeError(f"Required input missing or empty: {required}")

print("PY_STEP=LOAD_FROZEN_FEATURES", flush=True)
candidate_rows = read_tsv(CANDIDATE_MANIFEST)
genes = [row["gene"].upper() for row in candidate_rows]
if len(genes) != len(set(genes)):
    raise RuntimeError("Candidate gene manifest contains duplicate genes")

promoter_rows = [row for row in read_tsv(PROMOTER_PROBES_FILE) if row.get("present_in_TCGA", "").lower() == "true"]
regulatory_rows = [row for row in read_tsv(REGULATORY_PROBES_FILE) if row.get("present_in_TCGA", "").lower() == "true"]
promoter_probes = [row["IlmnID"] for row in promoter_rows]
regulatory_probes = [row["IlmnID"] for row in regulatory_rows]

print("PY_STEP=FREEZE_TCGA_TARGET_MATRICES", flush=True)
rna_samples, rna_rows, rna_feature_count = read_wide_targets(RNA_FILE, set(genes), upper_ids=True)
meth_samples, regulatory_matrix, meth_feature_count = read_wide_targets(
    METH_FILE, set(regulatory_probes), upper_ids=False
)
promoter_matrix = {
    probe: regulatory_matrix[probe]
    for probe in promoter_probes
    if probe in regulatory_matrix
}

missing_genes = sorted(set(genes) - set(rna_rows))
missing_promoter = sorted(set(promoter_probes) - set(promoter_matrix))
missing_regulatory = sorted(set(regulatory_probes) - set(regulatory_matrix))
if missing_genes or missing_promoter or missing_regulatory:
    raise RuntimeError(
        f"Frozen features missing: genes={missing_genes}; promoter={missing_promoter}; regulatory={missing_regulatory}"
    )

write_matrix(RUN_DIR / "tcga_expression_8genes.tsv.gz", rna_samples, rna_rows, genes)
write_matrix(
    RUN_DIR / "tcga_methylation_acox2_promoter_3probes.tsv.gz",
    meth_samples,
    promoter_matrix,
    promoter_probes,
)
write_matrix(
    RUN_DIR / "tcga_methylation_acox2_5prime_5probes.tsv.gz",
    meth_samples,
    regulatory_matrix,
    regulatory_probes,
)

print("PY_STEP=FREEZE_TCGA_CLINICAL_AND_SURVIVAL", flush=True)
clinical_file, resolved_fields = discover_clinical_file()
clinical = read_clinical(clinical_file, resolved_fields)
cdr_header, cdr_rows = read_cdr(CDR_FILE)
cdr_by_patient = {row["bcr_patient_barcode"].upper(): row for row in cdr_rows}

rna_set = set(rna_samples)
meth_set = set(meth_samples)
all_samples = sorted(rna_set | meth_set)
cdr_fields = ["OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time"]
sample_manifest: list[dict[str, object]] = []
for sid in all_samples:
    code = sample_code(sid)
    pid = patient_id(sid)
    phenotype = clinical.get(sid, {})
    cdr = cdr_by_patient.get(pid, {})
    row: dict[str, object] = {
        "sample_id": sid,
        "patient_id": pid,
        "sample_type_code": code,
        "barcode_group": barcode_group(code),
        "in_RNA": sid in rna_set,
        "in_methylation": sid in meth_set,
        "in_RNA_methylation_exact": sid in rna_set and sid in meth_set,
        "clinical_exact_match": sid in clinical,
        "CDR_patient_match": pid in cdr_by_patient,
        "sample_type": phenotype.get("sample_type", "NA"),
        "PAM50": phenotype.get("PAM50", "NA"),
        "AJCC_stage": phenotype.get("AJCC_stage", "NA"),
        "ER_status": phenotype.get("ER_status", "NA"),
        "PR_status": phenotype.get("PR_status", "NA"),
        "HER2_status": phenotype.get("HER2_status", "NA"),
    }
    for field in cdr_fields:
        row[f"CDR_{field}"] = cdr.get(field, "NA")
    sample_manifest.append(row)

sample_fields = [
    "sample_id",
    "patient_id",
    "sample_type_code",
    "barcode_group",
    "in_RNA",
    "in_methylation",
    "in_RNA_methylation_exact",
    "clinical_exact_match",
    "CDR_patient_match",
    "sample_type",
    "PAM50",
    "AJCC_stage",
    "ER_status",
    "PR_status",
    "HER2_status",
    *[f"CDR_{field}" for field in cdr_fields],
]
write_tsv(RUN_DIR / "tcga_sample_manifest.tsv", sample_manifest, sample_fields)
write_tsv(RUN_DIR / "tcga_cdr_brca.tsv", cdr_rows, cdr_header)

joint_samples = sorted(rna_set & meth_set)
joint_primary = [sid for sid in joint_samples if sample_code(sid) == "01"]
joint_normal = [sid for sid in joint_samples if sample_code(sid) == "11"]
write_tsv(
    RUN_DIR / "tcga_joint_primary_samples.tsv",
    ({"sample_id": sid, "patient_id": patient_id(sid)} for sid in joint_primary),
    ["sample_id", "patient_id"],
)
write_tsv(
    RUN_DIR / "tcga_joint_normal_samples.tsv",
    ({"sample_id": sid, "patient_id": patient_id(sid)} for sid in joint_normal),
    ["sample_id", "patient_id"],
)

print("PY_STEP=FREEZE_GEO_PROBE_LEVEL_INPUTS", flush=True)
availability = read_tsv(GEO_AVAILABILITY)
geo_frozen = read_tsv(GEO_SAMPLE_SHEET)
pair_rows = read_tsv(GEO_PAIR_SHEET)
geo_source_files = {"GSE70947": GSE70947_FILE, "GSE65194": GSE65194_FILE}
geo_counts: dict[str, int] = {}
for dataset, matrix_file in geo_source_files.items():
    dataset_availability = [row for row in availability if row["dataset"] == dataset]
    target_probes = {
        probe
        for row in dataset_availability
        if row["availability"] == "PRESENT"
        for probe in row["present_probes"].split(";")
        if probe and probe != "NA"
    }
    geo_samples, geo_rows, _ = extract_geo_matrix(matrix_file, target_probes)
    if set(geo_rows) != target_probes:
        raise RuntimeError(f"{dataset} frozen target probe mismatch")
    write_matrix(
        RUN_DIR / f"{dataset}_candidate_probe_matrix.tsv.gz",
        geo_samples,
        geo_rows,
        sorted(target_probes),
    )
    included_rows = [
        row for row in geo_frozen if row["dataset"] == dataset and row["analysis_decision"] == "include"
    ]
    included_ids = {row["sample"] for row in included_rows}
    if not included_ids.issubset(set(geo_samples)):
        raise RuntimeError(f"{dataset} frozen included samples absent from expression matrix")
    write_tsv(
        RUN_DIR / f"{dataset}_included_sample_sheet.tsv",
        included_rows,
        list(included_rows[0].keys()),
    )
    geo_counts[dataset] = len(included_rows)

valid_pairs = sum(valid_pair_flag(row.get("valid_pair", "")) for row in pair_rows)

print("PY_STEP=WRITE_PROVENANCE_AND_VALIDATE", flush=True)
rna_primary = sum(sample_code(sid) == "01" for sid in rna_samples)
rna_normal = sum(sample_code(sid) == "11" for sid in rna_samples)
clinical_rna_exact = sum(sid in clinical for sid in rna_samples)
observed = {
    "RNA_SAMPLES": len(rna_samples),
    "RNA_GENES": len(rna_rows),
    "RNA_PRIMARY": rna_primary,
    "RNA_NORMAL": rna_normal,
    "METH_SAMPLES": len(meth_samples),
    "PROMOTER_PROBES": len(promoter_matrix),
    "REGULATORY_PROBES": len(regulatory_matrix),
    "RNA_METH_EXACT": len(joint_samples),
    "JOINT_PRIMARY": len(joint_primary),
    "JOINT_NORMAL": len(joint_normal),
    "CLINICAL_RNA_EXACT": clinical_rna_exact,
    "CDR_BRCA_ROWS": len(cdr_rows),
    "GEO_INCLUDED_TISSUES": sum(geo_counts.values()),
    "GSE70947_INCLUDED": geo_counts.get("GSE70947", 0),
    "GSE65194_INCLUDED": geo_counts.get("GSE65194", 0),
    "GSE70947_VALID_PAIRS": valid_pairs,
}

source_paths = [
    RNA_FILE,
    METH_FILE,
    clinical_file,
    CDR_FILE,
    GSE70947_FILE,
    GSE65194_FILE,
    CANDIDATE_MANIFEST,
    GEO_AVAILABILITY,
    GEO_SAMPLE_SHEET,
    GEO_PAIR_SHEET,
    PROMOTER_PROBES_FILE,
    REGULATORY_PROBES_FILE,
]
source_manifest = []
for path in source_paths:
    stat = path.stat()
    source_manifest.append(
        {
            "path": str(path),
            "size_bytes": stat.st_size,
            "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
            "sha256": sha256(path) if stat.st_size <= 50 * 1024 * 1024 else "NOT_COMPUTED_GT_50MB",
        }
    )
write_tsv(
    RUN_DIR / "source_file_manifest.tsv",
    source_manifest,
    ["path", "size_bytes", "mtime_utc", "sha256"],
)

summary = [
    "FREEZE_SCOPE=TCGA_GEO_ANALYSIS_READY_INPUTS",
    f"CLINICAL_FILE={clinical_file}",
    "CLINICAL_FIELDS=" + ";".join(f"{key}:{value}" for key, value in sorted(resolved_fields.items())),
    f"RNA_SOURCE_FEATURES={rna_feature_count}",
    f"METHYLATION_SOURCE_FEATURES={meth_feature_count}",
]
all_pass = True
for metric, expected_value in EXPECTED.items():
    value = observed[metric]
    status = "PASS" if value == expected_value else "FAIL"
    all_pass = all_pass and status == "PASS"
    summary.append(f"{metric}={value};EXPECTED={expected_value};STATUS={status}")
summary.extend(
    [
        "TCGA_GENE_COLLAPSE=NOT_APPLICABLE_GENE_LEVEL_SOURCE",
        "GEO_PROBE_COLLAPSE=DEFERRED_UNTIL_PROBE_CONCORDANCE_AUDIT",
        "ACOX2_PROMOTER_SCORE=DEFERRED_UNTIL_3PROBE_CONCORDANCE_AUDIT",
        "RAW_DATA_MODIFIED=NO",
        "ANALYSIS_READY_FREEZE_PASS" if all_pass else "ANALYSIS_READY_FREEZE_FAIL",
    ]
)
(RUN_DIR / "analysis_ready_freeze_summary.txt").write_text(
    "\n".join(summary) + "\n", encoding="utf-8"
)

output_manifest = []
for path in sorted(RUN_DIR.iterdir()):
    if path.is_file() and path.name not in {"analysis_ready_freeze.log", "output_file_manifest.tsv"}:
        output_manifest.append(
            {"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)}
        )
write_tsv(
    RUN_DIR / "output_file_manifest.tsv",
    output_manifest,
    ["file", "size_bytes", "sha256"],
)

print("\n".join(summary))
if not all_pass:
    raise SystemExit(2)
PY

export PROJECT_ROOT RUN_DIR GEO_FREEZE_RUN CANDIDATE_AUDIT_RUN PROMOTER_MAP_RUN

{
  echo "RUN_DIR=${RUN_DIR}"
  echo "GEO_FREEZE_RUN=${GEO_FREEZE_RUN}"
  echo "CANDIDATE_AUDIT_RUN=${CANDIDATE_AUDIT_RUN}"
  echo "PROMOTER_MAP_RUN=${PROMOTER_MAP_RUN}"
  echo "STEP=PYTHON_SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
  echo "SYNTAX_CHECK=PASS"
  echo "STEP=BUILD_ANALYSIS_READY_FREEZE"
  conda run --no-capture-output -n batch_bench python "${PY_SCRIPT}"
  tail -n 1 "${RUN_DIR}/analysis_ready_freeze_summary.txt" | grep -qx 'ANALYSIS_READY_FREEZE_PASS'
  printf 'LATEST_ANALYSIS_READY_FREEZE=%s\n' "${RUN_DIR}" > "${READY_ROOT}/latest_analysis_ready_freeze.txt"
  echo "ANALYSIS_LOGIC_PASS"
} 2>&1 | tee "${LOG_FILE}"
