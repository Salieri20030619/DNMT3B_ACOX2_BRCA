#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${AUDIT_ROOT}/candidate_gene_audit_${STAMP}"
PY_SCRIPT="${RUN_DIR}/candidate_gene_audit.py"
LOG_FILE="${RUN_DIR}/candidate_gene_audit.log"
mkdir -p "${RUN_DIR}"

cat > "${PY_SCRIPT}" <<'PY'
#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import math
import os
import re
from pathlib import Path


PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
RUN_DIR = Path(os.environ["RUN_DIR"])
RAW_DIR = PROJECT_ROOT / "00_raw_data"

FILES = {
    "tcga": RAW_DIR / "TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz",
    "gpl13607": RAW_DIR / "GEO/platforms/GPL13607_old_annotations.txt.gz",
    "gpl570": RAW_DIR / "GEO/platforms/GPL570.annot.gz",
    "gse70947": RAW_DIR / "GEO/GSE70947/GSE70947_series_matrix.txt.gz",
    "gse65194": RAW_DIR / "GEO/GSE65194/GSE65194_series_matrix.txt.gz",
}

GENES = [
    {
        "rank": 0,
        "gene": "DNMT3B",
        "role": "regulator",
        "kd_log2fc": "NA",
        "kd_fold_change": "NA",
        "kd_fdr": "NA",
        "expected_dnmt3b_corr": "NA",
    },
    {
        "rank": 1,
        "gene": "RGS16",
        "role": "candidate_metabolism_related",
        "kd_log2fc": 1.793,
        "kd_fold_change": 3.47,
        "kd_fdr": 0.00334,
        "expected_dnmt3b_corr": "negative",
    },
    {
        "rank": 2,
        "gene": "GPRC5A",
        "role": "candidate_metabolism_related",
        "kd_log2fc": -1.538,
        "kd_fold_change": 0.344,
        "kd_fdr": 0.00334,
        "expected_dnmt3b_corr": "positive",
    },
    {
        "rank": 3,
        "gene": "NPNT",
        "role": "candidate_metabolism_related",
        "kd_log2fc": -1.192,
        "kd_fold_change": 0.438,
        "kd_fdr": 0.01455,
        "expected_dnmt3b_corr": "positive",
    },
    {
        "rank": 4,
        "gene": "AKR1C2",
        "role": "candidate_metabolism_related",
        "kd_log2fc": -1.224,
        "kd_fold_change": 0.428,
        "kd_fdr": 0.03246,
        "expected_dnmt3b_corr": "positive",
    },
    {
        "rank": 5,
        "gene": "ADCY5",
        "role": "candidate_metabolism_related",
        "kd_log2fc": 1.137,
        "kd_fold_change": 2.20,
        "kd_fdr": 0.01954,
        "expected_dnmt3b_corr": "negative",
    },
    {
        "rank": 6,
        "gene": "FGFR4",
        "role": "candidate_metabolism_related",
        "kd_log2fc": 1.392,
        "kd_fold_change": 2.63,
        "kd_fdr": 0.01954,
        "expected_dnmt3b_corr": "negative",
    },
    {
        "rank": 7,
        "gene": "ACOX2",
        "role": "methylation_target",
        "kd_log2fc": "NA",
        "kd_fold_change": "NA",
        "kd_fdr": "NA",
        "expected_dnmt3b_corr": "to_be_estimated",
    },
]
GENE_NAMES = [row["gene"] for row in GENES]


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value


def valid_number(value: str) -> bool:
    value = strip_quotes(value)
    if value.lower() in {"", "na", "nan", "null", "--"}:
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def write_tsv(name: str, rows: list[dict], fields: list[str]) -> None:
    with (RUN_DIR / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


for label, path in FILES.items():
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Required file missing or empty: {label}={path}")

print("PY_STEP=WRITE_CANDIDATE_MANIFEST", flush=True)
manifest_fields = [
    "rank",
    "gene",
    "role",
    "kd_log2fc",
    "kd_fold_change",
    "kd_fdr",
    "expected_dnmt3b_corr",
]
write_tsv("candidate_gene_manifest.tsv", GENES, manifest_fields)

print("PY_STEP=AUDIT_TCGA_EXPRESSION", flush=True)
tcga_rows: dict[str, dict] = {}
with gzip.open(FILES["tcga"], "rt", encoding="utf-8", errors="replace") as handle:
    header = [strip_quotes(value) for value in handle.readline().rstrip("\r\n").split("\t")]
    sample_count = len(header) - 1
    for line in handle:
        fields = line.rstrip("\r\n").split("\t")
        gene = strip_quotes(fields[0]).upper()
        if gene not in GENE_NAMES:
            continue
        values = fields[1:]
        tcga_rows[gene] = {
            "gene": gene,
            "row_found": True,
            "sample_values": len(values),
            "valid_numeric_values": sum(valid_number(value) for value in values),
            "missing_values": sum(not valid_number(value) for value in values),
        }

tcga_audit = []
for gene in GENE_NAMES:
    tcga_audit.append(
        tcga_rows.get(
            gene,
            {
                "gene": gene,
                "row_found": False,
                "sample_values": 0,
                "valid_numeric_values": 0,
                "missing_values": sample_count,
            },
        )
    )
write_tsv(
    "tcga_candidate_gene_availability.tsv",
    tcga_audit,
    ["gene", "row_found", "sample_values", "valid_numeric_values", "missing_values"],
)


def platform_hits(path: Path, platform: str) -> list[dict]:
    patterns = {
        gene: re.compile(rf"(?<![A-Z0-9]){re.escape(gene)}(?![A-Z0-9])")
        for gene in GENE_NAMES
    }
    hits: dict[tuple[str, str, str], dict] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("!", "#")) or "\t" not in line:
                continue
            fields = line.rstrip("\r\n").split("\t")
            probe = strip_quotes(fields[0])
            upper = line.upper()
            for gene, pattern in patterns.items():
                if pattern.search(upper):
                    key = (platform, gene, probe)
                    hits[key] = {
                        "platform": platform,
                        "gene": gene,
                        "probe_id": probe,
                        "annotation_record": line.rstrip("\r\n"),
                    }
    return list(hits.values())


print("PY_STEP=AUDIT_GEO_PLATFORM_ANNOTATIONS", flush=True)
hits = platform_hits(FILES["gpl13607"], "GPL13607")
hits += platform_hits(FILES["gpl570"], "GPL570")
write_tsv(
    "geo_candidate_probe_annotations.tsv",
    hits,
    ["platform", "gene", "probe_id", "annotation_record"],
)


def matrix_probes(path: Path, targets: set[str]) -> tuple[int, set[str]]:
    features = 0
    present: set[str] = set()
    in_table = False
    header_seen = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            if not header_seen:
                header_seen = True
                continue
            probe = strip_quotes(line.split("\t", 1)[0])
            features += 1
            if probe in targets:
                present.add(probe)
    return features, present


probe_sets: dict[tuple[str, str], set[str]] = {}
for row in hits:
    probe_sets.setdefault((row["platform"], row["gene"]), set()).add(row["probe_id"])

targets_13607 = {probe for (platform, _), probes in probe_sets.items() if platform == "GPL13607" for probe in probes}
targets_570 = {probe for (platform, _), probes in probe_sets.items() if platform == "GPL570" for probe in probes}

print("PY_STEP=AUDIT_GEO_MATRIX_PRESENCE", flush=True)
features_70947, present_70947 = matrix_probes(FILES["gse70947"], targets_13607)
features_65194, present_65194 = matrix_probes(FILES["gse65194"], targets_570)

geo_rows = []
for dataset, platform, present in (
    ("GSE70947", "GPL13607", present_70947),
    ("GSE65194", "GPL570", present_65194),
):
    for gene in GENE_NAMES:
        annotated = probe_sets.get((platform, gene), set())
        found = annotated & present
        geo_rows.append(
            {
                "dataset": dataset,
                "platform": platform,
                "gene": gene,
                "annotated_probe_count": len(annotated),
                "present_probe_count": len(found),
                "present_probes": ";".join(sorted(found)) if found else "NA",
                "availability": "PRESENT" if found else "ABSENT",
            }
        )
write_tsv(
    "geo_candidate_gene_availability.tsv",
    geo_rows,
    [
        "dataset",
        "platform",
        "gene",
        "annotated_probe_count",
        "present_probe_count",
        "present_probes",
        "availability",
    ],
)

tcga_present = sum(bool(row["row_found"]) for row in tcga_audit)
gse70947_present = sum(
    row["availability"] == "PRESENT" for row in geo_rows if row["dataset"] == "GSE70947"
)
gse65194_present = sum(
    row["availability"] == "PRESENT" for row in geo_rows if row["dataset"] == "GSE65194"
)
all_tcga_valid = tcga_present == len(GENE_NAMES) and all(
    row["valid_numeric_values"] == sample_count for row in tcga_audit
)

summary = [
    "AUDIT_SCOPE=DNMT3B_ACOX2_AND_SIX_CANDIDATES",
    f"TCGA_SAMPLES={sample_count}",
    f"TCGA_REQUIRED_GENES={len(GENE_NAMES)}",
    f"TCGA_PRESENT_GENES={tcga_present}",
    f"GSE70947_FEATURES={features_70947}",
    f"GSE70947_PRESENT_GENES={gse70947_present}",
    f"GSE65194_FEATURES={features_65194}",
    f"GSE65194_PRESENT_GENES={gse65194_present}",
    "KD_DIRECTION_USE=HYPOTHESIS_ONLY_NOT_FILTER",
    "RAW_DATA_MODIFIED=NO",
    "CANDIDATE_GENE_AUDIT_PASS" if all_tcga_valid else "CANDIDATE_GENE_AUDIT_FAIL",
]
(RUN_DIR / "candidate_gene_audit_summary.txt").write_text(
    "\n".join(summary) + "\n", encoding="utf-8"
)
print("\n".join(summary))
if not all_tcga_valid:
    raise SystemExit(2)
PY

export PROJECT_ROOT RUN_DIR

{
  echo "RUN_DIR=${RUN_DIR}"
  echo "STEP=PYTHON_SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
  echo "SYNTAX_CHECK=PASS"
  echo "STEP=CANDIDATE_GENE_AUDIT"
  conda run --no-capture-output -n batch_bench python "${PY_SCRIPT}"
  tail -n 1 "${RUN_DIR}/candidate_gene_audit_summary.txt" | grep -qx 'CANDIDATE_GENE_AUDIT_PASS'
  printf 'LATEST_CANDIDATE_GENE_AUDIT=%s\n' "${RUN_DIR}" > "${AUDIT_ROOT}/latest_candidate_gene_audit.txt"
  echo "ANALYSIS_LOGIC_PASS"
} 2>&1 | tee "${LOG_FILE}"
