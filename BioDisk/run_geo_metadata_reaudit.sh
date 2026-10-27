#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
LATEST_FILE="${AUDIT_ROOT}/latest_run.txt"

if [[ ! -s "${LATEST_FILE}" ]]; then
  echo "ERROR=latest_audit_pointer_missing:${LATEST_FILE}" >&2
  exit 1
fi

SOURCE_RUN="$(sed -n 's/^LATEST_RUN=//p' "${LATEST_FILE}")"
if [[ -z "${SOURCE_RUN}" || ! -d "${SOURCE_RUN}" ]]; then
  echo "ERROR=latest_audit_run_invalid:${SOURCE_RUN}" >&2
  exit 1
fi

SOURCE_TSV="${SOURCE_RUN}/geo_sample_metadata_audit.tsv"
if [[ ! -s "${SOURCE_TSV}" ]]; then
  echo "ERROR=geo_metadata_missing:${SOURCE_TSV}" >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${AUDIT_ROOT}/geo_reaudit_${STAMP}"
PY_SCRIPT="${RUN_DIR}/geo_metadata_reaudit.py"
LOG_FILE="${RUN_DIR}/geo_metadata_reaudit.log"
mkdir -p "${RUN_DIR}"

cat > "${PY_SCRIPT}" <<'PY'
#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


SOURCE_TSV = Path(os.environ["SOURCE_TSV"])
RUN_DIR = Path(os.environ["RUN_DIR"])


def compact_pair_key(title: str) -> str:
    value = title.lower()
    value = re.sub(
        r"adenocarcinoma|adjacent|healthy|normal|tumou?r|cancer|breast|tissue|sample|paired|case",
        "",
        value,
    )
    value = re.sub(r"[^a-z0-9]", "", value)
    value = re.sub(r"^[nt]([0-9]+)$", r"\1", value)
    value = re.sub(r"^([0-9]+)[nt]$", r"\1", value)
    return value


def classify(dataset: str, title: str, metadata: str) -> tuple[str, str]:
    text = f"{title} | {metadata}".lower()
    is_cell_line = bool(re.search(r"cell[ _-]?line|cultured cells", text))
    if dataset == "GSE70947":
        if re.search(r"adjacent|normal|healthy|non[- _]?tumou?r", text):
            return "Normal", "Tissue"
        if re.search(r"adenocarcinoma|tumou?r|cancer", text):
            return "Tumor", "Tissue"
        return "Unresolved", "Unresolved"
    if is_cell_line:
        return "Cell_line", "Cell_line"
    if re.search(r"normal breast|normal tissue|healthy breast", text):
        return "Normal", "Tissue"
    if re.search(r"non[- _]?tnbc|non[- _]?triple", text):
        return "Non_TNBC", "Tissue"
    if re.search(r"tnbc|triple[- _]?negative", text):
        return "TNBC", "Tissue"
    if re.search(r"luminal|her2|tumou?r|cancer|carcinoma", text):
        return "Non_TNBC", "Tissue"
    return "Unresolved", "Unresolved"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


with SOURCE_TSV.open("r", encoding="utf-8", newline="") as handle:
    source_rows = list(csv.DictReader(handle, delimiter="\t"))

required = {"dataset", "sample", "title", "combined_metadata"}
if not source_rows or not required.issubset(source_rows[0]):
    raise RuntimeError("Required GEO metadata columns are absent")

rows: list[dict] = []
for source in source_rows:
    dataset = source["dataset"]
    group, material = classify(dataset, source["title"], source["combined_metadata"])
    rows.append(
        {
            "dataset": dataset,
            "sample": source["sample"],
            "group": group,
            "material": material,
            "pair_key": compact_pair_key(source["title"]) if dataset == "GSE70947" else "",
            "title": source["title"],
            "combined_metadata": source["combined_metadata"],
        }
    )

write_tsv(
    RUN_DIR / "geo_curated_sample_sheet.tsv",
    rows,
    ["dataset", "sample", "group", "material", "pair_key", "title", "combined_metadata"],
)

counts = Counter((row["dataset"], row["group"]) for row in rows)
count_rows = [
    {"dataset": dataset, "group": group, "n": n}
    for (dataset, group), n in sorted(counts.items())
]
write_tsv(RUN_DIR / "geo_corrected_group_counts.tsv", count_rows, ["dataset", "group", "n"])

pairs: dict[str, set[str]] = defaultdict(set)
pair_samples: dict[str, list[str]] = defaultdict(list)
for row in rows:
    if row["dataset"] == "GSE70947" and row["pair_key"]:
        pairs[row["pair_key"]].add(row["group"])
        pair_samples[row["pair_key"]].append(row["sample"])

valid_pair_keys = {
    key
    for key, groups in pairs.items()
    if groups == {"Tumor", "Normal"} and len(pair_samples[key]) == 2
}
pair_rows = [
    {
        "pair_key": key,
        "groups": ";".join(sorted(pairs[key])),
        "samples": ";".join(pair_samples[key]),
        "valid_pair": key in valid_pair_keys,
    }
    for key in sorted(pairs)
]
write_tsv(
    RUN_DIR / "gse70947_pair_audit.tsv",
    pair_rows,
    ["pair_key", "groups", "samples", "valid_pair"],
)

observed = {
    "GSE70947_Tumor": counts[("GSE70947", "Tumor")],
    "GSE70947_Normal": counts[("GSE70947", "Normal")],
    "GSE70947_valid_pairs": len(valid_pair_keys),
    "GSE65194_TNBC_tissue": counts[("GSE65194", "TNBC")],
    "GSE65194_Non_TNBC_tissue": counts[("GSE65194", "Non_TNBC")],
    "GSE65194_Normal_tissue": counts[("GSE65194", "Normal")],
    "GSE65194_Cell_line": counts[("GSE65194", "Cell_line")],
    "GSE65194_Unresolved": counts[("GSE65194", "Unresolved")],
}
expected = {
    "GSE70947_Tumor": 148,
    "GSE70947_Normal": 148,
    "GSE70947_valid_pairs": 148,
    "GSE65194_TNBC_tissue": 41,
    "GSE65194_Non_TNBC_tissue": 89,
    "GSE65194_Normal_tissue": 11,
    "GSE65194_Cell_line": 37,
    "GSE65194_Unresolved": 0,
}

summary = [f"SOURCE_AUDIT_RUN={SOURCE_TSV.parent}"]
all_pass = True
for metric, expected_value in expected.items():
    value = observed[metric]
    status = "PASS" if value == expected_value else "FAIL"
    all_pass = all_pass and status == "PASS"
    summary.append(f"{metric}={value};EXPECTED={expected_value};STATUS={status}")
summary.append("RAW_DATA_MODIFIED=NO")
summary.append("GEO_GROUP_REAUDIT_PASS" if all_pass else "GEO_GROUP_REAUDIT_FAIL")

(RUN_DIR / "geo_reaudit_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))
if not all_pass:
    print("MISMATCHED_EXAMPLES_FILE=geo_curated_sample_sheet.tsv")
    raise SystemExit(2)
PY

export SOURCE_TSV RUN_DIR

{
  echo "RUN_DIR=${RUN_DIR}"
  echo "SOURCE_TSV=${SOURCE_TSV}"
  echo "STEP=SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
  echo "SYNTAX_CHECK=PASS"
  echo "STEP=GEO_METADATA_REAUDIT"
  conda run --no-capture-output -n batch_bench python "${PY_SCRIPT}"
  tail -n 1 "${RUN_DIR}/geo_reaudit_summary.txt" | grep -qx 'GEO_GROUP_REAUDIT_PASS'
  printf 'LATEST_GEO_REAUDIT=%s\n' "${RUN_DIR}" > "${AUDIT_ROOT}/latest_geo_reaudit.txt"
  echo "ANALYSIS_LOGIC_PASS"
} 2>&1 | tee "${LOG_FILE}"
