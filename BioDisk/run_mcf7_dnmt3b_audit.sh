#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
RESULT_ROOT="${PROJECT_ROOT}/04_results/4_1_expression_clinical"
RAW_MATRIX="${PROJECT_ROOT}/00_raw_data/GEO/GSE65194/GSE65194_series_matrix.txt.gz"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${RESULT_ROOT}/mcf7_audit_${STAMP}"
PY_SCRIPT="${RUN_DIR}/mcf7_dnmt3b_audit.py"
LOG_FILE="${RUN_DIR}/mcf7_dnmt3b_audit.log"

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
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


RUN_DIR = Path(os.environ["RUN_DIR"])
RAW_MATRIX = Path(os.environ["RAW_MATRIX"])
GEO_FREEZE_RUN = Path(os.environ["GEO_FREEZE_RUN"])
CANDIDATE_AUDIT_RUN = Path(os.environ["CANDIDATE_AUDIT_RUN"])
EXPECTED_CELL_LINES = int(os.environ.get("EXP_CELL_LINES", "14"))

SAMPLE_SHEET = GEO_FREEZE_RUN / "geo_frozen_sample_sheet.tsv"
AVAILABILITY = CANDIDATE_AUDIT_RUN / "geo_candidate_gene_availability.tsv"
MCF7_PATTERN = re.compile(r"(?<![A-Z0-9])MCF[-_ ]?7(?![A-Z0-9])", re.IGNORECASE)


def strip_quotes(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        return text[1:-1]
    return text


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_target_probes(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    rows = frame.loc[
        (frame["dataset"] == "GSE65194")
        & (frame["gene"] == "DNMT3B")
        & (frame["availability"] == "PRESENT")
    ]
    if rows.shape[0] != 1:
        raise RuntimeError(f"Expected one GSE65194 DNMT3B availability row; found {rows.shape[0]}")
    return [probe for probe in rows.iloc[0]["present_probes"].split(";") if probe and probe != "NA"]


def extract_probes(path: Path, targets: set[str]) -> tuple[list[str], dict[str, list[float]]]:
    samples: list[str] = []
    rows: dict[str, list[float]] = {}
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            fields = [strip_quotes(item) for item in line.split("\t")]
            if not samples:
                samples = fields[1:]
                continue
            probe = fields[0]
            if probe not in targets:
                continue
            values = [float(value) if value.lower() not in {"", "na", "nan", "null", "--"} else math.nan for value in fields[1:]]
            if len(values) != len(samples):
                raise RuntimeError(f"Value count mismatch for probe {probe}")
            rows[probe] = values
    return samples, rows


for path in [RAW_MATRIX, SAMPLE_SHEET, AVAILABILITY]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Required input missing or empty: {path}")

print("PY_STEP=LOCATE_DNMT3B_PROBE", flush=True)
probes = read_target_probes(AVAILABILITY)
samples, probe_rows = extract_probes(RAW_MATRIX, set(probes))
if set(probe_rows) != set(probes):
    raise RuntimeError(f"DNMT3B probe mismatch: expected={probes}; found={sorted(probe_rows)}")

matrix = pd.DataFrame(probe_rows, index=samples)
if len(probes) == 1:
    expression = matrix[probes[0]]
    aggregation = "SINGLE_PROBE"
else:
    expression = matrix.mean(axis=1, skipna=True)
    aggregation = "MEAN_OF_AVAILABLE_PROBES"

print("PY_STEP=IDENTIFY_CELL_LINES_AND_MCF7", flush=True)
metadata = pd.read_csv(SAMPLE_SHEET, sep="\t", dtype=str, keep_default_na=False)
cell_lines = metadata.loc[
    (metadata["dataset"] == "GSE65194")
    & ((metadata["material"] == "Cell_line") | (metadata["analysis_decision"] == "exclude_cell_line"))
].copy()
cell_lines["expression"] = cell_lines["sample"].map(expression)
cell_lines["search_text"] = (
    cell_lines["title"].fillna("") + " | " + cell_lines["combined_metadata"].fillna("")
)
cell_lines["is_MCF7"] = cell_lines["search_text"].map(lambda text: bool(MCF7_PATTERN.search(text)))

if cell_lines["sample"].duplicated().any():
    raise RuntimeError("Duplicated GSE65194 cell-line sample identifiers")
if cell_lines["expression"].isna().any():
    missing = cell_lines.loc[cell_lines["expression"].isna(), "sample"].tolist()
    raise RuntimeError(f"Missing DNMT3B expression for cell-line samples: {missing}")

cell_lines["rank_desc"] = cell_lines["expression"].rank(method="min", ascending=False).astype(int)
cell_lines["percentile"] = cell_lines["expression"].rank(method="average", pct=True) * 100
mean_value = float(cell_lines["expression"].mean())
sd_value = float(cell_lines["expression"].std(ddof=1))
cell_lines["z_score"] = (cell_lines["expression"] - mean_value) / sd_value if sd_value > 0 else math.nan
cell_lines = cell_lines.sort_values(["rank_desc", "sample"]).reset_index(drop=True)
cell_lines["plot_label"] = [f"CL{index:02d}" for index in range(1, len(cell_lines) + 1)]
mcf7_indices = cell_lines.index[cell_lines["is_MCF7"]].tolist()
for sequence, index in enumerate(mcf7_indices, start=1):
    cell_lines.loc[index, "plot_label"] = "MCF7" if len(mcf7_indices) == 1 else f"MCF7-{sequence}"

mcf7 = cell_lines.loc[cell_lines["is_MCF7"]].copy()
if mcf7.empty:
    decision = "MCF7_NOT_IDENTIFIED"
    mcf7_median = mcf7_percentile = mcf7_rank = math.nan
elif float(mcf7["percentile"].median()) >= 75:
    decision = "MCF7_HIGH_WITHIN_CELL_LINES"
    mcf7_median = float(mcf7["expression"].median())
    mcf7_percentile = float(mcf7["percentile"].median())
    mcf7_rank = int(mcf7["rank_desc"].min())
elif float(mcf7["percentile"].median()) <= 25:
    decision = "MCF7_LOW_WITHIN_CELL_LINES"
    mcf7_median = float(mcf7["expression"].median())
    mcf7_percentile = float(mcf7["percentile"].median())
    mcf7_rank = int(mcf7["rank_desc"].min())
else:
    decision = "MCF7_INTERMEDIATE_WITHIN_CELL_LINES"
    mcf7_median = float(mcf7["expression"].median())
    mcf7_percentile = float(mcf7["percentile"].median())
    mcf7_rank = int(mcf7["rank_desc"].min())

output_fields = [
    "plot_label",
    "sample",
    "title",
    "is_MCF7",
    "expression",
    "rank_desc",
    "percentile",
    "z_score",
    "combined_metadata",
]
cell_lines.to_csv(
    RUN_DIR / "gse65194_cell_line_dnmt3b.tsv",
    sep="\t",
    index=False,
    columns=output_fields,
)
cell_lines[["plot_label", "sample", "title", "is_MCF7"]].to_csv(
    RUN_DIR / "gse65194_cell_line_label_map.tsv", sep="\t", index=False
)

print("PY_STEP=GENERATE_CELL_LINE_FIGURE", flush=True)
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)
plot_data = cell_lines.sort_values("expression", ascending=True)
colors = ["#D62728" if value else "#6B8EAD" for value in plot_data["is_MCF7"]]
fig, ax = plt.subplots(figsize=(7.0, 5.6), constrained_layout=True)
ax.barh(plot_data["plot_label"], plot_data["expression"], color=colors, edgecolor="#444444", linewidth=0.5)
ax.axvline(float(cell_lines["expression"].median()), color="#333333", linestyle="--", linewidth=1.0)
ax.set_xlabel("DNMT3B")
ax.set_ylabel("")
ax.set_title("Cell lines", fontsize=12, fontweight="bold")
ax.grid(axis="x", color="#DDDDDD", linewidth=0.7)
ax.set_axisbelow(True)
ax.legend(
    handles=[Patch(facecolor="#D62728", label="MCF7"), Patch(facecolor="#6B8EAD", label="Other")],
    frameon=False,
    loc="lower right",
)
fig.savefig(RUN_DIR / "Fig4_1_MCF7_DNMT3B.pdf", bbox_inches="tight")
fig.savefig(RUN_DIR / "Fig4_1_MCF7_DNMT3B.png", dpi=400, bbox_inches="tight")
plt.close(fig)

caption = """Figure 4.1 supplementary cell-line audit. DNMT3B expression among the 14 GSE65194 cell-line samples. MCF7 samples are highlighted in red; other cell-line samples are shown in blue. The dashed line indicates the median across all cell-line samples. Short labels CL01-CL14 are mapped to original sample titles in gse65194_cell_line_label_map.tsv. Expression was derived from the frozen GSE65194 DNMT3B probe. This comparison is descriptive because individual cell lines are not biological replicates. Abbreviations: CL, cell line; MCF7, Michigan Cancer Foundation-7 breast cancer cell line.
"""
(RUN_DIR / "Fig4_1_MCF7_caption.txt").write_text(caption, encoding="utf-8")

count_pass = len(cell_lines) == EXPECTED_CELL_LINES
probe_pass = len(probes) >= 1
summary = [
    "AUDIT_SCOPE=GSE65194_MCF7_DNMT3B_EXPRESSION",
    f"DNMT3B_PROBES={';'.join(probes)}",
    f"PROBE_AGGREGATION={aggregation}",
    f"CELL_LINE_SAMPLES={len(cell_lines)};EXPECTED={EXPECTED_CELL_LINES};STATUS={'PASS' if count_pass else 'FAIL'}",
    f"MCF7_SAMPLES={len(mcf7)}",
    f"CELL_LINE_MEDIAN={float(cell_lines['expression'].median()):.8g}",
    f"MCF7_MEDIAN_EXPRESSION={mcf7_median:.8g}",
    f"MCF7_MEDIAN_PERCENTILE={mcf7_percentile:.8g}",
    f"MCF7_BEST_RANK={mcf7_rank}",
    f"MCF7_EXPRESSION_DECISION={decision}",
    "HIGH_DEFINITION=PERCENTILE_GE_75_WITHIN_CELL_LINE_PANEL",
    "INFERENCE=DESCRIPTIVE_NOT_REPLICATE_BASED",
    "RAW_AND_FROZEN_INPUTS_MODIFIED=NO",
    "MCF7_TARGETED_AUDIT_PASS" if count_pass and probe_pass else "MCF7_TARGETED_AUDIT_FAIL",
]
(RUN_DIR / "mcf7_dnmt3b_audit_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

manifest_rows = []
for path in sorted(RUN_DIR.iterdir()):
    if path.is_file() and path.name not in {"mcf7_dnmt3b_audit.log", "output_file_manifest.tsv"}:
        manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
with (RUN_DIR / "output_file_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "sha256"], delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

print("\n".join(summary))
if not count_pass or not probe_pass:
    raise SystemExit(2)
PY

export RUN_DIR RAW_MATRIX GEO_FREEZE_RUN CANDIDATE_AUDIT_RUN

{
  echo "RUN_DIR=${RUN_DIR}"
  echo "GEO_FREEZE_RUN=${GEO_FREEZE_RUN}"
  echo "CANDIDATE_AUDIT_RUN=${CANDIDATE_AUDIT_RUN}"
  echo "STEP=PYTHON_DEPENDENCY_CHECK"
  conda run -n batch_bench python -c 'import matplotlib, numpy, pandas'
  echo "DEPENDENCY_CHECK=PASS"
  echo "STEP=PYTHON_SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
  echo "SYNTAX_CHECK=PASS"
  echo "STEP=MCF7_TARGETED_AUDIT"
  conda run --no-capture-output -n batch_bench python "${PY_SCRIPT}"
  tail -n 1 "${RUN_DIR}/mcf7_dnmt3b_audit_summary.txt" | grep -qx 'MCF7_TARGETED_AUDIT_PASS'
  printf 'LATEST_MCF7_AUDIT=%s\n' "${RUN_DIR}" > "${RESULT_ROOT}/latest_mcf7_audit.txt"
  echo "ANALYSIS_LOGIC_PASS"
} 2>&1 | tee "${LOG_FILE}"
