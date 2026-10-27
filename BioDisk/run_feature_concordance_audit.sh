#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
READY_ROOT="${PROJECT_ROOT}/02_analysis_ready"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
FEATURE_ROOT="${PROJECT_ROOT}/03_feature_audit"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${FEATURE_ROOT}/concordance_${STAMP}"
PY_SCRIPT="${RUN_DIR}/feature_concordance_audit.py"
LOG_FILE="${RUN_DIR}/feature_concordance_audit.log"

pointer_value() {
  local pointer_file="$1"
  local pointer_key="$2"
  if [[ ! -s "${pointer_file}" ]]; then
    echo "ERROR=pointer_missing:${pointer_file}" >&2
    return 1
  fi
  sed -n "s/^${pointer_key}=//p" "${pointer_file}" | tail -n 1
}

READY_RUN="$(pointer_value "${READY_ROOT}/latest_analysis_ready_freeze.txt" "LATEST_ANALYSIS_READY_FREEZE")"
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
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(os.environ["PROJECT_ROOT"])
RUN_DIR = Path(os.environ["RUN_DIR"])
READY_RUN = Path(os.environ["READY_RUN"])
CANDIDATE_AUDIT_RUN = Path(os.environ["CANDIDATE_AUDIT_RUN"])
MIN_CORR_N = int(os.environ.get("MIN_CORR_N", "20"))

EXPECTED = {
    "ACOX2_PRIMARY_SAMPLES": int(os.environ.get("EXP_ACOX2_PRIMARY_SAMPLES", "783")),
    "ACOX2_PROMOTER_PROBES": int(os.environ.get("EXP_ACOX2_PROMOTER_PROBES", "3")),
    "GSE70947_INCLUDED": int(os.environ.get("EXP_GSE70947_INCLUDED", "296")),
    "GSE65194_INCLUDED": int(os.environ.get("EXP_GSE65194_INCLUDED", "164")),
    "GSE70947_GENES": int(os.environ.get("EXP_GSE70947_GENES", "7")),
    "GSE65194_GENES": int(os.environ.get("EXP_GSE65194_GENES", "8")),
}

PROMOTER_MATRIX = READY_RUN / "tcga_methylation_acox2_promoter_3probes.tsv.gz"
PRIMARY_SAMPLES = READY_RUN / "tcga_joint_primary_samples.tsv"
SAMPLE_MANIFEST = READY_RUN / "tcga_sample_manifest.tsv"
GEO_AVAILABILITY = CANDIDATE_AUDIT_RUN / "geo_candidate_gene_availability.tsv"

DATASETS = {
    "GSE70947": {
        "matrix": READY_RUN / "GSE70947_candidate_probe_matrix.tsv.gz",
        "samples": READY_RUN / "GSE70947_included_sample_sheet.tsv",
    },
    "GSE65194": {
        "matrix": READY_RUN / "GSE65194_candidate_probe_matrix.tsv.gz",
        "samples": READY_RUN / "GSE65194_included_sample_sheet.tsv",
    },
}


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    kwargs = {"mode": "rt", "encoding": "utf-8-sig", "newline": ""} if path.suffix == ".gz" else {
        "mode": "r",
        "encoding": "utf-8-sig",
        "newline": "",
    }
    with opener(path, **kwargs) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    if frame.shape[1] < 2 or frame.columns[0] != "feature_id":
        raise RuntimeError(f"Invalid frozen matrix: {path}")
    if frame["feature_id"].duplicated().any():
        raise RuntimeError(f"Duplicated feature identifiers: {path}")
    frame = frame.set_index("feature_id")
    return frame.apply(pd.to_numeric, errors="coerce")


def write_matrix(path: Path, frame: pd.DataFrame) -> None:
    output = frame.copy()
    output.index.name = "feature_id"
    output.to_csv(path, sep="\t", compression="gzip", na_rep="NA")


def finite_values(series: pd.Series) -> np.ndarray:
    values = series.to_numpy(dtype=float)
    return values[np.isfinite(values)]


def probe_qc_rows(dataset: str, gene: str, matrix: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for probe, values in matrix.iterrows():
        finite = finite_values(values)
        n_valid = int(finite.size)
        n_total = int(values.size)
        if n_valid:
            q25, median, q75 = np.quantile(finite, [0.25, 0.50, 0.75])
            mad = float(np.median(np.abs(finite - median)))
            minimum = float(np.min(finite))
            maximum = float(np.max(finite))
            variance = float(np.var(finite, ddof=1)) if n_valid > 1 else 0.0
        else:
            q25 = median = q75 = mad = minimum = maximum = variance = math.nan
        rows.append(
            {
                "dataset": dataset,
                "gene": gene,
                "probe_id": probe,
                "n_total": n_total,
                "n_valid": n_valid,
                "missing_rate": 1.0 - n_valid / n_total if n_total else math.nan,
                "min": minimum,
                "q25": q25,
                "median": median,
                "q75": q75,
                "max": maximum,
                "IQR": q75 - q25 if n_valid else math.nan,
                "MAD": mad,
                "variance": variance,
            }
        )
    return rows


def pairwise_spearman(dataset: str, gene: str, matrix: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    correlation = matrix.T.corr(method="spearman", min_periods=MIN_CORR_N)
    rows: list[dict] = []
    for probe_a, probe_b in combinations(matrix.index, 2):
        pair = matrix.loc[[probe_a, probe_b]].T.dropna()
        rho = correlation.loc[probe_a, probe_b]
        rows.append(
            {
                "dataset": dataset,
                "gene": gene,
                "probe_a": probe_a,
                "probe_b": probe_b,
                "n_complete": pair.shape[0],
                "spearman_rho": rho,
            }
        )
    return correlation, rows


def safe_median(values: list[float]) -> float:
    finite = [value for value in values if np.isfinite(value)]
    return float(np.median(finite)) if finite else math.nan


def select_representative(qc_rows: list[dict]) -> str:
    eligible = [row for row in qc_rows if row["missing_rate"] <= 0.20 and row["n_valid"] >= MIN_CORR_N]
    if not eligible:
        raise RuntimeError("No probe passed missingness and sample-count requirements")
    eligible.sort(key=lambda row: (-float(row["IQR"]), -float(row["MAD"]), row["probe_id"]))
    return str(eligible[0]["probe_id"])


def zmedian(matrix: pd.DataFrame) -> pd.Series:
    centered = matrix.sub(matrix.mean(axis=1), axis=0)
    scale = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = centered.div(scale, axis=0)
    return z.median(axis=0, skipna=True)


def cronbach_alpha(matrix: pd.DataFrame) -> float:
    complete = matrix.T.dropna()
    k = complete.shape[1]
    if k < 2 or complete.shape[0] < MIN_CORR_N:
        return math.nan
    item_variance = complete.var(axis=0, ddof=1).sum()
    total_variance = complete.sum(axis=1).var(ddof=1)
    if not np.isfinite(total_variance) or total_variance <= 0:
        return math.nan
    return float(k / (k - 1) * (1 - item_variance / total_variance))


def promoter_pca(matrix: pd.DataFrame) -> tuple[float, dict[str, float], int]:
    complete = matrix.T.dropna()
    if complete.shape[0] < MIN_CORR_N:
        return math.nan, {probe: math.nan for probe in matrix.index}, complete.shape[0]
    standard_deviation = complete.std(axis=0, ddof=1)
    if (standard_deviation <= 0).any():
        return math.nan, {probe: math.nan for probe in matrix.index}, complete.shape[0]
    scaled = (complete - complete.mean(axis=0)) / standard_deviation
    _, singular_values, right_vectors = np.linalg.svd(scaled.to_numpy(), full_matrices=False)
    eigenvalues = singular_values**2
    explained = float(eigenvalues[0] / eigenvalues.sum())
    loadings = {probe: float(right_vectors[0, index]) for index, probe in enumerate(complete.columns)}
    return explained, loadings, complete.shape[0]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


required = [PROMOTER_MATRIX, PRIMARY_SAMPLES, SAMPLE_MANIFEST, GEO_AVAILABILITY]
for values in DATASETS.values():
    required.extend([values["matrix"], values["samples"]])
for path in required:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Required input missing or empty: {path}")

print("PY_STEP=ACOX2_PROMOTER_CONCORDANCE", flush=True)
promoter_all = read_matrix(PROMOTER_MATRIX)
primary_ids = [row["sample_id"] for row in read_rows(PRIMARY_SAMPLES)]
missing_primary = sorted(set(primary_ids) - set(promoter_all.columns))
if missing_primary:
    raise RuntimeError(f"Primary samples absent from promoter matrix: {len(missing_primary)}")
promoter_primary = promoter_all.loc[:, primary_ids]

promoter_qc = probe_qc_rows("TCGA_BRCA_PRIMARY", "ACOX2", promoter_primary)
promoter_corr, promoter_pairs = pairwise_spearman(
    "TCGA_BRCA_PRIMARY", "ACOX2", promoter_primary
)
promoter_rhos = [float(row["spearman_rho"]) for row in promoter_pairs if pd.notna(row["spearman_rho"])]
promoter_median_rho = safe_median(promoter_rhos)
promoter_min_rho = min(promoter_rhos) if promoter_rhos else math.nan
promoter_alpha = cronbach_alpha(promoter_primary)
pc1_explained, pc1_loadings, pca_complete_n = promoter_pca(promoter_primary)
loading_values = [value for value in pc1_loadings.values() if np.isfinite(value)]
loadings_same_sign = bool(loading_values) and (
    all(value > 0 for value in loading_values) or all(value < 0 for value in loading_values)
)
max_missing = max(float(row["missing_rate"]) for row in promoter_qc)

promoter_composite_ok = (
    max_missing <= 0.05
    and len(promoter_rhos) == 3
    and promoter_min_rho >= 0.30
    and promoter_median_rho >= 0.50
    and np.isfinite(pc1_explained)
    and pc1_explained >= 0.60
    and loadings_same_sign
)
promoter_decision = "USE_MEAN_BETA_PRIMARY_PLUS_PER_PROBE" if promoter_composite_ok else "USE_PER_PROBE_PRIMARY_ONLY"

write_rows(
    RUN_DIR / "acox2_promoter_probe_qc.tsv",
    promoter_qc,
    ["dataset", "gene", "probe_id", "n_total", "n_valid", "missing_rate", "min", "q25", "median", "q75", "max", "IQR", "MAD", "variance"],
)
write_rows(
    RUN_DIR / "acox2_promoter_pairwise_spearman.tsv",
    promoter_pairs,
    ["dataset", "gene", "probe_a", "probe_b", "n_complete", "spearman_rho"],
)
loadings_rows = [
    {"probe_id": probe, "PC1_loading": loading}
    for probe, loading in pc1_loadings.items()
]
write_rows(RUN_DIR / "acox2_promoter_pc1_loadings.tsv", loadings_rows, ["probe_id", "PC1_loading"])

valid_counts = promoter_all.notna().sum(axis=0)
mean_beta = promoter_all.mean(axis=0, skipna=True).where(valid_counts >= 2)
composite_rows = [
    {
        "sample_id": sample,
        "promoter_mean_beta": mean_beta.loc[sample],
        "valid_probe_count": int(valid_counts.loc[sample]),
        "recommended_primary_use": promoter_composite_ok,
    }
    for sample in promoter_all.columns
]
write_rows(
    RUN_DIR / "acox2_promoter_candidate_composite.tsv",
    composite_rows,
    ["sample_id", "promoter_mean_beta", "valid_probe_count", "recommended_primary_use"],
)

print("PY_STEP=GEO_PROBE_CONCORDANCE_AND_COLLAPSE", flush=True)
availability = read_rows(GEO_AVAILABILITY)
all_geo_qc: list[dict] = []
all_geo_pairs: list[dict] = []
all_geo_decisions: list[dict] = []
observed_geo_samples: dict[str, int] = {}
observed_geo_genes: dict[str, int] = {}

for dataset, paths in DATASETS.items():
    matrix = read_matrix(paths["matrix"])
    sample_rows = read_rows(paths["samples"])
    included_samples = [row["sample"] for row in sample_rows]
    absent = sorted(set(included_samples) - set(matrix.columns))
    if absent:
        raise RuntimeError(f"{dataset} included samples absent from matrix: {len(absent)}")
    matrix = matrix.loc[:, included_samples]
    mapping_rows = [
        row for row in availability if row["dataset"] == dataset and row["availability"] == "PRESENT"
    ]
    primary_gene_rows: dict[str, pd.Series] = {}
    sensitivity_rows: dict[str, pd.Series] = {}
    for mapping in mapping_rows:
        gene = mapping["gene"]
        probes = [probe for probe in mapping["present_probes"].split(";") if probe and probe != "NA"]
        missing_probes = sorted(set(probes) - set(matrix.index))
        if missing_probes:
            raise RuntimeError(f"{dataset} {gene} missing frozen probes: {missing_probes}")
        gene_matrix = matrix.loc[probes]
        qc_rows = probe_qc_rows(dataset, gene, gene_matrix)
        all_geo_qc.extend(qc_rows)
        _, pair_rows = pairwise_spearman(dataset, gene, gene_matrix)
        all_geo_pairs.extend(pair_rows)
        rhos = [float(row["spearman_rho"]) for row in pair_rows if pd.notna(row["spearman_rho"])]
        median_rho = safe_median(rhos)
        min_rho = min(rhos) if rhos else math.nan
        representative = select_representative(qc_rows)
        if len(probes) == 1:
            method = "SINGLE_PROBE"
            primary = gene_matrix.loc[representative]
        elif len(rhos) == len(pair_rows) and min_rho >= 0.50 and median_rho >= 0.70:
            method = "MEAN_RAW_CONCORDANT_PROBES"
            primary = gene_matrix.mean(axis=0, skipna=True)
        else:
            method = "MAX_IQR_REPRESENTATIVE_PROBE"
            primary = gene_matrix.loc[representative]
        primary_gene_rows[gene] = primary
        sensitivity_rows[gene] = zmedian(gene_matrix)
        all_geo_decisions.append(
            {
                "dataset": dataset,
                "gene": gene,
                "probe_count": len(probes),
                "valid_pairwise_correlations": len(rhos),
                "median_pairwise_rho": median_rho,
                "minimum_pairwise_rho": min_rho,
                "primary_method": method,
                "representative_probe": representative,
                "sensitivity_method": "MEDIAN_PROBE_ZSCORE",
            }
        )

    primary_frame = pd.DataFrame(primary_gene_rows).T.loc[:, included_samples]
    sensitivity_frame = pd.DataFrame(sensitivity_rows).T.loc[:, included_samples]
    write_matrix(RUN_DIR / f"{dataset}_gene_expression_primary.tsv.gz", primary_frame)
    write_matrix(RUN_DIR / f"{dataset}_gene_expression_sensitivity_zmedian.tsv.gz", sensitivity_frame)
    observed_geo_samples[dataset] = len(included_samples)
    observed_geo_genes[dataset] = primary_frame.shape[0]

write_rows(
    RUN_DIR / "geo_probe_qc.tsv",
    all_geo_qc,
    ["dataset", "gene", "probe_id", "n_total", "n_valid", "missing_rate", "min", "q25", "median", "q75", "max", "IQR", "MAD", "variance"],
)
write_rows(
    RUN_DIR / "geo_probe_pairwise_spearman.tsv",
    all_geo_pairs,
    ["dataset", "gene", "probe_a", "probe_b", "n_complete", "spearman_rho"],
)
write_rows(
    RUN_DIR / "geo_probe_collapse_decisions.tsv",
    all_geo_decisions,
    ["dataset", "gene", "probe_count", "valid_pairwise_correlations", "median_pairwise_rho", "minimum_pairwise_rho", "primary_method", "representative_probe", "sensitivity_method"],
)

observed = {
    "ACOX2_PRIMARY_SAMPLES": promoter_primary.shape[1],
    "ACOX2_PROMOTER_PROBES": promoter_primary.shape[0],
    "GSE70947_INCLUDED": observed_geo_samples.get("GSE70947", 0),
    "GSE65194_INCLUDED": observed_geo_samples.get("GSE65194", 0),
    "GSE70947_GENES": observed_geo_genes.get("GSE70947", 0),
    "GSE65194_GENES": observed_geo_genes.get("GSE65194", 0),
}

summary = [
    "AUDIT_SCOPE=ACOX2_PROMOTER_AND_GEO_PROBE_CONCORDANCE",
    f"MIN_CORRELATION_COMPLETE_N={MIN_CORR_N}",
]
all_pass = True
for metric, expected in EXPECTED.items():
    value = observed[metric]
    status = "PASS" if value == expected else "FAIL"
    all_pass = all_pass and status == "PASS"
    summary.append(f"{metric}={value};EXPECTED={expected};STATUS={status}")
summary.extend(
    [
        f"ACOX2_MAX_PROBE_MISSING_RATE={max_missing:.8g}",
        f"ACOX2_MIN_PAIRWISE_RHO={promoter_min_rho:.8g}",
        f"ACOX2_MEDIAN_PAIRWISE_RHO={promoter_median_rho:.8g}",
        f"ACOX2_CRONBACH_ALPHA={promoter_alpha:.8g}",
        f"ACOX2_PC1_EXPLAINED_VARIANCE={pc1_explained:.8g}",
        f"ACOX2_PCA_COMPLETE_SAMPLES={pca_complete_n}",
        f"ACOX2_PC1_LOADINGS_SAME_SIGN={str(loadings_same_sign).upper()}",
        f"ACOX2_PROMOTER_COMPOSITE_DECISION={promoter_decision}",
        "ACOX2_COMPOSITE_RULE=MISS_LE_0.05;MIN_RHO_GE_0.30;MEDIAN_RHO_GE_0.50;PC1_GE_0.60;LOADINGS_SAME_SIGN",
    ]
)
for row in all_geo_decisions:
    summary.append(
        "GEO_DECISION="
        + ";".join(
            [
                row["dataset"],
                row["gene"],
                str(row["probe_count"]),
                row["primary_method"],
                row["representative_probe"],
            ]
        )
    )
summary.extend(
    [
        "RAW_AND_FROZEN_INPUTS_MODIFIED=NO",
        "FEATURE_COLLAPSE_AUDIT_PASS" if all_pass else "FEATURE_COLLAPSE_AUDIT_FAIL",
    ]
)
(RUN_DIR / "feature_concordance_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

output_manifest = []
for path in sorted(RUN_DIR.iterdir()):
    if path.is_file() and path.name not in {"feature_concordance_audit.log", "output_file_manifest.tsv"}:
        output_manifest.append(
            {"file": path.name, "size_bytes": path.stat().st_size, "sha256": file_sha256(path)}
        )
write_rows(
    RUN_DIR / "output_file_manifest.tsv",
    output_manifest,
    ["file", "size_bytes", "sha256"],
)

print("\n".join(summary))
if not all_pass:
    raise SystemExit(2)
PY

export PROJECT_ROOT RUN_DIR READY_RUN CANDIDATE_AUDIT_RUN

{
  echo "RUN_DIR=${RUN_DIR}"
  echo "READY_RUN=${READY_RUN}"
  echo "CANDIDATE_AUDIT_RUN=${CANDIDATE_AUDIT_RUN}"
  echo "STEP=PYTHON_DEPENDENCY_CHECK"
  conda run -n batch_bench python -c 'import numpy, pandas'
  echo "DEPENDENCY_CHECK=PASS"
  echo "STEP=PYTHON_SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
  echo "SYNTAX_CHECK=PASS"
  echo "STEP=FEATURE_CONCORDANCE_AUDIT"
  conda run --no-capture-output -n batch_bench python "${PY_SCRIPT}"
  tail -n 1 "${RUN_DIR}/feature_concordance_summary.txt" | grep -qx 'FEATURE_COLLAPSE_AUDIT_PASS'
  printf 'LATEST_FEATURE_CONCORDANCE_AUDIT=%s\n' "${RUN_DIR}" > "${FEATURE_ROOT}/latest_feature_concordance_audit.txt"
  echo "ANALYSIS_LOGIC_PASS"
} 2>&1 | tee "${LOG_FILE}"
