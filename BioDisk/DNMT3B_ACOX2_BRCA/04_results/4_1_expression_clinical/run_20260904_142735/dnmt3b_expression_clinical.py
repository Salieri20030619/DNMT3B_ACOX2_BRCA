#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import math
import os
import re
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


RUN_DIR = Path(os.environ["RUN_DIR"])
READY_RUN = Path(os.environ["READY_RUN"])
FEATURE_RUN = Path(os.environ["FEATURE_RUN"])

TCGA_MATRIX = READY_RUN / "tcga_expression_8genes.tsv.gz"
TCGA_META = READY_RUN / "tcga_sample_manifest.tsv"
GSE70947_MATRIX = FEATURE_RUN / "GSE70947_gene_expression_primary.tsv.gz"
GSE65194_MATRIX = FEATURE_RUN / "GSE65194_gene_expression_primary.tsv.gz"
GSE70947_META = READY_RUN / "GSE70947_included_sample_sheet.tsv"
GSE65194_META = READY_RUN / "GSE65194_included_sample_sheet.tsv"

EXPECTED = {
    "TCGA_PRIMARY": int(os.environ.get("EXP_TCGA_PRIMARY", "1097")),
    "TCGA_NORMAL": int(os.environ.get("EXP_TCGA_NORMAL", "114")),
    "TCGA_PAIRS": int(os.environ.get("EXP_TCGA_PAIRS", "114")),
    "PAM50_FIELD_VALID": int(os.environ.get("EXP_PAM50_FIELD_VALID", "522")),
    "STAGE_FIELD_VALID": int(os.environ.get("EXP_STAGE_FIELD_VALID", "779")),
    "ER_FIELD_VALID": int(os.environ.get("EXP_ER_FIELD_VALID", "782")),
    "PR_FIELD_VALID": int(os.environ.get("EXP_PR_FIELD_VALID", "781")),
    "HER2_FIELD_VALID": int(os.environ.get("EXP_HER2_FIELD_VALID", "776")),
    "GSE70947_PAIRS": int(os.environ.get("EXP_GSE70947_PAIRS", "148")),
    "GSE65194_INCLUDED": int(os.environ.get("EXP_GSE65194_INCLUDED", "164")),
    "GSE65194_NORMAL": int(os.environ.get("EXP_GSE65194_NORMAL", "11")),
    "GSE65194_TNBC": int(os.environ.get("EXP_GSE65194_TNBC", "55")),
    "GSE65194_HER2": int(os.environ.get("EXP_GSE65194_HER2", "39")),
    "GSE65194_LUMA": int(os.environ.get("EXP_GSE65194_LUMA", "29")),
    "GSE65194_LUMB": int(os.environ.get("EXP_GSE65194_LUMB", "30")),
}

PALETTE = {
    "Normal": "#4C78A8",
    "Primary": "#E45756",
    "Basal": "#D62728",
    "HER2": "#E377C2",
    "LumA": "#59A14F",
    "LumB": "#F28E2B",
    "Stg I": "#9ECAE1",
    "Stg II": "#6BAED6",
    "Stg III": "#3182BD",
    "Stg IV": "#08519C",
    "Negative": "#4C78A8",
    "Positive": "#E45756",
    "TNBC": "#D62728",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gene(path: Path, gene: str) -> pd.Series:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if frame.columns[0] != "feature_id":
        raise RuntimeError(f"Invalid matrix header: {path}")
    hits = frame.loc[frame["feature_id"].astype(str).str.upper() == gene.upper()]
    if hits.shape[0] != 1:
        raise RuntimeError(f"Expected one {gene} row in {path}; found {hits.shape[0]}")
    values = pd.to_numeric(hits.iloc[0, 1:], errors="coerce")
    values.index = frame.columns[1:]
    values.name = gene
    return values


def clean_text(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "null", "--", "[not available]", "[not applicable]"}:
        return None
    return text


def normalize_pam50(value: object) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    compact = re.sub(r"[^a-z0-9]", "", text.lower())
    if "basal" in compact:
        return "Basal"
    if "her2" in compact:
        return "HER2"
    if compact in {"luma", "luminala"}:
        return "LumA"
    if compact in {"lumb", "luminalb"}:
        return "LumB"
    if "normal" in compact:
        return "Normal"
    return None


def normalize_stage(value: object) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    upper = text.upper().replace("STAGE", "").strip()
    if re.search(r"\bIV\b|^4", upper):
        return "Stg IV"
    if re.search(r"\bIII[A-C]?\b|^3", upper):
        return "Stg III"
    if re.search(r"\bII[A-C]?\b|^2", upper):
        return "Stg II"
    if re.search(r"\bI[A-C]?\b|^1", upper):
        return "Stg I"
    return None


def normalize_receptor(value: object) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    lower = text.lower()
    if "positive" in lower or lower in {"pos", "+"}:
        return "Positive"
    if "negative" in lower or lower in {"neg", "-"}:
        return "Negative"
    return None


def fdr_bh(values: list[float]) -> list[float]:
    result = np.full(len(values), np.nan, dtype=float)
    valid = np.array([index for index, value in enumerate(values) if np.isfinite(value)], dtype=int)
    if valid.size == 0:
        return result.tolist()
    p = np.array([values[index] for index in valid], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    result[valid] = adjusted[inverse]
    return result.tolist()


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size == 0 or y.size == 0:
        return math.nan
    ranks = stats.rankdata(np.concatenate([x, y]))
    u = ranks[: x.size].sum() - x.size * (x.size + 1) / 2
    return float(2 * u / (x.size * y.size) - 1)


def paired_rank_biserial(diff: np.ndarray) -> float:
    values = diff[np.isfinite(diff) & (diff != 0)]
    if values.size == 0:
        return math.nan
    ranks = stats.rankdata(np.abs(values))
    positive = ranks[values > 0].sum()
    negative = ranks[values < 0].sum()
    return float((positive - negative) / (positive + negative))


def group_summary(data: pd.DataFrame, dataset: str, analysis: str) -> list[dict]:
    rows: list[dict] = []
    for group, subset in data.groupby("group", observed=True):
        values = subset["expression"].dropna().to_numpy(dtype=float)
        rows.append(
            {
                "dataset": dataset,
                "analysis": analysis,
                "group": group,
                "n": values.size,
                "mean": np.mean(values) if values.size else math.nan,
                "sd": np.std(values, ddof=1) if values.size > 1 else math.nan,
                "median": np.median(values) if values.size else math.nan,
                "q25": np.quantile(values, 0.25) if values.size else math.nan,
                "q75": np.quantile(values, 0.75) if values.size else math.nan,
            }
        )
    return rows


def pairwise_tests(data: pd.DataFrame, family: str, order: list[str]) -> list[dict]:
    rows: list[dict] = []
    available = [group for group in order if group in set(data["group"].dropna())]
    for group1, group2 in combinations(available, 2):
        x = data.loc[data["group"] == group1, "expression"].dropna().to_numpy(dtype=float)
        y = data.loc[data["group"] == group2, "expression"].dropna().to_numpy(dtype=float)
        if x.size < 2 or y.size < 2:
            continue
        test = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
        rows.append(
            {
                "family": family,
                "group1": group1,
                "group2": group2,
                "n1": x.size,
                "n2": y.size,
                "test": "Mann-Whitney U",
                "statistic": float(test.statistic),
                "median_diff_group1_minus_group2": float(np.median(x) - np.median(y)),
                "cliffs_delta": cliffs_delta(x, y),
                "p_value": float(test.pvalue),
            }
        )
    adjusted = fdr_bh([row["p_value"] for row in rows])
    for row, fdr in zip(rows, adjusted):
        row["FDR_BH_within_family"] = fdr
    return rows


def p_label(value: float) -> str:
    if not np.isfinite(value):
        return "P=NA"
    if value < 0.001:
        return f"P={value:.1e}"
    return f"P={value:.3f}"


def violin_box(ax: plt.Axes, data: pd.DataFrame, order: list[str], title: str, subtitle: str) -> None:
    present = [group for group in order if group in set(data["group"].dropna())]
    colors = {group: PALETTE.get(group, "#888888") for group in present}
    sns.violinplot(
        data=data,
        x="group",
        y="expression",
        order=present,
        palette=colors,
        inner=None,
        cut=0,
        linewidth=0.7,
        saturation=0.85,
        ax=ax,
        hue="group",
        hue_order=present,
        legend=False,
    )
    sns.boxplot(
        data=data,
        x="group",
        y="expression",
        order=present,
        width=0.22,
        showfliers=False,
        color="white",
        linecolor="#333333",
        linewidth=0.8,
        ax=ax,
    )
    ax.set_title(title, fontsize=11, fontweight="bold", pad=16)
    ax.text(0.5, 1.01, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("DNMT3B")
    ax.tick_params(axis="x", rotation=0)


def paired_plot(ax: plt.Axes, pairs: pd.DataFrame, title: str, subtitle: str) -> None:
    for _, row in pairs.iterrows():
        ax.plot([0, 1], [row["Normal"], row["Primary"]], color="#888888", alpha=0.16, linewidth=0.55)
    ax.scatter(np.zeros(len(pairs)), pairs["Normal"], color=PALETTE["Normal"], s=9, alpha=0.55, edgecolors="none")
    ax.scatter(np.ones(len(pairs)), pairs["Primary"], color=PALETTE["Primary"], s=9, alpha=0.55, edgecolors="none")
    ax.set_xticks([0, 1], ["Normal", "Primary"])
    ax.set_title(title, fontsize=11, fontweight="bold", pad=16)
    ax.text(0.5, 1.01, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("DNMT3B")


for path in [TCGA_MATRIX, TCGA_META, GSE70947_MATRIX, GSE65194_MATRIX, GSE70947_META, GSE65194_META]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Required input missing or empty: {path}")

print("PY_STEP=LOAD_AND_ALIGN_DATA", flush=True)
tcga_expr = read_gene(TCGA_MATRIX, "DNMT3B")
tcga_meta = pd.read_csv(TCGA_META, sep="\t", dtype=str).set_index("sample_id")
tcga = tcga_meta.join(tcga_expr.rename("expression"), how="inner")
tcga.index.name = "sample_id"
tcga["expression"] = pd.to_numeric(tcga["expression"], errors="coerce")

tn = tcga.loc[tcga["barcode_group"].isin(["Primary_Tumor", "Solid_Tissue_Normal"])].copy()
tn["group"] = tn["barcode_group"].map({"Primary_Tumor": "Primary", "Solid_Tissue_Normal": "Normal"})
tn = tn.dropna(subset=["expression"])

pair_source = tn.reset_index()[["sample_id", "patient_id", "group", "expression"]]
pair_counts = pair_source.groupby(["patient_id", "group"]).size().unstack(fill_value=0)
valid_patient_ids = pair_counts.index[(pair_counts.get("Primary", 0) == 1) & (pair_counts.get("Normal", 0) == 1)]
tcga_pairs = (
    pair_source.loc[pair_source["patient_id"].isin(valid_patient_ids)]
    .pivot(index="patient_id", columns="group", values="expression")
    .dropna(subset=["Primary", "Normal"])
)

primary = tcga.loc[tcga["barcode_group"] == "Primary_Tumor"].copy()
primary["pam50_group"] = primary["PAM50"].map(normalize_pam50)
primary["stage_group"] = primary["AJCC_stage"].map(normalize_stage)
for receptor in ["ER", "PR", "HER2"]:
    primary[f"{receptor}_group"] = primary[f"{receptor}_status"].map(normalize_receptor)

pam50 = primary.dropna(subset=["expression", "pam50_group"]).rename(columns={"pam50_group": "group"})
stage = primary.dropna(subset=["expression", "stage_group"]).rename(columns={"stage_group": "group"})
receptor_data = {
    receptor: primary.dropna(subset=["expression", f"{receptor}_group"]).rename(columns={f"{receptor}_group": "group"})
    for receptor in ["ER", "PR", "HER2"]
}

g709_expr = read_gene(GSE70947_MATRIX, "DNMT3B")
g709_meta = pd.read_csv(GSE70947_META, sep="\t", dtype=str).set_index("sample")
g709 = g709_meta.join(g709_expr.rename("expression"), how="inner")
g709.index.name = "sample"
g709["expression"] = pd.to_numeric(g709["expression"], errors="coerce")
g709["group"] = g709["disease_group"].map({"Tumor": "Primary", "Normal": "Normal"})
g709_pair_source = g709.reset_index()[["sample", "pair_key", "group", "expression"]].dropna()
g709_pair_counts = g709_pair_source.groupby(["pair_key", "group"]).size().unstack(fill_value=0)
g709_valid_keys = g709_pair_counts.index[(g709_pair_counts.get("Primary", 0) == 1) & (g709_pair_counts.get("Normal", 0) == 1)]
g709_pairs = (
    g709_pair_source.loc[g709_pair_source["pair_key"].isin(g709_valid_keys)]
    .pivot(index="pair_key", columns="group", values="expression")
    .dropna(subset=["Primary", "Normal"])
)

g651_expr = read_gene(GSE65194_MATRIX, "DNMT3B")
g651_meta = pd.read_csv(GSE65194_META, sep="\t", dtype=str).set_index("sample")
g651 = g651_meta.join(g651_expr.rename("expression"), how="inner")
g651.index.name = "sample"
g651["expression"] = pd.to_numeric(g651["expression"], errors="coerce")
g651["group"] = np.where(g651["disease_group"] == "Normal", "Normal", g651["subtype"])
g651 = g651.dropna(subset=["expression", "group"])

print("PY_STEP=STATISTICAL_TESTS", flush=True)
global_tests: list[dict] = []

def add_unpaired(name: str, data: pd.DataFrame, group1: str, group2: str) -> None:
    x = data.loc[data["group"] == group1, "expression"].dropna().to_numpy(dtype=float)
    y = data.loc[data["group"] == group2, "expression"].dropna().to_numpy(dtype=float)
    test = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
    global_tests.append(
        {
            "analysis": name,
            "test": "Mann-Whitney U",
            "n_total": x.size + y.size,
            "group_n": f"{group1}:{x.size};{group2}:{y.size}",
            "statistic": float(test.statistic),
            "effect_name": "Cliffs_delta_group1_vs_group2",
            "effect": cliffs_delta(x, y),
            "p_value": float(test.pvalue),
        }
    )


def add_paired(name: str, pairs: pd.DataFrame) -> None:
    diff = pairs["Primary"].to_numpy(dtype=float) - pairs["Normal"].to_numpy(dtype=float)
    test = stats.wilcoxon(diff, alternative="two-sided", zero_method="wilcox", method="auto")
    global_tests.append(
        {
            "analysis": name,
            "test": "Wilcoxon signed-rank",
            "n_total": len(diff),
            "group_n": f"pairs:{len(diff)}",
            "statistic": float(test.statistic),
            "effect_name": "Matched_rank_biserial_Primary_minus_Normal",
            "effect": paired_rank_biserial(diff),
            "p_value": float(test.pvalue),
        }
    )


def add_kruskal(name: str, data: pd.DataFrame, order: list[str]) -> None:
    groups = [data.loc[data["group"] == group, "expression"].dropna().to_numpy(dtype=float) for group in order]
    groups = [values for values in groups if values.size >= 2]
    test = stats.kruskal(*groups)
    global_tests.append(
        {
            "analysis": name,
            "test": "Kruskal-Wallis",
            "n_total": sum(len(values) for values in groups),
            "group_n": ";".join(f"{group}:{sum(data['group'] == group)}" for group in order if sum(data["group"] == group) > 0),
            "statistic": float(test.statistic),
            "effect_name": "epsilon_squared",
            "effect": max(0.0, float((test.statistic - len(groups) + 1) / (sum(len(values) for values in groups) - len(groups)))),
            "p_value": float(test.pvalue),
        }
    )


add_unpaired("TCGA_tumor_vs_normal", tn, "Primary", "Normal")
add_paired("TCGA_paired_tumor_normal", tcga_pairs)
add_kruskal("TCGA_PAM50", pam50, ["LumA", "LumB", "HER2", "Basal", "Normal"])
add_kruskal("TCGA_stage", stage, ["Stg I", "Stg II", "Stg III", "Stg IV"])

stage_codes = stage["group"].map({"Stg I": 1, "Stg II": 2, "Stg III": 3, "Stg IV": 4})
stage_trend = stats.spearmanr(stage_codes, stage["expression"], nan_policy="omit")
global_tests.append(
    {
        "analysis": "TCGA_stage_trend",
        "test": "Spearman trend",
        "n_total": int(stage_codes.notna().sum()),
        "group_n": "ordinal_stage_1_to_4",
        "statistic": float(stage_trend.statistic),
        "effect_name": "Spearman_rho",
        "effect": float(stage_trend.statistic),
        "p_value": float(stage_trend.pvalue),
    }
)
for receptor, data in receptor_data.items():
    add_unpaired(f"TCGA_{receptor}_status", data, "Positive", "Negative")
add_paired("GSE70947_paired_tumor_normal", g709_pairs)
add_kruskal("GSE65194_subtype", g651, ["Normal", "LumA", "LumB", "HER2", "TNBC"])

global_fdr = fdr_bh([row["p_value"] for row in global_tests])
for row, fdr in zip(global_tests, global_fdr):
    row["FDR_BH_all_global_tests"] = fdr

pairwise_rows: list[dict] = []
pairwise_rows.extend(pairwise_tests(pam50, "TCGA_PAM50", ["LumA", "LumB", "HER2", "Basal", "Normal"]))
pairwise_rows.extend(pairwise_tests(stage, "TCGA_stage", ["Stg I", "Stg II", "Stg III", "Stg IV"]))
pairwise_rows.extend(pairwise_tests(g651, "GSE65194_subtype", ["Normal", "LumA", "LumB", "HER2", "TNBC"]))

summary_rows: list[dict] = []
summary_rows.extend(group_summary(tn, "TCGA", "tumor_normal"))
summary_rows.extend(group_summary(pam50, "TCGA", "PAM50"))
summary_rows.extend(group_summary(stage, "TCGA", "stage"))
for receptor, data in receptor_data.items():
    summary_rows.extend(group_summary(data, "TCGA", f"{receptor}_status"))
summary_rows.extend(group_summary(g709, "GSE70947", "tumor_normal"))
summary_rows.extend(group_summary(g651, "GSE65194", "subtype"))

pd.DataFrame(global_tests).to_csv(RUN_DIR / "4_1_global_tests.tsv", sep="\t", index=False, na_rep="NA")
pd.DataFrame(pairwise_rows).to_csv(RUN_DIR / "4_1_pairwise_tests.tsv", sep="\t", index=False, na_rep="NA")
pd.DataFrame(summary_rows).to_csv(RUN_DIR / "4_1_group_summary.tsv", sep="\t", index=False, na_rep="NA")

analysis_data = []
for dataset, analysis, frame, id_name in [
    ("TCGA", "tumor_normal", tn.reset_index(), "sample_id"),
    ("TCGA", "PAM50", pam50.reset_index(), "sample_id"),
    ("TCGA", "stage", stage.reset_index(), "sample_id"),
    ("GSE70947", "tumor_normal", g709.reset_index(), "sample"),
    ("GSE65194", "subtype", g651.reset_index(), "sample"),
]:
    for _, row in frame.iterrows():
        analysis_data.append(
            {
                "dataset": dataset,
                "analysis": analysis,
                "sample_id": row[id_name],
                "patient_or_pair_id": row.get("patient_id", row.get("pair_key", "NA")),
                "group": row["group"],
                "DNMT3B_expression": row["expression"],
            }
        )
pd.DataFrame(analysis_data).to_csv(
    RUN_DIR / "4_1_analysis_data.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA"
)

print("PY_STEP=GENERATE_FIGURES", flush=True)
sns.set_theme(style="whitegrid", context="paper")
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

test_map = {row["analysis"]: row for row in global_tests}
fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), constrained_layout=True)
violin_box(
    axes[0, 0],
    tn,
    ["Normal", "Primary"],
    "TCGA T/N",
    p_label(test_map["TCGA_tumor_vs_normal"]["p_value"]),
)
paired_plot(
    axes[0, 1],
    tcga_pairs,
    "TCGA Pair",
    p_label(test_map["TCGA_paired_tumor_normal"]["p_value"]),
)
violin_box(
    axes[0, 2],
    pam50,
    ["LumA", "LumB", "HER2", "Basal", "Normal"],
    "PAM50",
    p_label(test_map["TCGA_PAM50"]["p_value"]),
)
violin_box(
    axes[1, 0],
    stage,
    ["Stg I", "Stg II", "Stg III", "Stg IV"],
    "Stage",
    f"rho={stage_trend.statistic:.2f}; {p_label(stage_trend.pvalue)}",
)
paired_plot(
    axes[1, 1],
    g709_pairs,
    "GEO Pair",
    p_label(test_map["GSE70947_paired_tumor_normal"]["p_value"]),
)
violin_box(
    axes[1, 2],
    g651,
    ["Normal", "LumA", "LumB", "HER2", "TNBC"],
    "GEO Subt.",
    p_label(test_map["GSE65194_subtype"]["p_value"]),
)
for label, ax in zip(list("ABCDEF"), axes.flat):
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")
fig.savefig(RUN_DIR / "Fig4_1_DNMT3B_main.pdf", bbox_inches="tight")
fig.savefig(RUN_DIR / "Fig4_1_DNMT3B_main.png", dpi=400, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.0), constrained_layout=True)
for label, receptor, ax in zip(list("ABC"), ["ER", "PR", "HER2"], axes):
    violin_box(
        ax,
        receptor_data[receptor],
        ["Negative", "Positive"],
        receptor,
        p_label(test_map[f"TCGA_{receptor}_status"]["p_value"]),
    )
    ax.text(-0.14, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")
fig.savefig(RUN_DIR / "Fig4_1_DNMT3B_receptor.pdf", bbox_inches="tight")
fig.savefig(RUN_DIR / "Fig4_1_DNMT3B_receptor.png", dpi=400, bbox_inches="tight")
plt.close(fig)

print("PY_STEP=VALIDATE_AND_SUMMARIZE", flush=True)
raw_valid = lambda series: int(series.map(clean_text).notna().sum())
g651_counts = g651["group"].value_counts().to_dict()
observed = {
    "TCGA_PRIMARY": int(sum(tn["group"] == "Primary")),
    "TCGA_NORMAL": int(sum(tn["group"] == "Normal")),
    "TCGA_PAIRS": tcga_pairs.shape[0],
    "PAM50_FIELD_VALID": raw_valid(primary["PAM50"]),
    "STAGE_FIELD_VALID": raw_valid(primary["AJCC_stage"]),
    "ER_FIELD_VALID": raw_valid(primary["ER_status"]),
    "PR_FIELD_VALID": raw_valid(primary["PR_status"]),
    "HER2_FIELD_VALID": raw_valid(primary["HER2_status"]),
    "GSE70947_PAIRS": g709_pairs.shape[0],
    "GSE65194_INCLUDED": g651.shape[0],
    "GSE65194_NORMAL": int(g651_counts.get("Normal", 0)),
    "GSE65194_TNBC": int(g651_counts.get("TNBC", 0)),
    "GSE65194_HER2": int(g651_counts.get("HER2", 0)),
    "GSE65194_LUMA": int(g651_counts.get("LumA", 0)),
    "GSE65194_LUMB": int(g651_counts.get("LumB", 0)),
}

audit_rows = []
all_pass = True
for metric, expected in EXPECTED.items():
    value = observed[metric]
    status = "PASS" if value == expected else "FAIL"
    all_pass = all_pass and status == "PASS"
    audit_rows.append({"metric": metric, "value": value, "expected": expected, "status": status})
pd.DataFrame(audit_rows).to_csv(RUN_DIR / "4_1_sample_audit.tsv", sep="\t", index=False)

caption = """Figure 4.1. DNMT3B expression and clinicopathological associations in breast cancer.
(A) TCGA-BRCA primary tumor versus solid-tissue normal samples. (B) Paired TCGA-BRCA tumor-normal comparison. (C) TCGA PAM50 subtypes. (D) DNMT3B expression across AJCC stages. (E) Paired external validation in GSE70947. (F) External molecular-subtype validation in GSE65194. The receptor-status figure shows ER, PR and HER2 comparisons in TCGA primary tumors. Violin plots show distributions and embedded boxes show medians and interquartile ranges. Paired panels connect samples from the same patient. Two-group comparisons used Mann-Whitney U tests or Wilcoxon signed-rank tests for paired data; multigroup comparisons used Kruskal-Wallis tests with BH-adjusted pairwise tests. Abbreviations: T/N, tumor/normal; Pair, paired comparison; Subt., subtype; Stg, AJCC stage; LumA, luminal A; LumB, luminal B; TNBC, triple-negative breast cancer; ER, estrogen receptor; PR, progesterone receptor; HER2, human epidermal growth factor receptor 2; rho, Spearman correlation coefficient.
"""
(RUN_DIR / "Fig4_1_caption.txt").write_text(caption, encoding="utf-8")

summary = [
    "ANALYSIS_SCOPE=4_1_DNMT3B_EXPRESSION_AND_CLINICOPATHOLOGY",
]
for row in audit_rows:
    summary.append(f"{row['metric']}={row['value']};EXPECTED={row['expected']};STATUS={row['status']}")
summary.extend(
    [
        f"PAM50_ANALYSIS_N={pam50.shape[0]}",
        f"STAGE_ANALYSIS_N={stage.shape[0]}",
        f"ER_BINARY_ANALYSIS_N={receptor_data['ER'].shape[0]}",
        f"PR_BINARY_ANALYSIS_N={receptor_data['PR'].shape[0]}",
        f"HER2_BINARY_ANALYSIS_N={receptor_data['HER2'].shape[0]}",
    ]
)
for row in global_tests:
    summary.append(
        "GLOBAL_TEST="
        + ";".join(
            [
                row["analysis"],
                row["test"],
                f"N={row['n_total']}",
                f"effect={row['effect']:.8g}",
                f"P={row['p_value']:.8g}",
                f"FDR={row['FDR_BH_all_global_tests']:.8g}",
            ]
        )
    )
summary.extend(
    [
        "MULTIPLE_TESTING=BH_WITHIN_PAIRWISE_FAMILY_AND_ACROSS_GLOBAL_TESTS",
        "RAW_AND_FROZEN_INPUTS_MODIFIED=NO",
        "FOUR_ONE_ANALYSIS_PASS" if all_pass else "FOUR_ONE_ANALYSIS_FAIL",
    ]
)
(RUN_DIR / "4_1_analysis_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

output_rows = []
for path in sorted(RUN_DIR.iterdir()):
    if path.is_file() and path.name not in {"dnmt3b_expression_clinical.log", "output_file_manifest.tsv"}:
        output_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": file_sha256(path)})
pd.DataFrame(output_rows).to_csv(RUN_DIR / "output_file_manifest.tsv", sep="\t", index=False)

print("\n".join(summary))
if not all_pass:
    raise SystemExit(2)
