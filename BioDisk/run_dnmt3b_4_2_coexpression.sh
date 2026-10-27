#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
OUT_ROOT="${PROJECT}/02_analysis/4_2_coexpression"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_ROOT}/run_${STAMP}"
PY_SCRIPT="${RUN_DIR}/analyze_4_2_coexpression.py"
LOG="${RUN_DIR}/4_2_coexpression.log"

mkdir -p "${RUN_DIR}"
exec > >(tee "${LOG}") 2>&1

echo "RUN_DIR=${RUN_DIR}"
echo "STEP=WRITE_PYTHON_ANALYSIS"

cat > "${PY_SCRIPT}" <<'PYTHON'
#!/usr/bin/env python3
from __future__ import annotations

import itertools
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu, spearmanr


GENES = ["DNMT3B", "RGS16", "GPRC5A", "NPNT", "AKR1C2", "ADCY5", "FGFR4", "ACOX2"]
CANDIDATES = ["RGS16", "GPRC5A", "NPNT", "AKR1C2", "ADCY5", "FGFR4"]
TARGET = "ACOX2"
EXPECTED_SIGN = {
    "RGS16": "negative",
    "GPRC5A": "positive",
    "NPNT": "positive",
    "AKR1C2": "positive",
    "ADCY5": "negative",
    "FGFR4": "negative",
    "ACOX2": "not_predefined",
}


def fail(message: str) -> None:
    raise RuntimeError(message)


def bh_adjust(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or not np.isfinite(p).all():
        fail("BH adjustment received invalid p-values")
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out


def canon_sample(value: object) -> str:
    text = str(value).strip().upper().replace(".", "-")
    match = re.match(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2})", text)
    return match.group(1) if match else text


def locate_expression(project: Path) -> Path:
    candidates = [
        path for path in project.rglob("tcga_expression_8genes.tsv.gz")
        if "4_2_coexpression" not in str(path)
    ]
    if not candidates:
        fail("tcga_expression_8genes.tsv.gz not found under project")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def load_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if frame.empty:
        fail(f"empty expression file: {path}")
    first_col = frame.columns[0]
    first_values = frame[first_col].astype(str).str.upper().str.strip()
    if set(GENES).issubset(set(first_values)):
        subset = frame.loc[first_values.isin(GENES)].copy()
        subset[first_col] = first_values.loc[subset.index]
        if subset[first_col].duplicated().any():
            fail("duplicated gene rows in expression input")
        matrix = subset.set_index(first_col).loc[GENES].T
        matrix.index = [canon_sample(value) for value in matrix.index]
    elif set(GENES).issubset(frame.columns):
        sample_options = ["sample_id", "sampleID", "sample", "Sample"]
        sample_col = next((col for col in sample_options if col in frame.columns), first_col)
        matrix = frame.set_index(sample_col)[GENES].copy()
        matrix.index = [canon_sample(value) for value in matrix.index]
    else:
        observed = sorted(set(first_values).intersection(GENES))
        fail(f"required genes not found; observed={observed}")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if matrix.index.duplicated().any():
        matrix = matrix.groupby(level=0).median()
    matrix.index.name = "sample_id"
    return matrix


def cliffs_delta(high: np.ndarray, low: np.ndarray) -> float:
    statistic = mannwhitneyu(high, low, alternative="two-sided").statistic
    return 2.0 * statistic / (len(high) * len(low)) - 1.0


def rho_class(value: float) -> str:
    magnitude = abs(value)
    if magnitude < 0.10:
        return "negligible"
    if magnitude < 0.30:
        return "weak"
    if magnitude < 0.50:
        return "moderate"
    return "strong"


def correlation_tables(primary: pd.DataFrame):
    all_rows = []
    for gene1, gene2 in itertools.combinations(GENES, 2):
        pair = primary[[gene1, gene2]].dropna()
        if len(pair) < 20:
            fail(f"insufficient complete samples for {gene1} vs {gene2}: {len(pair)}")
        rho, p_value = spearmanr(pair[gene1], pair[gene2])
        if not np.isfinite(rho) or not np.isfinite(p_value):
            fail(f"invalid Spearman result for {gene1} vs {gene2}")
        all_rows.append({
            "source": gene1,
            "target": gene2,
            "n_complete": len(pair),
            "spearman_rho": rho,
            "p_value": p_value,
            "sign": "positive" if rho > 0 else "negative" if rho < 0 else "zero",
            "strength": rho_class(rho),
        })
    all_pairs = pd.DataFrame(all_rows)
    all_pairs["FDR_BH_all_28_pairs"] = bh_adjust(all_pairs["p_value"])
    all_pairs["abs_rho"] = all_pairs["spearman_rho"].abs()
    all_pairs["primary_network_edge"] = (
        (all_pairs["FDR_BH_all_28_pairs"] < 0.05)
        & (all_pairs["abs_rho"] >= 0.10)
    )

    dnmt = all_pairs.loc[
        (all_pairs["source"] == "DNMT3B") | (all_pairs["target"] == "DNMT3B")
    ].copy()
    dnmt["gene"] = np.where(dnmt["source"] == "DNMT3B", dnmt["target"], dnmt["source"])
    dnmt = dnmt[[
        "gene", "n_complete", "spearman_rho", "p_value", "sign", "strength"
    ]].copy()
    dnmt["FDR_BH_DNMT3B_7_tests"] = bh_adjust(dnmt["p_value"])
    dnmt["role"] = np.where(dnmt["gene"].isin(CANDIDATES), "candidate", "target")
    dnmt["KD_expected_sign"] = dnmt["gene"].map(EXPECTED_SIGN)
    observed = np.where(dnmt["spearman_rho"] > 0, "positive", "negative")
    dnmt["KD_direction_concordance"] = np.where(
        dnmt["KD_expected_sign"] == "not_predefined",
        "not_assessed",
        np.where(observed == dnmt["KD_expected_sign"], "concordant", "discordant"),
    )
    dnmt["candidate_order"] = dnmt["gene"].map({gene: i for i, gene in enumerate(CANDIDATES + [TARGET])})
    dnmt = dnmt.sort_values("candidate_order").drop(columns="candidate_order").reset_index(drop=True)
    return all_pairs, dnmt


def stratified_tests(primary: pd.DataFrame):
    cutoff = float(primary["DNMT3B"].median())
    group = np.where(primary["DNMT3B"] <= cutoff, "Low", "High")
    grouped = primary.copy()
    grouped.insert(0, "DNMT3B_group", group)
    rows = []
    for gene in CANDIDATES + [TARGET]:
        high = grouped.loc[grouped["DNMT3B_group"] == "High", gene].dropna().to_numpy(dtype=float)
        low = grouped.loc[grouped["DNMT3B_group"] == "Low", gene].dropna().to_numpy(dtype=float)
        if len(high) < 20 or len(low) < 20:
            fail(f"insufficient high/low samples for {gene}")
        statistic, p_value = mannwhitneyu(high, low, alternative="two-sided")
        rows.append({
            "gene": gene,
            "role": "candidate" if gene in CANDIDATES else "target",
            "n_high": len(high),
            "n_low": len(low),
            "median_high": np.median(high),
            "median_low": np.median(low),
            "median_diff_high_minus_low": np.median(high) - np.median(low),
            "cliffs_delta_high_vs_low": cliffs_delta(high, low),
            "mann_whitney_U": statistic,
            "p_value": p_value,
        })
    tests = pd.DataFrame(rows)
    tests["FDR_BH_7_high_low_tests"] = bh_adjust(tests["p_value"])
    return grouped, tests, cutoff


def gene_summary(primary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in GENES:
        values = primary[gene].dropna().to_numpy(dtype=float)
        rows.append({
            "gene": gene,
            "role": "anchor" if gene == "DNMT3B" else "candidate" if gene in CANDIDATES else "target",
            "n": len(values),
            "missing": int(primary[gene].isna().sum()),
            "mean": np.mean(values),
            "sd": np.std(values, ddof=1),
            "median": np.median(values),
            "q25": np.quantile(values, 0.25),
            "q75": np.quantile(values, 0.75),
        })
    return pd.DataFrame(rows)


def correlation_matrices(all_pairs: pd.DataFrame):
    rho = pd.DataFrame(np.eye(len(GENES)), index=GENES, columns=GENES)
    fdr = pd.DataFrame(np.zeros((len(GENES), len(GENES))), index=GENES, columns=GENES)
    for row in all_pairs.itertuples(index=False):
        rho.loc[row.source, row.target] = row.spearman_rho
        rho.loc[row.target, row.source] = row.spearman_rho
        fdr.loc[row.source, row.target] = row.FDR_BH_all_28_pairs
        fdr.loc[row.target, row.source] = row.FDR_BH_all_28_pairs
    return rho, fdr


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6, alpha=0.8)


def draw_main_figure(
    dnmt: pd.DataFrame,
    tests: pd.DataFrame,
    rho_matrix: pd.DataFrame,
    fdr_matrix: pd.DataFrame,
    run_dir: Path,
) -> None:
    fig = plt.figure(figsize=(12.5, 10.0), dpi=220)
    grid = fig.add_gridspec(2, 2, wspace=0.30, hspace=0.35)

    ax1 = fig.add_subplot(grid[0, 0])
    ordered = dnmt.sort_values("spearman_rho")
    colors = ["#D95F59" if value > 0 else "#4477AA" for value in ordered["spearman_rho"]]
    ax1.barh(ordered["gene"], ordered["spearman_rho"], color=colors, alpha=0.90)
    ax1.axvline(0, color="#333333", linewidth=0.9)
    for index, row in enumerate(ordered.itertuples(index=False)):
        marker = "*" if row.FDR_BH_DNMT3B_7_tests < 0.05 else "ns"
        offset = 0.012 if row.spearman_rho >= 0 else -0.012
        align = "left" if row.spearman_rho >= 0 else "right"
        ax1.text(row.spearman_rho + offset, index, f"{row.spearman_rho:.2f} {marker}", va="center", ha=align, fontsize=8)
    ax1.set_xlabel("Spearman rho")
    ax1.set_title("A  DNMT3B correlation", loc="left", fontweight="bold")
    limit = max(0.35, float(ordered["spearman_rho"].abs().max()) + 0.12)
    ax1.set_xlim(-limit, limit)
    style_axis(ax1)

    ax2 = fig.add_subplot(grid[0, 1])
    ordered_tests = tests.sort_values("cliffs_delta_high_vs_low")
    colors2 = ["#D95F59" if value > 0 else "#4477AA" for value in ordered_tests["cliffs_delta_high_vs_low"]]
    ax2.barh(ordered_tests["gene"], ordered_tests["cliffs_delta_high_vs_low"], color=colors2, alpha=0.90)
    ax2.axvline(0, color="#333333", linewidth=0.9)
    for index, row in enumerate(ordered_tests.itertuples(index=False)):
        marker = "*" if row.FDR_BH_7_high_low_tests < 0.05 else "ns"
        value = row.cliffs_delta_high_vs_low
        offset = 0.012 if value >= 0 else -0.012
        align = "left" if value >= 0 else "right"
        ax2.text(value + offset, index, f"{value:.2f} {marker}", va="center", ha=align, fontsize=8)
    ax2.set_xlabel("Cliff's delta")
    ax2.set_title("B  High vs low", loc="left", fontweight="bold")
    limit2 = max(0.35, float(ordered_tests["cliffs_delta_high_vs_low"].abs().max()) + 0.12)
    ax2.set_xlim(-limit2, limit2)
    style_axis(ax2)

    ax3 = fig.add_subplot(grid[1, 0])
    outer = CANDIDATES + [TARGET]
    angles = np.linspace(0, 2 * np.pi, len(outer), endpoint=False) + np.pi / 2
    positions = {"DNMT3B": (0.0, 0.0)}
    positions.update({gene: (math.cos(angle), math.sin(angle)) for gene, angle in zip(outer, angles)})
    for row in dnmt.itertuples(index=False):
        x2, y2 = positions[row.gene]
        color = "#D95F59" if row.spearman_rho > 0 else "#4477AA"
        significant = row.FDR_BH_DNMT3B_7_tests < 0.05
        ax3.plot(
            [0, x2], [0, y2], color=color,
            linewidth=0.8 + 6.0 * abs(row.spearman_rho),
            alpha=0.90 if significant else 0.22,
            linestyle="-" if significant else "--",
            zorder=1,
        )
        ax3.text(x2 * 0.52, y2 * 0.52, f"{row.spearman_rho:.2f}", fontsize=7, ha="center", va="center", backgroundcolor="white")
    ax3.scatter([0], [0], s=1450, color="#6A3D9A", edgecolor="white", linewidth=1.5, zorder=3)
    ax3.text(0, 0, "DNMT3B", color="white", fontweight="bold", ha="center", va="center", fontsize=7.5, zorder=4)
    for gene in outer:
        x, y = positions[gene]
        color = "#2A9D8F" if gene in CANDIDATES else "#777777"
        ax3.scatter([x], [y], s=1250, color=color, edgecolor="white", linewidth=1.3, zorder=3)
        ax3.text(x, y, gene, color="white", fontweight="bold", ha="center", va="center", fontsize=7.2, zorder=4)
    ax3.set_xlim(-1.35, 1.35)
    ax3.set_ylim(-1.25, 1.25)
    ax3.set_aspect("equal")
    ax3.axis("off")
    ax3.set_title("C  DNMT3B network", loc="left", fontweight="bold")
    ax3.legend(handles=[
        Line2D([0], [0], color="#D95F59", lw=3, label="Positive"),
        Line2D([0], [0], color="#4477AA", lw=3, label="Negative"),
    ], loc="lower left", frameon=False, fontsize=8)

    ax4 = fig.add_subplot(grid[1, 1])
    image = ax4.imshow(rho_matrix.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax4.set_xticks(range(len(GENES)))
    ax4.set_yticks(range(len(GENES)))
    ax4.set_xticklabels(GENES, rotation=45, ha="right", fontsize=8)
    ax4.set_yticklabels(GENES, fontsize=8)
    for i in range(len(GENES)):
        for j in range(len(GENES)):
            value = rho_matrix.iloc[i, j]
            star = "*" if i != j and fdr_matrix.iloc[i, j] < 0.05 else ""
            color = "white" if abs(value) > 0.55 else "black"
            ax4.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=6.5, color=color)
    ax4.set_title("D  Correlation matrix", loc="left", fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax4, fraction=0.046, pad=0.04)
    colorbar.set_label("Spearman rho")

    fig.suptitle("DNMT3B coexpression in TCGA-BRCA", fontsize=16, fontweight="bold", y=0.99)
    fig.savefig(run_dir / "Fig4_2_DNMT3B_coexpression_main.png", dpi=400, bbox_inches="tight")
    fig.savefig(run_dir / "Fig4_2_DNMT3B_coexpression_main.pdf", bbox_inches="tight")
    plt.close(fig)


def draw_stratified_heatmap(grouped: pd.DataFrame, rho_matrix: pd.DataFrame, run_dir: Path) -> None:
    genes = CANDIDATES + [TARGET]
    distance = 1.0 - rho_matrix.loc[genes, genes].to_numpy()
    np.fill_diagonal(distance, 0.0)
    distance = np.clip((distance + distance.T) / 2.0, 0.0, 2.0)
    order = leaves_list(linkage(squareform(distance, checks=False), method="average"))
    ordered_genes = [genes[index] for index in order]

    sorted_frame = grouped.sort_values(["DNMT3B_group", "DNMT3B"], ascending=[True, True]).copy()
    sorted_frame["DNMT3B_group"] = pd.Categorical(
        sorted_frame["DNMT3B_group"], categories=["Low", "High"], ordered=True
    )
    sorted_frame = sorted_frame.sort_values(["DNMT3B_group", "DNMT3B"])
    values = sorted_frame[ordered_genes].T
    row_sd = values.std(axis=1, ddof=1).replace(0, np.nan)
    zscore = values.sub(values.mean(axis=1), axis=0).div(row_sd, axis=0).clip(-3, 3)
    if zscore.isna().all(axis=1).any():
        fail("constant gene encountered in stratified heatmap")
    low_n = int((sorted_frame["DNMT3B_group"] == "Low").sum())

    fig = plt.figure(figsize=(12.0, 4.2), dpi=220)
    grid = fig.add_gridspec(2, 1, height_ratios=[0.18, 4.0], hspace=0.04)
    ax_top = fig.add_subplot(grid[0, 0])
    group_colors = np.zeros((1, len(sorted_frame), 3))
    group_colors[:, :low_n, :] = np.array([68, 119, 170]) / 255.0
    group_colors[:, low_n:, :] = np.array([217, 95, 89]) / 255.0
    ax_top.imshow(group_colors, aspect="auto")
    ax_top.set_xticks([])
    ax_top.set_yticks([0])
    ax_top.set_yticklabels(["Group"], fontsize=8)
    ax_top.spines[:].set_visible(False)

    ax = fig.add_subplot(grid[1, 0])
    image = ax.imshow(zscore.to_numpy(), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto", interpolation="nearest")
    ax.axvline(low_n - 0.5, color="white", linewidth=1.4)
    ax.set_yticks(range(len(ordered_genes)))
    ax.set_yticklabels(ordered_genes)
    ax.set_xticks([low_n / 2, low_n + (len(sorted_frame) - low_n) / 2])
    ax.set_xticklabels([f"Low (n={low_n})", f"High (n={len(sorted_frame)-low_n})"])
    ax.set_xlabel("DNMT3B group")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.018, pad=0.015)
    colorbar.set_label("Row z-score")
    fig.suptitle("Candidate expression by DNMT3B group", fontweight="bold", fontsize=14, y=0.98)
    fig.subplots_adjust(top=0.88)
    fig.savefig(run_dir / "Fig4_2_DNMT3B_stratified_heatmap.png", dpi=400, bbox_inches="tight")
    fig.savefig(run_dir / "Fig4_2_DNMT3B_stratified_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if len(sys.argv) != 3:
        fail("usage: analyze_4_2_coexpression.py PROJECT RUN_DIR")
    project = Path(sys.argv[1])
    run_dir = Path(sys.argv[2])
    expression_path = locate_expression(project)
    expression = load_expression(expression_path)
    primary = expression.loc[[sample.endswith("-01") for sample in expression.index]].copy()
    if len(primary) != 1097:
        fail(f"primary sample count mismatch: observed={len(primary)} expected=1097")
    missing_rate = primary.isna().mean()
    if (missing_rate > 0.05).any():
        fail(f"gene missing rate exceeds 5%: {missing_rate.to_dict()}")

    all_pairs, dnmt = correlation_tables(primary)
    grouped, high_low, cutoff = stratified_tests(primary)
    summaries = gene_summary(primary)
    rho_matrix, fdr_matrix = correlation_matrices(all_pairs)

    primary.reset_index().to_csv(
        run_dir / "TCGA_primary_expression_8genes.tsv.gz", sep="\t", index=False
    )
    grouped.reset_index().to_csv(
        run_dir / "TCGA_primary_expression_8genes_grouped.tsv.gz", sep="\t", index=False
    )
    summaries.to_csv(run_dir / "TCGA_8gene_expression_summary.tsv", sep="\t", index=False)
    dnmt.to_csv(run_dir / "DNMT3B_target_spearman.tsv", sep="\t", index=False)
    all_pairs.to_csv(run_dir / "all_pairwise_spearman.tsv", sep="\t", index=False)
    high_low.to_csv(run_dir / "DNMT3B_high_low_tests.tsv", sep="\t", index=False)
    rho_matrix.to_csv(run_dir / "spearman_rho_matrix.tsv", sep="\t", index_label="gene")
    fdr_matrix.to_csv(run_dir / "spearman_FDR_matrix.tsv", sep="\t", index_label="gene")

    nodes = summaries.rename(columns={"gene": "id"}).copy()
    nodes["label"] = nodes["id"]
    nodes["KD_expected_sign"] = nodes["id"].map(EXPECTED_SIGN).fillna("anchor")
    nodes.to_csv(run_dir / "cytoscape_nodes.tsv", sep="\t", index=False)
    edges = all_pairs.rename(columns={
        "spearman_rho": "weight",
        "FDR_BH_all_28_pairs": "FDR",
    }).copy()
    edges["interaction"] = "coexpression"
    edges.to_csv(run_dir / "cytoscape_edges.tsv", sep="\t", index=False)

    draw_main_figure(dnmt, high_low, rho_matrix, fdr_matrix, run_dir)
    draw_stratified_heatmap(grouped, rho_matrix, run_dir)

    summary_path = run_dir / "4_2_coexpression_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("ANALYSIS_SCOPE=4_2_DNMT3B_COEXPRESSION_NETWORK\n")
        handle.write(f"EXPRESSION_FILE={expression_path}\n")
        handle.write(f"TCGA_PRIMARY_SAMPLES={len(primary)};EXPECTED=1097;STATUS=PASS\n")
        handle.write(f"GENES={','.join(GENES)}\n")
        handle.write(f"CORE_CANDIDATES={','.join(CANDIDATES)}\n")
        handle.write("ACOX2_ROLE=PREDEFINED_TARGET_REFERENCE_NOT_CORE_SIX\n")
        handle.write(f"DNMT3B_MEDIAN_CUTOFF={cutoff:.8g}\n")
        handle.write(f"DNMT3B_LOW_N={(grouped['DNMT3B_group']=='Low').sum()}\n")
        handle.write(f"DNMT3B_HIGH_N={(grouped['DNMT3B_group']=='High').sum()}\n")
        for row in dnmt.itertuples(index=False):
            handle.write(
                f"DNMT3B_COR={row.gene};N={row.n_complete};rho={row.spearman_rho:.8g};"
                f"P={row.p_value:.8g};FDR={row.FDR_BH_DNMT3B_7_tests:.8g};"
                f"strength={row.strength};KD_concordance={row.KD_direction_concordance}\n"
            )
        for row in high_low.itertuples(index=False):
            handle.write(
                f"HIGH_LOW={row.gene};N_high={row.n_high};N_low={row.n_low};"
                f"delta={row.cliffs_delta_high_vs_low:.8g};P={row.p_value:.8g};"
                f"FDR={row.FDR_BH_7_high_low_tests:.8g}\n"
            )
        handle.write(f"PRIMARY_NETWORK_EDGES={int(all_pairs['primary_network_edge'].sum())}\n")
        handle.write("PRIMARY_EDGE_RULE=FDR_LT_0.05_AND_ABS_RHO_GE_0.10\n")
        handle.write("KD_DIRECTION_USE=HYPOTHESIS_ONLY_NOT_FILTER\n")
        handle.write("MULTIPLE_TESTING=BH_WITHIN_DNMT3B_CORRELATIONS_ALL_PAIRS_AND_HIGH_LOW_FAMILIES\n")
        handle.write("RAW_AND_FROZEN_INPUTS_MODIFIED=NO\n")
        handle.write("FOUR_TWO_COEXPRESSION_PASS\n")

    print("ANALYSIS_SCOPE=4_2_DNMT3B_COEXPRESSION_NETWORK")
    print(f"EXPRESSION_FILE={expression_path}")
    print(f"TCGA_PRIMARY_SAMPLES={len(primary)};EXPECTED=1097;STATUS=PASS")
    print(f"GENES={','.join(GENES)}")
    print(f"DNMT3B_MEDIAN_CUTOFF={cutoff:.8g}")
    for row in dnmt.itertuples(index=False):
        print(
            f"DNMT3B_COR={row.gene};N={row.n_complete};rho={row.spearman_rho:.8g};"
            f"P={row.p_value:.8g};FDR={row.FDR_BH_DNMT3B_7_tests:.8g};"
            f"strength={row.strength};KD_concordance={row.KD_direction_concordance}"
        )
    for row in high_low.itertuples(index=False):
        print(
            f"HIGH_LOW={row.gene};N_high={row.n_high};N_low={row.n_low};"
            f"delta={row.cliffs_delta_high_vs_low:.8g};P={row.p_value:.8g};"
            f"FDR={row.FDR_BH_7_high_low_tests:.8g}"
        )
    print(f"PRIMARY_NETWORK_EDGES={int(all_pairs['primary_network_edge'].sum())}")
    print("RAW_AND_FROZEN_INPUTS_MODIFIED=NO")
    print("FOUR_TWO_COEXPRESSION_PASS")


if __name__ == "__main__":
    main()
PYTHON

echo "STEP=BASH_SYNTAX_CHECK"
bash -n "$0"
echo "BASH_SYNTAX_CHECK=PASS"

echo "STEP=PYTHON_DEPENDENCY_CHECK"
conda run -n batch_bench python -c 'import matplotlib, numpy, pandas, scipy; print("PYTHON_DEPENDENCY_CHECK=PASS")'

echo "STEP=PYTHON_SYNTAX_CHECK"
conda run -n batch_bench python -m py_compile "${PY_SCRIPT}"
echo "PYTHON_SYNTAX_CHECK=PASS"
echo "SYNTAX_CHECK=PASS"

echo "STEP=COEXPRESSION_ANALYSIS"
conda run -n batch_bench python "${PY_SCRIPT}" "${PROJECT}" "${RUN_DIR}"

for required in \
    "${RUN_DIR}/TCGA_primary_expression_8genes.tsv.gz" \
    "${RUN_DIR}/TCGA_primary_expression_8genes_grouped.tsv.gz" \
    "${RUN_DIR}/TCGA_8gene_expression_summary.tsv" \
    "${RUN_DIR}/DNMT3B_target_spearman.tsv" \
    "${RUN_DIR}/all_pairwise_spearman.tsv" \
    "${RUN_DIR}/DNMT3B_high_low_tests.tsv" \
    "${RUN_DIR}/spearman_rho_matrix.tsv" \
    "${RUN_DIR}/spearman_FDR_matrix.tsv" \
    "${RUN_DIR}/cytoscape_nodes.tsv" \
    "${RUN_DIR}/cytoscape_edges.tsv" \
    "${RUN_DIR}/Fig4_2_DNMT3B_coexpression_main.png" \
    "${RUN_DIR}/Fig4_2_DNMT3B_coexpression_main.pdf" \
    "${RUN_DIR}/Fig4_2_DNMT3B_stratified_heatmap.png" \
    "${RUN_DIR}/Fig4_2_DNMT3B_stratified_heatmap.pdf" \
    "${RUN_DIR}/4_2_coexpression_summary.txt"; do
    if [[ ! -s "${required}" ]]; then
        echo "ERROR=required_output_missing_or_empty:${required}" >&2
        exit 1
    fi
done

if ! grep -qx 'FOUR_TWO_COEXPRESSION_PASS' "${RUN_DIR}/4_2_coexpression_summary.txt"; then
    echo "ERROR=success_marker_missing" >&2
    exit 1
fi

ln -sfn "${RUN_DIR}" "${OUT_ROOT}/latest"
echo "ANALYSIS_LOGIC_PASS"
echo "LATEST_LINK=${OUT_ROOT}/latest"
echo "EXPECTED_SUCCESS_MARKER=FOUR_TWO_COEXPRESSION_PASS"
