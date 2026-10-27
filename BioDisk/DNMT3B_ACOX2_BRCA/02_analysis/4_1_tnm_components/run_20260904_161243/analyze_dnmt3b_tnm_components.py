#!/usr/bin/env python3
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr


def fail(message: str) -> None:
    raise RuntimeError(message)


def bh_adjust(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)
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
        if "4_1_tnm_components" not in str(path)
    ]
    if not candidates:
        fail("tcga_expression_8genes.tsv.gz not found under project")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def load_dnmt3b_expression(path: Path) -> pd.Series:
    frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    if frame.empty:
        fail(f"empty expression file: {path}")
    first_col = frame.columns[0]
    first_values = frame[first_col].astype(str).str.upper().str.strip()
    if (first_values == "DNMT3B").any():
        row = frame.loc[first_values == "DNMT3B"].iloc[0]
        series = pd.to_numeric(row.drop(labels=[first_col]), errors="coerce")
        series.index = [canon_sample(value) for value in series.index]
    elif "DNMT3B" in frame.columns:
        sample_candidates = ["sample_id", "sampleID", "sample", "Sample"]
        sample_col = next((col for col in sample_candidates if col in frame.columns), first_col)
        series = pd.Series(
            pd.to_numeric(frame["DNMT3B"], errors="coerce").to_numpy(),
            index=[canon_sample(value) for value in frame[sample_col]],
            name="DNMT3B",
        )
    else:
        fail(f"DNMT3B not found in expression file: {path}")
    series = series.dropna()
    if series.index.duplicated().any():
        series = series.groupby(level=0).median()
    series.name = "DNMT3B"
    return series


def resolve_sample_column(columns) -> str:
    preferred = ["sampleID", "sample_id", "sample", "Sample", "submitter_id"]
    lower_map = {str(col).lower(): str(col) for col in columns}
    for candidate in preferred:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    fail("sample identifier column not found in clinical table")


def normalize_column_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def resolve_component_column(columns, component: str) -> str:
    exact = {
        "T": [
            "AJCC_Tumor_Pathologic_pT_nature2012",
            "AJCC_Tumor_Pathologic_pT",
            "pathologic_T",
            "ajcc_pathologic_t",
            "pathologic_t_stage",
        ],
        "N": [
            "AJCC_Nodes_Pathologic_pN_nature2012",
            "AJCC_Nodes_Pathologic_pN",
            "pathologic_N",
            "ajcc_pathologic_n",
            "pathologic_n_stage",
        ],
        "M": [
            "AJCC_Metastasis_Pathologic_pM_nature2012",
            "AJCC_Metastasis_Pathologic_pM",
            "pathologic_M",
            "ajcc_pathologic_m",
            "pathologic_m_stage",
        ],
    }
    by_norm = {normalize_column_name(str(col)): str(col) for col in columns}
    for candidate in exact[component]:
        normalized = normalize_column_name(candidate)
        if normalized in by_norm:
            return by_norm[normalized]

    component_lower = component.lower()
    scored = []
    for col in columns:
        original = str(col)
        normalized = normalize_column_name(original)
        score = 0
        if "pathologic" in normalized:
            score += 5
        if "ajcc" in normalized:
            score += 3
        if f"p{component_lower}" in normalized:
            score += 6
        if component == "T" and "tumor" in normalized:
            score += 3
        if component == "N" and ("node" in normalized or "nodes" in normalized):
            score += 3
        if component == "M" and "metast" in normalized:
            score += 3
        if normalized.endswith(component_lower):
            score += 2
        if score >= 8:
            scored.append((score, original))
    if not scored:
        related = [str(col) for col in columns if any(
            token in str(col).lower() for token in ["ajcc", "pathologic", "tumor", "node", "metast"]
        )]
        fail(f"p{component} field not resolved; related columns={related[:40]}")
    scored.sort(key=lambda item: (-item[0], len(item[1]), item[1]))
    if len(scored) > 1 and scored[0][0] == scored[1][0]:
        fail(f"ambiguous p{component} fields: {scored[:5]}")
    return scored[0][1]


def parse_stage(value: object, component: str):
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper().replace(" ", "")
    if text in {"", "NA", "N/A", "NAN", "NOTAVAILABLE", "NOTREPORTED", "UNKNOWN", "[NOTAVAILABLE]"}:
        return np.nan
    text = text.replace("PATHOLOGIC", "").replace("STAGE", "")
    text = re.sub(r"^[CPYR]+(?=[TNM])", "", text)
    match = re.search(fr"{component}([0-4])", text)
    if not match:
        return np.nan
    level = int(match.group(1))
    allowed = {"T": {1, 2, 3, 4}, "N": {0, 1, 2, 3}, "M": {0, 1}}[component]
    return level if level in allowed else np.nan


def cliffs_delta(values1: np.ndarray, values2: np.ndarray) -> float:
    statistic = mannwhitneyu(values1, values2, alternative="two-sided").statistic
    return 2.0 * statistic / (len(values1) * len(values2)) - 1.0


def summarize_component(frame: pd.DataFrame, component: str, order: list[int]):
    label_col = f"p{component}"
    valid = frame.dropna(subset=[label_col, "DNMT3B"]).copy()
    valid[label_col] = valid[label_col].astype(int)
    valid = valid.loc[valid[label_col].isin(order)].copy()
    groups = [valid.loc[valid[label_col] == level, "DNMT3B"].to_numpy(dtype=float) for level in order]
    if any(len(values) < 2 for values in groups):
        fail(f"p{component} has a group with fewer than 2 samples: {[len(values) for values in groups]}")

    summary_rows = []
    for level, values in zip(order, groups):
        summary_rows.append({
            "component": f"p{component}",
            "group": f"{component}{level}",
            "ordinal": level,
            "n": len(values),
            "mean": np.mean(values),
            "sd": np.std(values, ddof=1),
            "median": np.median(values),
            "q25": np.quantile(values, 0.25),
            "q75": np.quantile(values, 0.75),
        })

    if len(order) == 2:
        omnibus_stat, omnibus_p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
        omnibus_test = "Mann-Whitney U"
        omnibus_effect_name = "Cliffs_delta_first_vs_second"
        omnibus_effect = cliffs_delta(groups[0], groups[1])
    else:
        omnibus_stat, omnibus_p = kruskal(*groups)
        omnibus_test = "Kruskal-Wallis"
        omnibus_effect_name = "epsilon_squared"
        n_total = sum(len(values) for values in groups)
        omnibus_effect = max((omnibus_stat - len(groups) + 1) / (n_total - len(groups)), 0.0)

    trend_rho, trend_p = spearmanr(valid[label_col], valid["DNMT3B"], nan_policy="omit")
    global_rows = [
        {
            "component": f"p{component}",
            "analysis": "omnibus",
            "test": omnibus_test,
            "n": len(valid),
            "statistic": omnibus_stat,
            "effect_name": omnibus_effect_name,
            "effect": omnibus_effect,
            "p_value": omnibus_p,
        },
        {
            "component": f"p{component}",
            "analysis": "ordinal_trend",
            "test": "Spearman trend",
            "n": len(valid),
            "statistic": trend_rho,
            "effect_name": "Spearman_rho",
            "effect": trend_rho,
            "p_value": trend_p,
        },
    ]

    pairwise_rows = []
    for index1, index2 in itertools.combinations(range(len(order)), 2):
        level1, level2 = order[index1], order[index2]
        values1, values2 = groups[index1], groups[index2]
        statistic, p_value = mannwhitneyu(values1, values2, alternative="two-sided")
        pairwise_rows.append({
            "component": f"p{component}",
            "group1": f"{component}{level1}",
            "group2": f"{component}{level2}",
            "n1": len(values1),
            "n2": len(values2),
            "test": "Mann-Whitney U",
            "statistic": statistic,
            "median_diff_group1_minus_group2": np.median(values1) - np.median(values2),
            "cliffs_delta": cliffs_delta(values1, values2),
            "p_value": p_value,
        })
    adjusted = bh_adjust([row["p_value"] for row in pairwise_rows])
    for row, fdr in zip(pairwise_rows, adjusted):
        row["FDR_BH_within_component"] = fdr
    return valid, summary_rows, global_rows, pairwise_rows


def add_panel(ax, valid, component: str, order: list[int], trend_fdr: float, rho: float) -> None:
    col = f"p{component}"
    data = [valid.loc[valid[col] == level, "DNMT3B"].to_numpy(dtype=float) for level in order]
    colors = {
        "T": ["#A8D5E5", "#79B9D4", "#468EBD", "#155C9B"],
        "N": ["#B7D7A8", "#82BD74", "#4C9B5F", "#216B43"],
        "M": ["#7AA6C2", "#D95F59"],
    }[component]
    parts = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("#333333")
        body.set_alpha(0.88)
    ax.boxplot(
        data,
        widths=0.22,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.4},
        boxprops={"facecolor": "white", "edgecolor": "#444444", "linewidth": 1.0},
        whiskerprops={"color": "#444444", "linewidth": 1.0},
        capprops={"color": "#444444", "linewidth": 1.0},
    )
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f"{component}{level}" for level in order])
    ax.set_title(f"p{component}\nrho={rho:.2f}; FDR={trend_fdr:.2g}", fontweight="bold", fontsize=11)
    ax.set_ylabel("DNMT3B")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    if len(sys.argv) != 3:
        fail("usage: analyze_dnmt3b_tnm_components.py PROJECT RUN_DIR")
    project = Path(sys.argv[1])
    run_dir = Path(sys.argv[2])
    clinical_path = project / "00_raw_data/TCGA_BRCA/TCGA.BRCA.clinical.tsv"
    if not clinical_path.is_file():
        fail(f"clinical file not found: {clinical_path}")
    expression_path = locate_expression(project)

    expression = load_dnmt3b_expression(expression_path)
    primary = expression.loc[[index.endswith("-01") for index in expression.index]].copy()
    if not 1050 <= len(primary) <= 1120:
        fail(f"unexpected primary-tumor expression count: {len(primary)}")

    clinical = pd.read_csv(clinical_path, sep="\t", low_memory=False, dtype=str)
    sample_col = resolve_sample_column(clinical.columns)
    t_col = resolve_component_column(clinical.columns, "T")
    n_col = resolve_component_column(clinical.columns, "N")
    m_col = resolve_component_column(clinical.columns, "M")
    clinical = clinical[[sample_col, t_col, n_col, m_col]].copy()
    clinical["sample_key"] = clinical[sample_col].map(canon_sample)
    clinical = clinical.drop_duplicates("sample_key", keep="first")

    analysis = pd.DataFrame({"sample_key": primary.index, "DNMT3B": primary.to_numpy(dtype=float)})
    analysis = analysis.merge(clinical, on="sample_key", how="left", validate="one_to_one")
    analysis["pT"] = analysis[t_col].map(lambda value: parse_stage(value, "T"))
    analysis["pN"] = analysis[n_col].map(lambda value: parse_stage(value, "N"))
    analysis["pM"] = analysis[m_col].map(lambda value: parse_stage(value, "M"))

    orders = {"T": [1, 2, 3, 4], "N": [0, 1, 2, 3], "M": [0, 1]}
    valid_frames = {}
    summaries = []
    global_tests = []
    pairwise_tests = []
    for component, order in orders.items():
        valid, component_summary, component_global, component_pairwise = summarize_component(
            analysis, component, order
        )
        valid_frames[component] = valid
        summaries.extend(component_summary)
        global_tests.extend(component_global)
        pairwise_tests.extend(component_pairwise)

    summary = pd.DataFrame(summaries)
    global_frame = pd.DataFrame(global_tests)
    pairwise = pd.DataFrame(pairwise_tests)
    global_frame["FDR_BH_all_TNM_tests"] = bh_adjust(global_frame["p_value"])

    analysis.to_csv(run_dir / "TCGA_DNMT3B_TNM_sample_level.tsv.gz", sep="\t", index=False)
    summary.to_csv(run_dir / "TCGA_DNMT3B_TNM_group_summary.tsv", sep="\t", index=False)
    global_frame.to_csv(run_dir / "TCGA_DNMT3B_TNM_global_tests.tsv", sep="\t", index=False)
    pairwise.to_csv(run_dir / "TCGA_DNMT3B_TNM_pairwise_tests.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2), dpi=220)
    for ax, component in zip(axes, ["T", "N", "M"]):
        trend = global_frame.loc[
            (global_frame["component"] == f"p{component}")
            & (global_frame["analysis"] == "ordinal_trend")
        ].iloc[0]
        add_panel(
            ax,
            valid_frames[component],
            component,
            orders[component],
            trend["FDR_BH_all_TNM_tests"],
            trend["effect"],
        )
    fig.suptitle("DNMT3B vs pathologic TNM components", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(run_dir / "Fig4_1_DNMT3B_TNM_components.png", dpi=400, bbox_inches="tight")
    fig.savefig(run_dir / "Fig4_1_DNMT3B_TNM_components.pdf", bbox_inches="tight")
    plt.close(fig)

    with open(run_dir / "TNM_component_analysis_summary.txt", "w", encoding="utf-8") as handle:
        handle.write("ANALYSIS_SCOPE=4_1_DNMT3B_PATHOLOGIC_TNM_COMPONENTS\n")
        handle.write(f"EXPRESSION_FILE={expression_path}\n")
        handle.write(f"CLINICAL_FILE={clinical_path}\n")
        handle.write(f"CLINICAL_T_FIELD={t_col}\n")
        handle.write(f"CLINICAL_N_FIELD={n_col}\n")
        handle.write(f"CLINICAL_M_FIELD={m_col}\n")
        handle.write(f"PRIMARY_EXPRESSION_SAMPLES={len(primary)}\n")
        for component in ["T", "N", "M"]:
            handle.write(f"p{component}_VALID_SAMPLES={len(valid_frames[component])}\n")
        for row in global_frame.itertuples(index=False):
            handle.write(
                f"GLOBAL_TEST={row.component};{row.analysis};{row.test};N={row.n};"
                f"effect={row.effect:.8g};P={row.p_value:.8g};FDR={row.FDR_BH_all_TNM_tests:.8g}\n"
            )
        handle.write("MULTIPLE_TESTING=BH_ACROSS_SIX_TNM_GLOBAL_TESTS_AND_WITHIN_PAIRWISE_COMPONENTS\n")
        handle.write("RAW_AND_FROZEN_INPUTS_MODIFIED=NO\n")
        handle.write("TNM_COMPONENT_ANALYSIS_PASS\n")

    print("ANALYSIS_SCOPE=4_1_DNMT3B_PATHOLOGIC_TNM_COMPONENTS")
    print(f"EXPRESSION_FILE={expression_path}")
    print(f"CLINICAL_FILE={clinical_path}")
    print(f"CLINICAL_T_FIELD={t_col}")
    print(f"CLINICAL_N_FIELD={n_col}")
    print(f"CLINICAL_M_FIELD={m_col}")
    print(f"PRIMARY_EXPRESSION_SAMPLES={len(primary)}")
    for component in ["T", "N", "M"]:
        print(f"p{component}_VALID_SAMPLES={len(valid_frames[component])}")
    for row in global_frame.itertuples(index=False):
        print(
            f"GLOBAL_TEST={row.component};{row.analysis};{row.test};N={row.n};"
            f"effect={row.effect:.8g};P={row.p_value:.8g};FDR={row.FDR_BH_all_TNM_tests:.8g}"
        )
    print("RAW_AND_FROZEN_INPUTS_MODIFIED=NO")
    print("TNM_COMPONENT_ANALYSIS_PASS")


if __name__ == "__main__":
    main()
