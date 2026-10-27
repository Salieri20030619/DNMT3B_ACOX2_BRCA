import csv
import gzip
import itertools
import math
import os
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu


def fail(message):
    raise RuntimeError(message)


def parse_series_matrix(path):
    data = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = row[0]
            values = row[1:]
            if key == "!Sample_title":
                data["title"] = values
            elif key == "!Sample_geo_accession":
                data["gsm"] = values
            elif key == "!Sample_platform_id":
                data["platform"] = values
            elif key == "!Sample_characteristics_ch1":
                fields = {
                    "external_id": "scan-b external id:",
                    "age": "age at diagnosis:",
                    "er": "er status:",
                    "pgr": "pgr status:",
                    "her2": "her2 status:",
                    "pam50": "pam50 subtype:",
                }
                for field, prefix in fields.items():
                    if values and all(value.lower().startswith(prefix) for value in values):
                        data[field] = [value.split(":", 1)[1].strip() for value in values]
            elif key == "!series_matrix_table_begin":
                break

    required = ["title", "gsm", "platform", "external_id", "er", "pgr", "her2", "pam50"]
    missing = [field for field in required if field not in data]
    if missing:
        fail(f"missing metadata fields in {path}: {','.join(missing)}")

    n = len(data["title"])
    if any(len(data[field]) != n for field in required):
        fail(f"metadata length mismatch in {path}")

    records = []
    for index in range(n):
        records.append({field: data[field][index] for field in required})
    return records


def read_target_gene(path, target_gene):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        sample_names = header[1:]
        target_rows = []
        for row in reader:
            if row and row[0].strip().upper() == target_gene.upper():
                values = np.array([
                    float(value) if value not in ("", "NA", "NaN", "nan") else np.nan
                    for value in row[1:]
                ], dtype=float)
                target_rows.append(values)
        if len(target_rows) != 1:
            fail(f"expected one {target_gene} row, found {len(target_rows)}")
        if len(sample_names) != len(target_rows[0]):
            fail("expression header and value count differ")
        return sample_names, target_rows[0]


def bh_adjust(p_values):
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 1.0
    m = len(p_values)
    for rank, index in reversed(list(enumerate(order, start=1))):
        running = min(running, p_values[index] * m / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def clinical_subtype(row):
    er, pgr, her2 = row["er"], row["pgr"], row["her2"]
    if "NA" in (er, pgr, her2):
        return None
    if er == "0" and pgr == "0" and her2 == "0":
        return "TNBC"
    if her2 == "1":
        return "HER2"
    if her2 == "0" and (er == "1" or pgr == "1"):
        return "Luminal"
    return None


def summarize_groups(frame, family, group_col, order):
    rows = []
    arrays = []
    for group in order:
        values = frame.loc[frame[group_col] == group, "DNMT3B"].dropna().to_numpy(dtype=float)
        arrays.append(values)
        rows.append({
            "family": family,
            "group": group,
            "n": len(values),
            "mean": np.mean(values),
            "sd": np.std(values, ddof=1),
            "median": np.median(values),
            "q25": np.quantile(values, 0.25),
            "q75": np.quantile(values, 0.75),
        })
    statistic, p_value = kruskal(*arrays)
    global_row = {
        "family": family,
        "test": "Kruskal-Wallis",
        "n_total": sum(len(values) for values in arrays),
        "statistic": statistic,
        "p_value": p_value,
    }

    pairwise = []
    for group1, group2 in itertools.combinations(order, 2):
        values1 = frame.loc[frame[group_col] == group1, "DNMT3B"].dropna().to_numpy(dtype=float)
        values2 = frame.loc[frame[group_col] == group2, "DNMT3B"].dropna().to_numpy(dtype=float)
        statistic, p_value = mannwhitneyu(values1, values2, alternative="two-sided")
        delta = 2.0 * statistic / (len(values1) * len(values2)) - 1.0
        pairwise.append({
            "family": family,
            "group1": group1,
            "group2": group2,
            "n1": len(values1),
            "n2": len(values2),
            "test": "Mann-Whitney U",
            "statistic": statistic,
            "median_diff_group1_minus_group2": np.median(values1) - np.median(values2),
            "cliffs_delta": delta,
            "p_value": p_value,
        })
    adjusted = bh_adjust([row["p_value"] for row in pairwise])
    for row, fdr in zip(pairwise, adjusted):
        row["FDR_BH_within_family"] = fdr
    return rows, global_row, pairwise


def add_boxplot(axis, frame, group_col, order, colors, title, p_value):
    data = [frame.loc[frame[group_col] == group, "DNMT3B"].dropna().to_numpy() for group in order]
    parts = axis.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("#333333")
        body.set_alpha(0.85)
    box = axis.boxplot(
        data,
        widths=0.22,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.4},
        boxprops={"facecolor": "white", "edgecolor": "#444444", "linewidth": 1.0},
        whiskerprops={"color": "#444444", "linewidth": 1.0},
        capprops={"color": "#444444", "linewidth": 1.0},
    )
    axis.set_xticks(range(1, len(order) + 1))
    axis.set_xticklabels(order)
    axis.set_title(f"{title}\nP={p_value:.2e}", fontweight="bold", fontsize=12)
    axis.set_ylabel("DNMT3B")
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def main():
    if len(sys.argv) != 5:
        fail("usage: analyze.py MATRIX_HISEQ MATRIX_NEXTSEQ EXPR_FILE RUN_DIR")
    matrix_hiseq, matrix_nextseq, expression_file, run_dir = sys.argv[1:]

    metadata_records = parse_series_matrix(matrix_hiseq) + parse_series_matrix(matrix_nextseq)
    metadata = {row["title"]: row for row in metadata_records}
    if len(metadata_records) != 3409 or len(metadata) != 3409:
        fail(f"expected 3409 expression columns, found {len(metadata_records)} metadata rows and {len(metadata)} unique titles")

    sample_names, expression = read_target_gene(expression_file, "DNMT3B")
    if len(sample_names) != 3409:
        fail(f"expected 3409 expression columns, found {len(sample_names)}")
    missing_metadata = sorted(set(sample_names) - set(metadata))
    if missing_metadata:
        fail(f"expression samples without metadata: {missing_metadata[:5]}")

    technical_rows = []
    for sample_name, value in zip(sample_names, expression):
        row = dict(metadata[sample_name])
        row["sample_name"] = sample_name
        row["DNMT3B"] = value
        row["biological_key"] = ".".join(row["external_id"].split(".")[:3])
        technical_rows.append(row)
    technical = pd.DataFrame(technical_rows)

    replicate_mask = technical["sample_name"].str.lower().str.endswith("repl")
    frame = technical.loc[~replicate_mask].copy()
    frame["n_technical"] = 1
    if frame["biological_key"].duplicated().any():
        fail("primary sample titles do not resolve to unique biological keys")
    frame["clinical_subtype"] = frame.apply(clinical_subtype, axis=1)

    n_biological = len(frame)
    n_replicates = int(replicate_mask.sum())
    if n_biological != 3273 or n_replicates != 136:
        fail(f"replicate collapse mismatch: biological={n_biological}, technical_excess={n_replicates}")

    pam50_order = ["LumA", "LumB", "Her2", "Basal", "Normal"]
    clinical_order = ["Luminal", "HER2", "TNBC"]
    expected_pam50 = {"LumA": 1657, "LumB": 729, "Her2": 327, "Basal": 339, "Normal": 221}
    expected_clinical = {"Luminal": 2346, "HER2": 338, "TNBC": 143}
    observed_pam50 = frame["pam50"].value_counts().to_dict()
    observed_clinical = frame["clinical_subtype"].value_counts().to_dict()
    if observed_pam50 != expected_pam50:
        fail(f"PAM50 counts differ from expected: {observed_pam50}")
    if observed_clinical != expected_clinical:
        fail(f"clinical subtype counts differ from expected: {observed_clinical}")

    summary_rows = []
    global_rows = []
    pairwise_rows = []
    for family, group_col, order in [
        ("PAM50", "pam50", pam50_order),
        ("Clinical", "clinical_subtype", clinical_order),
    ]:
        summary, global_row, pairwise = summarize_groups(frame, family, group_col, order)
        summary_rows.extend(summary)
        global_rows.append(global_row)
        pairwise_rows.extend(pairwise)

    summary = pd.DataFrame(summary_rows)
    global_tests = pd.DataFrame(global_rows)
    pairwise = pd.DataFrame(pairwise_rows)

    frame.to_csv(os.path.join(run_dir, "GSE96058_DNMT3B_sample_level.tsv.gz"), sep="\t", index=False)
    summary.to_csv(os.path.join(run_dir, "GSE96058_DNMT3B_group_summary.tsv"), sep="\t", index=False)
    global_tests.to_csv(os.path.join(run_dir, "GSE96058_DNMT3B_global_tests.tsv"), sep="\t", index=False)
    pairwise.to_csv(os.path.join(run_dir, "GSE96058_DNMT3B_pairwise_tests.tsv"), sep="\t", index=False)

    clinical_summary = summary.loc[summary["family"] == "Clinical"].set_index("group")
    tnbc_median = clinical_summary.loc["TNBC", "median"]
    her2_median = clinical_summary.loc["HER2", "median"]
    luminal_median = clinical_summary.loc["Luminal", "median"]
    tnbc_vs_her2 = pairwise.loc[
        (pairwise["family"] == "Clinical")
        & (pairwise["group1"] == "HER2")
        & (pairwise["group2"] == "TNBC")
    ].iloc[0]
    clinical_global_p = global_tests.loc[global_tests["family"] == "Clinical", "p_value"].iloc[0]
    pam50_global_p = global_tests.loc[global_tests["family"] == "PAM50", "p_value"].iloc[0]

    decision = (
        tnbc_median > her2_median > luminal_median
        and tnbc_vs_her2["FDR_BH_within_family"] < 0.05
        and clinical_global_p < 0.001
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), dpi=200)
    add_boxplot(
        axes[0], frame, "pam50", pam50_order,
        ["#5B8FF9", "#61DDAA", "#F6BD16", "#E8684A", "#9C9C9C"],
        "PAM50", pam50_global_p,
    )
    add_boxplot(
        axes[1], frame.dropna(subset=["clinical_subtype"]), "clinical_subtype", clinical_order,
        ["#5B8FF9", "#F6BD16", "#E8684A"],
        "IHC", clinical_global_p,
    )
    fig.suptitle("GSE96058: DNMT3B expression", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(run_dir, "Fig4_1_GSE96058_DNMT3B_subtype.png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(run_dir, "Fig4_1_GSE96058_DNMT3B_subtype.pdf"), bbox_inches="tight")
    plt.close(fig)

    with open(os.path.join(run_dir, "GSE96058_analysis_summary.txt"), "w", encoding="utf-8") as handle:
        handle.write("ANALYSIS_SCOPE=GSE96058_DNMT3B_SUBTYPE_REPLACEMENT\n")
        handle.write(f"EXPRESSION_COLUMNS={len(technical)}\n")
        handle.write(f"BIOLOGICAL_SAMPLES={n_biological}\n")
        handle.write(f"TECHNICAL_REPLICATES_EXCLUDED={n_replicates}\n")
        for group in pam50_order:
            row = summary[(summary["family"] == "PAM50") & (summary["group"] == group)].iloc[0]
            handle.write(f"PAM50_{group}_N={int(row['n'])};MEDIAN={row['median']:.8g}\n")
        for group in clinical_order:
            row = summary[(summary["family"] == "Clinical") & (summary["group"] == group)].iloc[0]
            handle.write(f"CLINICAL_{group}_N={int(row['n'])};MEDIAN={row['median']:.8g}\n")
        handle.write(f"CLINICAL_GLOBAL_P={clinical_global_p:.8g}\n")
        handle.write(f"TNBC_VS_HER2_FDR={tnbc_vs_her2['FDR_BH_within_family']:.8g}\n")
        handle.write(f"TNBC_HIGHEST_DECISION={'PASS' if decision else 'FAIL'}\n")
        handle.write(
            "REPLACEMENT_DECISION="
            + ("USE_GSE96058_FOR_SUBTYPE_VALIDATION" if decision else "DO_NOT_REPLACE_WITH_GSE96058")
            + "\n"
        )

    print("GSE96058_EXPRESSION_COLUMNS=3409")
    print(f"GSE96058_BIOLOGICAL_SAMPLES={n_biological};EXPECTED=3273;STATUS=PASS")
    print(f"GSE96058_TECHNICAL_REPLICATES_EXCLUDED={n_replicates};EXPECTED=136;STATUS=PASS")
    for group in pam50_order:
        row = summary[(summary["family"] == "PAM50") & (summary["group"] == group)].iloc[0]
        print(f"PAM50_GROUP={group};N={int(row['n'])};MEDIAN={row['median']:.8g}")
    for group in clinical_order:
        row = summary[(summary["family"] == "Clinical") & (summary["group"] == group)].iloc[0]
        print(f"CLINICAL_GROUP={group};N={int(row['n'])};MEDIAN={row['median']:.8g}")
    print(f"CLINICAL_GLOBAL_TEST=Kruskal-Wallis;P={clinical_global_p:.8g}")
    print(
        "TNBC_VS_HER2="
        f"median_diff={tnbc_median - her2_median:.8g};"
        f"cliffs_delta={-tnbc_vs_her2['cliffs_delta']:.8g};"
        f"P={tnbc_vs_her2['p_value']:.8g};"
        f"FDR={tnbc_vs_her2['FDR_BH_within_family']:.8g}"
    )
    print(f"TNBC_HIGHEST_DECISION={'PASS' if decision else 'FAIL'}")
    print(
        "REPLACEMENT_DECISION="
        + ("USE_GSE96058_FOR_SUBTYPE_VALIDATION" if decision else "DO_NOT_REPLACE_WITH_GSE96058")
    )
    if not decision:
        fail("GSE96058 did not meet the predefined alternative-validation decision")


if __name__ == "__main__":
    main()
