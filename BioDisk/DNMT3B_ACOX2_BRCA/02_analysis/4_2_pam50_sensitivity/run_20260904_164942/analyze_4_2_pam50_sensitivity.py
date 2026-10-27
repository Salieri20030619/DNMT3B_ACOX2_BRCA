#!/usr/bin/env python3
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, rankdata, spearmanr, t as t_dist


GENES = [
    "DNMT3B", "RGS16", "GPRC5A", "NPNT",
    "AKR1C2", "ADCY5", "FGFR4", "ACOX2",
]
TEST_GENES = [
    "RGS16", "GPRC5A", "NPNT",
    "AKR1C2", "ADCY5", "FGFR4", "ACOX2",
]
CORE_GENES = [
    "RGS16", "GPRC5A", "NPNT",
    "AKR1C2", "ADCY5", "FGFR4",
]
SUBTYPES = ["LumA", "LumB", "HER2", "Basal"]
EXPECTED_COUNTS = {
    "LumA": 231,
    "LumB": 127,
    "HER2": 58,
    "Basal": 98,
    "Normal-like": 8,
}


def fail(message: str) -> None:
    raise RuntimeError(message)


def bh_adjust(values) -> np.ndarray:
    p = np.asarray(values, dtype=float)

    if (
        p.ndim != 1
        or len(p) == 0
        or not np.isfinite(p).all()
    ):
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
    match = re.match(
        r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2})",
        text,
    )
    return match.group(1) if match else text


def locate_expression(project: Path) -> Path:
    candidates = [
        path
        for path in project.rglob(
            "tcga_expression_8genes.tsv.gz"
        )
        if "4_2_pam50_sensitivity" not in str(path)
    ]

    if not candidates:
        fail("tcga_expression_8genes.tsv.gz not found")

    return max(
        candidates,
        key=lambda path: path.stat().st_mtime,
    )


def load_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
    )

    if frame.empty:
        fail(f"empty expression file: {path}")

    first_col = frame.columns[0]
    symbols = (
        frame[first_col]
        .astype(str)
        .str.upper()
        .str.strip()
    )

    if set(GENES).issubset(set(symbols)):
        subset = frame.loc[
            symbols.isin(GENES)
        ].copy()

        subset[first_col] = symbols.loc[subset.index]

        if subset[first_col].duplicated().any():
            fail("duplicated gene rows")

        matrix = (
            subset
            .set_index(first_col)
            .loc[GENES]
            .T
        )

        matrix.index = [
            canon_sample(value)
            for value in matrix.index
        ]

    elif set(GENES).issubset(frame.columns):
        sample_col = next(
            (
                col
                for col in [
                    "sample_id",
                    "sampleID",
                    "sample",
                    "Sample",
                ]
                if col in frame.columns
            ),
            first_col,
        )

        matrix = (
            frame
            .set_index(sample_col)[GENES]
            .copy()
        )

        matrix.index = [
            canon_sample(value)
            for value in matrix.index
        ]

    else:
        fail(
            "required eight genes not found "
            "in expression input"
        )

    matrix = matrix.apply(
        pd.to_numeric,
        errors="coerce",
    )

    if matrix.index.duplicated().any():
        matrix = matrix.groupby(
            level=0
        ).median()

    matrix.index.name = "sample_id"
    return matrix


def resolve_column(
    columns,
    preferred: list[str],
    tokens: list[str],
) -> str:
    lower = {
        str(col).lower(): str(col)
        for col in columns
    }

    for candidate in preferred:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    matches = []

    for col in columns:
        normalized = re.sub(
            r"[^a-z0-9]+",
            "",
            str(col).lower(),
        )

        if all(
            token in normalized
            for token in tokens
        ):
            matches.append(str(col))

    if len(matches) != 1:
        fail(
            "clinical column resolution failed: "
            f"preferred={preferred}; "
            f"matches={matches}"
        )

    return matches[0]


def normalize_subtype(
    value: object,
) -> str | None:
    if pd.isna(value):
        return None

    text = re.sub(
        r"[^A-Z0-9]+",
        "",
        str(value).upper(),
    )

    mapping = {
        "LUMA": "LumA",
        "LUMINALA": "LumA",
        "LUMB": "LumB",
        "LUMINALB": "LumB",
        "HER2": "HER2",
        "HER2ENRICHED": "HER2",
        "BASAL": "Basal",
        "BASALLIKE": "Basal",
        "NORMAL": "Normal-like",
        "NORMALLIKE": "Normal-like",
    }

    return mapping.get(text)


def load_clinical(
    project: Path,
) -> tuple[pd.DataFrame, Path, str, str]:
    path = (
        project
        / "00_raw_data"
        / "TCGA_BRCA"
        / "TCGA.BRCA.clinical.tsv"
    )

    if not path.is_file():
        fail(f"clinical file missing: {path}")

    frame = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    sample_col = resolve_column(
        frame.columns,
        [
            "sampleID",
            "sample_id",
            "sample",
            "Sample",
            "submitter_id",
        ],
        ["sample"],
    )

    pam50_col = resolve_column(
        frame.columns,
        [
            "PAM50_mRNA_nature2012",
            "PAM50",
            "pam50",
        ],
        ["pam50"],
    )

    out = frame[
        [sample_col, pam50_col]
    ].copy()

    out["sample_id"] = out[
        sample_col
    ].map(canon_sample)

    out["PAM50"] = out[
        pam50_col
    ].map(normalize_subtype)

    out["has_pam50"] = (
        out["PAM50"]
        .notna()
        .astype(int)
    )

    out = out.sort_values(
        ["sample_id", "has_pam50"],
        ascending=[True, False],
    )

    out = out.drop_duplicates(
        "sample_id",
        keep="first",
    )

    return (
        out[["sample_id", "PAM50"]],
        path,
        sample_col,
        pam50_col,
    )


def raw_correlations(
    primary: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for gene in TEST_GENES:
        pair = primary[
            ["DNMT3B", gene]
        ].dropna()

        rho, p_value = spearmanr(
            pair["DNMT3B"],
            pair[gene],
        )

        rows.append({
            "gene": gene,
            "n": len(pair),
            "raw_rho": rho,
            "raw_p": p_value,
        })

    out = pd.DataFrame(rows)
    out["raw_FDR"] = bh_adjust(
        out["raw_p"]
    )

    return out


def residualize(
    values: pd.Series,
    design: np.ndarray,
) -> np.ndarray:
    ranked = rankdata(
        values.to_numpy(dtype=float),
        method="average",
    )

    coefficients = np.linalg.lstsq(
        design,
        ranked,
        rcond=None,
    )[0]

    fitted = design @ coefficients
    return ranked - fitted


def adjusted_correlations(
    annotated: pd.DataFrame,
) -> pd.DataFrame:
    dummies = pd.get_dummies(
        annotated["PAM50"],
        drop_first=True,
        dtype=float,
    )

    design = np.column_stack([
        np.ones(len(annotated)),
        dummies.to_numpy(dtype=float),
    ])

    covariate_count = (
        design.shape[1] - 1
    )

    x_resid = residualize(
        annotated["DNMT3B"],
        design,
    )

    rows = []

    for gene in TEST_GENES:
        y_resid = residualize(
            annotated[gene],
            design,
        )

        rho = float(
            np.corrcoef(
                x_resid,
                y_resid,
            )[0, 1]
        )

        degrees_freedom = (
            len(annotated)
            - covariate_count
            - 2
        )

        denominator = max(
            1.0 - rho * rho,
            np.finfo(float).eps,
        )

        statistic = (
            rho
            * math.sqrt(
                degrees_freedom
                / denominator
            )
        )

        p_value = (
            2.0
            * t_dist.sf(
                abs(statistic),
                df=degrees_freedom,
            )
        )

        rows.append({
            "gene": gene,
            "n": len(annotated),
            "covariates": "PAM50",
            "dummy_variables": covariate_count,
            "degrees_freedom": degrees_freedom,
            "adjusted_rho": rho,
            "adjusted_p": p_value,
        })

    out = pd.DataFrame(rows)
    out["adjusted_FDR"] = bh_adjust(
        out["adjusted_p"]
    )

    return out


def within_subtype_correlations(
    annotated: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for subtype in SUBTYPES:
        subset = annotated.loc[
            annotated["PAM50"] == subtype
        ]

        for gene in TEST_GENES:
            pair = subset[
                ["DNMT3B", gene]
            ].dropna()

            rho, p_value = spearmanr(
                pair["DNMT3B"],
                pair[gene],
            )

            rows.append({
                "subtype": subtype,
                "gene": gene,
                "n": len(pair),
                "spearman_rho": rho,
                "p_value": p_value,
            })

    out = pd.DataFrame(rows)

    out["FDR_BH_all_28_tests"] = (
        bh_adjust(out["p_value"])
    )

    out[
        "FDR_BH_within_subtype"
    ] = np.nan

    for subtype in SUBTYPES:
        mask = (
            out["subtype"] == subtype
        )

        out.loc[
            mask,
            "FDR_BH_within_subtype",
        ] = bh_adjust(
            out.loc[
                mask,
                "p_value",
            ]
        )

    return out


def subtype_meta_analysis(
    within: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for gene in TEST_GENES:
        subset = within.loc[
            within["gene"] == gene
        ].copy()

        rho = np.clip(
            subset["spearman_rho"]
            .to_numpy(dtype=float),
            -0.999999,
            0.999999,
        )

        weights = (
            subset["n"]
            .to_numpy(dtype=float)
            - 3.0
        )

        fisher_z = np.arctanh(rho)

        pooled_z = float(
            np.sum(
                weights * fisher_z
            )
            / np.sum(weights)
        )

        standard_error = (
            1.0
            / math.sqrt(
                float(np.sum(weights))
            )
        )

        z_statistic = (
            pooled_z
            / standard_error
        )

        p_value = (
            2.0
            * norm.sf(
                abs(z_statistic)
            )
        )

        pooled_rho = math.tanh(
            pooled_z
        )

        ci_low = math.tanh(
            pooled_z
            - 1.96 * standard_error
        )

        ci_high = math.tanh(
            pooled_z
            + 1.96 * standard_error
        )

        q_statistic = float(
            np.sum(
                weights
                * (
                    fisher_z
                    - pooled_z
                ) ** 2
            )
        )

        q_df = len(subset) - 1
        q_p = float(
            chi2.sf(
                q_statistic,
                q_df,
            )
        )

        if q_statistic > 0:
            i2 = max(
                (
                    q_statistic
                    - q_df
                )
                / q_statistic,
                0.0,
            ) * 100.0
        else:
            i2 = 0.0

        rows.append({
            "gene": gene,
            "subtypes": len(subset),
            "n_total": int(
                subset["n"].sum()
            ),
            "meta_rho": pooled_rho,
            "CI95_low": ci_low,
            "CI95_high": ci_high,
            "meta_p": p_value,
            "heterogeneity_Q": q_statistic,
            "heterogeneity_df": q_df,
            "heterogeneity_p": q_p,
            "I2_percent": i2,
        })

    out = pd.DataFrame(rows)
    out["meta_FDR"] = bh_adjust(
        out["meta_p"]
    )

    return out


def classify_results(
    combined: pd.DataFrame,
) -> pd.DataFrame:
    out = combined.copy()

    out["direction_stable"] = (
        np.sign(out["raw_rho"])
        == np.sign(
            out["adjusted_rho"]
        )
    )

    out["absolute_change"] = (
        out["adjusted_rho"].abs()
        - out["raw_rho"].abs()
    )

    labels = []

    for row in out.itertuples(
        index=False
    ):
        if (
            np.sign(row.raw_rho)
            != np.sign(
                row.adjusted_rho
            )
        ):
            label = (
                "direction_reversed"
            )

        elif (
            row.adjusted_FDR < 0.05
            and abs(
                row.adjusted_rho
            ) >= 0.10
        ):
            label = (
                "retained_after_adjustment"
            )

        elif (
            row.raw_FDR < 0.05
            and abs(
                row.raw_rho
            ) >= 0.10
            and abs(
                row.adjusted_rho
            ) < 0.10
        ):
            label = (
                "attenuated_below_0.10"
            )

        elif (
            row.raw_FDR < 0.05
            and row.adjusted_FDR >= 0.05
        ):
            label = (
                "not_significant_after_adjustment"
            )

        else:
            label = (
                "effect_below_0.10"
            )

        labels.append(label)

    out[
        "adjustment_decision"
    ] = labels

    out["role"] = np.where(
        out["gene"].isin(
            CORE_GENES
        ),
        "core_candidate",
        "target_reference",
    )

    return out


def draw_figure(
    combined: pd.DataFrame,
    within: pd.DataFrame,
    meta: pd.DataFrame,
    run_dir: Path,
) -> None:
    columns = [
        "Raw",
        "Adj",
        "LumA",
        "LumB",
        "HER2",
        "Basal",
        "Meta",
    ]

    effects = pd.DataFrame(
        index=TEST_GENES,
        columns=columns,
        dtype=float,
    )

    fdrs = pd.DataFrame(
        index=TEST_GENES,
        columns=columns,
        dtype=float,
    )

    indexed = combined.set_index(
        "gene"
    )

    meta_indexed = meta.set_index(
        "gene"
    )

    effects["Raw"] = indexed[
        "raw_rho"
    ]

    effects["Adj"] = indexed[
        "adjusted_rho"
    ]

    effects["Meta"] = meta_indexed[
        "meta_rho"
    ]

    fdrs["Raw"] = indexed[
        "raw_FDR"
    ]

    fdrs["Adj"] = indexed[
        "adjusted_FDR"
    ]

    fdrs["Meta"] = meta_indexed[
        "meta_FDR"
    ]

    for subtype in SUBTYPES:
        subset = (
            within.loc[
                within["subtype"]
                == subtype
            ]
            .set_index("gene")
        )

        effects[subtype] = subset[
            "spearman_rho"
        ]

        fdrs[subtype] = subset[
            "FDR_BH_all_28_tests"
        ]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.5, 5.5),
        dpi=220,
        gridspec_kw={
            "width_ratios": [
                1.35,
                1.0,
            ]
        },
    )

    ax1, ax2 = axes

    image = ax1.imshow(
        effects.to_numpy(),
        cmap="RdBu_r",
        vmin=-0.5,
        vmax=0.5,
        aspect="auto",
    )

    ax1.set_xticks(
        range(len(columns))
    )

    ax1.set_xticklabels(
        columns
    )

    ax1.set_yticks(
        range(len(TEST_GENES))
    )

    ax1.set_yticklabels(
        TEST_GENES
    )

    ax1.axhline(
        5.5,
        color="black",
        linewidth=1.0,
    )

    for i in range(
        len(TEST_GENES)
    ):
        for j in range(
            len(columns)
        ):
            value = effects.iloc[
                i,
                j,
            ]

            star = (
                "*"
                if fdrs.iloc[
                    i,
                    j,
                ] < 0.05
                else ""
            )

            color = (
                "white"
                if abs(value) >= 0.35
                else "black"
            )

            ax1.text(
                j,
                i,
                f"{value:.2f}{star}",
                ha="center",
                va="center",
                fontsize=8,
                color=color,
            )

    ax1.set_title(
        "A  PAM50 check",
        loc="left",
        fontweight="bold",
    )

    colorbar = fig.colorbar(
        image,
        ax=ax1,
        fraction=0.046,
        pad=0.04,
    )

    colorbar.set_label(
        "Spearman rho"
    )

    y = np.arange(
        len(TEST_GENES)
    )

    bar_height = 0.34

    ax2.barh(
        y + bar_height / 2,
        indexed.loc[
            TEST_GENES,
            "raw_rho",
        ],
        height=bar_height,
        color="#888888",
        label="Raw",
    )

    ax2.barh(
        y - bar_height / 2,
        indexed.loc[
            TEST_GENES,
            "adjusted_rho",
        ],
        height=bar_height,
        color="#E76F51",
        label="Adj",
    )

    ax2.axvline(
        0,
        color="black",
        linewidth=0.8,
    )

    ax2.axvline(
        -0.10,
        color="#999999",
        linewidth=0.7,
        linestyle="--",
    )

    ax2.axvline(
        0.10,
        color="#999999",
        linewidth=0.7,
        linestyle="--",
    )

    ax2.axhline(
        5.5,
        color="black",
        linewidth=1.0,
    )

    ax2.set_yticks(y)

    ax2.set_yticklabels(
        TEST_GENES
    )

    ax2.invert_yaxis()

    ax2.set_xlabel(
        "Spearman rho"
    )

    ax2.set_title(
        "B  Raw vs Adj",
        loc="left",
        fontweight="bold",
    )

    ax2.legend(
        frameon=False,
        loc="lower right",
    )

    ax2.grid(
        axis="x",
        color="#DDDDDD",
        linewidth=0.6,
    )

    ax2.spines[
        "top"
    ].set_visible(False)

    ax2.spines[
        "right"
    ].set_visible(False)

    fig.suptitle(
        "DNMT3B PAM50 check",
        fontsize=15,
        fontweight="bold",
        y=0.99,
    )

    fig.tight_layout()

    fig.savefig(
        run_dir
        / "Fig4_2_DNMT3B_PAM50_sensitivity.png",
        dpi=400,
        bbox_inches="tight",
    )

    fig.savefig(
        run_dir
        / "Fig4_2_DNMT3B_PAM50_sensitivity.pdf",
        bbox_inches="tight",
    )

    plt.close(fig)


def main() -> None:
    if len(sys.argv) != 3:
        fail(
            "usage: "
            "analyze_4_2_pam50_sensitivity.py "
            "PROJECT RUN_DIR"
        )

    project = Path(
        sys.argv[1]
    )

    run_dir = Path(
        sys.argv[2]
    )

    expression_path = (
        locate_expression(project)
    )

    expression = load_expression(
        expression_path
    )

    primary = expression.loc[
        [
            sample.endswith("-01")
            for sample
            in expression.index
        ]
    ].copy()

    if len(primary) != 1097:
        fail(
            "primary sample count "
            f"mismatch: {len(primary)}"
        )

    if (
        primary
        .isna()
        .mean()
        .max()
        > 0.05
    ):
        fail(
            "expression missing rate "
            "exceeds 5%"
        )

    (
        clinical,
        clinical_path,
        sample_col,
        pam50_col,
    ) = load_clinical(project)

    merged = (
        primary
        .reset_index()
        .merge(
            clinical,
            on="sample_id",
            how="left",
            validate="one_to_one",
        )
    )

    observed_counts = (
        merged["PAM50"]
        .value_counts()
        .to_dict()
    )

    for subtype, expected in (
        EXPECTED_COUNTS.items()
    ):
        observed = int(
            observed_counts.get(
                subtype,
                0,
            )
        )

        if observed != expected:
            fail(
                "PAM50 count mismatch "
                f"for {subtype}: "
                f"observed={observed}; "
                f"expected={expected}"
            )

    annotated = merged.loc[
        merged["PAM50"].isin(
            SUBTYPES
        )
    ].copy()

    if len(annotated) != 514:
        fail(
            "four-subtype analysis "
            f"count mismatch: "
            f"{len(annotated)}"
        )

    raw = raw_correlations(
        primary
    )

    adjusted = adjusted_correlations(
        annotated
    )

    within = (
        within_subtype_correlations(
            annotated
        )
    )

    meta = subtype_meta_analysis(
        within
    )

    combined = (
        raw
        .merge(
            adjusted,
            on="gene",
            validate="one_to_one",
        )
        .merge(
            meta,
            on="gene",
            validate="one_to_one",
        )
    )

    combined = classify_results(
        combined
    )

    merged.to_csv(
        run_dir
        / "TCGA_primary_8genes_PAM50.tsv.gz",
        sep="\t",
        index=False,
    )

    raw.to_csv(
        run_dir
        / "DNMT3B_raw_spearman.tsv",
        sep="\t",
        index=False,
    )

    adjusted.to_csv(
        run_dir
        / "DNMT3B_PAM50_adjusted_spearman.tsv",
        sep="\t",
        index=False,
    )

    within.to_csv(
        run_dir
        / "DNMT3B_within_PAM50_spearman.tsv",
        sep="\t",
        index=False,
    )

    meta.to_csv(
        run_dir
        / "DNMT3B_PAM50_fixed_effect_meta.tsv",
        sep="\t",
        index=False,
    )

    combined.to_csv(
        run_dir
        / "DNMT3B_PAM50_sensitivity_combined.tsv",
        sep="\t",
        index=False,
    )

    draw_figure(
        combined,
        within,
        meta,
        run_dir,
    )

    summary_path = (
        run_dir
        / "4_2_pam50_sensitivity_summary.txt"
    )

    with summary_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "ANALYSIS_SCOPE="
            "4_2_PAM50_ADJUSTED_AND_"
            "WITHIN_SUBTYPE_SENSITIVITY\n"
        )

        handle.write(
            f"EXPRESSION_FILE="
            f"{expression_path}\n"
        )

        handle.write(
            f"CLINICAL_FILE="
            f"{clinical_path}\n"
        )

        handle.write(
            f"CLINICAL_SAMPLE_FIELD="
            f"{sample_col}\n"
        )

        handle.write(
            f"PAM50_FIELD="
            f"{pam50_col}\n"
        )

        handle.write(
            "TCGA_PRIMARY_SAMPLES="
            "1097;EXPECTED=1097;"
            "STATUS=PASS\n"
        )

        for subtype in [
            "LumA",
            "LumB",
            "HER2",
            "Basal",
            "Normal-like",
        ]:
            label = subtype.replace(
                "-",
                "_",
            )

            count = observed_counts.get(
                subtype,
                0,
            )

            handle.write(
                f"PAM50_{label}_N="
                f"{count}\n"
            )

        handle.write(
            "PAM50_FOUR_SUBTYPE_"
            "ANALYSIS_N=514;"
            "EXPECTED=514;"
            "STATUS=PASS\n"
        )

        handle.write(
            "NORMAL_LIKE_EXCLUDED_N=8\n"
        )

        for row in combined.itertuples(
            index=False
        ):
            handle.write(
                f"SENSITIVITY={row.gene};"
                f"raw_rho={row.raw_rho:.8g};"
                f"raw_FDR={row.raw_FDR:.8g};"
                f"adj_rho="
                f"{row.adjusted_rho:.8g};"
                f"adj_FDR="
                f"{row.adjusted_FDR:.8g};"
                f"meta_rho="
                f"{row.meta_rho:.8g};"
                f"meta_FDR="
                f"{row.meta_FDR:.8g};"
                f"I2="
                f"{row.I2_percent:.6g};"
                f"decision="
                f"{row.adjustment_decision}\n"
            )

        handle.write(
            "ADJUSTED_METHOD="
            "RANK_RESIDUAL_PARTIAL_"
            "SPEARMAN_CONTROLLING_PAM50\n"
        )

        handle.write(
            "META_METHOD="
            "FIXED_EFFECT_FISHER_Z_"
            "ACROSS_LUMA_LUMB_HER2_BASAL\n"
        )

        handle.write(
            "SUBTYPE_MULTIPLE_TESTING="
            "BH_ALL_28_AND_WITHIN_"
            "EACH_SUBTYPE\n"
        )

        handle.write(
            "RAW_AND_FROZEN_INPUTS_"
            "MODIFIED=NO\n"
        )

        handle.write(
            "FOUR_TWO_PAM50_"
            "SENSITIVITY_PASS\n"
        )

    print(
        "ANALYSIS_SCOPE="
        "4_2_PAM50_ADJUSTED_AND_"
        "WITHIN_SUBTYPE_SENSITIVITY"
    )

    print(
        f"EXPRESSION_FILE="
        f"{expression_path}"
    )

    print(
        f"CLINICAL_FILE="
        f"{clinical_path}"
    )

    print(
        "TCGA_PRIMARY_SAMPLES="
        "1097;EXPECTED=1097;"
        "STATUS=PASS"
    )

    for subtype in [
        "LumA",
        "LumB",
        "HER2",
        "Basal",
        "Normal-like",
    ]:
        label = subtype.replace(
            "-",
            "_",
        )

        count = observed_counts.get(
            subtype,
            0,
        )

        print(
            f"PAM50_{label}_N="
            f"{count}"
        )

    print(
        "PAM50_FOUR_SUBTYPE_"
        "ANALYSIS_N=514;"
        "EXPECTED=514;"
        "STATUS=PASS"
    )

    for row in combined.itertuples(
        index=False
    ):
        print(
            f"SENSITIVITY={row.gene};"
            f"raw_rho={row.raw_rho:.8g};"
            f"raw_FDR={row.raw_FDR:.8g};"
            f"adj_rho="
            f"{row.adjusted_rho:.8g};"
            f"adj_FDR="
            f"{row.adjusted_FDR:.8g};"
            f"meta_rho="
            f"{row.meta_rho:.8g};"
            f"meta_FDR="
            f"{row.meta_FDR:.8g};"
            f"I2="
            f"{row.I2_percent:.6g};"
            f"decision="
            f"{row.adjustment_decision}"
        )

    print(
        "RAW_AND_FROZEN_INPUTS_"
        "MODIFIED=NO"
    )

    print(
        "FOUR_TWO_PAM50_"
        "SENSITIVITY_PASS"
    )


if __name__ == "__main__":
    main()
