#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


SOURCE_TSV = Path(os.environ["SOURCE_TSV"])
RUN_DIR = Path(os.environ["RUN_DIR"])


def pair_key(title: str) -> str:
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


def classify_gse70947(text: str) -> tuple[str, str, str, str]:
    if re.search(r"adjacent|normal|healthy|non[- _]?tumou?r", text):
        return "Tissue", "Normal", "Normal", "include"
    if re.search(r"adenocarcinoma|tumou?r|cancer", text):
        return "Tissue", "Tumor", "NA", "include"
    return "Unresolved", "Unresolved", "Unresolved", "exclude_unresolved"


def classify_gse65194(text: str) -> tuple[str, str, str, str]:
    if re.search(r"cell[ _-]?line|cultured cells", text):
        return "Cell_line", "Model", "Cell_line", "exclude_cell_line"
    if re.search(r"normal breast|normal tissue|healthy breast", text):
        return "Tissue", "Normal", "Normal", "include"
    if re.search(r"tnbc|triple[- _]?negative", text):
        return "Tissue", "Tumor", "TNBC", "include"
    if re.search(r"luminal[ _-]*a|\bluma\b", text):
        return "Tissue", "Tumor", "LumA", "include"
    if re.search(r"luminal[ _-]*b|\blumb\b", text):
        return "Tissue", "Tumor", "LumB", "include"
    if re.search(r"her[- _]?2|her2", text):
        return "Tissue", "Tumor", "HER2", "include"
    return "Unresolved", "Unresolved", "Unresolved", "exclude_unresolved"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


with SOURCE_TSV.open("r", encoding="utf-8", newline="") as handle:
    source_rows = list(csv.DictReader(handle, delimiter="\t"))

required = {"dataset", "sample", "title", "combined_metadata"}
if not source_rows or not required.issubset(source_rows[0]):
    raise RuntimeError("Required columns are absent from GEO metadata audit")

frozen: list[dict] = []
for row in source_rows:
    dataset = row["dataset"]
    text = f"{row['title']} | {row['combined_metadata']}".lower()
    if dataset == "GSE70947":
        material, disease_group, subtype, decision = classify_gse70947(text)
        key = pair_key(row["title"])
    elif dataset == "GSE65194":
        material, disease_group, subtype, decision = classify_gse65194(text)
        key = ""
    else:
        material, disease_group, subtype, decision = (
            "Unresolved",
            "Unresolved",
            "Unresolved",
            "exclude_unresolved",
        )
        key = ""
    frozen.append(
        {
            "dataset": dataset,
            "sample": row["sample"],
            "material": material,
            "disease_group": disease_group,
            "subtype": subtype,
            "pair_key": key,
            "analysis_decision": decision,
            "title": row["title"],
            "combined_metadata": row["combined_metadata"],
        }
    )

fields = [
    "dataset",
    "sample",
    "material",
    "disease_group",
    "subtype",
    "pair_key",
    "analysis_decision",
    "title",
    "combined_metadata",
]
write_tsv(RUN_DIR / "geo_frozen_sample_sheet.tsv", frozen, fields)

counts = Counter(
    (row["dataset"], row["material"], row["disease_group"], row["subtype"])
    for row in frozen
)
count_rows = [
    {
        "dataset": key[0],
        "material": key[1],
        "disease_group": key[2],
        "subtype": key[3],
        "n": value,
    }
    for key, value in sorted(counts.items())
]
write_tsv(
    RUN_DIR / "geo_frozen_group_counts.tsv",
    count_rows,
    ["dataset", "material", "disease_group", "subtype", "n"],
)

pairs: dict[str, list[dict]] = defaultdict(list)
for row in frozen:
    if row["dataset"] == "GSE70947" and row["pair_key"]:
        pairs[row["pair_key"]].append(row)

pair_rows: list[dict] = []
valid_pairs = 0
for key, rows in sorted(pairs.items()):
    groups = {row["disease_group"] for row in rows}
    valid = len(rows) == 2 and groups == {"Tumor", "Normal"}
    valid_pairs += int(valid)
    pair_rows.append(
        {
            "pair_key": key,
            "tumor_sample": next(
                (row["sample"] for row in rows if row["disease_group"] == "Tumor"), "NA"
            ),
            "normal_sample": next(
                (row["sample"] for row in rows if row["disease_group"] == "Normal"), "NA"
            ),
            "valid_pair": valid,
        }
    )
write_tsv(
    RUN_DIR / "gse70947_frozen_pairs.tsv",
    pair_rows,
    ["pair_key", "tumor_sample", "normal_sample", "valid_pair"],
)

observed = {
    "TOTAL_SAMPLES": len(frozen),
    "GSE70947_TUMOR": counts[("GSE70947", "Tissue", "Tumor", "NA")],
    "GSE70947_NORMAL": counts[("GSE70947", "Tissue", "Normal", "Normal")],
    "GSE70947_VALID_PAIRS": valid_pairs,
    "GSE65194_TNBC": counts[("GSE65194", "Tissue", "Tumor", "TNBC")],
    "GSE65194_HER2": counts[("GSE65194", "Tissue", "Tumor", "HER2")],
    "GSE65194_LUMA": counts[("GSE65194", "Tissue", "Tumor", "LumA")],
    "GSE65194_LUMB": counts[("GSE65194", "Tissue", "Tumor", "LumB")],
    "GSE65194_NORMAL": counts[("GSE65194", "Tissue", "Normal", "Normal")],
    "GSE65194_CELL_LINE": counts[("GSE65194", "Cell_line", "Model", "Cell_line")],
    "UNRESOLVED": sum(row["material"] == "Unresolved" for row in frozen),
    "INCLUDED_TISSUES": sum(row["analysis_decision"] == "include" for row in frozen),
}
expected = {
    "TOTAL_SAMPLES": 474,
    "GSE70947_TUMOR": 148,
    "GSE70947_NORMAL": 148,
    "GSE70947_VALID_PAIRS": 148,
    "GSE65194_TNBC": 55,
    "GSE65194_HER2": 39,
    "GSE65194_LUMA": 29,
    "GSE65194_LUMB": 30,
    "GSE65194_NORMAL": 11,
    "GSE65194_CELL_LINE": 14,
    "UNRESOLVED": 0,
    "INCLUDED_TISSUES": 460,
}

summary = [
    f"SOURCE_TSV={SOURCE_TSV}",
    "FREEZE_POLICY=ALL_PATIENT_TISSUES_INCLUDED;CELL_LINES_EXCLUDED",
]
all_pass = True
for metric, expected_value in expected.items():
    value = observed[metric]
    status = "PASS" if value == expected_value else "FAIL"
    all_pass = all_pass and status == "PASS"
    summary.append(f"{metric}={value};EXPECTED={expected_value};STATUS={status}")
summary.extend(
    [
        "RAW_DATA_MODIFIED=NO",
        "GEO_SAMPLE_FREEZE_PASS" if all_pass else "GEO_SAMPLE_FREEZE_FAIL",
    ]
)

(RUN_DIR / "geo_sample_freeze_summary.txt").write_text(
    "\n".join(summary) + "\n", encoding="utf-8"
)
print("\n".join(summary))
if not all_pass:
    raise SystemExit(2)
