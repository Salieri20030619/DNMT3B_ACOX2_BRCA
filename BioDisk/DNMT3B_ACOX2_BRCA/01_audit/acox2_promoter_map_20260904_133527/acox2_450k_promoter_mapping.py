#!/usr/bin/env python3
from __future__ import annotations

import csv
import codecs
import gzip
import hashlib
import os
from pathlib import Path


MANIFEST_FILE = Path(os.environ["MANIFEST_FILE"])
METH_FILE = Path(os.environ["METH_FILE"])
RUN_DIR = Path(os.environ["RUN_DIR"])
TARGET_GENE = "ACOX2"
TARGET_TRANSCRIPT_BASE = "NM_003500"
PROMOTER_GROUPS = {"TSS200", "TSS1500"}
REGULATORY_5P_GROUPS = {"TSS200", "TSS1500", "5'UTR", "1STEXON"}
MIN_MANIFEST_ROWS = int(os.environ.get("MIN_MANIFEST_ROWS", "480000"))
MIN_TCGA_PROBES = int(os.environ.get("MIN_TCGA_PROBES", "480000"))


def write_tsv(name: str, rows: list[dict], fields: list[str]) -> None:
    with (RUN_DIR / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_annotation(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";")]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tcga_probe_ids(path: Path) -> tuple[int, set[str]]:
    probes: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError("TCGA methylation matrix is empty")
        for line in handle:
            probe = line.split("\t", 1)[0].strip().strip('"')
            if probe:
                probes.add(probe)
    return len(probes), probes


def detect_manifest_encoding(path: Path) -> str:
    with path.open("rb") as handle:
        prefix = handle.read(8192)
    if prefix.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return "utf-16"
    if b"\x00" in prefix:
        even_nuls = prefix[0::2].count(0)
        odd_nuls = prefix[1::2].count(0)
        return "utf-16-le" if odd_nuls >= even_nuls else "utf-16-be"
    return "utf-8-sig"


def section_name(line: str) -> str:
    first_field = line.rstrip("\r\n").split(",", 1)[0]
    return first_field.strip().strip('"').lstrip("\ufeff")


def assay_reader(path: Path):
    encoding = detect_manifest_encoding(path)
    handle = path.open("r", encoding=encoding, errors="strict", newline="")
    assay_found = False
    for line in handle:
        if section_name(line) == "[Assay]":
            assay_found = True
            break
    if not assay_found:
        handle.close()
        raise RuntimeError("[Assay] section not found in Illumina manifest")
    reader = csv.DictReader(handle)
    required = {
        "IlmnID",
        "CHR",
        "MAPINFO",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group",
    }
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        fields = ",".join(reader.fieldnames or [])
        handle.close()
        raise RuntimeError(f"Required manifest columns absent; fields={fields}")
    return handle, reader, encoding


print("PY_STEP=LOAD_TCGA_PROBE_IDS", flush=True)
tcga_probe_count, tcga_probes = tcga_probe_ids(METH_FILE)

print("PY_STEP=PARSE_ILLUMINA_MANIFEST", flush=True)
manifest_rows = 0
acox2_by_probe: dict[str, dict] = {}
handle, reader, manifest_encoding = assay_reader(MANIFEST_FILE)
try:
    for row in reader:
        ilmn_id = (row.get("IlmnID") or "").strip()
        if not ilmn_id or ilmn_id.startswith("["):
            break
        manifest_rows += 1
        names = split_annotation(row.get("UCSC_RefGene_Name"))
        accessions = split_annotation(row.get("UCSC_RefGene_Accession"))
        groups = split_annotation(row.get("UCSC_RefGene_Group"))
        indices = [index for index, name in enumerate(names) if name.upper() == TARGET_GENE]
        if not indices:
            continue

        gene_groups: set[str] = set()
        gene_accessions: set[str] = set()
        transcript_groups: set[str] = set()
        for index in indices:
            group = groups[index].upper() if index < len(groups) and groups[index] else "NA"
            accession = accessions[index] if index < len(accessions) and accessions[index] else "NA"
            gene_groups.add(group)
            gene_accessions.add(accession)
            if accession.split(".", 1)[0].upper() == TARGET_TRANSCRIPT_BASE:
                transcript_groups.add(group)

        promoter_groups = gene_groups & PROMOTER_GROUPS
        transcript_promoter_groups = transcript_groups & PROMOTER_GROUPS
        regulatory_groups = gene_groups & REGULATORY_5P_GROUPS
        existing = acox2_by_probe.get(ilmn_id)
        record = {
            "IlmnID": ilmn_id,
            "CHR": row.get("CHR") or "NA",
            "MAPINFO": row.get("MAPINFO") or "NA",
            "Strand": row.get("Strand") or row.get("Forward_Sequence") or "NA",
            "ACOX2_accessions": ";".join(sorted(gene_accessions)),
            "ACOX2_groups": ";".join(sorted(gene_groups)),
            "promoter_groups": ";".join(sorted(promoter_groups)) or "NA",
            "regulatory_5prime_groups": ";".join(sorted(regulatory_groups)) or "NA",
            "NM_003500_groups": ";".join(sorted(transcript_groups)) or "NA",
            "NM_003500_promoter_groups": ";".join(sorted(transcript_promoter_groups)) or "NA",
            "Relation_to_UCSC_CpG_Island": row.get("Relation_to_UCSC_CpG_Island") or "NA",
            "UCSC_CpG_Islands_Name": row.get("UCSC_CpG_Islands_Name") or "NA",
            "present_in_TCGA": ilmn_id in tcga_probes,
        }
        if existing is None:
            acox2_by_probe[ilmn_id] = record
finally:
    handle.close()

all_rows = sorted(acox2_by_probe.values(), key=lambda row: row["IlmnID"])
promoter_rows = [row for row in all_rows if row["promoter_groups"] != "NA"]
regulatory_rows = [row for row in all_rows if row["regulatory_5prime_groups"] != "NA"]
transcript_promoter_rows = [
    row for row in all_rows if row["NM_003500_promoter_groups"] != "NA"
]

fields = [
    "IlmnID",
    "CHR",
    "MAPINFO",
    "Strand",
    "ACOX2_accessions",
    "ACOX2_groups",
    "promoter_groups",
    "regulatory_5prime_groups",
    "NM_003500_groups",
    "NM_003500_promoter_groups",
    "Relation_to_UCSC_CpG_Island",
    "UCSC_CpG_Islands_Name",
    "present_in_TCGA",
]
write_tsv("acox2_450k_all_gene_probes.tsv", all_rows, fields)
write_tsv("acox2_450k_promoter_TSS1500_TSS200.tsv", promoter_rows, fields)
write_tsv("acox2_450k_regulatory_5prime.tsv", regulatory_rows, fields)
write_tsv("acox2_450k_NM003500_promoter.tsv", transcript_promoter_rows, fields)

promoter_present = sum(row["present_in_TCGA"] for row in promoter_rows)
regulatory_present = sum(row["present_in_TCGA"] for row in regulatory_rows)
transcript_promoter_present = sum(row["present_in_TCGA"] for row in transcript_promoter_rows)
mapping_pass = (
    manifest_rows >= MIN_MANIFEST_ROWS
    and tcga_probe_count >= MIN_TCGA_PROBES
    and len(all_rows) > 0
    and len(promoter_rows) > 0
    and promoter_present > 0
)

summary = [
    "AUDIT_SCOPE=ACOX2_450K_PROMOTER_MAPPING",
    f"MANIFEST_FILE={MANIFEST_FILE}",
    f"MANIFEST_SHA256={file_sha256(MANIFEST_FILE)}",
    f"MANIFEST_ENCODING={manifest_encoding}",
    f"MANIFEST_ASSAY_ROWS={manifest_rows}",
    f"MIN_MANIFEST_ROWS_REQUIRED={MIN_MANIFEST_ROWS}",
    f"TCGA_UNIQUE_PROBES={tcga_probe_count}",
    f"MIN_TCGA_PROBES_REQUIRED={MIN_TCGA_PROBES}",
    f"ACOX2_ALL_GENE_PROBES={len(all_rows)}",
    f"ACOX2_PROMOTER_TSS1500_TSS200={len(promoter_rows)}",
    f"ACOX2_PROMOTER_PRESENT_TCGA={promoter_present}",
    f"ACOX2_REGULATORY_5PRIME={len(regulatory_rows)}",
    f"ACOX2_REGULATORY_5PRIME_PRESENT_TCGA={regulatory_present}",
    f"ACOX2_NM003500_PROMOTER={len(transcript_promoter_rows)}",
    f"ACOX2_NM003500_PROMOTER_PRESENT_TCGA={transcript_promoter_present}",
    "PRIMARY_PROMOTER_DEFINITION=TSS1500_OR_TSS200",
    "SENSITIVITY_DEFINITION=TSS1500_OR_TSS200_OR_5UTR_OR_1STEXON",
    "RAW_DATA_MODIFIED=NO_EXCEPT_MANIFEST_ADDITION",
    "ACOX2_PROMOTER_MAPPING_PASS" if mapping_pass else "ACOX2_PROMOTER_MAPPING_FAIL",
]
(RUN_DIR / "acox2_promoter_mapping_summary.txt").write_text(
    "\n".join(summary) + "\n", encoding="utf-8"
)
print("\n".join(summary))
if not mapping_pass:
    raise SystemExit(2)
