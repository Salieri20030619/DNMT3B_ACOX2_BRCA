#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
RAW_DIR="${PROJECT_ROOT}/00_raw_data"
MANIFEST_DIR="${PROJECT_ROOT}/00_manifest"
LOG_DIR="${PROJECT_ROOT}/00_logs"

mkdir -p \
  "${RAW_DIR}/TCGA_BRCA" \
  "${RAW_DIR}/GEO/GSE70947" \
  "${RAW_DIR}/GEO/GSE65194" \
  "${RAW_DIR}/GEO/platforms" \
  "${RAW_DIR}/reference" \
  "${MANIFEST_DIR}" \
  "${LOG_DIR}"

for required_command in wget gzip unzip grep zgrep sha256sum find sort xargs; do
  if ! command -v "${required_command}" >/dev/null 2>&1; then
    echo "ERROR=required_command_not_found:${required_command}" >&2
    exit 1
  fi
done

MANIFEST="${MANIFEST_DIR}/download_manifest.tsv"
cat > "${MANIFEST}" <<'MANIFEST_EOF'
source	dataset	filename	url	purpose
TCGA	TCGA_BRCA	TCGA.BRCA.HiSeqV2.tsv.gz	https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz	4.1-4.3 gene-level RNA-seq expression
TCGA	TCGA_BRCA	TCGA.BRCA.clinical.tsv	https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix	4.1 clinical stage receptor status and PAM50
TCGA	TCGA_BRCA	TCGA.BRCA.HumanMethylation450.tsv.gz	https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HumanMethylation450.gz	4.4 CpG-level methylation beta values
TCGA	TCGA_BRCA	TCGA-CDR-SupplementalTableS1.xlsx	https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81	4.3 curated OS PFI DSS and DFI endpoints
GEO	GSE70947	GSE70947_series_matrix.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70947/matrix/GSE70947_series_matrix.txt.gz	4.1 paired tumor-adjacent validation expression
GEO	GSE65194	GSE65194_series_matrix.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65194/matrix/GSE65194_series_matrix.txt.gz	4.1 TNBC non-TNBC and normal validation expression
GPL	GPL13607	GPL13607_old_annotations.txt.gz	https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13607/suppl/GPL13607_old_annotations.txt.gz	GSE70947 probe annotation
GPL	GPL570	GPL570.annot.gz	https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz	GSE65194 probe annotation
REF	ACOX2	NM_003500.4.fasta	https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_003500.4&rettype=fasta&retmode=text	ACOX2 RefSeq transcript sequence
REF	ACOX2	NM_003500.4.gb	https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_003500.4&rettype=gbwithparts&retmode=text	ACOX2 RefSeq transcript annotation
MANIFEST_EOF

download_one() {
  local url="$1"
  local dest="$2"
  local part="${dest}.part"

  if [[ -s "${dest}" ]]; then
    echo "DOWNLOAD_SKIP=${dest}"
    return 0
  fi

  echo "DOWNLOAD_START=${dest}"
  wget \
    --continue \
    --tries=20 \
    --timeout=60 \
    --read-timeout=120 \
    --retry-connrefused \
    --progress=dot:giga \
    --output-document="${part}" \
    "${url}"
  mv -f "${part}" "${dest}"
  echo "DOWNLOAD_DONE=${dest}"
}

while IFS=$'\t' read -r source dataset filename url purpose; do
  if [[ "${source}" == "source" ]]; then
    continue
  fi

  case "${source}" in
    TCGA)
      destination="${RAW_DIR}/TCGA_BRCA/${filename}"
      ;;
    GEO)
      destination="${RAW_DIR}/GEO/${dataset}/${filename}"
      ;;
    GPL)
      destination="${RAW_DIR}/GEO/platforms/${filename}"
      ;;
    REF)
      destination="${RAW_DIR}/reference/${filename}"
      ;;
    *)
      echo "ERROR=unknown_source:${source}" >&2
      exit 1
      ;;
  esac

  download_one "${url}" "${destination}"
done < "${MANIFEST}"

echo "VALIDATION_START"

while IFS= read -r gz_file; do
  gzip -t "${gz_file}"
  echo "GZIP_PASS=${gz_file}"
done < <(find "${RAW_DIR}" -type f -name '*.gz' | sort)

CDR_FILE="${RAW_DIR}/TCGA_BRCA/TCGA-CDR-SupplementalTableS1.xlsx"
unzip -tqq "${CDR_FILE}"
echo "XLSX_PASS=${CDR_FILE}"

FASTA_FILE="${RAW_DIR}/reference/NM_003500.4.fasta"
GB_FILE="${RAW_DIR}/reference/NM_003500.4.gb"
grep -q 'NM_003500.4' "${FASTA_FILE}"
grep -q 'NM_003500.4' "${GB_FILE}"
grep -q 'ACOX2' "${GB_FILE}"
echo "REFSEQ_PASS=NM_003500.4_ACOX2"

EXPR_FILE="${RAW_DIR}/TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz"
if ! zgrep -m 1 -E '^DNMT3B([[:space:]]|\|)' "${EXPR_FILE}" >/dev/null; then
  echo "ERROR=DNMT3B_not_found_in_TCGA_expression" >&2
  exit 1
fi
if ! zgrep -m 1 -E '^ACOX2([[:space:]]|\|)' "${EXPR_FILE}" >/dev/null; then
  echo "ERROR=ACOX2_not_found_in_TCGA_expression" >&2
  exit 1
fi
echo "TCGA_TARGET_ROWS_PASS=DNMT3B,ACOX2"

CLIN_FILE="${RAW_DIR}/TCGA_BRCA/TCGA.BRCA.clinical.tsv"
grep -m 1 -q '^sampleID' "${CLIN_FILE}"
grep -m 1 -q 'PAM50_mRNA_nature2012' "${CLIN_FILE}"
grep -m 1 -q 'AJCC_Stage_nature2012' "${CLIN_FILE}"
echo "TCGA_CLINICAL_FIELDS_PASS=sampleID,PAM50,AJCC_stage"

GPL13607_FILE="${RAW_DIR}/GEO/platforms/GPL13607_old_annotations.txt.gz"
GPL570_FILE="${RAW_DIR}/GEO/platforms/GPL570.annot.gz"
for platform_file in "${GPL13607_FILE}" "${GPL570_FILE}"; do
  target_hits="${MANIFEST_DIR}/$(basename "${platform_file}").target_hits.txt"
  zgrep -i -E 'DNMT3B|ACOX2' "${platform_file}" > "${target_hits}"
  if ! grep -q -i 'DNMT3B' "${target_hits}"; then
    echo "ERROR=DNMT3B_not_found_in_platform:${platform_file}" >&2
    exit 1
  fi
  if ! grep -q -i 'ACOX2' "${target_hits}"; then
    echo "ERROR=ACOX2_not_found_in_platform:${platform_file}" >&2
    exit 1
  fi
  echo "GEO_PLATFORM_TARGET_PASS=${platform_file}"
done

CHECKSUM_FILE="${MANIFEST_DIR}/sha256sums.txt"
find "${RAW_DIR}" -type f ! -name '*.part' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "${CHECKSUM_FILE}"

du -sh "${RAW_DIR}" | tee "${MANIFEST_DIR}/download_size.txt"
date -Is > "${MANIFEST_DIR}/download_completed_at.txt"

echo "DOWNLOAD_AND_AUDIT_PASS"
