#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA"
AUDIT_ROOT="${PROJECT_ROOT}/01_audit"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${AUDIT_ROOT}/run_${STAMP}"
LOG_FILE="${RUN_DIR}/audit.log"
PYTHON_SCRIPT="${RUN_DIR}/audit_dnmt3b_acox2_brca.py"
SOURCE_SCRIPT="${1:-/mnt/BioDisk/wangmingyuan/audit_dnmt3b_acox2_brca.py}"

mkdir -p "${RUN_DIR}"

if [[ ! -s "${SOURCE_SCRIPT}" ]]; then
  echo "ERROR=python_audit_script_missing:${SOURCE_SCRIPT}" >&2
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR=conda_not_found" >&2
  exit 1
fi

cp -p "${SOURCE_SCRIPT}" "${PYTHON_SCRIPT}"
export RUN_DIR
export PROJECT_ROOT

{
  echo "===== AUDIT START ====="
  echo "STAMP=${STAMP}"
  echo "HOST=$(hostname)"
  echo "PROJECT_ROOT=${PROJECT_ROOT}"
  echo "RUN_DIR=${RUN_DIR}"

  echo "STEP=PYTHON_SYNTAX_CHECK"
  conda run -n batch_bench python -m py_compile "${PYTHON_SCRIPT}"
  echo "PYTHON_COMPILE_PASS"

  echo "STEP=RUN_AUDIT"
  conda run --no-capture-output -n batch_bench python "${PYTHON_SCRIPT}"

  if ! grep -qx 'AUDIT_PASS' "${RUN_DIR}/audit_summary.txt"; then
    echo "ERROR=audit_success_marker_missing" >&2
    exit 1
  fi

  echo "LATEST_RUN=${RUN_DIR}" > "${AUDIT_ROOT}/latest_run.txt"
  echo "AUDIT_OUTPUT_FILES=$(find "${RUN_DIR}" -maxdepth 1 -type f | wc -l)"
  echo "AUDIT_AND_LOGIC_PASS"
  echo "===== AUDIT END ====="
} 2>&1 | tee "${LOG_FILE}"
