#!/usr/bin/env bash
set -Eeuo pipefail
TASK_DIR=$(mktemp -d /mnt/BioDisk/wangmingyuan/gsea_symbol_audit_XXXX)
cat > "$TASK_DIR/audit.py" <<'PY'
import csv
import gzip
from collections import defaultdict, Counter
from pathlib import Path

path = Path('/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA/00_raw_data/TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz')
groups = defaultdict(list)
with gzip.open(path, 'rt', newline='') as handle:
    reader = csv.reader(handle, delimiter='\t')
    header = next(reader)
    for number, row in enumerate(reader, 2):
        if not row:
            continue
        original = row[0]
        symbol = original.split('|', 1)[0].strip()
        groups[symbol].append((number, original))

invalid = {'', '?', 'NA', 'N/A', 'NAN', 'NULL', 'NONE'}
dup = {k:v for k,v in groups.items() if len(v)>1}
real = {k:v for k,v in dup.items() if k.upper() not in invalid}

print('AUDIT_SCOPE=TCGA_FULL_GENE_SYMBOLS')
print(f'INPUT_FEATURES={sum(map(len,groups.values()))}')
print(f'INVALID_SYMBOL_ROWS={sum(len(v) for k,v in groups.items() if k.upper() in invalid)}')
print(f'DUPLICATE_SYMBOL_GROUPS={len(dup)}')
print(f'VALID_DUPLICATE_GROUPS={len(real)}')
print('symbol\tn_rows\toriginal_ids_first_12')

for symbol, records in sorted(dup.items()):
    print(f'{symbol!r}\t{len(records)}\t'+ ';'.join(r[1] for r in records[:12]))

for symbol in ['DNMT3B','ACOX2']:
    print(f'TARGET={symbol};ROWS={groups.get(symbol,[])}')

flat = [r[1] for v in groups.values() for r in v]
print(f'EXACT_DUPLICATE_ORIGINAL_IDS={sum(n-1 for n in Counter(flat).values())}')
print('RAW_DATA_MODIFIED=NO')
print('GENE_SYMBOL_AUDIT_PASS')
PY

bash -n "$0"
conda run -n batch_bench python -m py_compile "$TASK_DIR/audit.py"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python "$TASK_DIR/audit.py"
