from pathlib import Path
import tarfile, io, json

root = Path('/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA')
src = root / '02_analysis/4_3_gsea/run_20260906_102834_vx9A'
gmt = root / '00_raw_data/reference/MSigDB_2025.1.Hs'
dest = root / '02_analysis/4_3_gsea/4_3_reestimate_inputs.tar.gz'

if not src.is_dir():
    raise SystemExit(f'SOURCE_NOT_FOUND={src}')

exts = {'.rnk', '.tsv', '.txt', '.csv', '.json',
        '.py', '.R', '.r', '.sh', '.gmt', '.log'}
items = []
for p in sorted(src.rglob('*')):
    if p.is_file() and not p.is_symlink():
        if p.suffix in exts or p.name.endswith(('.tsv.gz', '.rnk.gz', '.csv.gz')):
            items.append((p, 'original_run/' + str(p.relative_to(src))))

for name in ['h.all.v2025.1.Hs.symbols.gmt',
             'c2.cp.reactome.v2025.1.Hs.symbols.gmt']:
    p = gmt / name
    if not p.is_file():
        raise SystemExit(f'GMT_NOT_FOUND={p}')
    items.append((p, 'reference/' + name))

print('ANALYSIS_LOGIC_PASS=READ_ONLY_EXPORT')
manifest = [{'archive_path': a, 'bytes': p.stat().st_size} for p, a in items]
part = dest.with_name(dest.name + '.part')

with tarfile.open(part, 'w:gz') as t:
    for p, a in items:
        t.add(p, arcname=a, recursive=False)
    payload = json.dumps(manifest, ensure_ascii=False, indent=2).encode()
    info = tarfile.TarInfo('export_manifest.json')
    info.size = len(payload)
    t.addfile(info, io.BytesIO(payload))

with tarfile.open(part, 'r:gz') as t:
    for m in t:
        if m.isfile():
            h = t.extractfile(m)
            while h.read(1024 * 1024):
                pass

part.replace(dest)
print(f'FILES={len(items)}')
print(f'RETURN_FILE={dest}')
print(f'SIZE_MIB={dest.stat().st_size / 1024**2:.2f}')
print('RAW_DATA_MODIFIED=NO')
print('GSEA_REESTIMATE_EXPORT_PASS')
