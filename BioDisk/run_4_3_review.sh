#!/usr/bin/env bash
set -Eeuo pipefail
BASE=/mnt/BioDisk/wangmingyuan/DNMT3B_ACOX2_BRCA/02_analysis/4_3_gsea
SOURCE="$BASE/run_20260906_102834_vx9A"
RUN=$(mktemp -d "$BASE/review_$(date +%Y%m%d_%H%M%S)_XXXX")
export SOURCE RUN

cat > "$RUN/review.py" <<'PY'
import os, shutil, tarfile
from pathlib import Path
import pandas as pd
import numpy as np

src, out = Path(os.environ['SOURCE']), Path(os.environ['RUN'])
assert 'FOUR_THREE_GSEA_PASS' in (src/'4_3_gsea_summary.txt').read_text()

r = pd.read_csv(src/'gsea_all_results.tsv', sep='\t')
required = {
    'analysis', 'collection', 'Term', 'NES',
    'NOM p-val', 'FDR q-val', 'Lead_genes'
}
assert required <= set(r), f'Missing fields: {required-set(r)}'
assert not r.duplicated(['analysis', 'collection', 'Term']).any()

for col in ['NES', 'NOM p-val', 'FDR q-val']:
    r[col] = pd.to_numeric(r[col], errors='raise')
assert np.isfinite(r[['NES', 'NOM p-val', 'FDR q-val']]).all().all()
assert r['FDR q-val'].between(0, 1).all()

methods = ['HL1097', 'Rho1097', 'Raw514', 'Adj514']
assert set(r.analysis) == set(methods)
lines = []

def log(x):
    print(x, flush=True)
    lines.append(x)

def save(df, name):
    df.to_csv(out/name, sep='\t', index=False)

log('REVIEW_SCOPE=4_3_REACTOME_AND_LEADING_EDGE')
log(f'SOURCE_RUN={src}')
log('ANALYSIS_LOGIC_PASS')

wide = r.pivot(
    index=['collection', 'Term'],
    columns='analysis',
    values=['NES', 'NOM p-val', 'FDR q-val']
)
wide.columns = [f'{a}__{b}' for a, b in wide.columns]
wide = wide.reset_index()
wide['NES_change_514'] = wide['NES__Adj514'] - wide['NES__Raw514']

rawsig = wide['FDR q-val__Raw514'] < .05
adjsig = wide['FDR q-val__Adj514'] < .05
same = np.sign(wide['NES__Raw514']) == np.sign(wide['NES__Adj514'])

wide['PAM50_result'] = np.select(
    [
        rawsig & adjsig & same,
        rawsig & adjsig & ~same,
        rawsig & ~adjsig,
        ~rawsig & adjsig
    ],
    [
        'SIGNIFICANT_SAME_DIRECTION',
        'SIGNIFICANT_OPPOSITE_DIRECTION',
        'NOT_SIGNIFICANT_AFTER_ADJ',
        'SIGNIFICANT_AFTER_ADJ'
    ],
    default='NEITHER_SIGNIFICANT'
)
save(wide, 'all_pathway_comparison.tsv')

# 提取供复核的代谢相关通路；保留原始检验范围与FDR。
pattern = (
    'GLYCOL|GLUCOSE|PYRUVATE|LACTATE|OXIDATIVE_PHOSPHORYLATION|'
    'RESPIRATORY_ELECTRON|TCA_CYCLE|CITRIC_ACID|FATTY_ACID|LIPID|'
    'PEROXISOM|CHOLESTEROL|TRIGLYCERIDE|KETONE|AMINO_ACID|'
    'GLUTAMINE|GLUTAMATE|PENTOSE|NUCLEOTIDE|PURINE|PYRIMIDINE'
)
m = wide[
    (wide.collection == 'Reactome')
    & wide.Term.str.contains(pattern, regex=True)
].copy()
m = m.sort_values(['FDR q-val__HL1097', 'Term'])

save(m, 'reactome_metabolic_comparison.tsv')
save(m.head(20), 'reactome_metabolic_top20.tsv')

log(f'REACTOME_METABOLIC_PATHWAYS={len(m)}')
for method in methods:
    n = (m[f'FDR q-val__{method}'] < .05).sum()
    log(f'REACTOME_METABOLIC_SIGNIFICANT={method};N={n}')

for _, row in m.head(10).iterrows():
    log(
        f'REACTOME={row.Term};'
        f'HL_NES={row["NES__HL1097"]:.4g};'
        f'HL_FDR={row["FDR q-val__HL1097"]:.4g};'
        f'ADJ_NES={row["NES__Adj514"]:.4g};'
        f'ADJ_FDR={row["FDR q-val__Adj514"]:.4g}'
    )

# 提取四套排序中糖酵解通路的leading-edge成员。
h = r[
    (r.collection == 'Hallmark')
    & (r.Term == 'HALLMARK_GLYCOLYSIS')
].copy()
assert len(h) == 4
save(h, 'glycolysis_results.tsv')

edges = {}
for method in methods:
    value = h.loc[h.analysis == method, 'Lead_genes'].iloc[0]
    assert pd.notna(value), f'No leading-edge genes: {method}'
    edges[method] = {
        g.strip() for g in str(value).split(';') if g.strip()
    }
    log(f'GLYCOLYSIS_LEADING_EDGE={method};N={len(edges[method])}')

genes = sorted(set.union(*edges.values()))
e = pd.DataFrame({'gene': genes})
for method in methods:
    e[f'LE_{method}'] = e.gene.isin(edges[method])

e['unadjusted_support_n'] = e[
    ['LE_HL1097', 'LE_Rho1097', 'LE_Raw514']
].sum(axis=1)

g = pd.read_csv(src/'all_gene_statistics.tsv.gz', sep='\t')
assert g.gene.is_unique
cols = [
    'gene', 'mean_diff', 'Welch_t', 'Welch_P', 'Welch_FDR',
    'rho_all', 'rho_514', 'rho_adj514'
]
e = e.merge(g[cols], on='gene', how='left', validate='one_to_one')
assert e.rho_all.notna().all()

e = e.sort_values(
    ['unadjusted_support_n', 'rho_all'],
    ascending=[False, False]
)
save(e, 'glycolysis_leading_edge_genes.tsv')

shared = edges['HL1097'] & edges['Rho1097'] & edges['Raw514']
log('GLYCOLYSIS_SHARED_UNADJUSTED=' + ';'.join(sorted(shared)))
log('NOTE=Adj514_glycolysis_leading_edge_is_descriptive_pathway_FDR_not_significant')
log('NOTE=PKM_gene_level_expression_does_not_resolve_PKM2_isoform')
log('NOTE=Metabolic_subset_uses_keywords_not_an_exhaustive_ontology')
log('NOTE=NES_change_is_descriptive_not_a_formal_between_model_test')
log('NOTE=Original_GSEA_FDR_retained_no_subset_recalculation')
log('SOURCE_FILES_MODIFIED=NO')
log('FOUR_THREE_REVIEW_PASS')

(out/'4_3_review_summary.txt').write_text('\n'.join(lines) + '\n')

# 汇集报告所需结果和现有图件。
for name in [
    '4_3_gsea_summary.txt',
    'gsea_all_results.tsv',
    'geneset_overlap.tsv',
    'figure_pathway_labels.tsv',
    'figure_captions.txt',
    'software_versions.json',
    'source_manifest.tsv'
]:
    shutil.copy2(src/name, out/name)

for p in src.glob('Fig4_3_*'):
    if p.suffix in {'.png', '.pdf'}:
        shutil.copy2(p, out/p.name)

archive = out/'4_3_review_return.tar.gz'
with tarfile.open(archive, 'w:gz') as tar:
    for p in sorted(out.iterdir()):
        if p.is_file() and p != archive:
            tar.add(p, arcname=p.name)

print(f'RETURN_ARCHIVE={archive}', flush=True)
PY

bash -n "$0"
conda run -n batch_bench python -m py_compile "$RUN/review.py"
echo SYNTAX_CHECK=PASS
conda run --no-capture-output -n batch_bench python "$RUN/review.py"
ln -sfnT "$RUN" "$BASE/review_latest"
