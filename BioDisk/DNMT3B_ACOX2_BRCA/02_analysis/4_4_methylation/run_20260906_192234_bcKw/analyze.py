import os, json, hashlib, tarfile
from pathlib import Path
import importlib.metadata as im
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root, out = Path(os.environ['ROOT']), Path(os.environ['RUN'])
ready = Path(os.environ.get(
    'READY', str(root/'02_analysis_ready/freeze_20260904_135405')
))
files = {
    'RNA': ready/'tcga_expression_8genes.tsv.gz',
    'Meth': ready/'tcga_methylation_acox2_promoter_3probes.tsv.gz',
    'Primary': ready/'tcga_joint_primary_samples.tsv'
}
probes = ['cg02259384', 'cg13705284', 'cg16587010']

for p in files.values():
    assert p.is_file(), f'Missing input: {p}'

expr = pd.read_csv(files['RNA'], sep='\t', index_col=0)
meth = pd.read_csv(files['Meth'], sep='\t', index_col=0)
samples = pd.read_csv(files['Primary'], sep='\t', dtype=str)

for mat in [expr, meth]:
    assert mat.index.is_unique and mat.columns.is_unique, (
        'Duplicate matrix identifiers'
    )

assert {'sample_id', 'patient_id'}.issubset(samples.columns)
assert samples.sample_id.is_unique and samples.patient_id.is_unique
assert samples[['sample_id', 'patient_id']].notna().all().all()
assert len(samples) == 783, (
    f'Expected 783 primary samples; found {len(samples)}'
)
assert samples.sample_id.str.split('-').str[3].str[:2].eq('01').all()
assert samples.patient_id.eq(samples.sample_id.str[:12]).all()
assert {'DNMT3B', 'ACOX2'}.issubset(expr.index)
assert set(meth.index) == set(probes), (
    'Promoter probe set differs from frozen definition'
)

ids = samples.sample_id.tolist()
assert set(ids).issubset(expr.columns)
assert set(ids).issubset(meth.columns)

e = expr.loc[['DNMT3B', 'ACOX2'], ids].apply(
    pd.to_numeric, errors='raise'
).T
m = meth.loc[probes, ids].apply(
    pd.to_numeric, errors='raise'
).T

assert np.isfinite(e.to_numpy()).all(), 'Nonfinite target expression'
assert not np.isinf(m.to_numpy()).any()
valid = m.to_numpy()[np.isfinite(m.to_numpy())]
assert len(valid) > 0 and ((valid >= 0) & (valid <= 1)).all(), (
    'Beta values outside [0,1]'
)

d = samples.set_index('sample_id').join(e).join(m)
d['ValidCpGs'] = m.notna().sum(axis=1)
d['PromMean'] = m.mean(axis=1).where(d.ValidCpGs >= 2)
d['Complete3'] = d.ValidCpGs.eq(3)
assert d.PromMean.notna().sum() >= 20

print('ANALYSIS_SCOPE=4_4_DNMT3B_ACOX2_PROMOTER', flush=True)
print('PRIMARY_SAMPLES=783;PROMOTER_PROBES=3', flush=True)
print(
    'PROMOTER_DEFINITION=MANIFEST_TSS1500_OR_TSS200;'
    'OBSERVED_ALL_TSS1500', flush=True
)
print(
    'COMPOSITE=MEAN_BETA_MIN_2_OF_3;MISSING_IMPUTATION=NONE',
    flush=True
)
print('PAM50_ADJUSTMENT=NO;NORMAL_COMPARISON=NOT_INCLUDED', flush=True)
print(
    f'COMPOSITE_VALID_N={d.PromMean.notna().sum()};'
    f'COMPLETE_3PROBE_N={d.Complete3.sum()}', flush=True
)
print('ANALYSIS_LOGIC_PASS', flush=True)


def bh(p):
    a = np.asarray(p, float)
    ix = np.argsort(a)
    q = np.empty(len(a))
    q[ix] = np.minimum.accumulate(
        (a[ix] * len(a) / np.arange(1, len(a) + 1))[::-1]
    )[::-1]
    return np.minimum(q, 1)


rows = []


def corr(frame, x, y, family, label):
    z = frame[[x, y]].dropna()
    row = {
        'family': family, 'comparison': label,
        'x': x, 'y': y, 'n': len(z),
        'rho': np.nan, 'P': np.nan
    }
    if len(z) < 20 or z[x].nunique() < 2 or z[y].nunique() < 2:
        row['status'] = 'NOT_ESTIMABLE'
    else:
        r, p = spearmanr(z[x], z[y])
        row.update(rho=float(r), P=float(p), status='OK')
    rows.append(row)


# Primary family: promoter mean plus three individual CpGs.
for y in ['PromMean'] + probes:
    corr(
        d, 'DNMT3B', y,
        'Primary_DNMT_Meth', 'DNMT3B_vs_' + y
    )

# Supporting expression links.
corr(
    d, 'DNMT3B', 'ACOX2',
    'Supporting_Expression', 'DNMT3B_vs_ACOX2'
)
corr(
    d, 'PromMean', 'ACOX2',
    'Supporting_Expression', 'PromMean_vs_ACOX2'
)

for x in probes:
    corr(
        d, x, 'ACOX2',
        'Supporting_CpG_Expression', x + '_vs_ACOX2'
    )

# Same promoter definition, restricted to all three probes observed.
cc = d.loc[d.Complete3]
corr(
    cc, 'DNMT3B', 'PromMean',
    'Complete3_Sensitivity', 'DNMT3B_vs_PromMean_complete3'
)
corr(
    cc, 'PromMean', 'ACOX2',
    'Complete3_Sensitivity', 'PromMean_vs_ACOX2_complete3'
)

r = pd.DataFrame(rows)
r['FDR_BH'] = np.nan
for _, ix in r.groupby('family').groups.items():
    good = [i for i in ix if np.isfinite(r.loc[i, 'P'])]
    if good:
        r.loc[good, 'FDR_BH'] = bh(r.loc[good, 'P'])

r.to_csv(
    out/'correlations.tsv', sep='\t', index=False, na_rep='NA'
)
d.reset_index().to_csv(
    out/'analysis_samples.tsv', sep='\t', index=False, na_rep='NA'
)

qc = []
for p in probes + ['PromMean']:
    v = d[p].dropna()
    qc.append({
        'feature': p, 'n': len(v),
        'missing_n': d[p].isna().sum(),
        'missing_rate': d[p].isna().mean(),
        'min': v.min(), 'median': v.median(),
        'max': v.max(), 'SD': v.std()
    })
pd.DataFrame(qc).to_csv(
    out/'methylation_qc.tsv', sep='\t', index=False
)

(out/'sources.json').write_text(json.dumps({
    k: {
        'file': str(p),
        'sha256': hashlib.sha256(p.read_bytes()).hexdigest()
    }
    for k, p in files.items()
}, indent=2))
(out/'versions.json').write_text(json.dumps({
    p: im.version(p)
    for p in ['pandas', 'numpy', 'scipy', 'matplotlib']
}, indent=2))

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11
})


def scatter(ax, x, y, row):
    z = d[[x, y]].dropna()
    ax.scatter(
        z[x], z[y], s=12, alpha=.35,
        color='#327E88', edgecolors='none'
    )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if row.status == 'OK':
        ax.set_title(
            f'N={int(row.n)}; rho={row.rho:.2f}\n'
            f'P={row.P:.2g}; q={row.FDR_BH:.2g}',
            fontsize=10
        )
    else:
        ax.set_title('Not estim.', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)


fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, label in zip(axes, [
    'DNMT3B_vs_PromMean',
    'PromMean_vs_ACOX2',
    'DNMT3B_vs_ACOX2'
]):
    row = r.loc[r.comparison == label].iloc[0]
    scatter(ax, row.x, row.y, row)
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(
        out/f'Fig4_4_Main.{ext}', dpi=250, bbox_inches='tight'
    )
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, p in zip(axes, probes):
    row = r.loc[r.comparison == 'DNMT3B_vs_' + p].iloc[0]
    scatter(ax, 'DNMT3B', p, row)
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(
        out/f'Fig4_4_CpGs.{ext}', dpi=250, bbox_inches='tight'
    )
plt.close(fig)

lines = [
    'ANALYSIS_SCOPE=4_4_DNMT3B_ACOX2_PROMOTER',
    f'PRIMARY_N={len(d)}',
    f'COMPOSITE_VALID_N={d.PromMean.notna().sum()}',
    f'COMPLETE_3PROBE_N={int(d.Complete3.sum())}',
    'MAIN_TEST=DNMT3B_VS_PROMOTER_MEAN_BETA',
    'PRIMARY_FDR_FAMILY=MEAN_PLUS_3_CPG_TESTS',
    'SECONDARY_FDR_FAMILIES=2_EXPRESSION_LINKS;'
    '3_CPG_EXPRESSION_LINKS;2_COMPLETE3_CHECKS'
]
for _, z in r.iterrows():
    lines.append(
        f'COR={z.comparison};N={z.n};rho={z.rho:.8g};'
        f'P={z.P:.8g};FDR={z.FDR_BH:.8g};STATUS={z.status}'
    )
lines.extend([
    'CAUSAL_DIRECTION=NOT_INFERRED_FROM_CORRELATION',
    'RAW_AND_FROZEN_INPUTS_MODIFIED=NO'
])
if (r.status != 'OK').any():
    lines.append('MODEL_REVIEW_REQUIRED=TRUE')
lines.append('FOUR_FOUR_ANALYSIS_COMPLETE')

(out/'4_4_summary.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines), flush=True)

(out/'figure_captions.txt').write_text(
    'Fig4_4_Main: Spearman associations in primary tumors with '
    'matched RNA and 450K data. PromMean: arithmetic mean beta of '
    'cg02259384, cg13705284 and cg16587010, requiring at least '
    '2 valid probes. DNMT3B and ACOX2: expression on the frozen '
    'source scale. Each panel reports its pairwise complete N, '
    'rho, two-sided asymptotic P, and within-family BH q. '
    'The first panel belongs to the 4-test DNMT3B-methylation '
    'family; the other two panels belong to the 2-test '
    'expression-link family. '
    'Fig4_4_CpGs: DNMT3B versus each promoter CpG beta; q belongs '
    'to the same 4-test primary family as promoter mean. '
    'Complete3 sensitivity uses the identical 3-CpG mean '
    'restricted to complete samples. Probe labels are manifest '
    'identifiers. All probes are annotated TSS1500 for '
    'ACOX2/NM_003500 in the frozen manifest; no new transcript-'
    'version-specific TSS remapping was performed. '
    'No methylation direction or direct regulatory mechanism '
    'is assumed.\n'
)

archive = out/'4_4_return.tar.gz'
with tarfile.open(archive, 'w:gz') as tar:
    for p in sorted(out.iterdir()):
        if p.is_file() and p != archive and p.name != 'run.log':
            tar.add(p, arcname=p.name)

print('RESULT_DIR=' + str(out), flush=True)
