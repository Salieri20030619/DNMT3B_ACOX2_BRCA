import os, json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

root, out = Path(os.environ['ROOT']), Path(os.environ['RUN'])
expr = root/'00_raw_data/TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz'
ref = root/'00_raw_data/reference/MSigDB_2025.1.Hs'
files = [ref/'h.all.v2025.1.Hs.symbols.gmt',
         ref/'c2.cp.reactome.v2025.1.Hs.symbols.gmt']
terms = [
'HALLMARK_GLYCOLYSIS','HALLMARK_FATTY_ACID_METABOLISM',
'HALLMARK_CHOLESTEROL_HOMEOSTASIS','HALLMARK_BILE_ACID_METABOLISM',
'HALLMARK_PEROXISOME','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
'HALLMARK_XENOBIOTIC_METABOLISM',
'REACTOME_METABOLISM_OF_POLYAMINES',
'REACTOME_METABOLISM_OF_NUCLEOTIDES',
'REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES',
'REACTOME_METABOLISM_OF_FOLATE_AND_PTERINES',
'REACTOME_GLYCOSAMINOGLYCAN_METABOLISM',
'REACTOME_GLYCOSPHINGOLIPID_METABOLISM',
'REACTOME_GLUTATHIONE_CONJUGATION',
'REACTOME_CHOLESTEROL_BIOSYNTHESIS',
'REACTOME_PEROXISOMAL_LIPID_METABOLISM']
sets = {}
for f in files:
    with f.open() as h:
        for line in h:
            a = line.rstrip('\n').split('\t')
            if a[0] in terms:
                sets[a[0]] = set(a[2:])
assert set(sets) == set(terms), f'Missing gene sets: {set(terms)-set(sets)}'
(out/'network_definition.json').write_text(json.dumps({
    'terms':terms,
    'selection':'Fixed annotation sets, no DE or correlation filtering',
    'target':'DNMT3B added to each network',
    'subsampling':100,'fraction':0.8,'seed':20260907,
    'primary_weight':'abs(Spearman)^2',
    'candidate_rule':
    'percentile>=90 and >=80% subsamples percentile>=90; exploratory convention'
}, indent=2))

print('STEP=LOAD_EXPRESSION', flush=True)
x = pd.read_csv(expr, sep='\t', index_col=0)
x.index = x.index.astype(str).str.split('|', regex=False).str[0]
x = x.loc[~x.index.str.upper().isin(['?','NA','NAN','NONE',''])]
assert x.index.is_unique
cols = [c for c in x
        if len(c.split('-'))>=4 and c.split('-')[3][:2]=='01']
assert len(cols)==1097 and len({c[:12] for c in cols})==1097
x = x[cols].apply(pd.to_numeric, errors='coerce')
universe = set().union(*sets.values())
requested = sorted(universe|{'DNMT3B'})
audit = pd.DataFrame({'gene':requested})
audit['present'] = audit.gene.isin(x.index)
valid = x.notna().all(axis=1) & np.isfinite(x).all(axis=1) & (x.std(axis=1)>0)
x = x.loc[valid & x.index.isin(requested)]
audit['usable'] = audit.gene.isin(x.index)
audit.to_csv(out/'gene_availability.tsv', sep='\t', index=False)
assert 'DNMT3B' in x.index and len(x)>100
pd.DataFrame({'sample_id':cols}).to_csv(
    out/'samples.tsv',sep='\t',index=False)
lookup = {s:i for i,s in enumerate(x.index)}
networks = {'MetUnion':list(range(len(x)))}
coverage=[]
for name,genes in sets.items():
    ids = [lookup[g] for g in sorted(genes|{'DNMT3B'}) if g in lookup]
    coverage.append({
        'network':name,'annotated':len(genes),
        'usable_metabolic':len(genes & set(x.index)),
        'target_native_member':'DNMT3B' in genes})
    if len(ids)>=16:
        networks[name] = ids
pd.DataFrame(coverage).to_csv(out/'coverage.tsv',sep='\t',index=False)
print(f'PRIMARY_N={len(cols)};UNION_GENES={len(x)};NETWORKS={len(networks)}',
      flush=True)
print('ANALYSIS_LOGIC_PASS',flush=True)
a = x.to_numpy(dtype=float)

def corr(a):
    r=rankdata(a,axis=1)
    r-=r.mean(axis=1,keepdims=True)
    norms=np.sqrt((r*r).sum(axis=1,keepdims=True))
    if (norms==0).any():
        raise RuntimeError('Constant feature in subsample; no silent imputation')
    r/=norms
    c=np.clip(r@r.T,-1,1)
    np.fill_diagonal(c,0)
    return c

def measures(c):
    ab=np.abs(c)
    return {
        'AbsRho2':(ab**2).sum(axis=1),
        'AbsRho1':ab.sum(axis=1),
        'AbsRho6':(ab**6).sum(axis=1),
        'Degree020':(ab>=.2).sum(axis=1),
        'Degree030':(ab>=.3).sum(axis=1)}

def percentile(v,j):
    return float((rankdata(v,method='average')[j]-1)/(len(v)-1)*100)

c=corr(a)
summary=[]
ranks=[]
for name,ids in networks.items():
    genes=x.index[ids]
    j=list(genes).index('DNMT3B')
    cc=c[np.ix_(ids,ids)]
    for method,v in measures(cc).items():
        summary.append({
            'network':name,'method':method,'n_genes':len(ids),
            'target_connectivity':v[j],
            'target_percentile':percentile(v,j),
            'target_rank':float(rankdata(-v)[j])})
        for g,score,r in zip(genes,v,rankdata(-v)):
            ranks.append({
                'network':name,'method':method,'gene':g,
                'connectivity':score,'rank':r})
    if name=='MetUnion':
        pd.DataFrame({'gene':genes,'rho':cc[j]}).query(
            'gene != "DNMT3B"').to_csv(
            out/'dnmt3b_union_correlations.tsv',sep='\t',index=False)
pd.DataFrame(ranks).to_csv(
    out/'all_node_rankings.tsv.gz',sep='\t',index=False,compression='gzip')
pd.DataFrame(summary).to_csv(
    out/'centrality_sensitivity.tsv',sep='\t',index=False)

print('STEP=SUBSAMPLING_100',flush=True)
rng=np.random.default_rng(20260907)
boot=[]
constant_audit=[]
for b in range(100):
    idx=rng.choice(a.shape[1],int(.8*a.shape[1]),replace=False)
    sub=a[:,idx]
    variable=np.ptp(sub,axis=1)>0
    assert variable[lookup['DNMT3B']], 'DNMT3B constant in subsample'
    retained=np.flatnonzero(variable)
    remap={int(old):new for new,old in enumerate(retained)}
    cb=corr(sub[variable])
    dropped=list(x.index[~variable])
    constant_audit.append({
        'iteration':b+1,'dropped_n':len(dropped),
        'dropped_genes':';'.join(dropped)})
    for name,ids in networks.items():
        available=[old for old in ids if old in remap]
        assert len(available)>=16, f'Too few variable genes: {name}'
        reduced=[remap[old] for old in available]
        j=available.index(lookup['DNMT3B'])
        v=(cb[np.ix_(reduced,reduced)]**2).sum(axis=1)
        boot.append({
            'iteration':b+1,'network':name,
            'percentile':percentile(v,j),
            'network_original_n':len(ids),
            'network_variable_n':len(available),
            'network_constant_n':len(ids)-len(available)})
    if (b+1)%10==0:
        print(f'SUBSAMPLING={b+1}/100;CONSTANT_GENES={len(dropped)}',flush=True)
pd.DataFrame(constant_audit).to_csv(
    out/'subsample_constant_genes.tsv',sep='\t',index=False)
bt=pd.DataFrame(boot)
assert bt.groupby('network').size().eq(100).all()
bt.to_csv(out/'subsampling.tsv',sep='\t',index=False)
result=pd.DataFrame(summary).query('method=="AbsRho2"').copy()
for i,r in result.iterrows():
    vals=bt.loc[bt.network==r.network,'percentile']
    freq=float((vals>=90).mean())
    result.loc[i,'top10_frequency']=freq
    result.loc[i,'percentile_p025']=vals.quantile(.025)
    result.loc[i,'percentile_p975']=vals.quantile(.975)
    result.loc[i,'decision']=(
        'STABLE_CANDIDATE'
        if r.target_percentile>=90 and freq>=.8
        else 'HUB_NOT_SUPPORTED_BY_RULE')
result.to_csv(out/'hub_summary.tsv',sep='\t',index=False)
with (out/'summary.txt').open('w') as h:
    h.write('SCOPE=FIXED_METABOLIC_UNION_AND_PATHWAY_NETWORKS\nPRIMARY_N=1097\n')
    h.write(result.to_string(index=False)+'\n')
    h.write('INTERPRETATION=COEXPRESSION_CONNECTIVITY;NOT_CAUSAL_REGULATION\n')
    h.write('RAW_DATA_MODIFIED=NO\nHUB_TEST_PASS\n')
source=[]
for f in [expr,*files]:
    sha=hashlib.sha256()
    with f.open('rb') as h:
        for chunk in iter(lambda:h.read(1024*1024),b''):
            sha.update(chunk)
    source.append({'path':str(f),'sha256':sha.hexdigest()})
pd.DataFrame(source).to_csv(out/'sources.tsv',sep='\t',index=False)
print(result.to_string(index=False),flush=True)
print('HUB_TEST_PASS',flush=True)
