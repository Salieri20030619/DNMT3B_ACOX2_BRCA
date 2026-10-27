import os, sys, re, json, hashlib, importlib.metadata
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp

root, out = Path(os.environ['PROJECT']), Path(os.environ['RUN'])
raw = root/'00_raw_data/TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz'
clinical = root/'00_raw_data/TCGA_BRCA/TCGA.BRCA.clinical.tsv'
summary = []
def log(s):
    summary.append(s)
    print(s, flush=True)
def require(ok, message):
    if not ok:
        raise RuntimeError(message)
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()
def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    ans = np.empty(len(p)); ans[order] = np.minimum(q,1)
    return ans

log('ANALYSIS_SCOPE=4_3_DNMT3B_GSEA')
log('STEP=LOAD_FULL_TRANSCRIPTOME')
x = pd.read_csv(raw, sep='\t', index_col=0)
log(f'SOURCE_FEATURES={len(x)};SOURCE_SAMPLES={x.shape[1]}')
require(x.shape == (20530,1218), 'Unexpected source dimensions')
require(x.columns.is_unique, 'Duplicate RNA sample IDs')
x.index = x.index.astype(str).str.split('|').str[0].str.strip()
require(x.index.is_unique, 'Duplicate gene symbols: resolve mapping before analysis')
ids = [s for s in x.columns if re.match(r'^TCGA-[^-]+-[^-]+-01(?:[A-Z]|$|-)',s)]
require(len(ids)==1097, f'Primary samples={len(ids)}, expected 1097')
require(len(set(s[:12] for s in ids))==1097, 'Multiple primary aliquots per patient')
x = x[ids].apply(pd.to_numeric, errors='coerce')
a = x.to_numpy(float)
valid = np.isfinite(a).all(axis=1) & (np.nanstd(a,axis=1)>0)
pd.DataFrame({'gene':x.index,'retained':valid}).to_csv(out/'gene_qc.tsv',sep='\t',index=False)
x = x.loc[valid]
require('DNMT3B' in x.index, 'DNMT3B absent')
dn = x.loc['DNMT3B'].to_numpy(float)
cut = float(np.median(dn)); high = dn > cut
require(abs(cut-6.8534)<1e-5, f'Expression scale mismatch: cutoff={cut}')
require(high.sum()==548 and (~high).sum()==549, 'High/low group mismatch')
log(f'TCGA_PRIMARY=1097;HIGH={high.sum()};LOW={(~high).sum()};CUTOFF={cut}')
log(f'FINITE_VARIABLE_GENES={len(x)};TRANSFORMATION=NONE_SOURCE_SCALE')

c = pd.read_csv(clinical,sep='\t',dtype=str)
require({'sampleID','PAM50_mRNA_nature2012'} <= set(c), 'Clinical fields missing')
require(c.sampleID.is_unique, 'Duplicate clinical sample IDs')
sub = c.set_index('sampleID')['PAM50_mRNA_nature2012'].reindex(ids)
sub = sub.str.strip().replace({'Luminal A':'LumA','Luminal B':'LumB',
                              'HER2-enriched':'HER2','Basal-like':'Basal'})
mask = sub.isin(['LumA','LumB','HER2','Basal']).to_numpy()
require(mask.sum()==514, f'PAM50 four-subtype samples={mask.sum()}, expected 514')
meta = pd.DataFrame({'sample_id':ids,'DNMT3B':dn,'group':np.where(high,'High','Low'),
                     'PAM50':sub.to_numpy(),'sensitivity_included':mask})
meta.to_csv(out/'sample_sheet.tsv',sep='\t',index=False)
log('PAM50_COUNTS='+json.dumps(sub[mask].value_counts().to_dict()))

# Continuous correlations and partial rank correlations.
def rank_corr(values, phenotype, design=None):
    y = stats.rankdata(values,axis=1)
    z = stats.rankdata(phenotype)
    if design is None:
        y -= y.mean(axis=1,keepdims=True); z -= z.mean()
    else:
        q, _ = np.linalg.qr(design,mode='reduced')
        y -= (y@q)@q.T; z -= q@(q.T@z)
    den = np.linalg.norm(y,axis=1)*np.linalg.norm(z)
    ans = np.full(len(values),np.nan)
    np.divide(y@z,den,out=ans,where=den>1e-10)
    return np.clip(ans,-1,1)

a = x.to_numpy(float)
t,p = stats.ttest_ind(a[:,high],a[:,~high],axis=1,equal_var=False)
diff = a[:,high].mean(axis=1)-a[:,~high].mean(axis=1)
rho = rank_corr(a,dn)
design = np.column_stack([np.ones(mask.sum()),pd.get_dummies(sub[mask],drop_first=True).to_numpy(float)])
rho_sub = rank_corr(a[:,mask],dn[mask])
rho_adj = rank_corr(a[:,mask],dn[mask],design)
genes = x.index
universe = (genes!='DNMT3B') & np.isfinite(t) & np.isfinite(rho) & np.isfinite(rho_sub) & np.isfinite(rho_adj)
# Keep one common gene universe for all four rankings.
full = pd.DataFrame({'gene':genes,'mean_diff':diff,'Welch_t':t,'Welch_P':p,
                     'rho_all':rho,'rho_514':rho_sub,'rho_adj514':rho_adj,'GSEA_include':universe})
full['Welch_FDR'] = bh(np.nan_to_num(p,nan=1))
full.to_csv(out/'all_gene_statistics.tsv.gz',sep='\t',index=False)
rankings = {'HL1097':t,'Rho1097':rho,'Raw514':rho_sub,'Adj514':rho_adj}
log(f'COMMON_RANKED_GENES={universe.sum()};DNMT3B_EXCLUDED_FROM_GSEA=YES')
require(universe.sum()>15000,'Too few ranked genes')
log('ANALYSIS_LOGIC_PASS')

results, overlap, curves = [], [], []
for label, scores in rankings.items():
    rnk = pd.DataFrame({'gene':genes[universe],'score':scores[universe]})
    rnk = rnk.sort_values(['score','gene'],ascending=[False,True])
    rnk.to_csv(out/f'{label}.rnk',sep='\t',index=False,header=False)
    log(f'RANK={label};TIED_SCORE_FRACTION={rnk.score.duplicated().mean():.6f}')
    for collection, filename in [('Hallmark','h.all.v2025.1.Hs.symbols.gmt'),
                                  ('Reactome','c2.cp.reactome.v2025.1.Hs.symbols.gmt')]:
        path = root/'00_raw_data/reference/MSigDB_2025.1.Hs'/filename
        sets = {}
        for line in path.read_text().splitlines():
            fields = line.split('\t')
            require(len(fields)>=3,'Invalid GMT content')
            sets[fields[0]] = list(dict.fromkeys(fields[2:]))
        require((len(sets)==50 if collection=='Hallmark' else len(sets)>1000),'Unexpected GMT contents')
        present = set(rnk.gene)
        for term,members in sets.items():
            n = len(set(members)&present)
            overlap.append([label,collection,term,len(members),n,15<=n<=500])
        log(f'STEP=GSEA;RANK={label};COLLECTION={collection}')
        job = out/f'{label}_{collection}'
        pre = gp.prerank(rnk=rnk,gene_sets=sets,outdir=str(job),min_size=15,max_size=500,
                         permutation_num=1000,weight=1,ascending=False,threads=4,
                         seed=20260905,no_plot=True,verbose=False)
        tab = pre.res2d.copy()
        tab.insert(0,'analysis',label); tab.insert(1,'collection',collection)
        for col in ['ES','NES','NOM p-val','FDR q-val','FWER p-val']:
            tab[col] = pd.to_numeric(tab[col],errors='raise')
        require(len(tab)>0 and np.isfinite(tab.NES).all(),'Invalid GSEA output')
        results.append(tab)
        key = 'HALLMARK_GLYCOLYSIS'
        if collection=='Hallmark' and key in pre.results:
            res = pre.results[key]
            curves.append((label,res))
            log(f'GLYCOLYSIS={label};NES={res["nes"]:.5g};P={res["pval"]:.5g};FDR={res["fdr"]:.5g}')
        log(f'GSEA_DONE={label}_{collection};PATHWAYS={len(tab)}')

allres = pd.concat(results,ignore_index=True)
allres.to_csv(out/'gsea_all_results.tsv',sep='\t',index=False)
focus = allres.Term.str.contains('GLYCOL|GLUCOSE|PYRUVATE|LACTATE|OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|LIPID|PEROXISOM',case=False,regex=True)
allres.loc[focus].to_csv(out/'gsea_metabolic_results.tsv',sep='\t',index=False)
pd.DataFrame(overlap,columns=['analysis','collection','term','set_size','overlap','eligible']).to_csv(out/'geneset_overlap.tsv',sep='\t',index=False)

# All Hallmark pathways, without significance-based selection.
h = allres[allres.collection=='Hallmark']
terms = sorted(h.Term.unique())
order = list(rankings)
nes = h.pivot(index='Term',columns='analysis',values='NES').reindex(index=terms,columns=order)
fdr = h.pivot(index='Term',columns='analysis',values='FDR q-val').reindex(index=terms,columns=order)
labels = [f'H{i:02d}' for i in range(1,len(terms)+1)]
pd.DataFrame({'label':labels,'pathway':terms}).to_csv(out/'figure_pathway_labels.tsv',sep='\t',index=False)
fig,ax = plt.subplots(figsize=(7,14))
lim = max(1,float(np.nanmax(np.abs(nes.to_numpy()))))
im = ax.imshow(nes.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
ax.set_xticks(range(4),order); ax.set_yticks(range(len(terms)),labels)
ax.set_title('Hallmark')
for i in range(len(terms)):
    for j in range(4):
        if np.isfinite(nes.iloc[i,j]):
            ax.text(j,i,f'{nes.iloc[i,j]:.1f}'+('*' if fdr.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=7)
fig.colorbar(im,ax=ax,label='NES'); fig.tight_layout()
for ext in ['png','pdf']:
    fig.savefig(out/f'Fig4_3_Hallmark.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
fig,axes = plt.subplots(2,2,figsize=(10,7))
for ax,(label,res) in zip(axes.ravel(),curves):
    ax.plot(res['RES'],color='#2b7a78',lw=1.5)
    ax.axhline(0,color='gray',lw=.7)
    ax.set_title(f'{label}\nNES={res["nes"]:.2f}; FDR={res["fdr"]:.3g}')
    ax.set_xlabel('Gene rank'); ax.set_ylabel('Run. ES')
    lo,hi=ax.get_ylim()
    ax.vlines(res['hits'],lo,lo+(hi-lo)*.035,color='black',lw=.3)
fig.suptitle('Glycolysis'); fig.tight_layout()
for ext in ['png','pdf']:
    fig.savefig(out/f'Fig4_3_Glycolysis.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
(out/'figure_captions.txt').write_text(
    'Hallmark: all tested Hallmark pathways; H01-H50 mapped in figure_pathway_labels.tsv.\n'
    'HL1097: high-versus-low Welch t ranking in 1097 primary tumors.\n'
    'Rho1097: continuous Spearman ranking in 1097 primary tumors.\n'
    'Raw514 and Adj514: unadjusted and PAM50-adjusted rank correlations in the SAME 514 tumors.\n'
    'NES: normalized enrichment score; Run. ES: running enrichment score.\n'
    'Positive NES: enrichment toward higher DNMT3B/group High. Stars: GSEA FDR<0.05.\n'
    'GSEA FDR is the native permutation-derived q-value, within each collection and ranking.\n'
    'Permutation count=1000, seed=20260905; nominal P=0 means below permutation resolution.\n'
    'Ranks use the common finite gene universe; DNMT3B is removed.\n',encoding='utf-8')
versions = {p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib','gseapy']}
(out/'software_versions.json').write_text(json.dumps(versions,indent=2))
files = [raw,clinical]+list((root/'00_raw_data/reference/MSigDB_2025.1.Hs').glob('*.gmt'))
pd.DataFrame([{'path':str(p),'sha256':sha(p)} for p in files]).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
log('GSEA_PERMUTATIONS=1000;MIN_SIZE=15;MAX_SIZE=500;SEED=20260905')
log('MULTIPLE_TESTING=NATIVE_GSEA_FDR_WITHIN_COLLECTION_AND_RANKING')
log('RAW_EXPRESSION_MODIFIED=NO')
log('FOUR_THREE_GSEA_PASS')
(out/'4_3_gsea_summary.txt').write_text('\n'.join(summary)+'\n')
