import os,json,hashlib,tarfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

root=Path(os.environ['ROOT']);out=Path(os.environ['RUN'])
source=root/'00_raw_data/TCGA_BRCA/TCGA.BRCA.HiSeqV2.tsv.gz'
genes='ACER2 ACOX2 ADCY5 AKR1C2 ALG9 ALOXE3 ASS1 ATP5F1E ATP6V0A4 B3GNT4 CERCAM CHSY1 COX6C DNMT3B EOGT ERAP1 FAXDC2 FGFR4 GALC GALNT10 GPRC5A GSTM5 HK2 KMT2B LDHA MGST1 MT-ATP6 MT-CO1 MT-CO2 MT-CO3 MT-CYB MT-ND1 MT-ND2 MT-ND4 MT-ND4L MT-ND5 MT-ND6 NDUFA3 NPNT PDE6B PIGT PIP5K1B PKM PNPO PRDM2 PRDM4 PRDM6 PRDX4 QPRT RGS16 SETD5 ST6GALNAC4 SYNJ1 TH UGT1A6 UQCRHL'.split()
genes=list(dict.fromkeys(genes))

def bh(p):
    p=np.asarray(p,float)
    q=np.full(len(p),np.nan)
    ok=np.flatnonzero(np.isfinite(p))
    order=ok[np.argsort(p[ok])]
    n=len(order)
    if n:
        q[order]=np.minimum(
            1,np.minimum.accumulate(
                (p[order]*n/np.arange(1,n+1))[::-1]
            )[::-1]
        )
    return q

print('STEP=LOAD_TCGA',flush=True)
x=pd.read_csv(source,sep='\t',index_col=0)
x.index=x.index.astype(str).str.split('|',regex=False).str[0].str.strip()
x=x.loc[~x.index.str.upper().isin(['?','','NA','NAN','NONE'])]
assert x.index.is_unique,'Duplicate valid gene symbols'

samples=[
    s for s in x
    if len(s.split('-'))>=4 and s.split('-')[3][:2]=='01'
]
assert len(samples)==1097,f'Unexpected primary N={len(samples)}'
assert len(set('-'.join(s.split('-')[:3]) for s in samples))==1097
x=x[samples].apply(pd.to_numeric,errors='coerce')
x=x.replace([np.inf,-np.inf],np.nan)
assert 'DNMT3B' in x.index
ref=x.loc['DNMT3B'].to_numpy(float)
assert np.isfinite(ref).all() and np.std(ref)>0
print('ANALYSIS_LOGIC_PASS;PRIMARY_N=1097',flush=True)

def cor(a,b):
    ok=np.isfinite(a)&np.isfinite(b)
    n=int(ok.sum())
    if n<20 or np.std(a[ok])==0 or np.std(b[ok])==0:
        return n,np.nan,np.nan
    r,p=spearmanr(a[ok],b[ok])
    return n,float(r),float(p)

rows=[]
for i,(g,v) in enumerate(x.iterrows()):
    if g=='DNMT3B':
        continue
    n,r,p=cor(ref,v.to_numpy(float))
    rows.append([g,n,r,p])
    if i%5000==0:
        print(f'GENES_DONE={i}',flush=True)

a=pd.DataFrame(rows,columns=['gene','n','rho','P'])
a['FDR_genome']=bh(a.P)
a.to_csv(out/'tcga_dnmt3b_all_genes.tsv',sep='\t',index=False)

av=pd.DataFrame({
    'gene':genes,
    'present':[g in x.index for g in genes]
})
av.to_csv(out/'gene_availability.tsv',sep='\t',index=False)

b=a[a.gene.isin(genes)].copy()
b['FDR_screen_pool']=bh(b.P)
b.to_csv(out/'tcga_dnmt3b_screen_pool.tsv',sep='\t',index=False)

present=[g for g in genes if g in x.index]
x.loc[present].to_csv(
    out/'tcga_screen_expression.tsv.gz',
    sep='\t',compression='gzip'
)

rows=[]
for i,g in enumerate(present):
    for h in present[i+1:]:
        n,r,p=cor(x.loc[g].to_numpy(float),x.loc[h].to_numpy(float))
        rows.append([g,h,n,r,p])
e=pd.DataFrame(rows,columns=['gene1','gene2','n','rho','P'])
e['FDR_all_pairs']=bh(e.P)
e.to_csv(out/'tcga_screen_all_pairs.tsv',sep='\t',index=False)

pd.DataFrame({
    'sample_id':samples,
    'patient_id':['-'.join(s.split('-')[:3]) for s in samples]
}).to_csv(out/'sample_sheet.tsv',sep='\t',index=False)

hashobj=hashlib.sha256()
with source.open('rb') as f:
    for chunk in iter(lambda:f.read(1048576),b''):
        hashobj.update(chunk)
(out/'source.json').write_text(json.dumps({
    'source':str(source),
    'sha256':hashobj.hexdigest(),
    'normalization':'SOURCE_SCALE_UNCHANGED',
    'scope':'4.2 exploratory rescreen; no 4.3 or 4.4 changes'
},indent=2))

summary=(
    f'PRIMARY_N=1097\n'
    f'SCREEN_REQUESTED={len(genes)}\n'
    f'SCREEN_PRESENT={len(present)}\n'
    f'GENOME_TESTED={int(a.P.notna().sum())}\n'
    'POOL_ROLE=ANNOTATION_PENDING_CANDIDATES_PLUS_CONTROLS\n'
    'RAW_DATA_MODIFIED=NO\n'
    'FOUR_TWO_RESCREEN_TCGA_PASS\n'
)
(out/'summary.txt').write_text(summary)
print(summary,flush=True)

archive=out/'4_2_rescreen_tcga_return.tar.gz'
with tarfile.open(archive,'w:gz') as tf:
    for f in sorted(out.iterdir()):
        if f.is_file() and f!=archive and f.name!='run.log':
            tf.add(f,arcname=f.name)
print('RETURN_ARCHIVE='+str(archive),flush=True)
