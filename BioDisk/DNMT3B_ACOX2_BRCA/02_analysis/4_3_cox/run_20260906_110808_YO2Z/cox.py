import os, re, json, warnings, hashlib, tarfile
from pathlib import Path
import importlib.metadata as im
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

root,out=Path(os.environ['ROOT']),Path(os.environ['RUN'])
sources={
 'expression':root/'02_analysis/4_3_gsea/run_20260906_102834_vx9A/sample_sheet.tsv',
 'cdr':root/'02_analysis_ready/freeze_20260904_135405/tcga_cdr_brca.tsv',
 'clinical':root/'00_raw_data/TCGA_BRCA/TCGA.BRCA.clinical.tsv'
}
notes=[]

def log(s):
    print(s,flush=True)
    notes.append(s)

def save(d,n):
    d.to_csv(out/n,sep='\t',index=False)

def col(d,options):
    for k in options:
        if k in d:
            return k
    raise ValueError(f'Missing {options}; available={list(d.columns)}')

def num(s):
    return pd.to_numeric(s,errors='coerce')

def bh(p):
    p=np.asarray(p,float)
    order=np.argsort(p)
    z=np.minimum.accumulate(
        (p[order]*len(p)/np.arange(1,len(p)+1))[::-1]
    )[::-1]
    ans=np.empty(len(p))
    ans[order]=np.minimum(1,z)
    return ans

a=pd.read_csv(sources['expression'],sep='\t',dtype={'sample_id':str})
b=pd.read_csv(sources['cdr'],sep='\t',dtype=str)
c=pd.read_csv(sources['clinical'],sep='\t',dtype=str)

assert len(a)==1097 and a.sample_id.is_unique
assert a.sample_id.str.match(r'^TCGA-[^-]+-[^-]+-01').all()
a['patient_id']=a.sample_id.str[:12]
assert a.patient_id.is_unique

patient=col(b,['bcr_patient_barcode','patient_id','submitter_id'])
age=col(b,['age_at_initial_pathologic_diagnosis','age_at_diagnosis','age'])
assert age=='age_at_initial_pathologic_diagnosis', (
    'Age unit requires explicit verification'
)

b['patient_id']=b[patient].str.strip().str.upper()
assert b.patient_id.is_unique
b['age_years']=num(b[age])
assert b.age_years.dropna().between(0,120).all(), (
    'Unexpected CDR age unit/range'
)

for ep in ['OS','PFI']:
    assert ep in b and ep+'.time' in b, f'Missing CDR endpoint {ep}'
    b[ep]=num(b[ep])
    b[ep+'.time']=num(b[ep+'.time'])

stagefield=col(c,['AJCC_Stage_nature2012'])
assert c.sampleID.is_unique
stageraw=c.set_index('sampleID')[stagefield]
a['stage_raw']=a.sample_id.map(stageraw)

def stage(v):
    s=str(v).strip().upper()
    s=re.sub(r'^(STAGE|STG)\s*','',s).replace(' ','')
    m=re.fullmatch(r'(IV|III|II|I)[ABC]?',s)
    return m.group(1) if m else np.nan

a['stage']=a.stage_raw.map(stage)
save(
    a.groupby(['stage_raw','stage'],dropna=False)
     .size().reset_index(name='n'),
    'stage_mapping.tsv'
)

d=a.merge(
    b[['patient_id','age_years','OS','OS.time','PFI','PFI.time']],
    on='patient_id',how='left',validate='one_to_one',indicator=True
)
assert (d['_merge']=='both').all(), 'RNA patients missing in CDR'
d['DNMT3B']=num(d.DNMT3B)
assert np.isfinite(d.DNMT3B).all()

sd=d.DNMT3B.std(ddof=1)
mu=d.DNMT3B.mean()
assert sd>0
d['DNMT_SD']=(d.DNMT3B-mu)/sd
d['Age10']=d.age_years/10

log('ANALYSIS_SCOPE=4_3_DNMT3B_COX_OS_PRIMARY_PFI_SECONDARY')
log(f'PATIENTS={len(d)};EXPRESSION_MEAN={mu:.8g};EXPRESSION_SD={sd:.8g}')
log(f'AGE_SOURCE=CDR:{age};STAGE_SOURCE={stagefield}')
log('MODEL=CONTINUOUS_DNMT3B_PER_SD_PLUS_AGE_PER_10YEARS_PLUS_CATEGORICAL_STAGE')
log('STAGE_REFERENCE=I;PAM50_INCLUDED=NO;TIME_UNIT=DAYS;PENALIZER=0')
log('ANALYSIS_LOGIC_PASS')

rows=[]
audits=[]
phrows=[]
statuses=[]

for ep in ['OS','PFI']:
    valid=(
        d[ep].isin([0,1])
        & np.isfinite(d[ep+'.time'])
        & (d[ep+'.time']>0)
    )
    d[ep+'_eligible']=valid
    v=d.loc[valid].copy()
    cc=v.dropna(subset=['age_years','stage']).copy()
    assert 'I' in set(cc.stage), 'Stage I reference absent'

    log(
        f'ENDPOINT={ep};VALID_N={len(v)};EVENTS={int(v[ep].sum())};'
        f'CC_N={len(cc)};CC_EVENTS={int(cc[ep].sum())};'
        f'EXCLUDED_TIME_OR_STATUS={(~valid).sum()}'
    )

    for st,part in cc.groupby('stage'):
        audits.append({
            'endpoint':ep,'stage':st,
            'n':len(part),'events':int(part[ep].sum())
        })

    stagevars=[]
    for st in ['II','III','IV']:
        if st in set(cc.stage):
            name={'II':'StageII','III':'StageIII','IV':'StageIV'}[st]
            cc[name]=(cc.stage==st).astype(float)
            stagevars.append(name)

    models=[
        ('UniAll',v,['DNMT_SD']),
        ('UniCC',cc,['DNMT_SD']),
        ('AgeCC',cc,['Age10']),
        ('StageCC',cc,stagevars),
        ('Multi',cc,['DNMT_SD','Age10']+stagevars)
    ]

    for name,frame,covs in models:
        fitdata=frame[[ep+'.time',ep]+covs].rename(
            columns={ep+'.time':'time',ep:'event'}
        ).astype(float)

        n,events=len(fitdata),int(fitdata.event.sum())
        status={
            'endpoint':ep,'model':name,'n':n,'events':events,
            'parameters':len(covs),
            'events_per_parameter':events/max(1,len(covs))
        }

        try:
            assert covs and events>0 and n>len(covs)+1, (
                'Insufficient model information'
            )
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                fit=CoxPHFitter(penalizer=0)
                fit.fit(fitdata,duration_col='time',event_col='event')
                ss=fit.summary

            warningtext=' | '.join(str(w.message) for w in caught)
            status['warnings']=warningtext
            status['status']='REVIEW_WARNING' if warningtext else 'OK'
            status['model_LR_P']=fit.log_likelihood_ratio_test().p_value
            assert np.isfinite(ss[['coef','se(coef)','p']]).all().all(), (
                'Nonfinite estimate'
            )

            for term,z in ss.iterrows():
                rows.append({
                    'endpoint':ep,'model':name,'term':term,
                    'n':n,'events':events,
                    'HR':z['exp(coef)'],
                    'CI_low':z['exp(coef) lower 95%'],
                    'CI_high':z['exp(coef) upper 95%'],
                    'coef':z['coef'],'SE':z['se(coef)'],
                    'P':z['p'],'status':status['status']
                })

            try:
                ph=proportional_hazard_test(
                    fit,fitdata,time_transform='rank'
                ).summary
                for term,z in ph.iterrows():
                    phrows.append({
                        'endpoint':ep,'model':name,
                        'term':term,'PH_P':z['p']
                    })
            except Exception as exc:
                status['PH_error']=str(exc)

            log(
                f'COX={ep};MODEL={name};N={n};EVENTS={events};'
                f'STATUS={status["status"]}'
            )
        except Exception as exc:
            status['status']='FAILED'
            status['error']=str(exc)
            log(f'COX={ep};MODEL={name};STATUS=FAILED;ERROR={exc}')

        statuses.append(status)

save(d,'survival_sample_audit.tsv')
save(pd.DataFrame(audits),'stage_event_counts.tsv')
save(pd.DataFrame(statuses),'model_status.tsv')
save(pd.DataFrame(phrows),'proportional_hazards_tests.tsv')

r=pd.DataFrame(rows)
assert len(r)>0, 'No models fitted'
save(r,'cox_all_coefficients.tsv')

focus=r[r.term=='DNMT_SD'].copy()
focus['FDR_endpoints']=np.nan
for model,ix in focus.groupby('model').groups.items():
    focus.loc[ix,'FDR_endpoints']=bh(focus.loc[ix,'P'])
save(focus,'cox_dnmt3b_results.tsv')

for _,z in focus.iterrows():
    log(
        f'DNMT3B_COX={z.endpoint};MODEL={z.model};HR={z.HR:.5g};'
        f'CI={z.CI_low:.5g},{z.CI_high:.5g};'
        f'P={z.P:.5g};FDR={z.FDR_endpoints:.5g}'
    )

plot=focus[
    np.isfinite(focus[['HR','CI_low','CI_high']]).all(axis=1)
    & (focus.CI_low>0)
].copy()

if len(plot):
    fig,ax=plt.subplots(figsize=(8,4.5))
    for i,(_,z) in enumerate(plot.iterrows()):
        color='#d96559' if z.endpoint=='OS' else '#497b9c'
        ax.plot([z.CI_low,z.CI_high],[i,i],color=color,lw=1.8)
        ax.scatter(z.HR,i,color=color,s=35)
    ax.set_yticks(range(len(plot)),plot.endpoint+'/'+plot.model)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.axvline(1,color='gray',ls='--')
    ax.set_xlabel('HR (95% CI)')
    ax.set_title('DNMT3B Cox')
    fig.tight_layout()
    for ext in ['png','pdf']:
        fig.savefig(
            out/f'Fig4_3_Cox.{ext}',dpi=300,bbox_inches='tight'
        )
    plt.close(fig)

(out/'methods_and_captions.txt').write_text(
    'OS: overall survival (primary). PFI: progression-free interval (secondary). CDR times in days.\n'
    'DNMT_SD: DNMT3B per 1 SD, standardized once across all 1097 primary tumors.\n'
    'UniAll: univariable on all endpoint-valid patients. UniCC: same model on complete cases.\n'
    'Multi: complete-case Cox adjusted for age per 10 years and categorical stage (I reference).\n'
    'AgeCC and StageCC: univariable clinical covariate models on complete cases.\n'
    'HR: hazard ratio; CI: confidence interval; plotted intervals are 95%.\n'
    'P: two-sided Wald test; FDR_endpoints: BH across OS and PFI separately within each DNMT3B model.\n'
    'Stage global significance: likelihood-ratio P in StageCC model_status.tsv.\n'
    'PH: rank-time Schoenfeld residual test per coefficient; PH_P<0.05 flags model review.\n'
    'No event imputation, no time-zero replacement, no optimized expression cutoff, no automatic penalization.\n'
)

(out/'software_versions.json').write_text(json.dumps(
    {p:im.version(p) for p in [
        'lifelines','pandas','numpy','scipy','matplotlib'
    ]},indent=2
))
save(pd.DataFrame([
    {
        'source':k,'path':str(p),
        'sha256':hashlib.sha256(p.read_bytes()).hexdigest()
    }
    for k,p in sources.items()
]),'source_manifest.tsv')

review=any(
    s['status']!='OK' or 'PH_error' in s for s in statuses
)
if phrows:
    review=review or any(
        p['PH_P']<.05 for p in phrows if p['model']=='Multi'
    )
log('COX_MODEL_REVIEW_REQUIRED='+str(review).upper())
log('RAW_AND_FROZEN_INPUTS_MODIFIED=NO')
log('FOUR_THREE_COX_COMPLETE')

(out/'4_3_cox_summary.txt').write_text('\n'.join(notes)+'\n')
archive=out/'4_3_cox_return.tar.gz'
with tarfile.open(archive,'w:gz') as t:
    for p in sorted(out.iterdir()):
        if p.is_file() and p!=archive and p.name!='run.log':
            t.add(p,arcname=p.name)
print(f'RETURN_ARCHIVE={archive}',flush=True)
