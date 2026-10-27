import os
from pathlib import Path

source=Path(os.environ['SOURCE_RUN'])/'hub.py'
out=Path(os.environ['RUN'])
s=source.read_text()
start=s.index('for b in range(100):\n')
end=s.index('bt=pd.DataFrame(boot)',start)
assert s.count('for b in range(100):\n')==1
assert 'cb=corr(a[:,idx])' in s[start:end]

block='''constant_audit=[]
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
    out/'subsample_constant_genes.tsv',sep='\\t',index=False)
'''
s=s[:start]+block+s[end:]
s=s.replace(
    "bt=pd.DataFrame(boot)",
    "bt=pd.DataFrame(boot)\nassert bt.groupby('network').size().eq(100).all()",
    1)
(out/'hub.py').write_text(s)
(out/'repair_policy.txt').write_text(
    'Full-cohort universe unchanged. Per subsample, remove constant genes before ranking.\n'
    'Percentile denominator uses variable nodes in that subsample.\n'
    'DNMT3B must remain variable; every network must have 100 valid iterations.\n'
    'No zero substitution for undefined correlations. Seed and thresholds unchanged.\n')
print('SOURCE_SCRIPT='+str(source))
print('GENERATED_SCRIPT='+str(out/'hub.py'))
