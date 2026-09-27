"""P10: registered pregnancy-level twin PE panel transport, GSE272342."""
import io, re, tarfile, itertools, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo, stats

GSE='GSE272342'; PDF=Path('results/pe_twin_p10.csv'); SD=Path('results/pe_twin_p10_samples.csv')
_, ann, meta=geo.parse_series_matrix(geo.download_matrices(GSE)[0])
ann=ann[ann.Sample_characteristics_ch1_2.eq('singleton or_twin: twin')].copy()
ann['pregnancy']=ann.Sample_title.str.extract(r', (\d+)[AB]')[0]
assert len(ann)==32 and ann.pregnancy.notna().all() and ann.pregnancy.nunique()==16
ann['label']=ann.Sample_characteristics_ch1_1.map({'clinical group: preeclampsia':'case','clinical group: normotensive':'control'})
assert ann.label.value_counts().to_dict()=={'control':18,'case':14}
assert ann.groupby('pregnancy').label.nunique().eq(1).all()
assert ann.groupby('pregnancy').size().eq(2).all()
used=set(pd.read_csv('results/gsm_to_study.csv').accession)
assert set(ann.index)&used in (set(), set(ann.index)), 'Partial prior inclusion or unexpected accession overlap'; assert all(pd.read_csv('results/gsm_to_study.csv').set_index('accession').loc[gsm,'study_accession']=='GSE272342' for gsm in set(ann.index)&used)
parts={}
with tarfile.open('data/raw/rnaseq/GSE272342/GSE272342_RAW.tar') as tar:
    members={re.match(r'(GSM\d+)_',m.name).group(1):m for m in tar if re.match(r'GSM\d+_',m.name)}
    assert set(ann.index)<=set(members)
    for gsm in ann.index:
        with tar.extractfile(members[gsm]) as fh:
            counts=pd.read_csv(fh,sep='\t',header=None,names=['gene','count'],compression='gzip').set_index('gene')['count']
        assert counts.index.is_unique and counts.notna().all() and (counts>=0).all()
        parts[gsm]=np.log2(1e6*counts/counts.sum()+1)
x=pd.DataFrame(parts);assert not x.isna().any().any()
preg=pd.DataFrame({pid:x[rows.index].mean(axis=1) for pid,rows in ann.groupby('pregnancy')})
lab=ann.groupby('pregnancy').label.first().loc[preg.columns]
assert lab.value_counts().to_dict()=={'control':9,'case':7}
g,v=stats.hedges_g(preg.loc[:,lab=='case'].to_numpy(),preg.loc[:,lab=='control'].to_numpy())
tp=ttest_ind(preg.loc[:,lab=='case'].to_numpy(),preg.loc[:,lab=='control'].to_numpy(),axis=1,equal_var=False).pvalue
out=pd.DataFrame({'gene':preg.index,'g':g,'variance':v,'welch_p':tp}).set_index('gene')
d=pd.read_csv('results/meta_discovery/preeclampsia.csv.gz').set_index('gene')
tags=pd.read_csv('results/split_final.csv').query("disease=='preeclampsia' and split=='validation'").tag
z=pd.concat([pd.read_csv(f'results/series/preeclampsia__{t}.csv.gz').set_index('gene').add_suffix(f'_{i}') for i,t in enumerate(tags)],axis=1)
r=stats.dersimonian_laird(z.filter(like='g_').to_numpy(float),z.filter(like='v_').to_numpy(float))
vdf=pd.DataFrame({k:r[k] for k in ['mu','p','k']},index=z.index)
sel=d.join(vdf,lsuffix='_d',rsuffix='_v').dropna();sel=sel[(sel.q<.05)&(sel.k_d>=8)&(sel.k_v>=3)&(sel.p_v<.01)&(sel.mu_d*sel.mu_v>0)].sort_values(['p_d']).index.tolist()
assert sel==['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1'],sel
out.loc[sel].to_csv(PDF)
ann[['Sample_title','Sample_characteristics_ch1_1','pregnancy','label']].rename_axis('gsm').to_csv(SD)
print('P10 32 GSMs, 16 pregnancies; selected',sel)
print(out.loc[sel].to_string());print('negative signs',int(out.loc[sel].g.lt(0).sum()),'/6')
print('Eligible old screen size:',len(sel),'other eligible down genes:',sum(sel2 not in sel for sel2 in sel))
print('Registered matched-direction 10k six-gene null is degenerate because old screen has exactly six eligible genes; not performed, no empirical p.')

# Descriptive post-outcome uncertainty sensitivity, never a substitute for the registered null.
X=preg.loc[sel].to_numpy(float); labels=lab.eq('case').to_numpy()
obs_t=ttest_ind(X[:,labels],X[:,~labels],axis=1,equal_var=False).statistic
perms=np.asarray(list(itertools.combinations(range(16),7)),dtype=int)
mat=np.zeros((len(perms),16),dtype=bool);mat[np.arange(len(perms))[:,None],perms]=True
n1,n0=7,9
s1=X@mat.T;s0=X.sum(axis=1,keepdims=True)-s1
sq1=(X**2)@mat.T;sq0=(X**2).sum(axis=1,keepdims=True)-sq1
v1=(sq1-s1**2/n1)/(n1-1);v0=(sq0-s0**2/n0)/(n0-1)
t=(s1/n1-s0/n0)/np.sqrt(v1/n1+v0/n0)
maxnull=np.nanmax(abs(t),axis=0)
family_p={gene:float(np.mean(maxnull>=abs(tt))) for gene,tt in zip(sel,obs_t)}
loo={}
for gene in sel:
 signs=[];effects=[]
 for pid in preg.columns:
  rem=preg.loc[gene,preg.columns!=pid]; y=lab.loc[rem.index]
  a=rem[y=='case'].to_numpy()[None,:];b=rem[y=='control'].to_numpy()[None,:]
  gg,_=stats.hedges_g(a,b);signs.append(bool(gg[0]<0));effects.append(float(gg[0]))
 loo[gene]={'negative_loocv':int(sum(signs)),'total':16,'g_min':float(min(effects)),'g_max':float(max(effects))}
res={'cohort':GSE,'unit':'pregnancy','case':7,'control':9,'observed_welch_t':dict(zip(sel,map(float,obs_t))),
     'posthoc_exact_familywise_abs_welch_t_p':family_p,'n_label_permutations':len(perms),'loo':loo,
     'caveat':'Post-outcome descriptive sensitivity, no gestational-age or twin chorionicity adjustment; not a registered test or discovery.'}
Path('results/pe_twin_p10_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
print('Exact pregnancy-label permutation max-|Welch t| familywise p:',family_p)
print('Leave-one-pregnancy-out negative counts:',{k:v['negative_loocv'] for k,v in loo.items()})

# Secondary post-outcome gestational-age sensitivity. Delivery week.days is not decimal weeks.
from scipy.stats import norm
meta=ann.groupby('pregnancy').first().loc[preg.columns]
def week_days(v):
    raw=v.split(': ',1)[1]; pieces=raw.split('.')
    assert len(pieces) in (1,2) and (len(pieces)==1 or (len(pieces[1])==1 and int(pieces[1])<=6))
    return int(pieces[0])+(int(pieces[1])/7 if len(pieces)>1 else 0)
weeks=meta.Sample_characteristics_ch1_3.map(week_days).to_numpy(float)
y=lab.eq('case').to_numpy(float)
Xreg=np.column_stack([np.ones(16),y,weeks-weeks.mean()]);
beta=np.linalg.lstsq(Xreg,preg.loc[sel].T.to_numpy(float),rcond=None)[0]
err=preg.loc[sel].T.to_numpy(float)-Xreg@beta
se=np.sqrt(np.diag(np.linalg.inv(Xreg.T@Xreg))[1] * (err**2).sum(axis=0)/(16-3))
from scipy.stats import t as tdist
p=2*tdist.sf(abs(beta[1]/se),df=13)
res['posthoc_delivery_age_adjusted']={gene:{'case_coefficient_log2cpm':float(b),'t_p':float(pp)} for gene,b,pp in zip(sel,beta[1],p)}
res['case_mean_delivery_age_weeks']=float(weeks[labels].mean())
res['control_mean_delivery_age_weeks']=float(weeks[~labels].mean())
res['delivery_age_weeks']={str(k):float(v) for k,v in zip(preg.columns,weeks)}
res['gestational_caveat']='Delivery age is post-outcome and potentially affected by disease; adjustment is descriptive only, may induce bias, and is not matched first-trimester replication.'
Path('results/pe_twin_p10_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n')
print('Post-outcome delivery-age adjusted associations:',res['posthoc_delivery_age_adjusted'])
