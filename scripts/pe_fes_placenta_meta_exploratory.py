"""Post-selection, post-outcome FES placental descriptive meta-analysis.

This CANNOT establish a novel biomarker. Contrasts are biologically nonexchangeable
and were all searched after discovery/validation selected the gene.
"""
import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm,t
from statsmodels.stats.meta_analysis import combine_effects
sources=[('GSE272342','twin pregnancy mean, 7 PE/9 controls','results/pe_twin_p10.csv'),('GSE190971','term/preterm placenta, 7 PE/6 controls','results/pe_plac_p12.csv'),('GSE186257','severe PE/healthy placenta, 26/18','results/pe_plac_p13_GSE186257.csv'),('GSE114691','PE-only placenta, 20/21; incomplete historical ENST map','results/pe_plac_p13_GSE114691.csv')]
rows=[]
for g,setting,file in sources:
 x=pd.read_csv(file).set_index('gene').loc['FES'];assert np.isfinite(x.g) and np.isfinite(x.variance) and x.variance>0
 rows.append(dict(gse=g,setting=setting,g=float(x.g),variance=float(x.variance),se=float(np.sqrt(x.variance))))
d=pd.DataFrame(rows);d.to_csv('results/pe_fes_placenta_meta_cohorts.csv',index=False)
def meta(x):
 fit=combine_effects(x.g.to_numpy(),x.variance.to_numpy(),method_re='iterated',use_t=False)
 mu=float(fit.mean_effect_re);se=float(fit.sd_eff_w_re)
 tau=float(np.sqrt(max(0,fit.tau2)))
 out={'k':len(x),'mu':mu,'se':se,'tau2':float(fit.tau2),'Q':float(fit.q),'I2_raw':float(fit.i2),'normal_reference_two_sided_p':float(2*norm.sf(abs(mu/se))),'normal_95_CI':[mu-1.96*se,mu+1.96*se]}
 if len(x)>2:
  # Descriptive Hartung-Knapp t reference with k-1 df, not a preregistered test.
  hkse=float(fit.sd_eff_w_re_hksj);out['hksj_se']=hkse;out['hksj_t_p']=float(2*t.sf(abs(mu/hkse),df=len(x)-1));out['hksj_95_CI']=[mu-t.ppf(.975,len(x)-1)*hkse,mu+t.ppf(.975,len(x)-1)*hkse]
  out['normal_prediction_interval']=[mu-1.96*np.sqrt(tau*tau+se*se),mu+1.96*np.sqrt(tau*tau+se*se)]
 return out
res={'source':'four previously inspected independent GSE placental contrasts, post-selected FES','cohorts':rows,'all':meta(d),'leave_one_study_out':{g:meta(d[d.gse.ne(g)]) for g in d.gse},'non_twin':meta(d[d.gse.ne('GSE272342')]),'caveat':'Post-selected FES and cohorts are strongly confounded by gestational age, IUGR/SGA, tissue composition, platform and subtype; the pooled estimate is descriptive and cannot identify PE-specific effect or biomarker.'}
Path('results/pe_fes_placenta_meta_exploratory.json').write_text(json.dumps(res,indent=2)+'\n')
print(d.to_string(index=False));print('all',res['all']);print('non-twin',res['non_twin']);print('LOO',[(k,round(v['mu'],3),round(v['hksj_t_p'],4)) for k,v in res['leave_one_study_out'].items()])
