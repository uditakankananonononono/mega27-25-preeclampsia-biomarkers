"""P3 PE ST3GAL2 pooled-effect estimator sensitivity; descriptive after registered failure."""
import json
import numpy as np,pandas as pd
from scipy.stats import norm
from statsmodels.stats.meta_analysis import combine_effects
x=pd.read_csv('results/fresh_cohort_tests.csv')
x=x[(x.disease=='preeclampsia') & x.gse.isin(['GSE204835','GSE296973','GSE303840','GSE306864'])].sort_values('gse')
assert len(x)==4 and x.g.notna().all()
rows=[]
for name,z in [('all_four',x),('placenta_and_villus_only',x[x.gse!='GSE296973'])]:
 for method in ['dl','iterated']:
  fit=combine_effects(z.g.to_numpy(),z.se.to_numpy()**2,method_re=method,row_names=list(z.gse))
  mu=float(fit.mean_effect_re);se=float(fit.sd_eff_w_re)
  rows.append(dict(stratum=name,method=method,n_cohorts=len(z),cohorts=';'.join(z.gse),mu=mu,se=se,tau2=float(fit.tau2),normal_one_sided_p_down=float(norm.cdf(mu/se)),normal_CI95=[mu-1.96*se,mu+1.96*se],leave_one_out=False))
 if name=='all_four':
  for omit in x.gse:
   y=x[x.gse!=omit]
   fit=combine_effects(y.g.to_numpy(),y.se.to_numpy()**2,method_re='iterated')
   mu=float(fit.mean_effect_re);se=float(fit.sd_eff_w_re)
   rows.append(dict(stratum='leave_one_out',method='iterated',n_cohorts=3,omitted=omit,cohorts=';'.join(y.gse),mu=mu,se=se,tau2=float(fit.tau2),normal_one_sided_p_down=float(norm.cdf(mu/se)),normal_CI95=[mu-1.96*se,mu+1.96*se],leave_one_out=True))
pd.DataFrame(rows).to_csv('results/pe_p3_estimator_sensitivity.csv',index=False)
print(pd.DataFrame(rows)[['stratum','method','omitted','mu','se','normal_one_sided_p_down']].to_string(index=False))
