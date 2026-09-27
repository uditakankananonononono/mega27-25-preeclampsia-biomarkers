"""Post-outcome descriptive SGA sensitivity of frozen PE genes in GSE186257.

Does not redefine P13, does not provide a new validation cohort or gene claim.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
G='GSE186257'; genes=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
_,a,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);x=pd.read_csv('data/raw/rnaseq/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz',sep='\t').set_index('Unnamed: 0')
assert len(a)==44 and set(a.Sample_title)==set(x.columns)
x.index=x.index.str.split('$').str[-1]
assert all(sum(x.index==g)==1 for g in genes)
a['group']=a.Sample_characteristics_ch1_0.map({'disease: Severe preeclampsia':'PE','disease: Non-hypertensive Control':'control'})
a['sga']=a.Sample_characteristics_ch1_2.map({'sga: yes':'yes','sga: no':'no'})
a['sex']=a.Sample_characteristics_ch1_1.map({'gender: Male':'male','gender: Female':'female'})
assert a[['group','sga','sex']].notna().all().all()
assert pd.crosstab(a.group,a.sga).to_dict()=={'no':{'PE':1,'control':10},'yes':{'PE':25,'control':8}}
xx=np.log2(x+1); out=[]
for setting,ids in [('all',a.index),('SGA_only',a.index[a.sga.eq('yes')]),('female_only',a.index[a.sex.eq('female')]),('male_only',a.index[a.sex.eq('male')])]:
 b=a.loc[ids];case=b.index[b.group.eq('PE')];control=b.index[b.group.eq('control')]
 # Subgroups with <2 per arm cannot estimate variance and are explicitly excluded.
 if min(len(case),len(control))<2:continue
 for gene in genes:
  ca=xx.loc[gene,b.loc[case,'Sample_title']].to_numpy()[None,:]
  co=xx.loc[gene,b.loc[control,'Sample_title']].to_numpy()[None,:]
  gg,v=stats.hedges_g(ca,co)
  out.append(dict(setting=setting,gene=gene,n_PE=len(case),n_control=len(control),hedges_g=float(gg[0]),welch_p=float(ttest_ind(ca[0],co[0],equal_var=False).pvalue)))
d=pd.DataFrame(out);d.to_csv('results/pe_p13_sga_sensitivity.csv',index=False)
Path('results/pe_p13_sga_sensitivity.json').write_text(json.dumps({'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186257','SGA_table':pd.crosstab(a.group,a.sga).to_dict(),'sex_table':pd.crosstab(a.group,a.sex).to_dict(),'note':'Post-outcome descriptive within-study strata, not P13 registered primary, not independent replication. Only one PE non-SGA so non-SGA PE contrast not estimable. SGA is possibly downstream of PE; conditioning on it does not causally adjust.'},indent=2)+'\n')
print(pd.crosstab(a.group,a.sga));print(pd.crosstab(a.group,a.sex));print(d[d.gene.eq('FES')].to_string(index=False));print(d[d.setting.eq('SGA_only')].to_string(index=False))
