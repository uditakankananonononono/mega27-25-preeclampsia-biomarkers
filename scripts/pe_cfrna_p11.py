"""P11 registered exploratory first-trimester cfRNA panel, GSE192902."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
G='GSE192902'; PANEL=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
_,a,_=geo.parse_series_matrix(geo.download_matrices(G)[0])
assert len(a)==404 and a.Sample_title.is_unique and a.index.is_unique
assert a.assign(mother=a.Sample_title.str.extract(r'^(\d+)_')[0]).groupby('mother').Sample_characteristics_ch1_2.nunique().eq(1).all()
out=[];samples=[];notes={}
for key,label in [('Discovery','Discovery'),('Validation1','Validation 1'),('Validation2','Validation 2')]:
 full=pd.read_csv(f'data/raw/rnaseq/{G}/{G}_counts_{key}_postQC.csv.gz')
 aa=a[a.Sample_characteristics_ch1_2.eq('cohort: '+label)].copy()
 assert len(full.columns)-2==len(aa) and set(full.columns[2:])==set(aa.Sample_title)
 assert full.gene_num.is_unique and full.gene_name.notna().all()
 assert full.iloc[:,2:].notna().all().all() and (full.iloc[:,2:]>=0).all().all()
 assert all(z not in full.gene_name.value_counts().loc[lambda x:x>1] for z in PANEL)
 aa=aa[aa.Sample_characteristics_ch1_1.eq('sampling time group: ≤12 weeks gestation')].copy()
 aa['mother']=aa.Sample_title.str.extract(r'^(\d+)_')[0].values
 aa['sample_order']=aa.Sample_title.str.extract(r'^\d+_(\d+)')[0].astype(int).values
 aa=aa.sort_values(['mother','sample_order','Sample_title']).drop_duplicates('mother',keep='first')
 assert aa.mother.is_unique and len(aa)>10
 aa['label']=aa.Sample_characteristics_ch1_0.map({'disease: control':'control','disease: pre-eclampsia':'case','disease: severe pre-eclampsia':'case'})
 assert aa.label.notna().all()
 names=aa.Sample_title.tolist();x=full.set_index('gene_name')[names];x=np.log2(1e6*x/x.sum(axis=0)+1)
 names_case=aa.loc[aa.label.eq('case'),'Sample_title'];names_control=aa.loc[aa.label.eq('control'),'Sample_title']
 assert len(names_case)>=3 and len(names_control)>=9
 for gene in PANEL:
  if gene not in x.index:
   out.append({'cohort':key,'gene':gene,'n_case':len(names_case),'n_control':len(names_control),'g':np.nan,'welch_p':np.nan,'measured':False});continue
  one=x.loc[gene,names_case].to_numpy(dtype=float)[None,:];zero=x.loc[gene,names_control].to_numpy(dtype=float)[None,:]
  g,v=stats.hedges_g(one,zero);p=ttest_ind(one,zero,axis=1,equal_var=False).pvalue[0]
  out.append({'cohort':key,'gene':gene,'n_case':len(names_case),'n_control':len(names_control),'g':g[0],'welch_p':p,'measured':True})
 aa[['Sample_title','mother','label','Sample_characteristics_ch1_0']].rename_axis('gsm').assign(cohort=key).to_csv(f'results/pe_cfrna_p11_{key}_samples.csv')
 samples.extend(aa.index.tolist());notes[key]={'n_all_postqc':len(full.columns)-2,'n_geo':len(a[a.Sample_characteristics_ch1_2.eq('cohort: '+label)]),'n_early_unique_mothers':len(aa),'n_case':len(names_case),'n_control':len(names_control)}
r=pd.DataFrame(out);r.to_csv('results/pe_cfrna_p11.csv',index=False)
assert len(samples)==len(set(samples)) and not set(samples)&set(pd.read_csv('results/gsm_to_study.csv').accession)
notes['all_selected_gsm']=samples
Path('results/pe_cfrna_p11.json').write_text(json.dumps(notes,indent=2)+'\n')
print(notes.keys(),{k:v for k,v in notes.items() if k!='all_selected_gsm'});print(r.to_string(index=False));print('negative signs',r.groupby('cohort').g.apply(lambda z:int(z.lt(0).sum())).to_dict())
