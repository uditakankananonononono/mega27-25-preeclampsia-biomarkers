"""Exploratory P14, selected FES in reused HDP source; NOT independent validation."""
import csv, gzip, io, json, hashlib
from pathlib import Path
import requests, numpy as np, pandas as pd
from scipy.stats import ttest_ind
from ubiomark.stats import hedges_g
BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303463/'
COUNT=Path('data/geo/GSE303463/raw_counts.csv.gz')
assert hashlib.sha256(COUNT.read_bytes()).hexdigest()=='d9ba5a9d683236f35b87286a700e6506f1039b121fde279a8dc460b6f36e49c5'
rows=[]
for platform in ['GPL24676','GPL34281']:
 url=BASE+f'matrix/GSE303463-{platform}_series_matrix.txt.gz'
 response=requests.get(url,timeout=30);response.raise_for_status()
 s=gzip.decompress(response.content).decode().split('!series_matrix_table_begin')[0]
 lines={}
 for line in s.splitlines():
  if line.startswith('!Sample_'):
   fields=next(csv.reader(io.StringIO(line),delimiter='\t'))
   lines.setdefault(fields[0],[]).append(fields[1:])
 titles=lines['!Sample_title'][0]; n=len(titles)
 assert all(len(x)==n for value in lines.values() for x in value)
 characteristics=lines['!Sample_characteristics_ch1']
 fields={r[0].split(':',1)[0].lower():r for r in characteristics}
 for i,title in enumerate(titles):
  d={'title':title,'gsm':lines['!Sample_geo_accession'][0][i],'platform':platform}
  for key,v in fields.items():d[key]=v[i].split(':',1)[1].strip()
  d['biosample']=next((v[i].split('/')[-1] for v in lines.get('!Sample_relation',[]) if 'BioSample:' in v[i]),'')
  rows.append(d)
meta=pd.DataFrame(rows).set_index('title')
mat=pd.read_csv(COUNT,index_col=0)
assert mat.columns.is_unique and mat.index.is_unique and mat.shape==(20481,168)
assert len(meta)==168 and meta.index.is_unique and set(meta.index)==set(mat.columns)
assert meta.gsm.is_unique and meta.biosample.is_unique and meta.biosample.ne('').all()
meta=meta.loc[mat.columns];meta['gestation_weeks']=pd.to_numeric(meta['gestational age_(weeks)'],errors='coerce')
assert meta['disease'].isin(['Normal','PE','siPE','cHTN','gHTN']).all()
assert set(meta['placental pathology'])=={'nonMVM','MVM'}
assert (mat.to_numpy()>=0).all() and (mat.to_numpy()%1==0).all()
lib=mat.sum(axis=0)
assert (lib>0).all()
expr=np.log2(1+mat.loc['ENSG00000182511$FES'].to_numpy()/lib.to_numpy()*1e6)
meta['fes_log2_cpm1']=expr
meta.to_csv('results/pe_p14_sample_map.csv',index_label='title')
primary=meta[meta.disease.isin(['PE','Normal'])].copy()
def contrast(q):
 a=q[q.disease.eq('PE')].fes_log2_cpm1.to_numpy();b=q[q.disease.eq('Normal')].fes_log2_cpm1.to_numpy()
 out={'n_pe':len(a),'n_normal':len(b)}
 if min(len(a),len(b))<3:return {**out,'status':'insufficient groups'}
 out.update(mean_pe=float(np.mean(a)),mean_normal=float(np.mean(b)),g=float(hedges_g(a.reshape(1,-1),b.reshape(1,-1))[0].item()),welch_p=float(ttest_ind(a,b,equal_var=False).pvalue));return out
result={'scope':'selected FES exploratory reuse-risk/confounding diagnostic, not independent cohort','source':BASE,'count_sha256':hashlib.sha256(COUNT.read_bytes()).hexdigest(),'metadata_url':[BASE+f'matrix/GSE303463-{p}_series_matrix.txt.gz' for p in ['GPL24676','GPL34281']],'counts':meta.groupby(['platform','disease','placental pathology']).size().rename('n').reset_index().to_dict('records'),'all_pe_normal':contrast(primary),'platforms':{p:contrast(q) for p,q in primary.groupby('platform')},'gestational_34plus':contrast(primary[primary.gestation_weeks.ge(34)]),'same_platform_age_34to39':contrast(primary[(primary.platform=='GPL34281') & primary.gestation_weeks.between(34,39)]),'same_platform_nonMVM':contrast(primary[(primary.platform=='GPL34281') & primary['placental pathology'].eq('nonMVM')]),'same_platform_MVM':contrast(primary[(primary.platform=='GPL34281') & primary['placental pathology'].eq('MVM')]),'complete_case_age_by_disease':primary.groupby('disease').gestation_weeks.agg(['min','max','median','count']).to_dict('index'),'complete_case_mvm_by_disease':primary.groupby(['disease','placental pathology']).size().rename('n').reset_index().to_dict('records'),'person_source_mapping':'UNRESOLVED: parent GSE explicitly reuses GSE186257, GSE234729 and PRJNA1027377; no independent validation claim or accession tally.'}
for platform,q in primary.groupby('platform'):
 z=q.dropna(subset=['gestation_weeks'])
 X=pd.DataFrame({'const':1.,'PE':z.disease.eq('PE').astype(float),'gestational_week':z.gestation_weeks.astype(float),'MVM':z['placental pathology'].eq('MVM').astype(float)})
 if len(z)>X.shape[1] and np.linalg.matrix_rank(X.to_numpy())==X.shape[1]:
  import statsmodels.api as sm
  fit=sm.OLS(z.fes_log2_cpm1.to_numpy(),X).fit();result.setdefault('platform_adjustment',{})[platform]={'n':len(z),'pe_beta':float(fit.params['PE']),'pe_p':float(fit.pvalues['PE']),'pe_95_ci':[float(v) for v in fit.conf_int().loc['PE']],'note':'post-diagnosis gestation/pathology may mediate and is not causal control'}
 else:result.setdefault('platform_adjustment',{})[platform]={'n':len(z),'status':'insufficient full-rank design'}
Path('results/pe_p14_hdp_reuse_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('counts','metadata_url')},indent=2))
