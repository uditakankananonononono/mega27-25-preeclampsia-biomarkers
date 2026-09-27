"""Post-outcome exploratory P21 sensitivity: A/B normal samples and LOO cases."""
import pathlib,re,json,hashlib,requests
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,mannwhitneyu
import sys
sys.path.insert(0,'src');from ubiomark import stats
M=pd.read_csv('results/pe_p21_sample_map.csv');H=pd.read_csv('data/raw/hgnc/hgnc_complete_set.txt',sep='\t',dtype=str)
p=pathlib.Path('data/raw/rnaseq/GSE148241')
# Include the two A/B control candidate repeats in a separate descriptive sensitivity.
rows=M[M.tissue=='Placenta'];all_x={}
for r in rows.itertuples(index=False):
 f=p/(r.gsm+'_'+r.title+'.raw.counts.txt.gz')
 if not f.exists():
  response=requests.get(r.url,timeout=30);response.raise_for_status();f.write_bytes(response.content)
 x=pd.read_csv(f,sep='\t',compression='gzip');assert list(x.columns)==['ENSEMBL_ID','raw.counts'] and x.ENSEMBL_ID.is_unique and (x['raw.counts']>=0).all()
 all_x[r.gsm]=x.set_index('ENSEMBL_ID')['raw.counts'];
X=pd.concat(all_x,axis=1);assert not X.isna().any().any();L=X.sum(axis=0)
labels=rows.set_index('gsm').status
cases=labels[labels!='Normal'].index.tolist();normals=labels[labels=='Normal'].index.tolist();primary=M[(M.use)&(M.status=='Normal')].gsm.tolist()
assert (len(cases),len(normals),len(primary))==(9,32,30)
Y=np.log2(X.div(L,axis=1)*1e6+1)
records=[]
for gene in ['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']:
 eid=H.loc[(H.symbol==gene)&(H.status=='Approved'),'ensembl_gene_id'].iloc[0];assert eid in Y.index
 for excluded_case in ['none']+cases:
  aa=[g for g in cases if g!=excluded_case]
  for control_group,bb in [('A_only',primary),('include_B',normals)] if excluded_case=='none' else [('A_only',primary)]:
   a=Y.loc[eid,aa].to_numpy();b=Y.loc[eid,bb].to_numpy()
   g,var=stats.hedges_g(a.reshape(1,-1),b.reshape(1,-1));pr=ttest_ind(a,b,equal_var=False).pvalue
   records.append(dict(gene=gene,excluded_case=excluded_case,control_group=control_group,g=float(g[0]),welch_p=float(pr),n_case=len(a),n_control=len(b)))
R=pd.DataFrame(records);R.to_csv('results/pe_p21_postoutcome_sensitivity.csv',index=False)
lib=pd.DataFrame(dict(gsm=L.index,library_total=L.to_numpy(),status=labels.loc[L.index].to_numpy()));lib.to_csv('results/pe_p21_library_totals.csv',index=False)
meta=dict(library_total_case_min=int(L[cases].min()),library_total_case_max=int(L[cases].max()),library_total_control_min=int(L[normals].min()),library_total_control_max=int(L[normals].max()),library_total_mannwhitney_p=float(mannwhitneyu(L[cases],L[normals]).pvalue),loe_case_9=list(cases))
pathlib.Path('results/pe_p21_postoutcome_sensitivity.json').write_text(json.dumps(meta,indent=2)+'\n')
print(R.groupby('gene').agg(min_g=('g','min'),max_g=('g','max'),max_p=('welch_p','max')));print(R[R.control_group=='include_B'][['gene','g','welch_p']]);print(meta)
