"""P21 prespecified at-delivery EOSPE placenta diagnostic; not prediction or discovery."""
import concurrent.futures,hashlib,json,pathlib,re,requests,sys
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
sys.path.insert(0,'src');from ubiomark import stats
BASE='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
series=requests.get(BASE+'GSE148241&targ=self&form=text&view=full',timeout=20).text
ids=re.findall(r'^!Series_sample_id = (GSM\d+)',series,re.M);assert len(ids)==len(set(ids))==43

def metadata(g):
 t=requests.get(BASE+g+'&targ=self&form=text&view=full',timeout=15).text
 f=lambda k:re.findall(r'^!Sample_'+k+r' = (.*)$',t,re.M)
 title=f('title')[0].strip(); ch=[x.strip() for x in f('characteristics_ch1')]
 status=next(x.split(': ',1)[1] for x in ch if x.startswith('subject status: '));tissue=next(x.split(': ',1)[1] for x in ch if x.startswith('tissue: '))
 files=[x.strip().replace('ftp://','https://') for x in re.findall(r'^!Sample_supplementary_file_\d+ = (.*)$',t,re.M)]
 path=[x for x in files if x.endswith('.raw.counts.txt.gz')];assert len(path)==1,(g,path)
 return dict(gsm=g,title=title,status=status,tissue=tissue,url=path[0])
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool: allrows=list(pool.map(metadata,ids))
M=pd.DataFrame(allrows);assert M.gsm.is_unique and M.title.is_unique
assert M.tissue.value_counts().to_dict()=={'Placenta':41,'Umbilical Cord Blood':2}
assert M.status.value_counts().to_dict()=={'Normal':34,'early-onset severe preeclampsia (EOSPE)':9}
old=pd.read_csv('results/dataset_manifest.csv');assert 'GSE148241' not in set(old.accession) and not set(M.gsm)&set(old.accession)
M['use']=(M.tissue=='Placenta')&(~M.title.isin(['PHNPR27B','PHNPR38B']))
assert M.loc[M.use,'status'].value_counts().to_dict()=={'Normal':30,'early-onset severe preeclampsia (EOSPE)':9}
M.to_csv('results/pe_p21_sample_map.csv',index=False)
p=pathlib.Path('data/raw/rnaseq/GSE148241');p.mkdir(parents=True,exist_ok=True)

def download(row):
 fname=row.gsm+'_'+row.title+'.raw.counts.txt.gz'; local=p/fname
 if not local.exists():
  r=requests.get(row.url,timeout=30);r.raise_for_status();local.write_bytes(r.content)
 v=pd.read_csv(local,sep='\t',compression='gzip')
 assert list(v.columns)==['ENSEMBL_ID','raw.counts'] and v.ENSEMBL_ID.is_unique and len(v)>10000
 assert pd.api.types.is_integer_dtype(v['raw.counts']) and (v['raw.counts']>=0).all()
 return row.gsm,v.set_index('ENSEMBL_ID')['raw.counts'],hashlib.sha256(local.read_bytes()).hexdigest()
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool: samples=list(pool.map(download,M[M.use].itertuples(index=False)))
sets=[set(v.index) for _,v,_ in samples];assert all(s==sets[0] for s in sets),[(len(s),len(s^sets[0])) for s in sets]
X=pd.concat({g:v for g,v,_ in samples},axis=1);X.index.name='ENSEMBL_ID';assert (X.sum(axis=0)>0).all()
(pd.DataFrame([dict(gsm=g,sha256=h) for g,_,h in samples])).to_csv('results/pe_p21_raw_checksums.csv',index=False)
H=pd.read_csv('data/raw/hgnc/hgnc_complete_set.txt',sep='\t',dtype=str)
genes=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']; R=[]
case=M.loc[M.use&(M.status!='Normal'),'gsm'].tolist();control=M.loc[M.use&(M.status=='Normal'),'gsm'].tolist()
Y=np.log2(X.div(X.sum(axis=0),axis=1)*1e6+1)
for gene in genes:
 hit=H[(H.symbol==gene)&(H.status=='Approved')]
 if len(hit)!=1 or pd.isna(hit.ensembl_gene_id.iloc[0]) or hit.ensembl_gene_id.iloc[0] not in X.index:
  R.append(dict(gene=gene,measured=False,reason='absent_or_ambiguous_mapping'));continue
 eid=hit.ensembl_gene_id.iloc[0]; a=Y.loc[eid,case].to_numpy();b=Y.loc[eid,control].to_numpy()
 g,var=stats.hedges_g(a.reshape(1,-1),b.reshape(1,-1)); pval=ttest_ind(a,b,equal_var=False).pvalue
 R.append(dict(gene=gene,ensembl_id=eid,measured=True,g=float(g[0]),v=float(var[0]),welch_p=float(pval),down=bool(g[0]<0),n_case=9,n_control=30,mean_case_log2=float(np.mean(a)),mean_control_log2=float(np.mean(b))))
R=pd.DataFrame(R);R.to_csv('results/pe_p21_panel.csv',index=False)
measured=int(R.measured.sum());down=int(R.down.fillna(False).sum());hits=R.loc[(R.down==True)&(R.welch_p<.05/6),'gene'].tolist()
out=dict(gse='GSE148241',source=BASE+'GSE148241',n_case=9,n_control=30,genes_measured=measured,genes_down=down,same_direction_bonferroni=hits,descriptive_panel_support=bool(measured==6 and down==6 and hits),n_raw_files=39,exclusions=['GSM4458451/52 cord blood','PHNPR27B and PHNPR38B possible same-woman B samples'],interpretation='2020 at-delivery early-onset severe PE versus normal placenta, postselected panel, not early pregnancy prediction or independently nominated discovery')
pathlib.Path('results/pe_p21_result.json').write_text(json.dumps(out,indent=2)+'\n');print(R.to_string(index=False));print(json.dumps(out,indent=2))
