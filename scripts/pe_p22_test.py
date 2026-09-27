"""P22 fixed PE six-gene panel in preterm-birth control placenta, descriptive only."""
import concurrent.futures,hashlib,json,pathlib,re,requests,sys
import pandas as pd,numpy as np
from scipy.stats import ttest_ind
sys.path.insert(0,'src');from ubiomark import stats
URL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
s=requests.get(URL+'GSE177049&targ=self&form=text&view=full',timeout=20).text
ids=re.findall(r'^!Series_sample_id = (GSM\d+)',s,re.M);assert len(ids)==len(set(ids))==20

def get(g):
 t=requests.get(URL+g+'&targ=self&form=text&view=full',timeout=15).text
 f=lambda k:re.findall(r'^!Sample_'+k+r' = (.*)$',t,re.M)
 title=f('title')[0].strip();c=[x.strip() for x in f('characteristics_ch1')]
 disease=next(x.split(': ',1)[1] for x in c if x.startswith('disease: '))
 return dict(gsm=g,title=title,disease=disease,assay='RNA' if title.endswith('- RNA Sequencing') else 'miRNA',token=title.split(' - ')[0])
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:rows=list(ex.map(get,ids))
M=pd.DataFrame(rows);assert M.gsm.is_unique and len(M[M.assay=='RNA'])==10 and M[M.assay=='RNA'].token.is_unique
R=M[M.assay=='RNA'].copy();assert R.disease.value_counts().to_dict()=={'early-onset preeclampsia':5,'preterm birth':5}
assert all(x.startswith('PE') if y=='early-onset preeclampsia' else x.startswith('N') for x,y in zip(R.token,R.disease))
old=pd.read_csv('results/dataset_manifest.csv');assert 'GSE177049' not in set(old.accession) and not set(R.gsm)&set(old.accession)
M.to_csv('results/pe_p22_sample_map.csv',index=False)
p=pathlib.Path('data/raw/rnaseq/GSE177049/GSE177049_gene-miRNA-RNA_expression_data.xlsx');p.parent.mkdir(parents=True,exist_ok=True)
if not p.exists():
 import urllib.request
 urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE177nnn/GSE177049/suppl/GSE177049_gene-miRNA-RNA_expression_data.xlsx',p)
sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='c882ad47d8f64c9c84cf3b643418518ee0f4839087c987ecd3205096557847d1'
x=pd.read_excel(p,sheet_name='gene expression',header=16)
assert set(R.token)==set(x.columns[11:21]) and len(x.columns)==21,(x.columns.tolist())
assert {'Gene_Name','PE_FPKM','N_FPKM'}.issubset(x.columns)
values=x.loc[:,R.token].apply(pd.to_numeric,errors='raise')
assert np.isfinite(values.to_numpy()).all() and (values.to_numpy()>=0).all()
x.Gene_Name=x.Gene_Name.astype(str).str.strip();dup=set(x.loc[x.Gene_Name.duplicated(keep=False),'Gene_Name'])
genes=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1'];results=[]
for gene in genes:
 z=x[x.Gene_Name==gene]
 if len(z)!=1 or gene in dup:
  results.append(dict(gene=gene,measured=False,reason='absent_or_ambiguous_symbol'));continue
 v=z[R.token].iloc[0].astype(float)
 a=np.log2(v[R[R.disease=='early-onset preeclampsia'].token].to_numpy()+1);b=np.log2(v[R[R.disease=='preterm birth'].token].to_numpy()+1)
 g,var=stats.hedges_g(a.reshape(1,-1),b.reshape(1,-1));pval=ttest_ind(a,b,equal_var=False).pvalue
 results.append(dict(gene=gene,measured=True,g=float(g[0]),v=float(var[0]),welch_p=float(pval),down=bool(g[0]<0),n_case=5,n_control=5,mean_case_log2=float(np.mean(a)),mean_control_log2=float(np.mean(b))))
T=pd.DataFrame(results);T.to_csv('results/pe_p22_panel.csv',index=False)
measured=int(T.measured.sum());down=int(T.down.fillna(False).sum());hits=T.loc[(T.down==True)&(T.welch_p<.05/6),'gene'].tolist()
out=dict(gse='GSE177049',source=URL+'GSE177049',n_case=5,n_control=5,genes_measured=measured,genes_down=down,same_direction_bonferroni=hits,descriptive_panel_support=bool(measured==6 and down==6 and hits),workbook_sha256=sha,interpretation='Published small 5/5 at-delivery placenta preterm-birth controls, not early prospective prediction, new discovery or specificity')
pathlib.Path('results/pe_p22_result.json').write_text(json.dumps(out,indent=2)+'\n');print(T.to_string(index=False));print(json.dumps(out,indent=2))
