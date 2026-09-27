"""GSE143953 late-placenta transport of six frozen preeclampsia DOWN genes."""
import re,json,hashlib, pathlib
import pandas as pd,numpy as np,requests
from scipy.stats import ttest_ind
import sys
sys.path.insert(0,'src')
from ubiomark import stats
url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
s=requests.get(url+'GSE143953&targ=self&form=text&view=full',timeout=20).text
ids=re.findall(r'^!Series_sample_id = (GSM\d+)',s,re.M);assert len(ids)==len(set(ids))==8
rows=[]
for gsm in ids:
 t=requests.get(url+gsm+'&targ=self&form=text&view=full',timeout=20).text
 def field(k):
  m=re.search(r'^!Sample_'+k+r' = (.*)$',t,re.M);assert m,(gsm,k);return m.group(1).strip()
 title=field('title');ds=[v for v in re.findall(r'^!Sample_characteristics_ch1 = disease: (.*)$',t,re.M)];assert len(ds)==1
 disease=ds[0].strip();assert (title.startswith('ECL') and disease=='preeclampsia') or (title.startswith('NEG') and disease=='control')
 rows.append(dict(gsm=gsm,title=title,disease=disease))
M=pd.DataFrame(rows);assert M.title.is_unique and M.disease.value_counts().to_dict()=={'preeclampsia':4,'control':4}
old=pd.read_csv('results/dataset_manifest.csv');assert not set(M.gsm)&set(old.accession) and 'GSE143953' not in set(old.accession)
M.to_csv('results/pe_p19_sample_map.csv',index=False)
path=pathlib.Path('data/raw/rnaseq/GSE143953/GSE143953_Expression_Gene.xlsx');path.parent.mkdir(parents=True,exist_ok=True)
if not path.exists():
 import urllib.request
 urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143953/suppl/GSE143953_Expression_Gene.xlsx',path)
assert hashlib.sha256(path.read_bytes()).hexdigest()=='42e900944734b1b461e86ff28c94c0767fae95fb74d97c52f534ab3757fa1136'
x=pd.read_excel(path,sheet_name='gene_exp',header=8)
assert set(M.title)==set(x.columns[5:]) and len(x.columns)==13
vals=x.loc[:,M.title].apply(pd.to_numeric,errors='raise');assert (vals.to_numpy()>=0).all()
x=x.assign(Gene_Name=x.Gene_Name.astype(str).str.strip())
# Symbol occurs in more than one workbook row: no outcome-based row selection.
dup=set(x.loc[x.Gene_Name.duplicated(keep=False),'Gene_Name']);genes=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
records=[]
for gene in genes:
 z=x[x.Gene_Name==gene]
 if len(z)!=1 or gene in dup:records.append(dict(gene=gene,measured=False,reason='absent_or_ambiguous_symbol'));continue
 v=z[M.title].iloc[0].astype(float)
 if not np.isfinite(v).all() or (v<0).any():records.append(dict(gene=gene,measured=False,reason='nonfinite_or_negative'));continue
 # Already normalized FPKM. No re-normalization with library-size CPM.
 a=np.log2(v[M.loc[M.disease=='preeclampsia','title']].to_numpy()+1);b=np.log2(v[M.loc[M.disease=='control','title']].to_numpy()+1)
 g,var=stats.hedges_g(a.reshape(1,-1),b.reshape(1,-1));p=ttest_ind(a,b,equal_var=False).pvalue
 records.append(dict(gene=gene,measured=bool(np.isfinite(g[0])),g=float(g[0]),v=float(var[0]),welch_p=float(p),down=bool(g[0]<0),n_case=4,n_control=4))
R=pd.DataFrame(records);R.to_csv('results/pe_p19_panel.csv',index=False)
n=int(R.measured.sum());down=int(R.down.fillna(False).sum());bonf=R[(R.down==True)&(R.welch_p<.05/6)].gene.tolist()
out=dict(gse='GSE143953',source=url+'GSE143953',n_case=4,n_control=4,genes_measured=n,genes_down=down,same_direction_bonferroni=bonf,descriptive_panel_support=bool(n==6 and down==6 and bonf),workbook_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),interpretation='Late placenta, historical small series, different from P9 early chorionic villi; no novel biomarker or published comparator')
open('results/pe_p19_result.json','w').write(json.dumps(out,indent=2)+'\n');print(R.to_string(index=False));print(json.dumps(out,indent=2))
