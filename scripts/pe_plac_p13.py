"""Registered P13 severe-PE placental arm; other cohort held for ID mapping."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
G='GSE186257';PANEL=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
SRC=Path('data/raw/rnaseq/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz')
_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);assert len(ann)==44
assert ann.Sample_title.is_unique and ann.index.is_unique
ann['label']=ann.Sample_characteristics_ch1_0.map({'disease: Severe preeclampsia':'case','disease: Non-hypertensive Control':'control'})
assert ann.label.value_counts().to_dict()=={'case':26,'control':18}
prior=pd.read_csv('results/gsm_to_study.csv').set_index('accession')
assert all(prior.loc[gsm,'study_accession']==G for gsm in set(ann.index)&set(prior.index)), 'GSM reused by other study'
x=pd.read_csv(SRC,sep='\t').set_index('Unnamed: 0')
assert x.index.is_unique and set(x.columns)==set(ann.Sample_title) and x.notna().all().all() and (x>=0).all().all()
assert x.index.str.match(r'^ENSG\d+\$[^$]*$').all()
x.index=x.index.str.split('$').str[-1]
# Duplicate gene symbols in the original annotated ID field cannot be resolved by selecting a favorable row.
assert all(x.index.tolist().count(g)==1 for g in PANEL)
x=np.log2(x+1) # Already deposited as DESeq2 normalized counts, not raw counts or CPM.
a=ann.loc[ann.label.eq('case'),'Sample_title'];b=ann.loc[ann.label.eq('control'),'Sample_title']
rows=[]
for gene in PANEL:
 if gene not in x.index:rows.append({'gene':gene,'g':np.nan,'variance':np.nan,'welch_p':np.nan,'measured':False});continue
 ca=x.loc[gene,a].to_numpy()[None,:];co=x.loc[gene,b].to_numpy()[None,:]
 g,v=stats.hedges_g(ca,co)
 rows.append({'gene':gene,'g':float(g[0]),'variance':float(v[0]),'welch_p':float(ttest_ind(ca[0],co[0],equal_var=False).pvalue),'measured':True})
out=pd.DataFrame(rows);out.to_csv('results/pe_plac_p13_GSE186257.csv',index=False)
ann[['Sample_title','label','Sample_characteristics_ch1_2']].rename_axis('gsm').to_csv('results/pe_plac_p13_GSE186257_samples.csv')
res={'gse':G,'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186257/suppl/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz','sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'n_case':26,'n_control':18,'measured':int(out.measured.sum()),'down':int((out.g<0).sum()),'bonferroni_down_genes':out.loc[(out.g<0)&(out.welch_p<.05/6),'gene'].tolist(),'matrix_kind':'DESeq2-normalized filtered, transformed log2(value+1); not raw counts'}
res['single_study_panel_pass']=bool(res['measured']==6 and res['down']==6 and len(res['bonferroni_down_genes'])>=1)
Path('results/pe_plac_p13_GSE186257.json').write_text(json.dumps(res,indent=2)+'\n');print(out.to_string(index=False));print(res)
