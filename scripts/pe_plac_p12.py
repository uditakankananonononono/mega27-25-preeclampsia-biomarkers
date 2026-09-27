"""P12 preregistered PLAC-arm transport, one placenta per person in GSE190971."""
import hashlib,json,re
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
G='GSE190971'; PANEL=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
SRC=Path('data/raw/rnaseq/GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz')
_,ann,meta=geo.parse_series_matrix(geo.download_matrices(G)[0]);assert len(ann)==41
plac=ann[ann.Sample_characteristics_ch1_0.eq('tissue: Placenta')].copy()
assert len(plac)==13 and plac.Sample_title.is_unique and plac.index.is_unique
plac['person']=plac.Sample_title.str.extract(r'^mRNA_(V\d+)_PLAC$')[0].values
assert plac.person.notna().all() and plac.person.is_unique
plac['label']=plac.Sample_characteristics_ch1_1.map({'disease state (normal or pe): PE':'case','disease state (normal or pe): NP':'control'})
assert plac.label.value_counts().to_dict()=={'case':7,'control':6}
# Other arm metadata must agree on phenotype per participant; they never enter test.
allmeta=ann.copy();allmeta['person']=allmeta.Sample_title.str.extract(r'^mRNA_(V\d+)_(?:PLAC|10K|150K)$')[0].values
assert allmeta.person.notna().all() and allmeta.groupby('person').Sample_characteristics_ch1_1.nunique().eq(1).all()
used=set(pd.read_csv('results/gsm_to_study.csv').accession)
assert all(pd.read_csv('results/gsm_to_study.csv').set_index('accession').loc[gsm,'study_accession']==G for gsm in set(plac.index)&used), 'GSM reused by other study'
counts=pd.read_csv(SRC,sep='\t').set_index('Gene_Symbol')
assert counts.index.is_unique and counts.notna().all().all() and (counts>=0).all().all()
colmap={c:re.fullmatch(r'(V\d+)_PLAC_(PE|NORMAL)',c).groups() for c in counts.columns}
assert len(colmap)==13 and len(set(x[0] for x in colmap.values()))==13
assert set(p for p,_ in colmap.values())==set(plac.person)
for col,(person,group) in colmap.items():
 row=plac[plac.person.eq(person)].iloc[0]
 assert group=={'case':'PE','control':'NORMAL'}[row.label]
 assert row.Sample_title=='mRNA_'+person+'_PLAC'
plac['count_column']=plac.person.map({x[0]:c for c,x in colmap.items()})
assert plac.count_column.notna().all()
assert (counts.sum(axis=0)>0).all()
x=np.log2(1e6*counts/counts.sum(axis=0)+1)
case=plac.loc[plac.label.eq('case'),'count_column'];control=plac.loc[plac.label.eq('control'),'count_column']
rows=[]
for gene in PANEL:
 if gene not in x.index: rows.append({'gene':gene,'g':np.nan,'variance':np.nan,'welch_p':np.nan,'measured':False});continue
 a=x.loc[gene,case].to_numpy()[None,:];b=x.loc[gene,control].to_numpy()[None,:]
 g,v=stats.hedges_g(a,b)
 p=ttest_ind(a[0],b[0],equal_var=False).pvalue
 rows.append({'gene':gene,'g':float(g[0]),'variance':float(v[0]),'welch_p':float(p),'measured':True})
out=pd.DataFrame(rows);out.to_csv('results/pe_plac_p12.csv',index=False)
plac[['Sample_title','person','label','count_column']].rename_axis('gsm').to_csv('results/pe_plac_p12_samples.csv')
res={'gse':G,'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190971/suppl/GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz','sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'n_case':7,'n_control':6,'measured':int(out.measured.sum()),'down':int((out.g<0).sum()),'bonferroni_down_genes':out.loc[(out.g<0)&(out.welch_p<.05/6),'gene'].tolist()}
res['registered_panel_pass']=bool(res['measured']==6 and res['down']==6 and len(res['bonferroni_down_genes'])>=1)
res['independence_caveat']='PLAC is one arm in a 41-GSM three-arm series; 10K and 150K share participant tokens. Term placenta tissue is not first-trimester maternal screening.'
Path('results/pe_plac_p12.json').write_text(json.dumps(res,indent=2)+'\n')
print(out.to_string(index=False));print(res)
