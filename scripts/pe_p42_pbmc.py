"""P42 fixed six DOWN PE genes in independent late PBMC pure-PE vs CTL."""
import csv,hashlib,json,re,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p42';p=folder/'GSE284329_Transcriptomics_fcounts.txt.gz';expected='f3968cba1cb9fc1d4a76f3373f40a0a2be3ec8952015ac6bd2b6ab465f8033f1';assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
rows=list(csv.DictReader((root/'results/pe_p42_samples.csv').open()));assert len(rows)==17 and len({z['gsm'] for z in rows})==len({z['column'] for z in rows})==17 and {g:sum(z['group']==g for z in rows) for g in {'CTL','GDM','PE','PG'}}=={'CTL':5,'GDM':5,'PE':3,'PG':4}
for z in rows:assert hashlib.sha256((folder/'gsm'/(z['gsm']+'.txt')).read_bytes()).hexdigest()==z['sha256']
x=pd.read_csv(p,sep='\t',index_col=0);assert x.index.is_unique and x.index.str.fullmatch(r'ENSG\d+').all() and set(x.columns[5:])=={z['column'] for z in rows}
y=x[[z['column'] for z in rows if z['group'] in ('PE','CTL')]].copy();assert np.isfinite(y.to_numpy()).all() and (y.to_numpy()>=0).all() and np.equal(y.to_numpy(),np.floor(y.to_numpy())).all()
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];mapping=dict(zip(h.ensembl_gene_id,h.symbol));y['gene']=y.index.map(mapping);expr=y.dropna(subset=['gene']).groupby('gene').sum();assert expr.index.is_unique;lib=expr.sum(axis=0);assert (lib>0).all()
cpm=expr.div(lib,axis=1)*1e6;log=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1);case=[z['column'] for z in rows if z['group']=='PE'];ctrl=[z['column'] for z in rows if z['group']=='CTL'];g,v=stats.hedges_g(log[case].to_numpy(),log[ctrl].to_numpy());E=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna();panel=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1'];obs=[z for z in panel if z in E.index]
result=[]
for gene in panel:
 if gene in obs:
  pvalue=float(ttest_ind(log.loc[gene,case].to_numpy(),log.loc[gene,ctrl].to_numpy(),equal_var=False).pvalue)
  result.append(dict(gene=gene,status='measured',g=float(E.loc[gene,'g']),v=float(E.loc[gene,'v']),welch_p=pvalue,down=bool(E.loc[gene,'g']<0)))
 else:result.append(dict(gene=gene,status='not uniquely mapped or expression filtered'))
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284329',source_sha256=expected,n_primary_case=3,n_primary_control=5,n_excluded_GDM=5,n_excluded_PG=4,n_unique_mapped_genes=len(expr),n_effect_genes=len(E),n_fixed_measured=len(obs),n_fixed_down=sum(z.get('down',False) for z in result),n_corrected_down=sum(z.get('down',False) and z.get('welch_p',1)<.05/6 for z in result),registered_descriptive_support=bool(len(obs)==6 and all(z.get('down') for z in result) and any(z.get('down') and z.get('welch_p',1)<.05/6 for z in result)),ambiguous_ensembl_ids=len(amb),genes=result,limitation='Only 3 pure-PE vs 5 control PBMC persons, late-pregnancy specimen, old postselected panel, not early CVS prediction or published benchmark.')
(root/'results/pe_p42_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
