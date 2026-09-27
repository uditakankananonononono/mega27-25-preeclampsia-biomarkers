"""P44 frozen six down-direction PE genes in archived 6/6 placental FPKM."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p44';p=folder/'GSE279757_Processed_Data.xlsx';expected='05fc9d8d9235ef4a9c2ec2060790c5a078121e0a94f8006035bb3456c37c42ca';assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
rows=list(csv.DictReader((root/'results/pe_p44_samples.csv').open()));assert len(rows)==12 and len({z['gsm'] for z in rows})==len({z['column'] for z in rows})==12
for z in rows:assert hashlib.sha256((folder/'gsm'/(z['gsm']+'.txt')).read_bytes()).hexdigest()==z['sha256']
x=pd.read_excel(p,dtype={'ID':str});assert x.ID.str.fullmatch(r'ENSG\d+').all() and x.ID.is_unique and set(z['column'] for z in rows)==set(x.columns[2:14]);y=x.set_index('ID')[[z['column'] for z in rows]].copy();assert np.isfinite(y.to_numpy()).all() and (y.to_numpy()>=0).all()
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];mapping=dict(zip(h.ensembl_gene_id,h.symbol));y['gene']=y.index.map(mapping);mapped=y.dropna(subset=['gene']);dups=set(mapped.loc[mapped.gene.duplicated(keep=False),'gene']);expr=mapped[~mapped.gene.isin(dups)].set_index('gene');assert expr.index.is_unique;log=np.log2(expr.loc[(expr>1).mean(axis=1)>=.2]+1)
case=[z['column'] for z in rows if z['group']=='PE'];ctrl=[z['column'] for z in rows if z['group']=='N'];assert len(case)==len(ctrl)==6;g,v=stats.hedges_g(log[case].to_numpy(),log[ctrl].to_numpy());E=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna();panel=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1'];rec=[]
for gene in panel:
 if gene in E.index:rec.append(dict(gene=gene,status='measured',g=float(E.loc[gene,'g']),v=float(E.loc[gene,'v']),welch_p=float(ttest_ind(log.loc[gene,case],log.loc[gene,ctrl],equal_var=False).pvalue),down=bool(E.loc[gene,'g']<0)))
 else:rec.append(dict(gene=gene,status='not uniquely mapped or expression filtered'))
count=sum(z.get('down',False) for z in rec);corrected=sum(z.get('down',False) and z.get('welch_p',1)<.05/6 for z in rec)
out=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279757',source_sha256=expected,n_case=6,n_control=6,n_source_ensg_rows=len(x),n_unique_mapped_genes=len(expr),n_effect_genes=len(E),n_fixed_measured=sum(z['status']=='measured' for z in rec),n_down=count,n_corrected_down=corrected,registered_descriptive_support=bool(len(rec)==6 and all(z.get('down') for z in rec) and corrected>=1),ambiguous_ensembl_ids=len(amb),duplicate_gene_symbols_excluded=len(dups),genes=rec,limitation='Six samples per arm, archived placental FPKM, at-delivery stage and selected panel do not prove early prediction, clinical utility, novel mechanism or benchmark superiority.')
(root/'results/pe_p44_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
