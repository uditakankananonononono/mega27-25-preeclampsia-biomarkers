"""P24 registered first-trimester chorionic-villus six-gene PE panel transport."""
import csv,hashlib,json,sys,urllib.request
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,openpyxl
from scipy.stats import ttest_ind
from ubiomark.stats import hedges_g
root=Path(__file__).resolve().parents[1]
matrix=root/'data/geo/p24/GSE295760_series_matrix.txt.gz'
workbook=root/'data/geo/p24/GSE295760_CVS_mRNA_gene_counts.xlsx'
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295760/'
for path,url,expected in [
 (matrix,base+'matrix/GSE295760_series_matrix.txt.gz','db778a7e309477ed18777dfa11d9beddc985f7d5ee96177e095675a3791e8d11'),
 (workbook,base+'suppl/GSE295760_CVS_mRNA_gene_counts.xlsx','bb9bfde3c93ae486070306fe1cd602b7315a164f105f89b67bc087481121fcf4')]:
 if not path.exists():
  path.parent.mkdir(parents=True,exist_ok=True)
  urllib.request.urlretrieve(url,path)
 assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, f'Source checksum changed: {path}'
from ubiomark.geo import parse_series_matrix
_,ann,_=parse_series_matrix(str(matrix))
assert len(ann)==38
m=ann.iloc[:19];small=ann.iloc[19:]
assert m['Sample_library_strategy'].eq('RNA-Seq').all()
assert small['Sample_title'].str.contains('small RNAseq').all()
assert m.index.is_unique and small.index.is_unique and set(m.index).isdisjoint(small.index)
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())}
assert set(m.index).isdisjoint(prior)
titles={str(v).split(',')[-1].strip():i for i,v in m['Sample_title'].items()}
assert len(titles)==19
out={i:v.split(': ',1)[-1] for i,v in m['Sample_characteristics_ch1_1'].items()}
assert Counter(out.values())==Counter({'Preterm Preeclampsia':4,'Term Preeclampsia':6,'Normotensive control':9})
w=openpyxl.load_workbook(workbook,read_only=True,data_only=True).active
rows=w.values; header=next(rows)
assert header[0]=='geneid' and len(header)==20
assert {str(x) for x in header[1:]}==set(titles)
sample_order=[titles[str(x)] for x in header[1:]]
fixed=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
targets={x:[] for x in fixed}
totals=np.zeros(19,dtype=float)
count=0
for row in rows:
 if not row or not row[0]:continue
 vals=np.array(row[1:20],dtype=float)
 assert vals.shape==(19,) and np.isfinite(vals).all() and (vals>=0).all() and np.equal(vals,np.floor(vals)).all()
 totals+=vals;count+=1
 if str(row[0]) in targets:targets[str(row[0])].append(vals)
assert count>10000 and (totals>0).all()
case1=np.array([out[x]=='Preterm Preeclampsia' for x in sample_order]);control=np.array([out[x]=='Normotensive control' for x in sample_order]);allcase=np.array([out[x] in ('Preterm Preeclampsia','Term Preeclampsia') for x in sample_order]);assert (case1.sum(),control.sum(),allcase.sum())==(4,9,10)
conflicts=[]
for gsm,row in small.iterrows():
 token=str(row['Sample_title']).split('gestation, ',1)[1].split(', small RNAseq')[0]
 assert token in titles
 mrna_gsm=titles[token]
 mrna_label=out[mrna_gsm]
 small_label=str(row['Sample_characteristics_ch1_1']).split(': ',1)[-1]
 if mrna_label!=small_label:
  conflicts.append(dict(token=token,mrna_gsm=mrna_gsm,mrna_outcome=mrna_label,small_rna_gsm=gsm,small_rna_outcome=small_label))
assert len(conflicts)==4, f'Cross-modality labels changed: {conflicts}'
result=[]
for gene in fixed:
 if len(targets[gene])!=1:
  result.append(dict(gene=gene,status='unmeasurable',rows=len(targets[gene])));continue
 x=np.log2(targets[gene][0]/totals*1e6+1)
 r=dict(gene=gene,status='measured',rows=1)
 for label,cases in [('preterm',case1),('all_pe_secondary',allcase)]:
  g=float(hedges_g(x[cases,None].T,x[control,None].T)[0][0]);p=float(ttest_ind(x[cases],x[control],equal_var=False).pvalue)
  r.update({label+'_g':g,label+'_p':p,label+'_down':bool(g<0),label+'_bonferroni':bool(g<0 and p<.05/6)})
 result.append(r)
passrule=bool(all(r.get('preterm_down') for r in result) and any(r.get('preterm_bonferroni') for r in result))
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295760',matrix_sha256=hashlib.sha256(matrix.read_bytes()).hexdigest(),workbook_sha256=hashlib.sha256(workbook.read_bytes()).hexdigest(),n_genes=count,n_samples=19,source_summary_says='14 samples (4 preterm PE, 4 term PE, 6 controls)',actual_mrna_metadata=dict(Counter(out.values())),sample_library_totals=totals.tolist(),panel=result,registered_primary_rule_pass=None,calculated_rule_under_mrna_labels=passrule,interpretation='uninterpretable: four same-token cross-modality outcome-label conflicts',label_conflicts=conflicts)
(root/'results/pe_p24_result.json').write_text(json.dumps(rec,indent=2)+'\n')
with (root/'results/pe_p24_panel.csv').open('w',newline='') as f:
 field=['gene','status','rows','preterm_g','preterm_p','preterm_down','preterm_bonferroni','all_pe_secondary_g','all_pe_secondary_p','all_pe_secondary_down','all_pe_secondary_bonferroni'];w=csv.DictWriter(f,fieldnames=field);w.writeheader();w.writerows(result)
with (root/'results/pe_p24_sample_map.csv').open('w',newline='') as f:
 field=['gsm','token','outcome','library_total'];w=csv.DictWriter(f,fieldnames=field);w.writeheader();w.writerows([dict(gsm=sample_order[i],token=header[i+1],outcome=out[sample_order[i]],library_total=int(totals[i])) for i in range(19)])
print(json.dumps({k:v for k,v in rec.items() if k!='sample_library_totals'},indent=2))
