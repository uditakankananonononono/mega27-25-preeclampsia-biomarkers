"""Registered P23 archival GSE44711 fixed PE panel test; see preregistered_predictions.md."""
import csv, json, sys, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.stats import ttest_ind
from ubiomark import geo, stats

root=Path(__file__).resolve().parents[1]
path=root/'data/raw/matrix/GSE44711_series_matrix.txt.gz'
expr, ann, _ = geo.parse_series_matrix(str(path))
titles=ann['Sample_title'].astype(str)
case=titles.index[titles.str.fullmatch(r'ChorionicVilli_EOPET_[1-8] \(expression\)')].tolist()
control=titles.index[titles.str.fullmatch(r'ChorionicVilli_CONTROL_[1-8] \(expression\)')].tolist()
assert len(case)==len(control)==8 and len(set(case+control))==16 and len(ann)==16
assert set(case+control)=={f'GSM{i}' for i in range(1089229,1089245)}
assert ann['Sample_platform_id'].eq('GPL10558').all()
old={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())}
assert not old.intersection(case+control), 'GSM overlap with existing manifest'
p2s=geo.probe_to_symbol('GPL10558')
genes=geo.to_gene_level(expr,p2s)
assert set(case+control).issubset(set(genes.columns))
fixed=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1']
rows=[]
for gene in fixed:
 if gene not in genes.index: rows.append(dict(gene=gene,status='unmeasurable'));continue
 a=genes.loc[gene,case].to_numpy(float);b=genes.loc[gene,control].to_numpy(float)
 if not np.isfinite(a).all() or not np.isfinite(b).all(): rows.append(dict(gene=gene,status='nonfinite'));continue
 g=float(stats.hedges_g(a[:,None].T,b[:,None].T)[0][0])
 p=float(ttest_ind(a,b,equal_var=False).pvalue)
 rows.append(dict(gene=gene,status='measured',g=g,welch_p=p,pe_mean=float(a.mean()),control_mean=float(b.mean()),down=bool(g<0),bonferroni=bool(g<0 and p<.05/6)))
ages=[]
for id in case+control:
 age=float(ann.loc[id,'Sample_characteristics_ch1_2'].split(':',1)[1]);ages.append(age)
maprows=[dict(gsm=id,title=titles[id],label='PE' if id in case else 'control',gestational_weeks=float(ann.loc[id,'Sample_characteristics_ch1_2'].split(':',1)[1])) for id in case+control]
result=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44711',matrix_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),n_pe=8,n_control=8,pe_gestational_range=[min(ages[:8]),max(ages[:8])],control_gestational_range=[min(ages[8:]),max(ages[8:])],panel=rows,registered_rule_pass=bool(all(r.get('down') for r in rows) and any(r.get('bonferroni') for r in rows)))
with (root/'results/pe_p23_panel.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gene','status','g','welch_p','pe_mean','control_mean','down','bonferroni']);w.writeheader();w.writerows(rows)
with (root/'results/pe_p23_sample_map.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(maprows[0]));w.writeheader();w.writerows(maprows)
(root/'results/pe_p23_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
