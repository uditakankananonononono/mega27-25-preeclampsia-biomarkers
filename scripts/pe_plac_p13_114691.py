"""P13 PE-only placenta versus controls in GSE114691 transcript-count release."""
import hashlib,json,re
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
from ubiomark import geo,stats
G='GSE114691';PANEL=['FES','GPAT3','FURIN','GRAMD1A','ZNF467','LPGAT1'];F='data/raw/rnaseq/GSE114691_MasterCount_'
_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);assert len(ann)==79 and ann.Sample_title.is_unique
prior=pd.read_csv('results/gsm_to_study.csv').set_index('accession')
assert all(prior.loc[gsm,'study_accession']==G for gsm in set(ann.index)&set(prior.index)), 'GSM reused by other study'
groups={'Control':'ControlONLY','Preeclampsia only':'PEONLY'}
def canon(s):return re.sub('[^A-Z0-9]','',s.upper())
parts=[];samples=[]
for label,file in groups.items():
 aa=ann[ann.Sample_characteristics_ch1_1.eq('disease group: '+label)].copy()
 assert len(aa)=={'Control':21,'Preeclampsia only':20}[label]
 src=Path(F+file+'.txt.gz');x=pd.read_csv(src,sep='\t').set_index('Unnamed: 0')
 assert x.index.is_unique and x.index.str.match(r'^ENST\d+$').all()
 assert x.notna().all().all() and (x>=0).all().all()
 keyed={canon(t):gsm for gsm,t in aa.Sample_description_0.items()}
 cols={c:keyed[canon(c.removeprefix('X'))] for c in x.columns}
 assert len(cols)==len(x.columns)==len(aa) and set(cols.values())==set(aa.index)
 for c,gsm in cols.items():samples.append({'gsm':gsm,'sample_title':aa.loc[gsm,'Sample_title'],'description':aa.loc[gsm,'Sample_description_0'],'matrix_column':c,'label':'control' if label=='Control' else 'case'})
 parts.append(x.rename(columns=cols))
x=pd.concat(parts,axis=1)
assert x.columns.is_unique and x.index.is_unique and not x.isna().any().any()
# Current Ensembl transcript membership: a limited gene-level extraction from this historical transcript-count release.
# Record complete map and number of historical IDs observed; never pick a favorable single transcript.
transcript={}
for gene in PANEL:
 j=json.loads(Path('results/external/p13_ensembl_transcript_map.json').read_text())['genes'][gene]
 transcript[gene]=set(j['transcript_ids'])&set(x.index)
 assert transcript[gene]
assert sum(len(z) for z in transcript.values())==len(set().union(*transcript.values()))
expr=pd.DataFrame({gene:x.loc[sorted(ids)].sum(axis=0) for gene,ids in transcript.items()}).T
expr=np.log2(1e6*expr/x.sum(axis=0)+1) # total transcript-library denominator; gene sums partial if old isoforms unmatched
s=pd.DataFrame(samples).set_index('gsm');a=s.index[s.label.eq('case')];b=s.index[s.label.eq('control')]
assert len(a)==20 and len(b)==21
rows=[]
for gene in PANEL:
 ca=expr.loc[gene,a].to_numpy()[None,:];co=expr.loc[gene,b].to_numpy()[None,:]
 g,v=stats.hedges_g(ca,co)
 rows.append({'gene':gene,'n_matched_transcripts':len(transcript[gene]),'g':float(g[0]),'variance':float(v[0]),'welch_p':float(ttest_ind(ca[0],co[0],equal_var=False).pvalue)})
out=pd.DataFrame(rows);out.to_csv('results/pe_plac_p13_GSE114691.csv',index=False);s.to_csv('results/pe_plac_p13_GSE114691_samples.csv')
res={'gse':G,'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114691/suppl/','n_case':20,'n_control':21,'down':int((out.g<0).sum()),'bonferroni_down_genes':out.loc[(out.g<0)&(out.welch_p<.05/6),'gene'].tolist(),'source_sha256':{file:hashlib.sha256(Path(F+file+'.txt.gz').read_bytes()).hexdigest() for file in groups.values()},'mapping':'Current Ensembl lookup/symbol/homo_sapiens/G?expand=1 transcript membership, only transcript IDs found in old release; incomplete historical isoform mapping is a limitation. Summed gene transcript counts then log2 CPM using total transcript counts for each sample; one sample per GEO title.'}
res['single_study_panel_pass']=bool(res['down']==6 and len(res['bonferroni_down_genes'])>=1)
res['matched_transcript_ids']={g:sorted(ids) for g,ids in transcript.items()}
Path('results/pe_plac_p13_GSE114691.json').write_text(json.dumps(res,indent=2)+'\n')
print(out.to_string(index=False));print({k:v for k,v in res.items() if k!='matched_transcript_ids'})
