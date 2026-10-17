from pathlib import Path
import pandas as pd,json,hashlib
r=Path('manuscript_reanalysis');out=r/'shared_candidate_genes_460/expanded_two_methods';ref=r/'egwas_verified_460/source/reference'
d=pd.read_csv(r/'candidate_overlap_verified_460/candidate_membership_common_universe.csv').fillna('')
d['n_methods']=d.simple_RDA.eq('candidate').astype(int)+d.partial_RDA.eq('candidate').astype(int)+d.eGWAS_per_trait.eq(True).astype(int)
d=d[d.n_methods>=2].sort_values(['chr','position']).copy();assert len(d)==116 and not d.SNP.duplicated().any()
d['methods']=d.apply(lambda x:'; '.join(k for k,v in [('simple RDA',x.simple_RDA=='candidate'),('pRDA',x.partial_RDA=='candidate'),('eGWAS',x.eGWAS_per_trait==True)] if v),axis=1)
a=pd.read_csv(ref/'Vunguiculata_540_v1.2.annotation_info.txt',sep='\t').fillna('')
ann={g:{'Arabidopsis_best_hit':'; '.join(sorted(set(x['Best-hit-arabi-name'])-{''})),'annotation_description':'; '.join(sorted(set(x['Best-hit-arabi-defline'])-{''}))} for g,x in a.groupby('locusName')}
genes={}
for line in open(ref/'Vunguiculata_540_v1.2.gene.gff3'):
 z=line.rstrip().split('\t')
 if len(z)==9 and z[2]=='gene':
  attrs=dict(x.split('=',1) for x in z[8].split(';') if '=' in x);genes[attrs['Name']]=(z[0],int(z[3]),int(z[4]))
# Hand-reviewed primary-study evidence. Each statement is restricted to the tested species/gene.
refs={
 'AT5G45340':('CYP707A3','Arabidopsis mutant and overexpression experiments support ABA catabolism and dehydration/rehydration responses.','Experimental evidence for the Arabidopsis best-hit gene','Umezawa et al. (2006)','https://doi.org/10.1111/j.1365-313X.2006.02683.x'),
 'AT2G44745':('WRKY12','Arabidopsis WRKY12 represses secondary cell-wall formation in pith; this is developmental evidence, not a demonstrated drought-tolerance effect.','Experimental evidence for the Arabidopsis best-hit gene','Wang et al. (2010)','https://doi.org/10.1073/pnas.1016436107'),
 'AT2G30070':('KUP1 / KT1','Arabidopsis KUP1 mediates high- and low-affinity potassium uptake in yeast; uptake was inhibited by NaCl.','Heterologous transport assay of the Arabidopsis best-hit gene','Fu and Luan (1998)','https://doi.org/10.1105/tpc.10.1.63'),
 'AT2G22910':('NAGS1','Overexpression of tomato SlNAGS1 in Arabidopsis increased ornithine and salt/drought tolerance; this is family-level evidence, not a test of AT2G22910.','Related tomato gene expressed in Arabidopsis; family-level support only','Kalamaki et al. (2009)','https://doi.org/10.1093/jxb/erp072'),
 'AT5G25190':('ESE3','Arabidopsis ESE3 transcript was induced by ethylene precursor and salt; functional salt-tolerance experiments in this paper focused on ESE1, not ESE3.','Expression evidence for the Arabidopsis best-hit gene','Zhang et al. (2011)','https://doi.org/10.1104/pp.111.179028')}
source='Vunguiculata_540_v1.2.annotation_info.txt; Vunguiculata_540_v1.2.gene.gff3'
source_ref='Lonardi et al. (2019), https://doi.org/10.1111/tpj.14349 (reference genome and annotation framework; not functional validation)'
def annotation(g):
 z=ann.get(g,{'Arabidopsis_best_hit':'','annotation_description':''}).copy();z['gene']=g
 v=refs.get(z['Arabidopsis_best_hit']);z.update(dict(zip(['gene_designation','published_function_evidence','literature_evidence_level','function_reference','function_reference_URL'],v or ('','No gene-specific primary study verified in this targeted review.','Annotation only; targeted literature review incomplete','',''))))
 z['annotation_source']=source;z['annotation_reference']=source_ref
 if not z['annotation_description']:z['annotation_description']='Not supplied'
 return z
records=[];overview=[]
ld=pd.read_csv(r/'egwas_verified_460/go_enrichment/Marker_evidence.csv');ld=ld[ld.significant_SNP.isin(d.SNP)].copy()
bim=pd.read_csv(r/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim',sep=r'\s+',header=None,names=['chr','SNP','cm','bp','A1','A2']).set_index('SNP')
for row in d.to_dict('records'):
 direct=sorted(g for g,(c,lo,hi) in genes.items() if c==row['chr'] and lo<=row['position']<=hi)
 z=row.copy();z['direct_genes']='; '.join(direct) or 'None';z['n_direct_genes']=len(direct)
 for key in ['Arabidopsis_best_hit','annotation_description','gene_designation','published_function_evidence','function_reference','function_reference_URL']:
  z[key]=' | '.join(f'{g}: {annotation(g)[key] or "Not supplied"}' for g in direct) or 'No direct gene-body overlap'
 z['annotation_reference']=source_ref;z['LD_review_scope']='Saved eGWAS LD evidence available' if row['eGWAS_per_trait'] else 'LD not assessed for RDA-only candidates in this table'
 overview.append(z)
 for g in direct:records.append({**row,**annotation(g),'mapping_evidence':'direct_gene_body_overlap','supporting_marker':row['SNP'],'LD_r2':'','marker_distance_bp':0})
 for e in ld[(ld.significant_SNP==row['SNP'])&(ld.evidence!='direct_overlap')].drop_duplicates(['gene','supporting_marker']).to_dict('records'):
  c,lo,hi=genes[e['gene']];m=bim.loc[e['supporting_marker']]
  assert c==row['chr']==m['chr'] and lo<=m.bp<=hi and e['r2']>.2 and abs(m.bp-row['position'])<1000000
  records.append({**row,**annotation(e['gene']),'mapping_evidence':'LD_supported_gene_body_marker','supporting_marker':e['supporting_marker'],'LD_r2':e['r2'],'marker_distance_bp':abs(m.bp-row['position'])})
summary=pd.DataFrame(overview);detail=pd.DataFrame(records)
summary.to_csv(out/'Table_S7b_116_shared_SNPs.csv',index=False)
detail.to_csv(out/'Table_S7c_gene_and_LD_evidence.csv',index=False)
pd.DataFrame([dict(Arabidopsis_locus=k,**dict(zip(['designation','evidence','evidence_level','citation','URL'],v))) for k,v in refs.items()]).to_csv(out/'Table_S7d_function_references.csv',index=False)
counts={'SNPs_at_least_two_methods':len(d),'three_methods':int((d.n_methods==3).sum()),'simple_and_partial_only':int(((d.simple_RDA=='candidate')&(d.partial_RDA=='candidate')&(~d.eGWAS_per_trait.astype(bool))).sum()),'SNPs_with_direct_gene':int((summary.n_direct_genes>0).sum()),'distinct_direct_genes':detail[detail.mapping_evidence=='direct_gene_body_overlap'].gene.nunique(),'gene_marker_evidence_rows':len(detail),'direct_genes_with_targeted_primary_references':detail[(detail.mapping_evidence=='direct_gene_body_overlap')&(detail.function_reference!='')].gene.nunique()}
counts['input_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [r/'candidate_overlap_verified_460/candidate_membership_common_universe.csv',ref/'Vunguiculata_540_v1.2.annotation_info.txt',ref/'Vunguiculata_540_v1.2.gene.gff3',r/'egwas_verified_460/go_enrichment/Marker_evidence.csv']}
(out/'validation.json').write_text(json.dumps(counts,indent=2)+'\n')
(out/'README.txt').write_text('Table S7b: one row per SNP detected by at least two methods. Membership uses per-trait eGWAS significance; across-92 status is separate. RDA/pRDA use adopted 3 SD screening. RDA predictor labels are strongest correlations, not independent per-variable tests. All significant eGWAS traits are retained. Direct annotation uses inclusive gene-body overlap; it is not a coding-effect prediction.\nTable S7c: all direct SNP-gene links plus previously computed eGWAS LD-supported gene-marker links. LD was not newly calculated for RDA-only candidates. Multiple rows are evidence records, not independent genes.\nTable S7d: targeted primary-literature references; not an exhaustive literature review. Experimental species and family-level versus exact-best-hit evidence are explicit. Cowpea function and causal adaptive variants are unvalidated. All missing annotations are retained. Simple RDA and pRDA are related methods, so their overlap is not independent validation.\nAnnotation descriptions reproduce the supplied reference annotation; genome-source citation is not evidence for individual gene function. No nearest-gene assignment or new GO analysis was performed. Main Table 4 remains the 15 pRDA/eGWAS intragenic SNPs.\n')
print(json.dumps({k:v for k,v in counts.items() if k!='input_sha256'},indent=2))
# Compact manuscript-facing version inspired by Omondi et al. Table 3.
compact=[]
for z in overview:
 gs=z['direct_genes'].split('; ') if z['n_direct_genes'] else []
 env=[]
 if z['simple_RDA']=='candidate':env.append('RDA: '+z['simple_assigned_predictor'])
 if z['partial_RDA']=='candidate':env.append('pRDA: '+z['partial_assigned_predictor'])
 if z['eGWAS_per_trait']:env.append('eGWAS: '+z['per_trait_traits'].replace(';',', '))
 compact.append({'SNP ID':z['SNP'],'Detected by':z['methods'],'Environmental variables':'; '.join(env),'Cowpea gene / Arabidopsis best hit':'; '.join(g+' / '+(annotation(g)['gene_designation'] or annotation(g)['Arabidopsis_best_hit'] or 'No best hit supplied') for g in gs) or 'No direct annotated gene overlap','Annotated function':'; '.join(g+': '+annotation(g)['annotation_description'] for g in gs) or 'See separate LD evidence table where available','Published functional evidence':'; '.join(g+': '+annotation(g)['published_function_evidence'] for g in gs if annotation(g)['function_reference']) or 'No gene-specific primary reference verified in this targeted review','References':'; '.join(annotation(g)['function_reference']+' '+annotation(g)['function_reference_URL'] for g in gs if annotation(g)['function_reference']) or 'Reference annotation: Lonardi et al. (2019), https://doi.org/10.1111/tpj.14349'})
pd.DataFrame(compact).to_csv(out/'Table_S7b_publication_candidate_table.csv',index=False)
short={
'AT5G45340':('ABA catabolism and dehydration/rehydration responses (Arabidopsis).','1'),
'AT2G44745':('Repression of secondary cell-wall formation in pith (Arabidopsis).','2'),
'AT2G30070':('High- and low-affinity potassium transport; NaCl-sensitive uptake in yeast assays.','3'),
'AT2G22910':('N-acetylglutamate synthesis; related tomato SlNAGS1 increased salt/drought tolerance when overexpressed in Arabidopsis.','4'),
'AT5G25190':('Ethylene- and salt-induced gene expression (Arabidopsis).','5')}
pub=[]
for z in overview:
 gs=z['direct_genes'].split('; ') if z['n_direct_genes'] else []
 env=[]
 if z['simple_RDA']=='candidate':env.append('RDA: '+z['simple_assigned_predictor'])
 if z['partial_RDA']=='candidate':env.append('pRDA: '+z['partial_assigned_predictor'])
 if z['eGWAS_per_trait']:env.append('eGWAS: '+z['per_trait_traits'].replace(';',', '))
 funcs=[]
 for g in gs:
  a0=annotation(g);fun,num=short.get(a0['Arabidopsis_best_hit'],(a0['annotation_description'],'A'))
  funcs.append((g+': ' if len(gs)>1 else '')+fun+' ['+num+']')
 pub.append({'SNP ID':z['SNP'],'Detected by':z['methods'],'Environmental variables':'; '.join(env),'Gene / protein':'; '.join(g+' / '+(annotation(g)['gene_designation'] or annotation(g)['Arabidopsis_best_hit'] or 'No best hit supplied') for g in gs) or 'No direct annotated gene overlap','Function':'; '.join(funcs) or 'Unassigned; available LD links are listed separately.'})
pd.DataFrame(pub).to_csv(out/'Table_S7b_publication_candidate_table.csv',index=False)
references=pd.read_csv(out/'Table_S7d_function_references.csv');references.insert(0,'Reference_number',[1,2,3,4,5]);references.to_csv(out/'Table_S7d_function_references.csv',index=False)
references=pd.concat([references,pd.DataFrame([{'Reference_number':'A','Arabidopsis_locus':'All annotation-only rows','designation':'Reference annotation','evidence':'Supplied Vunguiculata v1.2 gene annotation and Arabidopsis best-hit descriptions. Genome paper supports the reference resource, not experimental validation of each function.','evidence_level':'Annotation source','citation':'Lonardi et al. (2019). The genome of cowpea (Vigna unguiculata [L.] Walp.). The Plant Journal 98:767-782.','URL':'https://doi.org/10.1111/tpj.14349'}])],ignore_index=True)
references.to_csv(out/'Table_S7d_function_references.csv',index=False)
