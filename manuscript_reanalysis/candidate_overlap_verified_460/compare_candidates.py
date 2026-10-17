"""Exact SNP overlaps on the common tested universe; run from project root."""
from pathlib import Path
import pandas as pd
import hashlib, json, itertools

root = Path('manuscript_reanalysis')
out = root/'candidate_overlap_verified_460'
rda = root/'rda_verified_460'
gwas = root/'egwas_verified_460'
sources = [rda/'input/marker_manifest.csv', rda/'results/climate_only_putative_outliers_3SD.csv', rda/'results/primary_adjusted_putative_outliers_3SD.csv', gwas/'gwas_summary/per_traitsignificant_associations.csv', gwas/'gwas_summary/complete_scansignificant_associations.csv', root/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim']
markers, simple, partial, hits, strict = [pd.read_csv(p) for p in sources[:5]]
bim = pd.read_csv(sources[5],sep=r'\s+',header=None,names=['chr','SNP','cm','position','A1','A2'])
assert len(bim)==41654 and not bim.SNP.duplicated().any()
rsamples=pd.read_csv(rda/'input/sample_manifest.csv')
gsamples=pd.read_csv(gwas/'input/sample_order.tsv',sep='\t')
assert len(rsamples)==len(gsamples)==460
assert list(rsamples.Accession)==list(gsamples.IID)
parsed=pd.read_csv(gwas/'post_gwas/parsed/per_trait/cowpea460_gwas_hitlst.txt',sep='\t')
assert hits.equals(parsed)
assert len(markers)==20027 and not markers.SNP.duplicated().any()
assert len(simple)==437 and len(partial)==763
assert set(markers.SNP)<=set(bim.SNP)
bm=bim.set_index('SNP')
for frame,idcol,chrcol,poscol,a1,a2 in [(markers,'SNP','PLINK_chromosome_code','position','A1','A2'),(simple,'SNP','PLINK_chromosome_code','position','A1','A2'),(partial,'SNP','PLINK_chromosome_code','position','A1','A2'),(hits,'rs','chr','ps','allele1','allele0')]:
 for row in frame.to_dict('records'):
  ref=bm.loc[row[idcol]]
  assert str(row[chrcol])==str(ref['chr']) and row[poscol]==ref.position
  assert {row[a1],row[a2]}=={ref.A1,ref.A2}
assert set(simple.SNP)|set(partial.SNP)<=set(markers.SNP)
assert len(hits)==239 and hits.rs.nunique()==98
assert len(strict)==34 and strict.rs.nunique()==17
assert set(zip(strict.ID,strict.rs))==set(zip(hits.loc[hits.p_lrt<.05/(41654*92),'ID'],hits.loc[hits.p_lrt<.05/(41654*92),'rs']))
assert (hits.p_lrt<.05/41654).all()
S=set(simple.SNP);P=set(partial.SNP);U=set(markers.SNP)
allids=S|P|set(hits.rs)
membership=bim[bim.SNP.isin(allids)].copy().drop(columns='cm')
membership['tested_by_RDA']=membership.SNP.isin(U)
membership['simple_RDA']=membership.SNP.map(lambda x:'candidate' if x in S else ('not_selected' if x in U else 'not_tested'))
membership['partial_RDA']=membership.SNP.map(lambda x:'candidate' if x in P else ('not_selected' if x in U else 'not_tested'))
for name,df in [('simple',simple),('partial',partial)]:
 for col in ['assigned_predictor','association_r','max_abs_loading_z','outlier_axes']:
  membership[name+'_'+col]=membership.SNP.map(df.set_index('SNP')[col])
for label,df in [('per_trait',hits),('across92',strict)]:
 membership['eGWAS_'+label]=membership.SNP.isin(df.rs)
 membership[label+'_traits']=membership.SNP.map(df.groupby('rs').ID.agg(lambda s:';'.join(sorted(set(s)))))
 membership[label+'_min_P']=membership.SNP.map(df.groupby('rs').p_lrt.min())
membership.to_csv(out/'candidate_membership_all.csv',index=False)
membership[membership.tested_by_RDA].to_csv(out/'candidate_membership_common_universe.csv',index=False)
membership[~membership.tested_by_RDA].to_csv(out/'eGWAS_hits_not_tested_by_RDA.csv',index=False)
membership[(membership.simple_RDA=='candidate')&(membership.partial_RDA=='candidate')&membership.eGWAS_per_trait].to_csv(out/'shared_all_three_per_trait.csv',index=False)
membership[(membership.partial_RDA=='candidate')&membership.eGWAS_per_trait].to_csv(out/'shared_partial_RDA_eGWAS.csv',index=False)
summary={}
for label,df in [('per_trait',hits),('across92',strict)]:
 G=set(df.rs);sets=[S,P,G&U]; rows=[]
 for bits in itertools.product([0,1],repeat=3):
  if not any(bits):continue
  ids=set(U)
  for b,s in zip(bits,sets):ids &= s if b else U-s
  rows.append(dict(zip(['simple_RDA','partial_RDA','eGWAS'],bits),n_SNPs=len(ids),SNPs=';'.join(sorted(ids))))
 pd.DataFrame(rows).to_csv(out/(label+'_exact_intersections.csv'),index=False)
 assert sum(x['n_SNPs'] for x in rows)==len(S|P|(G&U))
 summary[label]={'eGWAS_all_hits':len(G),'eGWAS_common_tested':len(G&U),'eGWAS_not_tested_by_RDA':len(G-U),'simple_and_eGWAS':len(S&G),'partial_and_eGWAS':len(P&G),'all_three':len(S&P&G),'union_common_tested':len(S|P|(G&U)),'union_all_hits':len(S|P|G),'exact_intersections':rows}
report={'common_tested_SNPs':len(U),'eGWAS_tested_SNPs':len(bim),'simple_candidates':len(S),'partial_candidates':len(P),'simple_partial_shared':len(S&P),'thresholds':{'RDA_sample_SD':3,'eGWAS_per_trait':.05/41654,'eGWAS_across92':.05/(41654*92)},'comparisons':summary,'checks':'SNP IDs, chromosome, position, allele sets, thresholds, unique counts and exact-intersection totals verified','input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'interpretation':'Exact typed-SNP overlap, not LD-region overlap or independent validation. Models use the same cohort, different predictor sets and different statistical criteria. No enrichment or overlap significance test performed.'}
(out/'overlap_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:{a:b for a,b in v.items() if a!='exact_intersections'} for k,v in summary.items()},indent=2))
