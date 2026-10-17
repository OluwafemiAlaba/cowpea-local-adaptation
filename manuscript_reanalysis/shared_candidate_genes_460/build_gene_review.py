from pathlib import Path
import pandas as pd
import json, hashlib, gzip

root=Path('manuscript_reanalysis');out=root/'shared_candidate_genes_460'
gw=root/'egwas_verified_460';ref=gw/'source/reference'
src=root/'candidate_overlap_verified_460/shared_partial_RDA_eGWAS.csv'
s=pd.read_csv(src);ids=set(s.SNP);assert len(ids)==22
e=pd.read_csv(gw/'go_enrichment/Marker_evidence.csv');e=e[e.significant_SNP.isin(ids)].copy()
a=pd.read_csv(gw/'go_enrichment/Candidate_genes.csv').fillna('').set_index('gene')
genes={}
for line in (ref/'Vunguiculata_540_v1.2.gene.gff3').open():
 if line.startswith('#'):continue
 z=line.rstrip().split('\t')
 if z[2]!='gene':continue
 attrs=dict(t.split('=',1) for t in z[8].split(';') if '=' in t)
 genes[attrs['Name']]=(z[0],int(z[3]),int(z[4]))
b=pd.read_csv(root/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim',sep=r'\s+',header=None,names=['chr','SNP','cm','bp','A1','A2']).set_index('SNP')
links={}
with gzip.open(gw/'post_gwas/tags/cowpea460_candidate_local_LD.ld.gz','rt') as f:
 h=next(f).split()
 for line in f:
  z=dict(zip(h,line.split()));x,y=z['SNP_A'],z['SNP_B']
  if x in ids or y in ids:links[frozenset([x,y])]=float(z['R2'])
for r in e.drop_duplicates(['significant_SNP','gene','supporting_marker']).to_dict('records'):
 x,y=r['significant_SNP'],r['supporting_marker'];c,start,end=genes[r['gene']]
 assert b.loc[y,'chr']==c and start<=b.loc[y,'bp']<=end
 if x!=y:
  assert b.loc[x,'chr']==b.loc[y,'chr'] and abs(b.loc[x,'bp']-b.loc[y,'bp'])<1000000
  assert links[frozenset([x,y])]>.2 and abs(links[frozenset([x,y])]-r['r2'])<1e-6
 else:assert r['evidence']=='direct_overlap'
for r in s.to_dict('records'):
 direct={g for g,(c,lo,hi) in genes.items() if c==r['chr'] and lo<=r['position']<=hi}
 observed=set(e[(e.significant_SNP==r['SNP'])&(e.supporting_marker==r['SNP'])].gene)
 assert direct==observed
hits=pd.read_csv(gw/'gwas_summary/per_traitsignificant_associations.csv')
rows=[]
for r in s.to_dict('records'):
 d=e[e.significant_SNP==r['SNP']];direct=sorted(set(d[d.evidence=='direct_overlap'].gene))
 minimum=hits[hits.rs==r['SNP']].sort_values(['p_lrt','ID']).iloc[0]
 r.update(direct_genes='; '.join(direct),Arabidopsis_best_hit='; '.join(a.loc[g,'Arabidopsis_best_hit'] or 'Not supplied' for g in direct),annotation_description='; '.join(a.loc[g,'Arabidopsis_description'] or 'Not supplied' for g in direct),n_direct_genes=len(direct),n_LD_only_genes=len(set(d.gene)-set(direct)),minimum_P_trait=minimum.ID)
 rows.append(r)
table=pd.DataFrame(rows)
table.to_csv(out/'all_22_shared_SNP_annotations.csv',index=False)
table[table.n_direct_genes>0].to_csv(out/'Table4_intragenic_candidates.csv',index=False)
detail=e.merge(a.reset_index()[['gene','Arabidopsis_best_hit','Arabidopsis_description','rice_best_hit']],on='gene',validate='many_to_one')
detail.to_csv(out/'all_22_gene_marker_evidence.csv',index=False)
summary=[]
for g,d in e.groupby('gene'):
 direct=sorted(set(d[d.evidence=='direct_overlap'].significant_SNP)); strict=d[d.across92_significant]
 summary.append({'gene':g,'direct_shared_SNPs':'; '.join(direct),'all_shared_SNPs':'; '.join(sorted(set(d.significant_SNP))),'evidence':'direct and LD' if direct and (d.evidence=='LD_linked_marker_overlap').any() else ('direct' if direct else 'LD only'),'across92_supported':len(strict)>0,'Arabidopsis_best_hit':a.loc[g,'Arabidopsis_best_hit'],'annotation_description':a.loc[g,'Arabidopsis_description']})
pd.DataFrame(summary).to_csv(out/'shared_candidate_gene_summary.csv',index=False)
tags=pd.read_csv(gw/'post_gwas/tags/per_trait/tag_membership.txt',sep='\t')
tags[tags.SNP.isin(ids)].to_csv(out/'shared_SNP_tag_membership.csv',index=False)
counts={'shared_SNPs':22,'intragenic_SNPs':int((table.n_direct_genes>0).sum()),'SNPs_without_direct_gene_overlap':int((table.n_direct_genes==0).sum()),'direct_genes':e[e.evidence=='direct_overlap'].gene.nunique(),'direct_genes_without_Arabidopsis_best_hit':sum(a.loc[g,'Arabidopsis_best_hit']=='' for g in set(e[e.evidence=='direct_overlap'].gene)),'direct_or_LD_genes':e.gene.nunique(),'LD_only_genes':len(set(e.gene)-set(e[e.evidence=='direct_overlap'].gene)),'across92_supported_genes':e[e.across92_significant].gene.nunique(),'checks':'All direct assignments independently matched to GFF gene bodies; every LD evidence pair checked against saved genotype LD; original threshold r2 > 0.20 and distance < 1 Mb retained. No function inferred from proximity alone.'}
paths=[src,gw/'go_enrichment/Marker_evidence.csv',gw/'go_enrichment/Candidate_genes.csv',ref/'Vunguiculata_540_v1.2.gene.gff3',ref/'Vunguiculata_540_v1.2.annotation_info.txt']
counts['input_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(out/'validation.json').write_text(json.dumps(counts,indent=2)+'\n')
print(json.dumps({k:v for k,v in counts.items() if k!='input_sha256'},indent=2))
