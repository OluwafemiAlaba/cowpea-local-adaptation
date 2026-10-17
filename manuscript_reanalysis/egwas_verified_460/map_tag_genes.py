"""Map all tag windows without selecting by gene function. 1-based inclusive coordinates."""
from pathlib import Path
import csv,gzip,json,hashlib
from collections import defaultdict
r=Path(__file__).resolve().parent; out=r/'post_gwas/genes'; out.mkdir(exist_ok=True)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ref=r/'source/reference'; ann={x['locusName']:x for x in read(ref/'Vunguiculata_540_v1.2.annotation_info.txt')}
length={z[0]:int(z[1]) for z in (l.split() for l in (ref/'Vunguiculata_540_v1.0.fa.fai').read_text().splitlines())}
genes=defaultdict(list)
for line in (ref/'Vunguiculata_540_v1.2.gene.gff3').open():
 if line.startswith('#'):continue
 z=line.rstrip().split('\t')
 if z[2]!='gene':continue
 attr=dict(x.split('=',1) for x in z[8].split(';') if '=' in x)
 genes[z[0]].append((int(z[3]),int(z[4]),z[6],attr['Name']))
for gs in genes.values():gs.sort()
bim=r.parent/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim'
loc={z[1]:(z[0],int(z[3])) for z in (l.split() for l in bim.read_text().splitlines())}
assert all(c in length and 1<=p<=length[c] for c,p in loc.values())
tags=read(r/'post_gwas/tags/per_trait/tag_SNPs_by_trait.txt')
strict={x['tag_SNP'] for x in read(r/'post_gwas/tags/complete_scan/tag_SNPs_by_trait.txt')}
bytag=defaultdict(list)
for t in tags:bytag[t['tag_SNP']].append(t)
ld=defaultdict(dict)
with gzip.open(r/'post_gwas/tags/cowpea460_candidate_local_LD.ld.gz','rt') as f:
 fields=next(f).split()
 for line in f:
  z=dict(zip(fields,line.split()));a=z['SNP_A'];b=z['SNP_B'];v=float(z['R2'])
  if a in bytag:ld[a][b]=v
rows=[]; summary=[]
for tag,ts in sorted(bytag.items()):
 c,p=loc[tag];ld[tag][tag]=1.;linked={s:v for s,v in ld[tag].items() if v>.2 and abs(loc[s][1]-p)<1000000}
 positions=[loc[s][1] for s in linked];lo,hi=min(positions),max(positions)
 members=set(s for t in ts for s in t['members'].split(','))
 window=[g for g in genes[c] if g[0]<=min(length[c],p+999999) and g[1]>=max(1,p-999999)]
 for a,b,strand,g in window:
  typed=[s for s in ld[tag] if a<=loc[s][1]<=b]
  linked_inside=[s for s in typed if s in linked]
  sig_inside=[s for s in members if a<=loc[s][1]<=b]
  info=ann.get(g,{})
  rows.append(dict(tag_SNP=tag,chromosome=c,tag_position=p,traits=','.join(t['trait'] for t in ts),complete_scan_tag=tag in strict,gene=g,start=a,end=b,strand=strand,distance_to_tag_bp=max(a-p,p-b,0),tag_overlaps_gene=a<=p<=b,significant_member_overlaps=','.join(sorted(sig_inside)),LD_linked_typed_markers_in_gene=','.join(sorted(linked_inside)),max_r2_typed_marker_in_gene=max((ld[tag][s] for s in typed),default=''),gene_overlaps_LD_marker_span=a<=hi and b>=lo,arabi_best_hit=info.get('Best-hit-arabi-name',''),arabi_description=info.get('Best-hit-arabi-defline',''),rice_best_hit=info.get('Best-hit-rice-name',''),rice_description=info.get('Best-hit-rice-defline',''),GO=info.get('GO','')))
 summary.append(dict(tag_SNP=tag,chromosome=c,position=p,traits=','.join(t['trait'] for t in ts),complete_scan_tag=tag in strict,min_p_lrt=min(float(t['p_lrt']) for t in ts),minor_allele_carriers=ts[0]['minor_allele_carriers'],local_typed_markers=len(ld[tag]),linked_typed_markers=len(linked),LD_marker_span_start=lo,LD_marker_span_end=hi,genes_in_fixed_window=len(window)))
write('all_genes_in_tag_windows.txt',rows);write('tag_region_summary.txt',summary)
write('genes_with_direct_marker_support.txt',[x for x in rows if x['tag_overlaps_gene'] or x['significant_member_overlaps'] or x['LD_linked_typed_markers_in_gene']])
meta={'status':'positional mapping complete; best hits are not verified orthologs','n_tags':len(bytag),'n_window_gene_pairs':len(rows),'n_distinct_window_genes':len({x['gene'] for x in rows}),'coordinate_system':'1-based inclusive','assembly':'Vunguiculata v1.0; annotation v1.2; historical VEP script references same assembly; all BIM coordinates within indexed sequence lengths. Verified genotype rebuild passed sequence-level allele/reference checks; see genotype_rebuild_460_20261014/audit/recovered_annotation_reference_check.json.','window':'tag +/-999999 bp, clipped to chromosome; all overlapping genes retained independently of annotation','LD_span':'descriptive envelope of typed markers with r2 >0.2, not a confidence interval or haplotype block','orthology':'Phytozome best hits only; functional transfer remains provisional','reference_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ref.iterdir() if p.is_file()}}
(out/'mapping_protocol.json').write_text(json.dumps(meta,indent=2)+'\n'); print(json.dumps(meta,indent=2))
