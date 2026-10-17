"""Fiscus-informed, per-trait lead-SNP grouping, with raw-genotype LD checks.

Uses r2 > 0.20 and distance <1 Mb. Greedily select smallest-P unassigned SNP
and its directly linked significant neighbours; no transitive chain merging.
This deterministic operational rule is explicit, not a claim that the paper
specified the exact same algorithm. Singletons and equal-P ties are retained.
"""
from pathlib import Path
from collections import defaultdict
import csv, gzip, hashlib, json, math, sys

r=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
out=r/'post_gwas/tags'
geno=r.parent/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001'
if not geno.with_suffix('.bed').exists():
    geno=r.parent/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001'
def read_tsv(p):
    with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
bim=[line.split() for line in geno.with_suffix('.bim').read_text().splitlines()]
idx={row[1]:i for i,row in enumerate(bim)}
loc={row[1]:(row[0],int(row[3])) for row in bim}
fam=[line.split() for line in geno.with_suffix('.fam').read_text().splitlines()]
assert len(fam)==460 and all(row[2:4]==['0','0'] for row in fam)
assert len(bim)==len(idx)==41654
raw=geno.with_suffix('.bed').read_bytes()
stride=(len(fam)+3)//4
assert raw[:3]==b'\x6c\x1b\x01' and len(raw)==3+stride*len(bim)
ids=(out/'per_trait_unique_snps.txt').read_text().splitlines()
assert len(ids)==len(set(ids)) and all(s in idx for s in ids)
G={s:[(2,None,1,0)[(raw[3+idx[s]*stride+j//4] >> (2*(j%4)))&3] for j in range(460)] for s in ids}
def calc_r2(a,b):
    xy=[(x,y) for x,y in zip(G[a],G[b]) if x is not None and y is not None]
    n=len(xy); sx=sum(x for x,y in xy); sy=sum(y for x,y in xy)
    vx=n*sum(x*x for x,y in xy)-sx*sx
    vy=n*sum(y*y for x,y in xy)-sy*sy
    cov=n*sum(x*y for x,y in xy)-sx*sy
    return (cov*cov/(vx*vy) if vx>0 and vy>0 else None),n
ldfile=out/'cowpea460_candidate_local_LD.ld.gz'
reported={}; focal_seen=set(); ld_rows=0; nonsignificant_neighbours=0
with gzip.open(ldfile,'rt') as f:
    fields=next(f).split()
    for line in f:
        z=dict(zip(fields,line.split())); a,b=z['SNP_A'],z['SNP_B']
        assert a in G and b in loc
        assert loc[a]==(z['CHR_A'],int(z['BP_A'])) and loc[b]==(z['CHR_B'],int(z['BP_B']))
        assert loc[a][0]==loc[b][0] and abs(loc[a][1]-loc[b][1])<=1000000
        focal_seen.add(a); ld_rows+=1
        if b in G and a!=b:
            key=tuple(sorted([a,b])); v=float(z['R2'])
            if key in reported: assert abs(reported[key]-v)<1e-12
            reported[key]=v
        elif b not in G: nonsignificant_neighbours+=1
pair_ld={}; pair_rows=[]; errors=[]
for i,a in enumerate(ids):
    for b in ids[i+1:]:
        if loc[a][0]!=loc[b][0] or abs(loc[a][1]-loc[b][1])>=1000000: continue
        key=tuple(sorted([a,b])); v,n=calc_r2(a,b)
        assert key in reported, f'Missing local candidate LD pair: {key}'
        if v is None:
            assert math.isnan(reported[key]); continue
        assert math.isfinite(reported[key]) and abs(v-reported[key])<=5.1e-7
        errors.append(abs(v-reported[key])); pair_ld[key]=v
        pair_rows.append(dict(SNP_A=a,SNP_B=b,chromosome=loc[a][0],distance_bp=abs(loc[a][1]-loc[b][1]),r2=v,jointly_called_plants=n,passes_r2_gt_0p2=v>.2))
write_tsv(out/'significant_SNP_pair_LD.txt',pair_rows,['SNP_A','SNP_B','chromosome','distance_bp','r2','jointly_called_plants','passes_r2_gt_0p2'])
def linked(a,b): return a==b or pair_ld.get(tuple(sorted([a,b])),-1)>.2
traits=(r/'input/trait_names.txt').read_text().splitlines()
summary={}
for level in ['per_trait','complete_scan']:
    hits=read_tsv(r/'post_gwas/parsed'/level/'cowpea460_gwas_hitlst.txt')
    assert len({(h['ID'],h['rs']) for h in hits})==len(hits)
    folder=out/level; folder.mkdir(exist_ok=True)
    tagrows=[]; members=[]; ties=[]
    for trait in traits:
        remaining={h['rs']:h for h in hits if h['ID']==trait}
        assigned=set(); count=0
        while remaining:
            lead=min(remaining,key=lambda s:(float(remaining[s]['p_lrt']),loc[s][0],loc[s][1],s))
            h=remaining[lead]; count+=1; group=f'{trait}:peak{count:03d}'
            group_ids=sorted([s for s in remaining if linked(lead,s)],key=lambda s:(loc[s][1],s))
            equal=[s for s in group_ids if float(remaining[s]['p_lrt'])==float(h['p_lrt'])]
            pmin=float(h['p_lrt'])
            for s in group_ids:
                assert s not in assigned and linked(lead,s)
                assert loc[s][0]==loc[lead][0] and abs(loc[s][1]-loc[lead][1])<1000000
                assert float(remaining[s]['p_lrt'])>=pmin
                assigned.add(s)
                members.append(dict(trait=trait,peak=group,tag_SNP=lead,SNP=s,chromosome=loc[s][0],position=loc[s][1],p_lrt=remaining[s]['p_lrt'],r2_with_tag=1.0 if s==lead else pair_ld[tuple(sorted([lead,s]))],distance_from_tag_bp=abs(loc[s][1]-loc[lead][1]),complete_scan_significant=remaining[s]['complete_scan_significant']))
            called=[v for v in G[lead] if v is not None]
            a1=sum(called)/(2*len(called)); carriers=sum(v>0 if a1<=.5 else v<2 for v in called)
            row=dict(trait=trait,peak=group,tag_SNP=lead,chromosome=loc[lead][0],position=loc[lead][1],p_lrt=h['p_lrt'],allele1=h['allele1'],allele0=h['allele0'],GEMMA_af=h['af'],beta=h['beta'],se=h['se'],n_significant_members=len(group_ids),significant_span_start=min(loc[s][1] for s in group_ids),significant_span_end=max(loc[s][1] for s in group_ids),members=','.join(group_ids),equal_P_members=','.join(equal),complete_scan_significant=h['complete_scan_significant'],minor_allele_carriers=carriers,carrier_flag='<10 carriers' if carriers<10 else 'No low-count flag')
            tagrows.append(row)
            for s in group_ids: del remaining[s]
        assert len(assigned)==sum(h['ID']==trait for h in hits)
    assert {(z['trait'],z['SNP']) for z in members}=={(h['ID'],h['rs']) for h in hits}
    assert len(members)==len(hits)
    for i,a in enumerate(tagrows):
        for b in tagrows[i+1:]:
            if a['trait']==b['trait']: assert not linked(a['tag_SNP'],b['tag_SNP'])
    unique=sorted({z['tag_SNP'] for z in tagrows})
    write_tsv(folder/'tag_SNPs_by_trait.txt',tagrows,list(tagrows[0]))
    write_tsv(folder/'tag_membership.txt',members,list(members[0]))
    (folder/'unique_tag_snps.txt').write_text('\n'.join(unique)+'\n')
    totals=[dict(trait=t,n_significant_SNPs=sum(h['ID']==t for h in hits),n_tags=sum(z['trait']==t for z in tagrows)) for t in traits]
    write_tsv(folder/'counts_by_trait.txt',totals,['trait','n_significant_SNPs','n_tags'])
    summary[level]={'input_associations':len(hits),'input_unique_SNPs':len({h['rs'] for h in hits}),'trait_specific_tags':len(tagrows),'unique_tag_SNPs':len(unique),'singleton_trait_specific_tags':sum(z['n_significant_members']==1 for z in tagrows),'tags_with_equal_P_ties':sum(',' in z['equal_P_members'] for z in tagrows),'unique_tags_under_10_carriers':len({z['tag_SNP'] for z in tagrows if z['minor_allele_carriers']<10})}
protocol={'reference':'Fiscus et al. 2024, G3, doi:10.1093/g3journal/jkae071','reference_url':'https://doi.org/10.1093/g3journal/jkae071','threshold_r2':0.2,'comparison':'strict >','max_distance_bp':1000000,'distance_comparison':'strict <','same_chromosome':True,'method':'Per trait and significance set, greedily choose minimum-P unassigned significant SNP and assign directly linked unassigned significant neighbours; retain singletons. No transitive merging.','ties':'Position then SNP ID; all equal-P members recorded.','r2_estimate':'Squared genotype allele-count correlation on jointly called plants; no imputation. Recomputed from BED for significant-marker pairs and checked against PLINK.','interpretation':'Operational peak representatives, not proven causal variants or globally independent loci. Significant-span bounds are not gene-search or haplotype-block boundaries.','decision_basis':'User requested tag selection, asked for best rule and directed review of Fiscus cowpea paper. Adopted published cutoff and lead-P principle; explicit greedy implementation as recommended in clarification.','results':summary}
(out/'tag_selection_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
audit={'status':'PASS','LD_rows':ld_rows,'focal_SNPs_in_LD_report':len(focal_seen),'candidate_pairs_verified_from_BED':len(pair_rows),'max_r2_absolute_difference':max(errors,default=0),'all_hit_rows_assigned_once':True,'singleton_retention':True,'within_trait_tags_not_linked_under_adopted_rule':True,'LD_sha256':hashlib.sha256(ldfile.read_bytes()).hexdigest(),'genotype_sha256':{ext:hashlib.sha256(geno.with_suffix(ext).read_bytes()).hexdigest() for ext in ['.bed','.bim','.fam']},'results':summary}
(r/'audit/tag_selection_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
