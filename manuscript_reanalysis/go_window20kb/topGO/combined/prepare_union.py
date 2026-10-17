from pathlib import Path
import csv,json,hashlib
root=Path(__file__).resolve().parents[2]
out=Path(__file__).resolve().parent
# Reuse the independently checked coordinate mapping without running historical writes.
f=root/'map_candidates.py'
ns={'__file__':str(f)}
exec(compile(f.read_text().split("sets = {'pRDA_3SD'")[0],str(f),'exec'),ns)
rows=ns['read'](root.parent/'rda_verified_460/results/climate_only_putative_outliers_3SD.csv')
simple={r['SNP'] for r in rows}
assert simple<=ns['universe_rda']
for r in rows: assert ns['loc'][r['SNP']]==(r['PLINK_chromosome_code'],int(r['position']))
sets={'simple_RDA':simple,'pRDA':ns['prda'],'eGWAS':ns['per']}
genes={k:{g for s in v for g in ns['mapping'][s]} for k,v in sets.items()}
fg=set.union(*genes.values());bg={g for s in ns['loc'] for g in ns['mapping'][s]}
assert fg<=bg
for name,v in [('candidate_genes',fg),('background_genes',bg)]:
 (out/'inputs'/f'combined_{name}.txt').write_text('\n'.join(sorted(v))+'\n')
with (out/'gene_provenance.csv').open('w') as h:
 w=csv.writer(h);w.writerow(['gene','simple_RDA','pRDA','eGWAS','candidate_SNPs'])
 for g in sorted(fg):
  snps=sorted(s for s in set.union(*sets.values()) if g in ns['mapping'][s])
  w.writerow([g,*[g in genes[k] for k in sets],';'.join(snps)])
summary={'SNPs_by_method':{k:len(v) for k,v in sets.items()},'genes_by_method':{k:len(v) for k,v in genes.items()},'union_genes':len(fg),'union_SNPs':len(set.union(*sets.values())),'background_genes':len(bg),'additional_simple_RDA_genes':len(genes['simple_RDA']-(genes['pRDA']|genes['eGWAS'])),'window_bp':20000,'background_SNPs':len(ns['loc'])}
(out/'union_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
