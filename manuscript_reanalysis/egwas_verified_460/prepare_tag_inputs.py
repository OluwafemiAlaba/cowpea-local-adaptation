"""Extract distinct significant SNP IDs, preserving both significance sets."""
from pathlib import Path
import csv, hashlib, json
r=Path(__file__).resolve().parent
out=r/'post_gwas/tags'
out.mkdir(exist_ok=True)
geno=r.parent/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001'
bim={x[1]:(x[0],int(x[3])) for x in (line.split() for line in geno.with_suffix('.bim').read_text().splitlines())}
fam=[line.split() for line in geno.with_suffix('.fam').read_text().splitlines()]
assert len(bim)==41654 and len(fam)==460
assert all(row[2:4]==['0','0'] for row in fam), 'PLINK LD uses founders; investigate nonfounder entries'
sets={}; sources={}
for level in ['per_trait','complete_scan']:
    p=r/'post_gwas/parsed'/level/'cowpea460_gwas_hitlst.txt'
    with p.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    ids=sorted({row['rs'] for row in rows})
    assert all(bim[row['rs']]==(row['chr'],int(row['ps'])) for row in rows)
    original=(p.parent/'cowpea460_gwas_hitlst_uniq_snps.txt').read_text().splitlines()
    assert ids==original
    (out/f'{level}_unique_snps.txt').write_text('\n'.join(ids)+'\n')
    sets[level]=ids
    sources[level]={'hitlist_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                    'associations':len(rows),'unique_SNPs':len(ids)}
assert set(sets['complete_scan'])<=set(sets['per_trait'])
report={'status':'PASS','n_plants':460,'all_founders':True,'n_QC_SNPs':len(bim),
 'sets':sources,'purpose':'Distinct IDs only; these are not yet LD-selected tags',
 'LD_preparation':'Measure same-chromosome LD within 1 Mb with all current QC SNPs; keep all r2 values. Adopted tag rule: r2 > 0.20, distance <1 Mb, per-trait minimum-P leads and directly linked neighbours.'}
(r/'audit/tag_inputs.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
