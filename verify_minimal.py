#!/usr/bin/env python3
"""Read-only checks; does not run models or install packages."""
from pathlib import Path
import csv,hashlib,json,math
root=Path(__file__).resolve().parent
data=root if (root/'manuscript_reanalysis/environment_460.csv').exists() else root.parent/'dryad'
with (root/'file_manifest.csv').open() as h:manifest=list(csv.DictReader(h))
for row in manifest:
 f=(root if row['package']=='github' else data)/row['relative_path']
 assert f.is_file(),str(f)
 assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'],str(f)
p=data/'manuscript_reanalysis'
with (p/'environment_460.csv').open() as h:env=list(csv.DictReader(h))
ids={r['accn'] for r in env};assert len(env)==len(ids)==460
g=p/'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001'
fam=[s.split() for s in g.with_suffix('.fam').read_text().splitlines()]
assert {r[1] for r in fam}==ids
assert len(g.with_suffix('.bim').read_text().splitlines())==41654
assert g.with_suffix('.bed').stat().st_size==3+115*41654
traits=(p/'egwas_verified_460/input/trait_names.txt').read_text().splitlines();assert len(set(traits))==len(traits)==92
assert all(math.isfinite(float(r[t])) for r in env for t in traits)
with (p/'egwas_verified_460/input/sample_order.tsv').open() as h:order=list(csv.DictReader(h,delimiter='\t'))
assert [r['IID'] for r in order]==[r[1] for r in fam]
print(json.dumps({'manifest_files_verified':len(manifest),'accessions':460,'SNPs':41654,'climate_variables':92,'sample_order':'verified by ID','full_model_rerun':False,'status':'PASS'},indent=2))
