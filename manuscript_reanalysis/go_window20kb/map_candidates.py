"""Map current pRDA/eGWAS SNPs and their tested universes to +/-20 kb genes.

Run from any directory with Python 3. Outputs are analysis intermediates.
No LD expansion or GO enrichment tests are performed.
"""
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
WINDOW = 20000
BIM = ROOT / 'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim'
GFF = ROOT / 'egwas_verified_460/source/reference/Vunguiculata_540_v1.2.gene.gff3'
ANN = ROOT / 'egwas_verified_460/source/reference/Vunguiculata_540_v1.2.annotation_info.txt'
LOAD = ROOT / 'rda_verified_460/results/primary_adjusted_SNP_loadings_scaling2.csv'
RDA = ROOT / 'rda_verified_460/results/primary_adjusted_putative_outliers_3SD.csv'
PER = ROOT / 'egwas_verified_460/post_gwas/parsed/per_trait/cowpea460_gwas_hitlst.txt'
STRICT = ROOT / 'egwas_verified_460/post_gwas/parsed/complete_scan/cowpea460_gwas_hitlst.txt'

def read(path, delimiter=','):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))

loc = {}
for line in BIM.read_text().splitlines():
    z = line.split()
    assert z[1] not in loc
    loc[z[1]] = (z[0], int(z[3]))
assert len(loc) == 41654
universe_rda = {r['SNP'] for r in read(LOAD)}
rda_rows = read(RDA)
prda = {r['SNP'] for r in rda_rows}
per_rows, strict_rows = read(PER, '\t'), read(STRICT, '\t')
per, strict = {r['rs'] for r in per_rows}, {r['rs'] for r in strict_rows}
assert len(universe_rda) == 20027 and len(prda) == 763
assert len(per) == 98 and len(strict) == 17
assert prda <= universe_rda <= loc.keys() and strict <= per <= loc.keys()
for r in rda_rows:
    assert loc[r['SNP']] == (r['PLINK_chromosome_code'], int(r['position']))
for r in per_rows + strict_rows:
    assert loc[r['rs']] == (r['chr'], int(r['ps']))

genes, by_chr = {}, defaultdict(list)
for line in GFF.open():
    if line.startswith('#'):
        continue
    z = line.rstrip().split('\t')
    if z[2] != 'gene':
        continue
    attrs = dict(a.split('=', 1) for a in z[8].split(';') if '=' in a)
    gene = attrs['Name']
    assert gene not in genes
    genes[gene] = dict(gene=gene, chromosome=z[0], start=int(z[3]), end=int(z[4]), strand=z[6])
    by_chr[z[0]].append(gene)
assert set(c for c, p in loc.values()) <= by_chr.keys()
annotations = defaultdict(lambda: defaultdict(set))
for r in read(ANN, '\t'):
    a = annotations[r['locusName']]
    a['GO'].update(x for x in r['GO'].split() if x.startswith('GO:'))
    for key in ('Best-hit-arabi-name', 'Best-hit-arabi-defline'):
        if r.get(key):
            a[key].add(r[key])

# SNP-centred inclusive windows intersect whole gene bodies, including introns.
mapping = {}
for snp, (chrom, pos) in loc.items():
    lo, hi = max(1, pos-WINDOW), pos+WINDOW
    mapping[snp] = sorted(g for g in by_chr[chrom]
                          if genes[g]['start'] <= hi and genes[g]['end'] >= lo)

# Independent gene-centred binary-search reconstruction checks all mappings.
chr_snps = defaultdict(list)
for snp, (chrom, pos) in loc.items():
    chr_snps[chrom].append((pos, snp))
for arr in chr_snps.values():
    arr.sort()
inverse = defaultdict(set)
for g, info in genes.items():
    arr = chr_snps[info['chromosome']]
    positions = [p for p, s in arr]
    first = bisect_left(positions, info['start']-WINDOW)
    last = bisect_right(positions, info['end']+WINDOW)
    for pos, snp in arr[first:last]:
        inverse[snp].add(g)
assert all(set(mapping[s]) == inverse[s] for s in loc)

sets = {'pRDA_3SD': (prda, universe_rda),
        'eGWAS_per_trait': (per, set(loc)),
        'eGWAS_across92': (strict, set(loc))}
summaries, payload, evidence = {}, {}, []
for name, (selected, tested) in sets.items():
    foreground = {g for s in selected for g in mapping[s]}
    background = {g for s in tested for g in mapping[s]}
    assert foreground <= background
    payload[name] = dict(candidate_SNPs=sorted(selected), tested_SNPs=sorted(tested),
                         candidate_genes=sorted(foreground), background_genes=sorted(background))
    summaries[name] = dict(candidate_SNPs=len(selected), tested_SNPs=len(tested),
        SNPs_with_genes=sum(bool(mapping[s]) for s in selected),
        SNPs_without_genes=sum(not mapping[s] for s in selected),
        unique_genes=len(foreground), genes_with_supplied_GO=sum(bool(annotations[g]['GO']) for g in foreground),
        background_genes=len(background), background_with_supplied_GO=sum(bool(annotations[g]['GO']) for g in background))
    for s in sorted(selected):
        chrom, pos = loc[s]
        for g in mapping[s]:
            info = genes[g]
            distance = max(info['start']-pos, pos-info['end'], 0)
            assert 0 <= distance <= WINDOW
            evidence.append(dict(analysis=name, SNP=s, chromosome=chrom, SNP_position=pos,
                window_start=max(1,pos-WINDOW), window_end=pos+WINDOW, **{k:v for k,v in info.items() if k!='chromosome'},
                distance_to_gene_body_bp=distance,
                evidence='intragenic' if distance==0 else 'within_20kb',
                supplied_GO=sorted(annotations[g]['GO']),
                Arabidopsis_best_hits=sorted(annotations[g]['Best-hit-arabi-name']),
                Arabidopsis_descriptions=sorted(annotations[g]['Best-hit-arabi-defline'])))

protocol = dict(status='Candidate and background mapping complete; GO enrichment not run',
    window_bp_each_side=WINDOW, coordinate_system='1-based inclusive',
    mapping='Whole annotated gene body overlaps [max(1,SNP-20000), SNP+20000]; same chromosome; no LD expansion',
    annotation='Vunguiculata_540_v1.2; genes counted once per analysis; supplied GO union across transcripts',
    cautions=['Proximity is not evidence of causality or LD.', 'Arabidopsis best hits are descriptive, not experimentally verified cowpea functions.',
              'GO annotation counts are prior to ontology/aspect reconciliation.', 'Prior LD-based enrichment results do not apply to these new sets.'],
    counts=summaries, validation='All 41,654 SNP mappings reproduced by independent gene-centred interval search',
    sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BIM,GFF,ANN,LOAD,RDA,PER,STRICT]})
(OUT/'mapping_protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
(OUT/'candidate_gene_evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')
(OUT/'candidate_and_background_sets.json').write_text(json.dumps(payload, indent=2)+'\n')
print(json.dumps(protocol, indent=2))
