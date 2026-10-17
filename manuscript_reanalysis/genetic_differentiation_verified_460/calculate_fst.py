#!/usr/bin/env python3
"""Array-marker WC84 FST, observed/expected heterozygosity, chromosome bootstrap.
Run from project root. Source genotypes and group memberships are not modified.
"""
from pathlib import Path
import csv, json, hashlib, itertools
import numpy as np
OUT=Path('manuscript_reanalysis/genetic_differentiation_verified_460')
PREFIX=OUT/'input/cowpea460_verified_geno005_maf001'

def write_csv(path, rows):
    rows=list(rows)
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def read_inputs():
    fam=[x.split() for x in PREFIX.with_suffix('.fam').read_text().splitlines()]
    bim=[x.split() for x in PREFIX.with_suffix('.bim').read_text().splitlines()]
    raw=np.fromfile(PREFIX.with_suffix('.bed'),dtype=np.uint8)
    assert raw[:3].tolist()==[108,27,1]
    n=len(fam); m=len(bim); stride=(n+3)//4
    assert len(raw)==3+m*stride
    packed=raw[3:].reshape(m,stride)
    codes=((packed[:,:,None] >> (2*np.arange(4))) & 3).reshape(m,stride*4)[:,:n]
    g=np.array([0,-1,1,2],dtype=np.int8)[codes] # A2 dosage; -1 missing
    q=list(csv.DictReader(open('manuscript_reanalysis/genetic_structure_verified_460/admixture/summary/K3_best_Q_assignments.csv')))
    q={x['Accession']:x for x in q}
    ids=[x[1] for x in fam]
    assert len(set(ids))==n==460 and set(ids)==set(q)
    assigned=np.array([q[x]['Assignment'] for x in ids])
    for x in ids:
        row=q[x];mx=max(float(row[f'Cluster{i}']) for i in [1,2,3])
        assert abs(mx-float(row['maxQ']))<1e-10
        assert (row['Assignment']=='Admixed') == (mx<.70)
    pruned=set((OUT/'input/structure_LD02.prune.in').read_text().split())
    masks={'all_QC_SNPs':np.ones(m,dtype=bool),'LD_pruned_SNPs':np.array([x[1] in pruned for x in bim])}
    assert m==41654 and masks['LD_pruned_SNPs'].sum()==len(pruned)==4194
    assert np.all((g>=0).mean(axis=1)>=.95)
    return g,ids,np.array([x[0] for x in bim]),bim,assigned,masks

def group_counts(g,idx):
    z=g[:,idx];n=(z>=0).sum(axis=1).astype(float)
    p=np.maximum(z,0).sum(axis=1)/(2*n)
    h=(z==1).sum(axis=1)/n
    return n,p,h

def wc_components(counts):
    # WC84 components for one of the two alleles; doubling cancels in FST.
    n,p,h=[np.column_stack(x) for x in zip(*counts)]
    r=n.shape[1];nt=n.sum(axis=1);nb=n.mean(axis=1)
    nc=(nt-(n*n).sum(axis=1)/nt)/(r-1)
    pb=(n*p).sum(axis=1)/nt
    s2=(n*(p-pb[:,None])**2).sum(axis=1)/(nb*(r-1))
    hb=(n*h).sum(axis=1)/nt
    a=nb/nc*(s2-(pb*(1-pb)-(r-1)/r*s2-hb/4)/(nb-1))
    b=nb/(nb-1)*(pb*(1-pb)-(r-1)/r*s2-(2*nb-1)/(4*nb)*hb)
    c=hb/2
    assert np.isfinite(a+b+c).all()
    return a,b,c

def main():
    g,ids,chrom,bim,assigned,masks=read_inputs()
    groups=['Cluster1','Cluster2','Cluster3']; idx={k:np.flatnonzero(assigned==k) for k in groups}
    assert [len(idx[k]) for k in groups]==[105,77,192]
    for k in groups:(OUT/f'{k}_ids.txt').write_text(''.join(ids[i]+'\n' for i in idx[k]))
    write_csv(OUT/'sample_membership.csv',[{'Accession':x,'Assignment':assigned[i],'included_cluster_FST':assigned[i]!='Admixed'} for i,x in enumerate(ids)])
    counts={k:group_counts(g,idx[k]) for k in groups}
    chromosomes=np.unique(chrom); assert len(chromosomes)==11
    rng=np.random.default_rng(460384)
    draws=rng.integers(0,len(chromosomes),size=(10000,len(chromosomes)))
    summary=[]; blocks=[]; loo=[]; diversity=[]; markerrows=[]
    for marker_set,mask in masks.items():
        for k in groups:
            n,p,h=counts[k];he=2*p*(1-p)
            diversity.append({'marker_set':marker_set,'cluster':k,'n_plants':len(idx[k]),'n_markers':int(mask.sum()),'mean_Ho':float(h[mask].mean()),'mean_He':float(he[mask].mean()),'mean_call_rate':float(n[mask].mean()/len(idx[k])),'polymorphic_markers':int(((p>0)&(p<1)&mask).sum())})
        for ks in list(itertools.combinations(groups,2))+[tuple(groups)]:
            name='_vs_'.join(ks) if len(ks)==2 else 'Global_3_clusters'
            a,b,c=wc_components([counts[k] for k in ks]);den=a+b+c
            A=np.array([a[mask&(chrom==ch)].sum() for ch in chromosomes])
            D=np.array([den[mask&(chrom==ch)].sum() for ch in chromosomes])
            fst=A.sum()/D.sum()
            boot=A[draws].sum(axis=1)/D[draws].sum(axis=1)
            ci=np.quantile(boot,[.025,.975])
            summary.append({'marker_set':marker_set,'comparison':name,'group1':ks[0],'group2':ks[1] if len(ks)==2 else 'all_three','n1':len(idx[ks[0]]) if len(ks)==2 else '', 'n2':len(idx[ks[1]]) if len(ks)==2 else '','markers_input':int(mask.sum()),'markers_positive_denominator':int((mask&(den>0)).sum()),'FST':float(fst),'CI_lower':float(ci[0]),'CI_upper':float(ci[1]),'bootstrap_replicates':len(draws),'chromosomes':len(chromosomes)})
            for j,ch in enumerate(chromosomes):
                blocks.append({'marker_set':marker_set,'comparison':name,'chromosome':ch,'n_markers':int((mask&(chrom==ch)).sum()),'numerator':float(A[j]),'denominator':float(D[j]),'FST':float(A[j]/D[j])})
                loo.append({'marker_set':marker_set,'comparison':name,'excluded_chromosome':ch,'FST':float((A.sum()-A[j])/(D.sum()-D[j]))})
            if marker_set=='all_QC_SNPs' and len(ks)==2:
                for i in range(len(bim)):
                    markerrows.append({'comparison':name,'chromosome':chrom[i],'position':int(bim[i][3]),'SNP':bim[i][1],'a':a[i],'b':b[i],'c':c[i],'FST':a[i]/den[i] if den[i]>0 else ''})
    write_csv(OUT/'FST_summary.csv',summary);write_csv(OUT/'marker_diversity.csv',diversity)
    write_csv(OUT/'chromosome_components.csv',blocks);write_csv(OUT/'leave_one_chromosome_out.csv',loo)
    write_csv(OUT/'WC84_marker_components.csv',markerrows)
    hashes={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted((OUT/'input').iterdir()) if x.is_file()}
    qpath=Path('manuscript_reanalysis/genetic_structure_verified_460/admixture/summary/K3_best_Q_assignments.csv')
    hashes[str(qpath)]=hashlib.sha256(qpath.read_bytes()).hexdigest()
    manifest={'cohort_n':460,'assigned_n':374,'admixed_excluded_from_cluster_FST':86,'groups':{k:len(idx[k]) for k in groups},'assignment_threshold':.70,'primary_SNPs':41654,'LD_pruned_sensitivity_SNPs':4194,'estimator':'Weir and Cockerham 1984 theta; ratio of sums of a/(a+b+c), not mean of per-SNP ratios','genotype_handling':'Diploid observed calls, no imputation, locus-specific sample sizes; preserve negative variance components','CI':'95% percentile chromosome-block bootstrap; 10000 draws of 11 entire chromosomes with replacement; fixed sampled plants and memberships','bootstrap_seed':460384,'diversity':'Mean per-marker Ho and He=2p(1-p) over common QC SNP set, including cluster-monomorphic markers; marker diversity, not genome-wide per-base pi','limitations':['Array ascertainment and MAF filtering affect these summaries.','Clusters selected from the same genotype data; FST describes this partition rather than independently validating discrete populations or adaptation.','Chromosome bootstrap measures genomic variation conditional on sampled accessions and memberships.','LD pruning changes marker weighting; sensitivity comparison is not a replacement estimator.'],'input_sha256':hashes,'numpy_version':np.__version__}
    (OUT/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(json.dumps(diversity,indent=2))

if __name__=='__main__':main()
