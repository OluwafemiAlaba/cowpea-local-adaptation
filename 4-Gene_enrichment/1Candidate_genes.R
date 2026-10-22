#!/usr/bin/env Rscript
# Pool pRDA and per-trait eGWAS genes within +/-20 kb; apply the same window to background SNPs.
source('config.R')
bim<-read.table(file.path(inputs,'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim'),
                col.names=c('chr','SNP','cm','position','A1','A2'))
p<-read.csv(result('rda','results/primary_adjusted_putative_outliers_3SD.csv'))
e<-read.delim(result('egwas','post_gwas/parsed/per_trait/cowpea460_gwas_hitlst.txt'))
stopifnot(nrow(bim)==41654,!anyDuplicated(bim$SNP),all(p$SNP %in% bim$SNP),all(e$rs %in% bim$SNP))
for (d in list(data.frame(SNP=p$SNP,chr=p$PLINK_chromosome_code,position=p$position),
               data.frame(SNP=e$rs,chr=e$chr,position=e$ps))) {
  b<-bim[match(d$SNP,bim$SNP),]
  stopifnot(all(as.character(d$chr)==as.character(b$chr)),all(d$position==b$position))
}
gff<-read.delim(input('egwas','source/reference/Vunguiculata_540_v1.2.gene.gff3'),
               header=FALSE,comment.char='#',quote='',fill=TRUE)
gff<-gff[gff$V3=='gene',]
genes<-data.frame(gene=sub('.*(?:^|;)Name=([^;]+).*','\\1',gff$V9,perl=TRUE),
                  chr=gff$V1,start=gff$V4,end=gff$V5)
stopifnot(!anyDuplicated(genes$gene))
window<-20000L; foreground<-unique(c(p$SNP,e$rs))
mapping<-lapply(seq_len(nrow(bim)),function(i) {
  b<-bim[i,]; genes$gene[genes$chr==b$chr & genes$start<=b$position+window & genes$end>=max(1,b$position-window)]
})
names(mapping)<-bim$SNP
bg<-sort(unique(unlist(mapping))); fg<-sort(unique(unlist(mapping[foreground])))
stopifnot(all(fg %in% bg))
writeLines(fg,output('go','inputs/without_simple_candidate_genes.txt'))
writeLines(bg,output('go','inputs/without_simple_background_genes.txt'))
rows<-lapply(foreground,function(s) data.frame(SNP=rep(s,length(mapping[[s]])),gene=mapping[[s]],
                    pRDA=rep(s %in% p$SNP,length(mapping[[s]])),eGWAS=rep(s %in% e$rs,length(mapping[[s]]))))
write.csv(do.call(rbind,rows),output('go','candidate_gene_links.csv'),row.names=FALSE)
cat(length(fg),'candidate genes;',length(bg),'background genes\n')
