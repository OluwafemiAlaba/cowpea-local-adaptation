#!/usr/bin/env Rscript
# Within each trait, choose the lowest-P lead and its directly linked significant SNPs.
# Same chromosome, distance <1 Mb and r2 >0.20; no transitive merging.
source('config.R')
Sys.setlocale('LC_COLLATE','C')
h <- read.delim(result('egwas','post_gwas/parsed/per_trait/cowpea460_gwas_hitlst.txt'))
bim <- read.table(file.path(inputs,'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim'),
                  col.names=c('chr','SNP','cm','position','A1','A2'))
fam <- read.table(file.path(inputs,'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.fam'))
stopifnot(nrow(fam)==460, all(fam[,3:4]==0), !anyDuplicated(bim$SNP))
ids <- sort(unique(h$rs)); b <- bim[match(ids,bim$SNP),]
raw <- read.table(result('egwas','post_gwas/tags/candidate_dosages.raw'),header=TRUE,check.names=FALSE)
stopifnot(identical(as.character(raw$IID),as.character(fam$V2)))
g <- as.matrix(raw[,paste(ids,b$A1,sep='_'),drop=FALSE]); colnames(g)<-ids
stopifnot(all(g[!is.na(g)] %in% 0:2))
r2 <- cor(g,use='pairwise.complete.obs')^2
ld <- read.table(gzfile(result('egwas','post_gwas/tags/cowpea460_candidate_local_LD.ld.gz')),header=TRUE)
z <- ld[ld$SNP_A %in% ids & ld$SNP_B %in% ids & ld$SNP_A!=ld$SNP_B,]
computed <- r2[cbind(match(z$SNP_A,ids),match(z$SNP_B,ids))]
stopifnot(all(is.na(computed)==is.na(z$R2)),max(abs(computed-z$R2),na.rm=TRUE)<5.1e-7)
near <- outer(b$chr,b$chr,'==') & abs(outer(b$position,b$position,'-'))<1e6
linked <- near & is.finite(r2) & r2>0.20; diag(linked)<-TRUE
dimnames(linked)<-list(ids,ids)
traits <- readLines(input('egwas','input/trait_names.txt'))
for (level in c('per_trait','complete_scan')) {
  hit <- if(level=='per_trait') h else h[h$complete_scan_significant,]
  tags <- members <- list()
  for (trait in traits) {
    remaining <- hit[hit$ID==trait,]; peak<-0L
    while(nrow(remaining)) {
      remaining<-remaining[order(remaining$p_lrt,remaining$chr,remaining$ps,remaining$rs),]
      lead<-remaining[1,]; peak<-peak+1L
      take<-linked[lead$rs,remaining$rs]
      block<-remaining[take,]; name<-sprintf('%s:peak%03d',trait,peak)
      tags[[length(tags)+1L]]<-data.frame(trait=trait,peak=name,tag_SNP=lead$rs,
        chromosome=lead$chr,position=lead$ps,p_lrt=lead$p_lrt,allele1=lead$allele1,
        allele0=lead$allele0,beta=lead$beta,se=lead$se,n_significant_members=nrow(block))
      members[[length(members)+1L]]<-data.frame(trait=trait,peak=name,tag_SNP=lead$rs,
        SNP=block$rs,r2_with_tag=r2[lead$rs,block$rs],distance_from_tag_bp=abs(block$ps-lead$ps))
      remaining<-remaining[!take,]
    }
  }
  t<-do.call(rbind,tags); m<-do.call(rbind,members)
  stopifnot(nrow(m)==nrow(hit), !anyDuplicated(m[c('trait','SNP')]))
  write_tsv(t,output('egwas',paste0('post_gwas/tags/',level,'/tag_SNPs_by_trait.txt')))
  write_tsv(m,output('egwas',paste0('post_gwas/tags/',level,'/tag_membership.txt')))
  writeLines(sort(unique(t$tag_SNP)),output('egwas',paste0('post_gwas/tags/',level,'/unique_tag_snps.txt')))
}
