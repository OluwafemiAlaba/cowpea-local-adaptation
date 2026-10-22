#!/usr/bin/env Rscript
# Summarize 92 GEMMA scans using likelihood-ratio P values.
source('config.R')
traits <- readLines(input('egwas', 'input/trait_names.txt'))
stopifnot(length(traits)==92, !anyDuplicated(traits))
expected <- scan(file.path(inputs, 'genotype_rebuild_460_20261014/cowpea460_verified_geno005_maf001.bim'),
                 what=list('', '', 0, 0, '', ''), quiet=TRUE)[[2]]
hits <- low <- vector('list', length(traits)); pve <- numeric(length(traits))
for (i in seq_along(traits)) {
  x <- read.delim(result('egwas', paste0('results/',traits[i],'.assoc.txt')), check.names=FALSE)
  stopifnot(nrow(x)==41654, !anyDuplicated(x$rs), setequal(x$rs,expected),
            all(is.finite(x$p_lrt)), all(x$p_lrt>0 & x$p_lrt<=1))
  x <- cbind(ID=traits[i], x)
  low[[i]] <- x[x$p_lrt<1e-4,]
  hits[[i]] <- x[x$p_lrt<0.05/41654,]
  log <- readLines(result('egwas', paste0('results/',traits[i],'.log.txt')))
  p <- grep('^## pve estimate in the null model =',log,value=TRUE)
  stopifnot(length(p)==1)
  pve[i] <- as.numeric(sub('.*= *','',p))
}
h <- do.call(rbind,hits); h$snps_tested <- 41654L
h$per_trait_cutoff <- 0.05/41654
h$complete_scan_cutoff <- 0.05/(41654*92)
h$complete_scan_significant <- h$p_lrt<h$complete_scan_cutoff
h <- h[order(match(h$ID,traits),h$chr,h$ps,h$rs),]
keep <- c('ID','rs','chr','ps','allele1','allele0','af','beta','se','p_lrt',
          'snps_tested','per_trait_cutoff','complete_scan_cutoff','complete_scan_significant')
h <- h[,keep]
h$significance_label <- ifelse(h$complete_scan_significant,'complete_scan_Bonferroni','per_trait_Bonferroni_only')
for (level in c('per_trait','complete_scan')) {
  z <- if(level=='per_trait') h else h[h$complete_scan_significant,]
  base <- paste0('post_gwas/parsed/',level,'/cowpea460_gwas_')
  write_tsv(z,output('egwas',paste0(base,'hitlst.txt')))
  writeLines(sort(unique(z$rs)),output('egwas',paste0(base,'hitlst_uniq_snps.txt')))
  counts <- data.frame(ID=traits,n_SNPs=as.integer(table(factor(z$ID,levels=traits))))
  write_tsv(counts,output('egwas',paste0(base,'counts_by_trait.txt')))
}
write_tsv(do.call(rbind,low),output('egwas','post_gwas/lowp.txt'))
write_tsv(data.frame(ID=traits,PVE=pve),output('egwas','post_gwas/PVE.txt'))
capture.output(sessionInfo(),file=output('egwas','post_gwas/sessionInfo.txt'))
