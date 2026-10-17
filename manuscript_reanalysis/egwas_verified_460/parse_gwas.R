#!/usr/bin/env Rscript
# Adapted from parse_gwas.R by cjfiscus (2022-02-14).
# Base R only. Run: Rscript parse_gwas.R /absolute/path/to/egwas_verified_460
# Uses the completed post-GWAS extracts; no new model fitting or LD selection.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- normalizePath(args[1],mustWork=TRUE)
options(stringsAsFactors=FALSE, digits=17)
Sys.setlocale('LC_COLLATE','C')
input <- file.path(root,'post_gwas')
out <- file.path(input,'parsed')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
traits <- readLines(file.path(root,'input/trait_names.txt'))
stopifnot(length(traits)==92,!anyDuplicated(traits))
tested <- read.table(file.path(input,'cowpea460_gwas_snps_tested.txt'),
                     col.names=c('snps_tested','ID'))
stopifnot(identical(tested$ID,traits),all(tested$snps_tested==41654))
low <- read.table(file.path(input,'cowpea460_gwas_snps_lowp.txt'),
                  col.names=c('ID','chr','rs','ps','n_miss','allele1','allele0',
                    'af','beta','se','logl_H1','l_remle','l_mle','p_wald','p_lrt','p_score'))
stopifnot(all(low$ID %in% traits),all(is.finite(low$p_lrt)),
          all(low$p_lrt>0 & low$p_lrt<1e-4),!anyDuplicated(low[c('ID','rs')]))
low$snps_tested <- tested$snps_tested[match(low$ID,tested$ID)]
low$per_trait_cutoff <- .05/low$snps_tested
low$complete_scan_cutoff <- .05/sum(tested$snps_tested)
low$complete_scan_significant <- low$p_lrt<low$complete_scan_cutoff

# Preserve original categories; reconcile PC IDs and correct display wording.
categories <- read.delim(file.path(root,'source/original_phenotype_categories.txt'))
categories$source_ID <- categories$ID
categories$source_PHENOTYPE <- categories$PHENOTYPE
pc <- grepl('^PC[0-9]+$',categories$ID)
categories$ID[pc] <- paste0('environment_',categories$ID[pc])
categories$PHENOTYPE[pc] <- paste('Environmental principal component',sub('PC','',categories$source_ID[pc]))
categories$PHENOTYPE <- sub('Vapor density','Vapour pressure deficit',categories$PHENOTYPE,fixed=TRUE)
monthly <- grepl('^(prec_|tmax_|tmin_|aridity|rhumav|vpdav)',categories$ID)
categories$PHENOTYPE[monthly] <- paste0(categories$PHENOTYPE[monthly],
    ' [season-aligned label; southern accessions shifted six months]')
stopifnot(!anyDuplicated(categories$ID),all(traits %in% categories$ID))
categories <- categories[match(traits,categories$ID),]
write_out <- function(x,path) write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write_out(categories,file.path(out,'cowpea460_phenotype_categories.txt'))
collapse_traits <- function(x) paste(traits[traits %in% x],collapse=',')
keep <- c('ID','rs','chr','ps','allele1','allele0','af','beta','se','p_lrt',
          'snps_tested','per_trait_cutoff','complete_scan_cutoff','complete_scan_significant')
levels <- c('per_trait','complete_scan')
summary <- list()
for(level in levels) {
    folder <- file.path(out,level)
    dir.create(folder,showWarnings=FALSE)
    hit <- low[if(level=='per_trait') low$p_lrt<low$per_trait_cutoff else low$complete_scan_significant,keep]
    hit <- hit[order(match(hit$ID,traits),hit$chr,hit$ps,hit$rs),]
    hit$significance_label <- ifelse(hit$complete_scan_significant,'complete_scan_Bonferroni','per_trait_Bonferroni_only')
    prefix <- file.path(folder,'cowpea460_gwas_')
    write_out(hit,paste0(prefix,'hitlst.txt'))
    loci <- unique(hit[c('rs','chr','ps')])
    stopifnot(!anyDuplicated(loci$rs))
    loci$Freq <- vapply(loci$rs,function(s)length(unique(hit$ID[hit$rs==s])),integer(1))
    loci$phenotypes <- vapply(loci$rs,function(s)collapse_traits(hit$ID[hit$rs==s]),character(1))
    loci$phenotype_categories <- vapply(loci$rs,function(s)paste(unique(categories$CATEGORY_1[categories$ID %in% hit$ID[hit$rs==s]]),collapse=','),character(1))
    loci$phenotype_categories2 <- vapply(loci$rs,function(s)paste(unique(categories$CATEGORY_2[categories$ID %in% hit$ID[hit$rs==s]]),collapse=','),character(1))
    loci <- loci[order(-loci$Freq,loci$chr,loci$ps,loci$rs),]
    write_out(loci[c('rs','Freq','phenotypes')],paste0(prefix,'hits_by_freq.txt'))
    write_out(loci,paste0(prefix,'hits_by_freq_cat.txt'))
    mtx <- table(factor(hit$rs,levels=sort(loci$rs)),factor(hit$ID,levels=traits))
    stopifnot(all(mtx %in% c(0,1)),sum(mtx)==nrow(hit))
    write.table(as.matrix(mtx),paste0(prefix,'hitmtx.txt'),sep='\t',quote=FALSE,col.names=NA)
    writeLines(sort(loci$rs),paste0(prefix,'hitlst_uniq_snps.txt'))
    counts <- categories[c('ID','PHENOTYPE','CATEGORY_1','CATEGORY_2')]
    counts$n_associations <- as.integer(table(factor(hit$ID,levels=traits)))
    write_out(counts,paste0(prefix,'counts_by_trait.txt'))
    catcounts <- data.frame(category=unique(categories$CATEGORY_1))
    catcounts$n_traits_tested <- vapply(catcounts$category,function(c)sum(categories$CATEGORY_1==c),integer(1))
    catcounts$n_associations <- vapply(catcounts$category,function(c)sum(hit$ID %in% categories$ID[categories$CATEGORY_1==c]),integer(1))
    catcounts$n_unique_SNPs <- vapply(catcounts$category,function(c)length(unique(hit$rs[hit$ID %in% categories$ID[categories$CATEGORY_1==c]])),integer(1))
    write_out(catcounts,paste0(prefix,'counts_by_category.txt'))
    summary[[level]] <- data.frame(threshold=level,n_associations=nrow(hit),n_unique_SNPs=nrow(loci),n_traits_with_hits=length(unique(hit$ID)))
}
write_out(do.call(rbind,summary),file.path(out,'summary_counts.txt'))
writeLines(c('92 separate climate-variable scans; environmental PCs excluded.',
 'Per trait: P < 0.05/41654. Complete scan: P < 0.05/(41654*92).',
 'Prespecified p_lrt used for all filtering; low-P extraction was P < 0.0001.',
 'per_trait/ retains all per-trait hits and labels those also passing the complete-scan cutoff.',
 'complete_scan/ is the stricter subset. These folders are not disjoint sets.',
 'Hit matrices include all 92 traits, including zero-hit columns; rows are SNPs.',
 'Freq is the number of associated traits, not an allele frequency.',
 'Category counts can overlap in SNP membership; do not sum unique SNPs across categories.',
 'SNP lists are associated markers, not LD-independent tag SNPs or validated causal genes.',
 'Monthly display labels reflect seasonal alignment, not a common calendar month across hemispheres.',
 'Environmental-PC category rows from the historical file are excluded from this scan.',
 'The removed duplicate precann is not reintroduced. Original category values are preserved.',
 'af, allele1 and allele0 are copied from GEMMA; af is not relabelled as MAF.',
 'Original parse_gwas.R is archived under source/post_gwas_scripts/.'),file.path(out,'README.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(do.call(rbind,summary),row.names=FALSE)
