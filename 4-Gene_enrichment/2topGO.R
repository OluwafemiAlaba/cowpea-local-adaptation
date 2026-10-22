#!/usr/bin/env Rscript
# topGO biological processes; elim Fisher P values are exploratory nominal evidence.
source('config.R')
suppressPackageStartupMessages(library(topGO))
dir.create(output('go','results'),recursive=TRUE,showWarnings=FALSE)
setwd(dirname(output('go','results')))
capture.output(sessionInfo(),file='results/sessionInfo.txt')
capture.output(AnnotationDbi::metadata(GO.db),file='results/GO_database_metadata.txt')
a <- read.delim(input('go','topGO/original/cowpea_go.txt'),check.names=FALSE)
a <- unique(a[,c('Gene Name','GO ID')]); a <- a[!is.na(a[[1]]) & grepl('^GO:[0-9]{7}$',a[[2]]),]
map <- lapply(split(a[[2]],a[[1]]),unique)
sets <- 'without_simple'
summary <- list()
for (s in sets) {
 bg <- unique(readLines(paste0('inputs/',s,'_background_genes.txt')))
 fg <- unique(readLines(paste0('inputs/',s,'_candidate_genes.txt')))
 stopifnot(all(fg %in% bg))
 selected <- setNames(as.integer(bg %in% fg),bg)
 go <- new('topGOdata',description=s,ontology='BP',allGenes=selected,
           geneSelectionFun=function(x)x==1,annot=annFUN.gene2GO,
           gene2GO=map[names(map)%in%bg],nodeSize=1)
 # Explicit historical elim default: child terms with P < .01 are eliminated upstream.
 elim <- runTest(go,algorithm='elim',statistic='fisher',cutOff=0.01)
 classic <- runTest(go,algorithm='classic',statistic='fisher')
 p <- score(elim); ids <- names(p)
 z <- termStat(go,ids); z$GO_ID <- rownames(z)
 z$Term <- unname(Term(GOTERM[ids]))
 z$elim_P <- as.numeric(p[ids]); z$classic_P <- as.numeric(score(classic)[ids])
 z$elim_BH_descriptive <- p.adjust(z$elim_P,'BH')
 z$classic_BH <- p.adjust(z$classic_P,'BH')
 z$nominal_elim_P_below_005 <- z$elim_P < .05
 genes <- genesInTerm(go,ids)
 z$candidate_genes <- vapply(genes[ids],function(x)paste(sort(intersect(x,fg)),collapse=';'),'')
 z$fold_enrichment <- ifelse(z$Expected>0,z$Significant/z$Expected,NA_real_)
 z <- z[order(z$elim_P,z$GO_ID),]
 write.csv(z,paste0('results/',s,'_ALL_BP_topGO.csv'),row.names=FALSE)
 write.csv(z[z$Significant>0,],paste0('results/',s,'_represented_BP.csv'),row.names=FALSE)
 saveRDS(list(data=go,elim=elim,classic=classic),paste0('results/',s,'_topGO.rds'))
 summary[[s]] <- data.frame(analysis=s,candidate_genes=length(fg),background_genes=length(bg),
 BP_feasible_background=numGenes(go),BP_feasible_candidates=numSigGenes(go),
 tested_terms=nrow(z),represented_terms=sum(z$Significant>0),
 nominal_elim_terms=sum(z$elim_P<.05),elim_BH_below_005=sum(z$elim_BH_descriptive<.05),
 classic_BH_below_005=sum(z$classic_BH<.05))
}
write.csv(do.call(rbind,summary),'results/without_simple_summary.csv',row.names=FALSE)
writeLines('COMPLETE','results/COMPLETE_without_simple')
