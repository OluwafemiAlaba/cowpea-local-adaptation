suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(pegas))
args <- commandArgs(TRUE)
stopifnot(length(args)==3)
vcf <- args[1]; keepfile <- args[2]; out <- args[3]
expected <- read.table(keepfile,stringsAsFactors=FALSE)[,2]
loci <- pegas::read.vcf(vcf)
genos <- pegas::loci2genind(loci,ploidy=2)
stopifnot(!anyDuplicated(indNames(genos)),setequal(indNames(genos),expected),nInd(genos)==length(expected))
# Scale by allele frequency and replace missing values with the mean,
# followed by ade4 PCA with no additional centring or scaling.
X <- adegenet::scaleGen(genos,NA.method='mean')
stopifnot(all(is.finite(X)))
fit <- ade4::dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=10)
scores <- data.frame(Accession=rownames(fit$li),fit$li,check.names=FALSE)
stopifnot(setequal(scores$Accession,expected))
write.csv(scores,paste0(out,'_scores.csv'),row.names=FALSE)
write.csv(data.frame(PC=seq_along(fit$eig),eigenvalue=fit$eig,percent=100*fit$eig/sum(fit$eig)),paste0(out,'_variance.csv'),row.names=FALSE)
saveRDS(fit,paste0(out,'_fit.rds'))
writeLines(c(paste('samples',nInd(genos)),paste('loci',nLoc(genos)),paste('scaled_allele_columns',ncol(X)),paste('input',vcf)),paste0(out,'_validation.txt'))
capture.output(sessionInfo(),file=paste0(out,'_sessionInfo.txt'))
