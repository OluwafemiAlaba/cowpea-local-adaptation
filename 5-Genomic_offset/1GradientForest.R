#!/usr/bin/env Rscript
# Full-data Gradient Forest fits for the three candidate sets and three recorded seeds.
source('config.R')
suppressPackageStartupMessages(library(gradientForest))
env<-read.csv(input('gf','comparison_input/environment_six.csv'),check.names=FALSE)
g<-read.csv(input('gf','comparison_input/genotypes_A1_dosage.csv'),check.names=FALSE)
stopifnot(identical(env$Accession,g$Accession),nrow(env)==460,ncol(env)==7)
x<-env[,-1,drop=FALSE]; y<-g[,-1,drop=FALSE]; train<-seq_len(nrow(x))
stopifnot(all(is.finite(as.matrix(x))))
sets<-lapply(c('partial_RDA','partial_RDA_eGWAS','all_three'),function(s)
  scan(input('gf','comparison_input',paste0(s,'_snps.txt')),what='',quiet=TRUE))
names(sets)<-c('partial_RDA','partial_RDA_eGWAS','all_three')
stopifnot(all(vapply(sets,function(s) !anyDuplicated(s) && all(s %in% names(y)),logical(1))))
curves<-list()
for(seed in c(460L,1460L,2460L)) {
  dest<-output('gf',paste0('comparison_results/seed',seed,'_FULL'))
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
for(s in names(sets)) {
  ys <- y[train,sets[[s]],drop=FALSE]
  # Mean-impute missing dosages and exclude constant SNPs for GF fitting.
  keep <- vapply(ys,function(z) length(unique(z[is.finite(z)]))>1,logical(1))
  dropped <- names(ys)[!keep]; ys <- ys[,keep,drop=FALSE]
  for(j in seq_len(ncol(ys))) ys[!is.finite(ys[,j]),j] <- mean(ys[,j],na.rm=TRUE)
  names(ys) <- make.names(names(ys))
  set.seed(seed)
  maxLevel <- floor(log2(length(train)*0.368/2))
  cat('GF',s,'FULL',seed,'\n')
  gf <- gradientForest::gradientForest(data=cbind(x[train,,drop=FALSE],ys),predictor.vars=names(x),response.vars=names(ys),ntree=500,mtry=2,maxLevel=maxLevel,corr.threshold=0.5,compact=TRUE,trace=FALSE)
  saveRDS(gf,file.path(dest,paste0(s,'_GF.rds')))
  # Save current transformed climates on the same model scale.
  tr <- predict(gf,x)
  write.csv(cbind(Accession=env$Accession,tr),file.path(dest,paste0(s,'_transformed_current.csv')),row.names=FALSE)
  writeLines(dropped,file.path(dest,paste0(s,'_training_constant_SNPs.txt')))
  curves[[paste(s,seed,sep='__')]] <- setNames(lapply(names(x),function(v) {z<-cumimp(gf,v);list(x=z$x,y=z$y)}),names(x))
  rm(gf); gc()
}

}
saveRDS(list(curves=curves),output('gf','landscape_inputs/curves_and_ranges.rds'))
capture.output(sessionInfo(),file=output('gf','sessionInfo.txt'))
