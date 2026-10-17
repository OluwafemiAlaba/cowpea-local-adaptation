args <- commandArgs(TRUE)
task <- as.integer(args[1])
stopifnot(!is.na(task), task >= 0, task < 18)
suppressPackageStartupMessages(library(gradientForest))
seed <- c(460L, 1460L, 2460L)[task %/% 6 + 1]
held <- c('FULL', paste0('F',1:5))[task %% 6 + 1]
dest <- file.path('comparison_results', paste0('seed',seed,'_',held))
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
env <- read.csv('comparison_input/environment_six.csv',check.names=FALSE)
g <- read.csv('comparison_input/genotypes_A1_dosage.csv',check.names=FALSE)
fold <- read.csv('comparison_input/folds.csv')
stopifnot(identical(env$Accession,g$Accession),identical(env$Accession,fold$Accession))
x <- env[,-1,drop=FALSE]; y <- g[,-1,drop=FALSE]
sets <- lapply(c('partial_RDA','partial_RDA_eGWAS','all_three'), function(s) scan(file.path('comparison_input',paste0(s,'_snps.txt')),what='',quiet=TRUE))
names(sets) <- c('partial_RDA','partial_RDA_eGWAS','all_three')
train <- if(held=='FULL') seq_len(nrow(x)) else which(fold$fold!=held)
test <- if(held=='FULL') integer() else which(fold$fold==held)
writeLines(capture.output(sessionInfo()), file.path(dest,'sessionInfo.txt'))

# Each SNP has identical predictor inputs and seed across membership sets, so
# fit it once per fold/seed and reuse its held-out score in all applicable sets.
if(length(test)) {
  rows <- vector('list',ncol(y))
  for(j in seq_len(ncol(y))) {
    observed_train <- train[is.finite(y[train,j])]
    observed_test <- test[is.finite(y[test,j])]
    mu <- mean(y[observed_train,j])
    pred <- rep(mu,length(observed_test))
    constant <- length(unique(y[observed_train,j]))<2
    if(!constant && length(observed_test)) {
      set.seed(seed+j)
      fit <- randomForest::randomForest(x=x[observed_train,,drop=FALSE],y=y[observed_train,j],ntree=500,mtry=2,nodesize=5)
      pred <- randomForest:::predict.randomForest(fit,x[observed_test,,drop=FALSE])
    }
    err <- (pred-y[observed_test,j])^2
    null <- (mu-y[observed_test,j])^2
    # Average errors within location first, then give locations equal weight.
    site_err <- tapply(err,fold$site_id[observed_test],mean)
    site_null <- tapply(null,fold$site_id[observed_test],mean)
    rows[[j]] <- data.frame(SNP=names(y)[j],seed=seed,fold=held,n_train=length(observed_train),n_test=length(observed_test),n_test_sites=length(site_err),training_constant=constant,mse=mean(site_err),baseline_mse=mean(site_null),skill=if(sum(site_null)>0) 1-sum(site_err)/sum(site_null) else NA_real_)
    if(j%%100==0) cat('RF',held,seed,j,'/',ncol(y),'\n')
  }
  scores <- do.call(rbind,rows)
  write.csv(scores,file.path(dest,'per_SNP_holdout_scores.csv'),row.names=FALSE)
  summary <- do.call(rbind,lapply(names(sets),function(s) {
    z<-scores[scores$SNP %in% sets[[s]],]
    data.frame(set=s,seed=seed,fold=held,SNPs=nrow(z),mean_mse=mean(z$mse),mean_baseline_mse=mean(z$baseline_mse),pooled_skill=1-sum(z$mse)/sum(z$baseline_mse),median_SNP_skill=median(z$skill,na.rm=TRUE),fraction_positive_skill=mean(z$skill>0,na.rm=TRUE),training_constant=sum(z$training_constant))
  }))
  write.csv(summary,file.path(dest,'set_holdout_scores.csv'),row.names=FALSE)
}

for(s in names(sets)) {
  ys <- y[train,sets[[s]],drop=FALSE]
  # GF requires complete responses: training-only mean imputation. This is not
  # used for held-out genotype scoring. Exclude training-constant responses.
  keep <- vapply(ys,function(z) length(unique(z[is.finite(z)]))>1,logical(1))
  dropped <- names(ys)[!keep]; ys <- ys[,keep,drop=FALSE]
  for(j in seq_len(ncol(ys))) ys[!is.finite(ys[,j]),j] <- mean(ys[,j],na.rm=TRUE)
  names(ys) <- make.names(names(ys))
  set.seed(seed)
  maxLevel <- floor(log2(length(train)*0.368/2))
  cat('GF',s,held,seed,'\n')
  gf <- gradientForest::gradientForest(data=cbind(x[train,,drop=FALSE],ys),predictor.vars=names(x),response.vars=names(ys),ntree=500,mtry=2,maxLevel=maxLevel,corr.threshold=0.5,compact=TRUE,trace=FALSE)
  saveRDS(gf,file.path(dest,paste0(s,'_GF.rds')))
  # Transform the same current-climate locations for cross-fit stability review;
  # these distances are descriptive, not future offsets or genotype predictions.
  tr <- predict(gf,x)
  write.csv(cbind(Accession=env$Accession,tr),file.path(dest,paste0(s,'_transformed_current.csv')),row.names=FALSE)
  writeLines(dropped,file.path(dest,paste0(s,'_training_constant_SNPs.txt')))
  rm(gf); gc()
}
writeLines('COMPLETE',file.path(dest,'COMPLETE'))
