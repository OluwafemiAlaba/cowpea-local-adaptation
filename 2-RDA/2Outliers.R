#!/usr/bin/env Rscript
# Putative outliers: absolute loading z score >3 on axes with permutation P <=0.05.
source('config.R')
suppressPackageStartupMessages(library(jsonlite))
dir.create(output('rda','results'),recursive=TRUE,showWarnings=FALSE)
setwd(dirname(output('rda', 'results')))
st<-readRDS(result('rda','results/analysis_state.rds'));Y<-st$Y_standardized
X<-st$environment_screened[,st$environment_selected,drop=FALSE]
C<-cbind(Intercept=1,as.matrix(st$structure),as.matrix(st$spatial_selected))
Yr<-qr.resid(qr(C),Y);Xr<-qr.resid(qr(C),as.matrix(X))
manifest<-read.csv(input('rda','input/marker_manifest.csv'));stopifnot(identical(manifest$SNP,colnames(Y)))
summary<-list()
for(prefix in c('primary_adjusted','climate_only')) {
 tests<-read.csv(result('rda','results',paste0(prefix,'_axis_tests.csv')),check.names=FALSE)
 axes<-tests$term[grepl('^RDA[0-9]+$',tests$term)&!is.na(tests$`Pr(>F)`)&tests$`Pr(>F)`<=.05]
 L<-as.matrix(read.csv(result('rda','results',paste0(prefix,'_SNP_loadings_scaling2.csv')),row.names=1,check.names=FALSE))
 stopifnot(identical(rownames(L),colnames(Y)))
 write.csv(data.frame(SNP=rownames(L),L,check.names=FALSE),file.path('results',paste0(prefix,'_SNP_loadings_scaling2.csv')),row.names=FALSE)
 if(!length(axes)) {
  write.csv(manifest[FALSE,],file.path('results',paste0(prefix,'_putative_outliers_3SD.csv')),row.names=FALSE)
  summary[[prefix]]<-list(axes_screened=axes,n_candidates=0,reason='No supported constrained axes');next
 }
 L<-L[,axes,drop=FALSE];Z<-scale(L);hits<-abs(Z)>3;selected<-rowSums(hits)>0
 correlations<-if(prefix=='primary_adjusted')cor(Yr,Xr)else cor(Y,X)
 best<-max.col(abs(correlations),ties.method='first')
 candidates<-data.frame(manifest[selected,,drop=FALSE],assigned_predictor=colnames(X)[best[selected]],association_r=correlations[cbind(which(selected),best[selected])],max_abs_loading_z=if(any(selected))apply(abs(Z[selected,,drop=FALSE]),1,max)else numeric(),outlier_axes=if(any(selected))apply(hits[selected,,drop=FALSE],1,function(z)paste(axes[z],collapse=';'))else character(),threshold_SD=rep(3,sum(selected)))
 write.csv(candidates,file.path('results',paste0(prefix,'_putative_outliers_3SD.csv')),row.names=FALSE)
 counts<-data.frame(predictor=colnames(X),n_candidates=as.integer(table(factor(candidates$assigned_predictor,levels=colnames(X)))))
 write.csv(counts,file.path('results',paste0(prefix,'_candidate_counts_by_predictor.csv')),row.names=FALSE)
 write.csv(data.frame(axis=axes,mean=colMeans(L),sd=apply(L,2,sd),lower=colMeans(L)-3*apply(L,2,sd),upper=colMeans(L)+3*apply(L,2,sd)),file.path('results',paste0(prefix,'_loading_cutoffs_3SD.csv')),row.names=FALSE)
 summary[[prefix]]<-list(axes_screened=axes,n_candidates=nrow(candidates),threshold_SD=3)
}
write_json(summary,'results/candidate_screening_summary.json',pretty=TRUE,auto_unbox=TRUE)
