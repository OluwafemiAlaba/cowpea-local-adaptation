#!/usr/bin/env Rscript
source('config.R')
work <- output('rda', 'results')
dir.create(work, recursive=TRUE, showWarnings=FALSE)
work <- dirname(work)
if (!file.exists(file.path(work, 'input'))) {
  if (!file.symlink(input('rda', 'input'), file.path(work, 'input')))
    stop('Could not link the deposited RDA inputs.')
}
setwd(work)
suppressPackageStartupMessages({library(vegan); library(adespatial); library(geosphere); library(data.table); library(jsonlite)})
options(stringsAsFactors=FALSE, width=150)
dir.create('results', showWarnings=FALSE)
set.seed(4602026)
logmsg <- function(...) {cat(format(Sys.time()), ..., '\n'); flush.console()}
savejson <- function(x,f) write_json(x,file.path('results',f),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
savecsv <- function(x,f) write.csv(x,file.path('results',f),row.names=FALSE)
namedcsv <- function(x,f) savecsv(data.frame(term=rownames(x),x,check.names=FALSE),f)
aligned <- function(f) {
  d <- as.data.frame(fread(file.path('input',f),check.names=FALSE))
  stopifnot(identical(d$Accession, ids)); rownames(d)<-d$Accession; d$Accession<-NULL; d
}
logmsg('Read observed genotype calls')
g <- as.data.frame(fread('input/genotypes_A1_observed.tsv',check.names=FALSE,na.strings='NA'))
ids <- g$Accession; rownames(g)<-ids; g$Accession<-NULL; G<-as.matrix(g); rm(g)
stopifnot(length(ids)==460,ncol(G)==20027,!anyDuplicated(ids),all(G[!is.na(G)] %in% 0:2))
pl<-as.data.frame(fread('input/plink_independent.raw',check.names=FALSE,na.strings='NA'))
markers<-as.data.frame(fread('input/marker_manifest.csv'))
stopifnot(identical(pl$IID,ids),identical(markers$SNP,colnames(G)),ncol(pl)==ncol(G)+6)
expected_names<-paste(markers$SNP,markers$A1,sep='_')
stopifnot(identical(colnames(pl)[-(1:6)],expected_names))
plg<-as.matrix(pl[,-(1:6)]);stopifnot(identical(unname(is.na(plg)),unname(is.na(G))),all(plg[!is.na(plg)]==G[!is.na(G)]))
savejson(list(n_accessions=nrow(G),n_markers=ncol(G),calls_compared=length(G),all_missingness_and_A1_dosages_match=TRUE), 'PLINK_export_validation.json')
rm(pl,plg)
E <- aligned('environment_aligned.csv'); P<-aligned('structure_PCs_aligned.csv')
M <- aligned('sample_manifest.csv'); stopifnot(ncol(E)==92,ncol(P)==4)
stopifnot(all(is.finite(as.matrix(E))),all(is.finite(as.matrix(P))))
missing <- colSums(is.na(G)); modes <- integer(ncol(G)); ties<-logical(ncol(G))
for(j in seq_len(ncol(G))) {
  counts<-tabulate(G[,j]+1,nbins=3); modes[j]<-which.max(counts)-1L
  ties[j]<-sum(counts==max(counts))>1
  G[is.na(G[,j]),j]<-modes[j]
}
savecsv(data.frame(SNP=colnames(G),missing_n=missing,imputation_mode_A1=modes,mode_tie=ties),'imputation_record.csv')
stopifnot(all(apply(G,2,sd)>0))
Y<-scale(G); rm(G)
stopifnot(abs(sum(Y^2)/(nrow(Y)-1)-ncol(Y))<1e-7)

# An orthogonal, lossless compression preserves all response sums of squares.
eg<-eigen(tcrossprod(Y),symmetric=TRUE)
positive<-eg$values>max(eg$values)*1e-12
Ysmall<-sweep(eg$vectors[,positive,drop=FALSE],2,sqrt(eg$values[positive]),'*')
rownames(Ysmall)<-ids
stopifnot(abs(sum(Ysmall^2)/sum(Y^2)-1)<1e-10)
logmsg('Response compression',ncol(Y),'SNPs to',ncol(Ysmall),'orthogonal dimensions')

# Remove the member of each highest-|r| pair with higher mean correlation.
retained<-colnames(E); corrlog<-data.frame(); viflog<-data.frame()
repeat {
  cr<-abs(cor(E[,retained,drop=FALSE]));diag(cr)<-NA_real_
  high<-max(cr,na.rm=TRUE);if(high<=.7) break
  ij<-which(cr==high,arr.ind=TRUE)[1,];av<-rowMeans(cr,na.rm=TRUE)
  gone<-ij[which.max(av[ij])]
  corrlog<-rbind(corrlog,data.frame(step=nrow(corrlog)+1,removed=retained[gone],correlated_with=retained[setdiff(ij,gone)],absolute_r=high))
  retained<-retained[-gone]
}
repeat {
  vv<-diag(solve(cor(E[,retained,drop=FALSE])))
  if(max(vv)<=5)break
  gone<-which.max(vv)
  viflog<-rbind(viflog,data.frame(step=nrow(viflog)+1,removed=retained[gone],VIF=vv[gone]))
  retained<-retained[-gone]
}
savecsv(corrlog,'climate_correlation_filter.csv');savecsv(viflog,'climate_VIF_filter.csv')
savecsv(data.frame(predictor=retained,VIF=diag(solve(cor(E[,retained,drop=FALSE])))),'climate_screened_predictors.csv')
X<-as.data.frame(scale(E[,retained,drop=FALSE]));stopifnot(qr(as.matrix(X))$rank==ncol(X))

# Build the spatial basis at distinct coordinates, then join scores to accessions.
site_ids<-sort(unique(M$site_id));first<-match(site_ids,M$site_id)
sites<-M[first,c('site_id','Longitude','Latitude')];rownames(sites)<-site_ids
D<-distm(as.matrix(sites[,c('Longitude','Latitude')]),fun=distGeo)/1000
dimnames(D)<-list(site_ids,site_ids)
threshold<-max(spantree(as.dist(D))$dist)
W<-1-(D/(4*threshold))^2;W[D>threshold*(1+1e-6)]<-0;diag(W)<-0
lw<-spdep::mat2listw(W,style='B',zero.policy=TRUE)
B<-scores.listw(lw,MEM.autocor='positive',store.listw=TRUE)
rownames(B)<-site_ids;colnames(B)<-paste0('MEM',seq_len(ncol(B)))
space_all<-as.data.frame(B[match(M$site_id,site_ids),,drop=FALSE]);rownames(space_all)<-ids
savecsv(data.frame(site_id=site_ids,B,check.names=FALSE),'spatial_MEMs_unique_coordinates.csv')
savejson(list(n_coordinate_locations=nrow(sites),distance='WGS84 ellipsoid geodesic km',MST_threshold_km=threshold,
              n_positive_MEM=ncol(B),weight='1 - (distance/(4*MST_threshold))^2 for neighbouring locations; diagonal zero',
              same_coordinate_handling='Basis calculated on distinct coordinates and joined to all accessions by site ID'), 'spatial_basis.json')

logmsg('Spatial global model and forward selection',ncol(space_all),'MEMs')
sp_global<-rda(Ysmall~.,data=space_all,scale=FALSE)
sp_test<-anova(sp_global,permutations=999)
namedcsv(as.data.frame(sp_test),'spatial_global_test.csv')
sp_cap<-RsquareAdj(sp_global)$adj.r.squared
if(is.finite(sp_cap) && sp_cap>0 && sp_test$`Pr(>F)`[1]<=.01) {
  sp_selected<-forward.sel(Ysmall,as.matrix(space_all),adjR2thresh=sp_cap,nperm=999,alpha=.01,Yscale=FALSE,Xscale=TRUE,R2more=0.001,verbose=FALSE)
  selS<-as.character(sp_selected$variables)
  # forward.sel can return column labels or integer indices depending on input.
  if(!all(selS %in% colnames(space_all))) selS<-colnames(space_all)[sp_selected$order]
  savecsv(sp_selected,'spatial_forward_selection.csv')
} else selS<-character()
stopifnot(all(selS %in% colnames(space_all)))
S<-space_all[,selS,drop=FALSE]
covars<-cbind(P,S);C<-as.matrix(covars)
stopifnot(qr(cbind(1,C))$rank==ncol(C)+1)
logmsg('Selected',ncol(S),'MEMs; climate screen retains',ncol(X),'predictors')

# Climate selection is conditional on geography AND the same four genetic PCs
null<-rda(Ysmall~1+Condition(C),data=X,scale=FALSE)
global<-rda(Ysmall~.+Condition(C),data=X,scale=FALSE)
global_test<-anova(global,permutations=9999)
namedcsv(as.data.frame(global_test),'adjusted_climate_global_test.csv')
climate_cap<-RsquareAdj(global)$adj.r.squared
if(is.finite(climate_cap) && climate_cap>0 && global_test$`Pr(>F)`[1]<=.01) {
  chosen<-ordiR2step(null,scope=formula(global),Pin=.01,permutations=999,R2scope=TRUE,trace=TRUE)
  selected<-intersect(attr(terms(chosen),'term.labels'),colnames(X))
  if(!is.null(chosen$anova))namedcsv(as.data.frame(chosen$anova),'climate_forward_selection.csv')
} else {chosen<-null;selected<-character()}
savecsv(data.frame(predictor=selected),'climate_final_predictors.csv')
saveRDS(list(ids=ids,sample_metadata=M,environment_all=E,environment_screened=X,environment_selected=selected,
             structure=P,spatial_selected=S,spatial_all=space_all,site_basis=B,Y_standardized=Y,Y_compressed=Ysmall,
             climate_global=global,climate_selected_compressed=chosen), 'results/analysis_state.rds')
if(!length(selected)) {
  savejson(list(status='No climate predictor passed conditional selection; do not declare candidate loci',
                global_adjusted_R2=climate_cap,global_test=as.data.frame(global_test)), 'run_status.json')
  capture.output(sessionInfo(),file='results/RDA_sessionInfo.txt')
  quit(status=0)
}

suppressPackageStartupMessages({library(vegan);library(jsonlite)})
parallel_cores <- max(1L, as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1')))
options(mc.cores=parallel_cores, width=150)
set.seed(4602027)
logmsg <- function(...) {cat(format(Sys.time()), ..., '\n');flush.console()}
savejson <- function(x,f) write_json(x,file.path('results',f),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
savecsv <- function(x,f) write.csv(x,file.path('results',f),row.names=FALSE)
namedcsv <- function(x,f) savecsv(data.frame(term=rownames(x),x,check.names=FALSE),f)
st<-readRDS('results/analysis_state.rds')
ids<-st$ids;M<-st$sample_metadata;E<-st$environment_all;X<-st$environment_screened
selected<-st$environment_selected;P<-st$structure;S<-st$spatial_selected;B<-st$site_basis
Y<-st$Y_standardized;Ysmall<-st$Y_compressed;C<-as.matrix(cbind(P,S))
retained<-colnames(X);selS<-colnames(S);site_ids<-sort(unique(M$site_id));first<-match(site_ids,M$site_id)
imp<-read.csv('results/imputation_record.csv');missing<-imp$missing_n;ties<-imp$mode_tie
stopifnot(nrow(Y)==460,ncol(Y)==20027,identical(rownames(Y),ids),length(site_ids)==348,
          identical(rownames(P),ids),identical(rownames(S),ids),all(is.finite(Y)),length(selected)>0)
logmsg('Validated selection checkpoint; resume final tests with',parallel_cores,'cores')
logmsg('Fit final adjusted and climate-only models; test axes')
Xs<-X[,selected,drop=FALSE]
primary<-rda(Y~.+Condition(C),data=Xs,scale=FALSE)
comparison<-rda(Y~.,data=Xs,scale=FALSE)
small_primary<-rda(Ysmall~.+Condition(C),data=Xs,scale=FALSE)
small_comparison<-rda(Ysmall~.,data=Xs,scale=FALSE)
stopifnot(max(abs(primary$CCA$eig-small_primary$CCA$eig))<1e-6,
          max(abs(comparison$CCA$eig-small_comparison$CCA$eig))<1e-6,
          abs(RsquareAdj(primary)$adj.r.squared-RsquareAdj(small_primary)$adj.r.squared)<1e-10)
saveRDS(primary,'results/primary_adjusted_RDA.rds');saveRDS(comparison,'results/climate_only_RDA.rds')
v<-vif.cca(primary);savecsv(data.frame(term=names(v),VIF=unname(v)),'primary_conditional_VIF.csv')
savejson(list(status='Fits verified; permutation tests pending',
              primary_R2=RsquareAdj(primary),climate_only_R2=RsquareAdj(comparison),
              primary_eigenvalues=primary$CCA$eig,climate_only_eigenvalues=comparison$CCA$eig,
              full_SNP_and_compressed_fits_agree=TRUE),'fit_checkpoint.json')
test_model<-function(model,prefix) {
  logmsg('Testing',prefix,'overall')
  a<-anova(model,permutations=9999,parallel=parallel_cores)
  namedcsv(as.data.frame(a),paste0(prefix,'_overall_test.csv'))
  logmsg('Testing',prefix,'axes')
  axes<-anova(model,by='axis',permutations=9999,parallel=parallel_cores)
  namedcsv(as.data.frame(axes),paste0(prefix,'_axis_tests.csv'))
  logmsg('Testing',prefix,'marginal terms')
  marginal<-anova(model,by='margin',permutations=9999,parallel=parallel_cores)
  namedcsv(as.data.frame(a),paste0(prefix,'_overall_test.csv'))
  namedcsv(as.data.frame(axes),paste0(prefix,'_axis_tests.csv'))
  namedcsv(as.data.frame(marginal),paste0(prefix,'_marginal_tests.csv'))
  list(overall=a,axes=axes,marginal=marginal)
}
primary_tests<-test_model(small_primary,'primary_adjusted')
comparison_tests<-test_model(small_comparison,'climate_only')

logmsg('Variation partitioning')
if(ncol(S)>0) {
  vp<-varpart(Ysmall,Xs,S,P,scale=FALSE)
  saveRDS(vp,'results/variation_partitioning.rds')
  namedcsv(vp$part$indfract,'variation_partitioning_individual_fractions.csv')
  namedcsv(vp$part$fract,'variation_partitioning_combined_fractions.csv')
  other<-list(pure_space=rda(Ysmall,as.matrix(S),as.matrix(cbind(Xs,P)),scale=FALSE),
              pure_structure=rda(Ysmall,as.matrix(P),as.matrix(cbind(Xs,S)),scale=FALSE))
  for(nm in names(other)) namedcsv(as.data.frame(anova(other[[nm]],permutations=9999,parallel=parallel_cores)),paste0(nm,'_test.csv'))
}

capture.output(sessionInfo(),file='results/RDA_sessionInfo.txt')

for(prefix in c('primary_adjusted','climate_only')) {
 m<-if(prefix=='primary_adjusted')primary else comparison
 L<-as.matrix(scores(m,display='species',choices=seq_len(m$CCA$rank),scaling=2))
 write.csv(data.frame(SNP=rownames(L),L,check.names=FALSE),file.path('results',paste0(prefix,'_SNP_loadings_scaling2.csv')),row.names=FALSE)
}
