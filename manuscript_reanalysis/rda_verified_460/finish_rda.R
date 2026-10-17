#!/usr/bin/env Rscript
# Resume final tests from the completed selection checkpoint, without reselection.
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
savejson(list(status='Fits verified; permutation and location checks pending',
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

# Equal-location check holds the primary predictor/space sets and SNP scaling
# fixed, averaging standardized dosages and genetic PCs within each coordinate.
# It is a robustness check, not independent replication or validation.
mean_by_site<-function(z) {a<-rowsum(as.matrix(z),M$site_id,reorder=TRUE);a/as.vector(table(factor(M$site_id,levels=rownames(a))))}
siteY<-mean_by_site(Ysmall); siteX<-Xs[first,,drop=FALSE]
siteC<-cbind(mean_by_site(P),as.matrix(B[,selS,drop=FALSE]))
stopifnot(identical(rownames(siteY),site_ids))
sitefit<-rda(siteY~.+Condition(siteC),data=siteX,scale=FALSE)
site_tests<-test_model(sitefit,'equal_location')
saveRDS(sitefit,'results/equal_location_check.rds')

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
summ<-function(m) list(n=nrow(m$CA$u),constrained_rank=m$CCA$rank,conditioning_rank=m$pCCA$rank,
                       R2=RsquareAdj(m),constrained_fraction_total=m$CCA$tot.chi/m$tot.chi,
                       eigenvalues=m$CCA$eig)
savejson(list(status='Models completed; candidate screening and manuscript figures require review',n_plants=460,n_locations=348,n_SNPs=ncol(Y),
              n_imputed_calls=sum(missing),n_mode_ties=sum(ties),climate_screened=retained,climate_selected=selected,spatial_selected=selS,
              primary=summ(primary),climate_only=summ(comparison),equal_location=summ(sitefit),
              response_scaling='All SNPs standardized once to variance 1 after per-SNP mode imputation; same scaling in every model',
              compression_validation='Final full-SNP and compressed-response model eigenvalues and adjusted R2 agree',
              permutations=list(selection=999,final_tests=9999,seed_spatial_selection=4602026,seed_climate_selection=4602028,seed_final_tests=4602027,parallel_cores=parallel_cores,scheme='Unrestricted reduced-model residual permutations'),
              caveats=c('Permutation tests after predictor and spatial selection are conditional descriptive tests, not independent validation.',
                        'Primary row permutations assume adequate exchangeability after conditioning; coordinate-level check reduces repeated-location weighting, not all spatial dependence.',
                        'Structure PCs derive from the same array cohort and may remove adaptive covariance shared with ancestry.',
                        'Array ascertainment, retained LD, provenance uncertainty, and imputation remain relevant.')), 'run_status.json')
capture.output(sessionInfo(),file='results/RDA_sessionInfo.txt')
logmsg('Complete')
