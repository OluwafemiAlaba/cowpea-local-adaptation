# Exact cached cumulative-importance transformation for predict.gradientForest
# with its existing default extrap=TRUE (linear end extrapolation).
transform_curves <- function(curves, data) {
 out <- matrix(NA_real_,nrow(data),length(gf_predictors),dimnames=list(NULL,gf_predictors))
 for(v in gf_predictors){
  ci<-curves[[v]];z<-data[,v]
  if(length(ci$x)==1L){out[,v]<-ci$y;next}
  val<-approx(ci$x,ci$y,xout=z,rule=1)$y
  xr<-range(ci$x);yr<-range(ci$y);dy<-diff(yr)
  hi<-z>xr[2];lo<-z<xr[1]
  val[hi]<-dy*((z[hi]-xr[1])/diff(xr)-1)+yr[2]
  val[lo]<-dy*(z[lo]-xr[1])/diff(xr)+yr[1]
  out[,v]<-val
 }
 stopifnot(all(is.finite(out)));out
}
range_count <- function(data,lo,hi) rowSums(sweep(data,2,lo,'<') | sweep(data,2,hi,'>'))
landscape_sets<-c('partial_RDA','partial_RDA_eGWAS','all_three')
landscape_seeds<-c(460,1460,2460)
landscape_layers<-c('partial_RDA_mean','partial_RDA_eGWAS_mean','all_three_mean',
 'partial_RDA_seed_min','partial_RDA_seed_max','current_n_outside','future_n_outside')
