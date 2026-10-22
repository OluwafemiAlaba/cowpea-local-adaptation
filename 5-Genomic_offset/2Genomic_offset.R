#!/usr/bin/env Rscript
# Site offsets from aligned current/future tables: no second monthly shift.
source('config.R')
predictors<-c('prec_08','prec_01','bio_15','bio_3','tmin_10','bio_19')
fold<-read.csv(input('gf','comparison_input/folds.csv'))
baseline<-read.csv(input('gf','baseline_audit/historical_raster_extraction.csv'))
env<-read.csv(input('gf','comparison_input/environment_six.csv'))
stopifnot(identical(fold$Accession,baseline$Accession),identical(env$Accession,baseline$Accession))
tolerance<-c(prec_08=0,prec_01=0,bio_15=.00051,bio_3=.000051,tmin_10=.000001,bio_19=0)
for(v in predictors) stopifnot(all(abs(env[[v]]-baseline[[v]])<=tolerance[v]))
first<-which(!duplicated(fold$site_id)); sites<-fold[first,c('site_id','Longitude','Latitude')]
current<-as.matrix(baseline[first,predictors]); stopifnot(nrow(current)==348,all(is.finite(current)))
cache<-readRDS(result('gf','landscape_inputs/curves_and_ranges.rds'))
# Matches predict.gradientForest's cumulative curves and linear end extrapolation.
transform<-function(curves, data) {
  z<-matrix(NA_real_,nrow(data),length(predictors),dimnames=list(NULL,predictors))
  for(v in predictors) {
    c<-curves[[v]]; a<-data[,v]
    if(length(c$x)==1L) {z[,v]<-c$y; next}
    z[,v]<-approx(c$x,c$y,xout=a,rule=1)$y
    xr<-range(c$x); yr<-range(c$y); dy<-diff(yr)
    hi<-a>xr[2]; lo<-a<xr[1]
    z[hi,v]<-dy*((a[hi]-xr[1])/diff(xr)-1)+yr[2]
    z[lo,v]<-dy*(a[lo]-xr[1])/diff(xr)+yr[1]
  }
  stopifnot(all(is.finite(z))); z
}
lo<-apply(current,2,min); hi<-apply(current,2,max)
rows<-list()
for(gcm in c('ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-HR'))
for(ssp in c('ssp126','ssp245','ssp585'))
for(period in c('2021-2040','2041-2060','2081-2100')) {
  f<-read.csv(input('gf','projection_results',gcm,paste(ssp,period,sep='_'),'future_extractions_all_sites.csv'))
  stopifnot(!anyDuplicated(f$site_id),setequal(f$site_id,sites$site_id))
  f<-f[match(sites$site_id,f$site_id),]
  stopifnot(all(f$Longitude==sites$Longitude),all(f$Latitude==sites$Latitude))
  future<-as.matrix(f[,predictors]);stopifnot(all(is.finite(future)))
  novel<-rowSums(sweep(future,2,lo,'<') | sweep(future,2,hi,'>'))
  for(key in names(cache$curves)) {
    ci<-cache$curves[[key]]; bits<-strsplit(key,'__',fixed=TRUE)[[1]]
    distance<-sqrt(rowSums((transform(ci,current)-transform(ci,future))^2))
    rows[[length(rows)+1L]]<-cbind(sites,climate_model=gcm,scenario=ssp,period=period,
             set=bits[1],seed=as.integer(bits[2]),offset=distance,n_predictors_outside=novel)
  }
}
d<-do.call(rbind,rows)
write.csv(d,output('gf','site_offsets_all_sets_seeds.csv'),row.names=FALSE)
means<-aggregate(offset~site_id+climate_model+scenario+period+set,d,mean)
write.csv(means,output('gf','site_offsets_seed_mean.csv'),row.names=FALSE)
