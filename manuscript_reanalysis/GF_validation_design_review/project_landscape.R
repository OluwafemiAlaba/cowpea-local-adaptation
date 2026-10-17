Sys.setenv(GDAL_CACHEMAX='128')
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(jsonlite))
source('projection_helpers.R');source('landscape_helpers.R')
task<-as.integer(commandArgs(TRUE)[1]);stopifnot(task %in% 0:26,file.exists('landscape_inputs/COMPLETE'))
grid<-expand.grid(period=c('2021-2040','2041-2060','2081-2100'),scenario=c('ssp126','ssp245','ssp585'),climate_model=c('ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-HR'),stringsAsFactors=FALSE)
e<-grid[task+1,];g<-e$climate_model;s<-e$scenario;p<-e$period
dest<-file.path('landscape_results',g,paste(s,p,sep='_'));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,'COMPLETE')))quit(status=0)
tmp<-file.path('landscape_tmp',paste0(task,'_',Sys.getpid()));dir.create(tmp,recursive=TRUE)
terraOptions(tempdir=tmp,memfrac=.1,memmax=2)
cache<-readRDS('landscape_inputs/curves_and_ranges.rds')
current<-rast('landscape_inputs/current_aligned.tif');ref<-current[[1]]
rd<-function(v)rast(file.path('projection_climate',paste0('wc2.1_30s_',v,'_',g,'_',s,'_',p,'_Africa_selected.tif')))
b<-rd('bioc');pr<-rd('prec');tm<-rd('tmin')
stopifnot(compareGeom(ref,b,stopOnError=FALSE),compareGeom(ref,pr,stopOnError=FALSE),compareGeom(ref,tm,stopOnError=FALSE))
inp<-c(current,b,pr,tm);readStart(inp)
out<-rast(ref,nlyrs=length(landscape_layers));names(out)<-landscape_layers
partial<-file.path(dest,'offsets.partial.tif');final<-file.path(dest,'offsets.tif')
writeStart(out,partial,overwrite=TRUE,wopt=list(datatype='FLT4S',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=YES')))
valid_cells<-0;current_only_missing<-0;maxima<-rep(-Inf,3)
for(row in seq(1,nrow(ref),by=64)){
 nr<-min(64,nrow(ref)-row+1);a<-readValues(inp,row=row,nrows=nr,mat=TRUE)
 south<-rep(yFromRow(ref,row:(row+nr-1))<0,each=ncol(ref))
 now<-a[,1:6,drop=FALSE];colnames(now)<-gf_predictors
 # Stack columns: current six + flag; BIO3/15/19; PREC1/2/7/8; TMIN4/10.
 future<-cbind(prec_08=ifelse(south,a[,12],a[,14]),prec_01=ifelse(south,a[,13],a[,11]),
  bio_15=a[,9],bio_3=a[,8],tmin_10=ifelse(south,a[,15],a[,16]),bio_19=a[,10])
 oknow<-rowSums(!is.finite(now))==0;ok<-oknow & rowSums(!is.finite(future))==0
 valid_cells<-valid_cells+sum(ok);current_only_missing<-current_only_missing+sum(oknow & !ok)
 vals<-matrix(NA_real_,nrow(a),7)
 if(any(ok)){
  now<-now[ok,,drop=FALSE];future<-future[ok,,drop=FALSE]
  for(j in seq_along(landscape_sets)){
   vals_sd<-vapply(landscape_seeds,function(sd){ci<-cache$curves[[paste(landscape_sets[j],sd,sep='__')]];sqrt(rowSums((transform_curves(ci,now)-transform_curves(ci,future))^2))},numeric(sum(ok)))
   vals[ok,j]<-rowMeans(vals_sd)
   if(j==1L){vals[ok,4]<-pmin(vals_sd[,1],vals_sd[,2],vals_sd[,3]);vals[ok,5]<-pmax(vals_sd[,1],vals_sd[,2],vals_sd[,3])}
  }
  vals[ok,6]<-a[ok,7];vals[ok,7]<-range_count(future,cache$lo,cache$hi)
  stopifnot(all(is.finite(vals[ok,])),all(vals[ok,1:5]>=0))
  maxima<-pmax(maxima,apply(vals[ok,1:3,drop=FALSE],2,max))
 }
 writeValues(out,vals,row,nr);gc(verbose=FALSE)
 if(row%%640==1)cat('Task',task,'rows',row,'/',nrow(ref),'valid cells',valid_cells,'\n')
}
writeStop(out);readStop(inp);stopifnot(file.rename(partial,final))
# Quantitative acceptance: native raster-cell offsets reproduce the previously
# accepted site offsets, with only float32 storage precision allowed.
site<-read.csv(file.path('projection_results',g,paste(s,p,sep='_'),'offsets_all_sets_seeds.csv'))
loc<-site[!duplicated(site$site_id),c('site_id','Longitude','Latitude')]
addresses<-read.csv('landscape_inputs/canonical_climate_cells.csv')
addresses<-addresses[match(loc$site_id,addresses$site_id),]
stopifnot(identical(loc$site_id,addresses$site_id))
v<-extract(rast(final),vect(addresses,geom=c('climate_cell_longitude','climate_cell_latitude'),crs='EPSG:4326'),ID=FALSE)
errors<-numeric(3)
for(j in 1:3){expected<-aggregate(offset~site_id,site[site$set==landscape_sets[j],],mean);z<-expected$offset[match(loc$site_id,expected$site_id)];errors[j]<-max(abs(v[,j]-z));stopifnot(all(is.finite(v[,j])),errors[j]<1e-7)}
write.csv(cbind(loc,v),file.path(dest,'raster_cell_site_check.csv'),row.names=FALSE)
write_json(list(status='PASS',climate_model=g,scenario=s,period=p,valid_cells=valid_cells,
 current_valid_future_missing_cells=current_only_missing,site_max_abs_errors=errors,maximum_offset_by_set=maxima,
 storage='Float32 output; calculations in double precision',layers=landscape_layers,
 flags='Counts of predictors outside sampled ranges for current and future climates separately; zero does not guarantee multivariate familiarity'),
 file.path(dest,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
# Display-only aggregation after native-cell predictions. It is not a refit or
# a prediction from averaged climates. Keep flags as any flagged native cell.
v<-rast(final);display<-c(aggregate(v[[1:5]],fact=12,fun='mean',na.rm=TRUE),aggregate(v[[6:7]],fact=12,fun='max',na.rm=TRUE))
names(display)<-landscape_layers
writeRaster(display,file.path(dest,'display_6arcmin.tif'),overwrite=TRUE,wopt=list(gdal='COMPRESS=DEFLATE'))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines('COMPLETE: native raster and site checks passed',file.path(dest,'COMPLETE'))
unlink(tmp,recursive=TRUE)
