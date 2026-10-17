Sys.setenv(GDAL_CACHEMAX='128')
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(gradientForest))
suppressPackageStartupMessages(library(jsonlite))
source('projection_helpers.R');source('landscape_helpers.R')
dir.create('landscape_inputs',showWarnings=FALSE)
tmp<-file.path('landscape_tmp',paste0('prepare_',Sys.getpid()));dir.create(tmp,recursive=TRUE)
terraOptions(tempdir=tmp,memfrac=.1,memmax=2)
m<-read.csv('historical_manifest.csv');m<-m[match(gf_predictors,m$predictor),]
# Expand portable home-relative raster paths before GDAL/terra reads.
project_root <- path.expand(Sys.getenv('COWPEA_ROOT', '~/data/cowpealandrace_460'))
for (field in c('path_north','path_south')) {
  m[[field]] <- vapply(m[[field]], function(p) {
    prefix <- '~/data/cowpealandrace_460/'
    if (startsWith(p, prefix)) file.path(project_root, substring(p, nchar(prefix)+1L)) else path.expand(p)
  }, character(1))
}
ref<-rast('projection_climate/wc2.1_30s_bioc_ACCESS-CM2_ssp126_2021-2040_Africa_selected.tif')[[1]]
stopifnot(abs(res(ref)[1]-1/120)<1e-10)
# VRT windows remain disk-backed; crop() on a list of global rasters can
# materialize several large crops before the streaming loop starts.
hist<-read.csv('baseline_audit/historical_raster_extraction.csv')
lo<-vapply(hist[gf_predictors],min,numeric(1));hi<-vapply(hist[gf_predictors],max,numeric(1))
if(!file.exists('landscape_inputs/current_aligned.tif')){
paths<-c(m$path_north,m$path_south)
layers<-lapply(seq_along(paths),function(i){
 vrt<-file.path('landscape_inputs',paste0('historical_window_',i,'.vrt'))
 status<-system2('gdal_translate',c('-q','-of','VRT','-projwin',xmin(ref),ymax(ref),xmax(ref),ymin(ref),shQuote(paths[i]),shQuote(vrt)))
 stopifnot(status==0,file.exists(vrt));r<-rast(vrt);stopifnot(compareGeom(ref,r,stopOnError=FALSE));r
})
inputs<-rast(layers);readStart(inputs)
output<-rast(ref,nlyrs=7);names(output)<-c(gf_predictors,'current_n_outside')
writeStart(output,'landscape_inputs/current_aligned.partial.tif',overwrite=TRUE,
 wopt=list(datatype='FLT8S',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=YES')))
for(row in seq(1,nrow(ref),by=64)){
 nr<-min(64,nrow(ref)-row+1);a<-readValues(inputs,row=row,nrows=nr,mat=TRUE)
 south<-rep(yFromRow(ref,row:(row+nr-1))<0,each=ncol(ref))
 current<-a[,1:6,drop=FALSE];current[south,]<-a[south,7:12,drop=FALSE]
 current<-sweep(sweep(current,2,m$multiplier,'*'),2,m$addend,'+')
 colnames(current)<-gf_predictors;ok<-rowSums(!is.finite(current))==0
 nout<-rep(NA_real_,nrow(a));nout[ok]<-range_count(current[ok,,drop=FALSE],lo,hi)
 writeValues(output,cbind(current,nout),row,nr);gc(verbose=FALSE)
 if(row%%640==1)cat('Baseline rows',row,'/',nrow(ref),'\n')
}
writeStop(output);readStop(inputs)
stopifnot(file.rename('landscape_inputs/current_aligned.partial.tif','landscape_inputs/current_aligned.tif'))
}
current<-rast('landscape_inputs/current_aligned.tif')
fold<-read.csv('comparison_input/folds.csv')
reference<-rast(m$path_north[m$predictor=='bio_3'])
cells<-cellFromXY(reference,as.matrix(fold[,c('Longitude','Latitude')]))
centers<-xyFromCell(reference,cells)
addresses<-cbind(fold[,c('Accession','site_id','Longitude','Latitude')],source_global_cell=cells,climate_cell_longitude=centers[,1],climate_cell_latitude=centers[,2])
write.csv(addresses,'landscape_inputs/canonical_climate_cells.csv',row.names=FALSE)
pts<-vect(addresses,geom=c('climate_cell_longitude','climate_cell_latitude'),crs='EPSG:4326')
stopifnot(compareGeom(current,ref,stopOnError=FALSE))
observed<-extract(current[[1:6]],pts,ID=FALSE)
stopifnot(identical(fold$Accession,hist$Accession),max(abs(as.matrix(observed)-as.matrix(hist[gf_predictors])))<1e-10)
# Cache curves, and prove cached evaluation matches the installed predictor,
# including extrapolated climates and all completed future site extractions.
futures<-list.files('projection_results',pattern='future_extractions_all_sites.csv',recursive=TRUE,full.names=TRUE)
test<-rbind(hist[gf_predictors],do.call(rbind,lapply(futures,function(f)read.csv(f)[gf_predictors])))
set.seed(460);random<-as.data.frame(lapply(seq_along(lo),function(j)runif(200,lo[j]-.5*(hi[j]-lo[j]),hi[j]+.5*(hi[j]-lo[j]))));names(random)<-gf_predictors
test<-rbind(test,random);curves<-list();errors<-list()
for(st in landscape_sets)for(sd in landscape_seeds){
 fit<-readRDS(file.path('comparison_results',paste0('seed',sd,'_FULL'),paste0(st,'_GF.rds')))
 ci<-setNames(lapply(gf_predictors,function(v){z<-cumimp(fit,v);list(x=z$x,y=z$y)}),gf_predictors)
 got<-transform_curves(ci,as.matrix(test));expected<-as.matrix(predict(fit,test))[,gf_predictors]
 err<-max(abs(got-expected));stopifnot(err<1e-12)
 key<-paste(st,sd,sep='__');curves[[key]]<-ci
 errors[[key]]<-data.frame(set=st,seed=sd,n_test=nrow(test),max_abs_error=err)
 rm(fit);gc()
}
saveRDS(list(curves=curves,lo=lo,hi=hi), 'landscape_inputs/curves_and_ranges.rds')
write.csv(do.call(rbind,errors),'landscape_inputs/cached_transform_validation.csv',row.names=FALSE)
world<-vect('landscape_inputs/ne_10m_admin_0_countries.geojson')
africa<-world[world$CONTINENT=='Africa',]
rasterize(africa,ref,field=1,background=0,filename='landscape_inputs/Africa_display_mask.tif',overwrite=TRUE,wopt=list(datatype='INT1U',gdal='COMPRESS=DEFLATE'))
write_json(list(status='PASS',grid=list(nrow=nrow(ref),ncol=ncol(ref),resolution=res(ref),extent=as.vector(ext(ref))),
 baseline_max_abs_error=max(abs(as.matrix(observed)-as.matrix(hist[gf_predictors]))),
 cached_prediction_max_error=max(vapply(errors,function(z)z$max_abs_error,numeric(1))),
 extrapolation='Default Gradient Forest linear extrapolation retained; flagged against sampled predictor ranges',
 mask='Natural Earth Africa continent display mask, not a cowpea cultivation or suitability mask'),
 'landscape_inputs/validation.json',pretty=TRUE,auto_unbox=TRUE)
writeLines('COMPLETE','landscape_inputs/COMPLETE')
unlink(tmp,recursive=TRUE)
