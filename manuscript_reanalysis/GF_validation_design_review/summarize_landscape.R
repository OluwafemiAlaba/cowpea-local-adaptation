Sys.setenv(GDAL_CACHEMAX='128')
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(jsonlite))
task<-as.integer(commandArgs(TRUE)[1]);stopifnot(task %in% 0:8)
grid<-expand.grid(period=c('2021-2040','2041-2060','2081-2100'),scenario=c('ssp126','ssp245','ssp585'),stringsAsFactors=FALSE)
e<-grid[task+1,];p<-e$period;s<-e$scenario
models<-c('ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-HR')
dirs<-file.path('landscape_results',models,paste(s,p,sep='_'))
stopifnot(all(file.exists(file.path(dirs,'COMPLETE'))))
dest<-file.path('landscape_ensemble',paste(s,p,sep='_'));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
tmp<-file.path('landscape_tmp',paste0('ensemble_',task,'_',Sys.getpid()));dir.create(tmp,recursive=TRUE)
terraOptions(tempdir=tmp,memfrac=.1,memmax=2)
rasters<-lapply(dirs,function(d)rast(file.path(d,'offsets.tif'))[[c(1,6,7)]])
stopifnot(compareGeom(rasters[[1]],rasters[[2]],stopOnError=FALSE),compareGeom(rasters[[1]],rasters[[3]],stopOnError=FALSE))
inp<-rast(rasters);ref<-inp[[1]];readStart(inp)
africa<-rast('landscape_inputs/Africa_display_mask.tif');readStart(africa)
out<-rast(ref,nlyrs=7);names(out)<-c('median_offset','gcm_min','gcm_max','gcm_range','current_outside','future_GCMs_outside','n_GCMs_complete')
writeStart(out,file.path(dest,'ensemble.partial.tif'),overwrite=TRUE,wopt=list(datatype='FLT4S',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=YES')))
counts<-c(africa_valid=0,current_outside=0,future_any_outside=0,future_all_outside=0,missing_model_cells=0)
for(row in seq(1,nrow(ref),by=64)){
 nr<-min(64,nrow(ref)-row+1);a<-readValues(inp,row=row,nrows=nr,mat=TRUE)
 mask<-readValues(africa,row=row,nrows=nr)
 v<-a[,c(1,4,7),drop=FALSE];ngcm<-rowSums(is.finite(v));ok<-ngcm==3
 vals<-matrix(NA_real_,nrow(a),7);vals[,7]<-ngcm
 if(any(ok)){
  mn<-pmin(v[ok,1],v[ok,2],v[ok,3]);mx<-pmax(v[ok,1],v[ok,2],v[ok,3])
  vals[ok,1]<-rowSums(v[ok,,drop=FALSE])-mn-mx;vals[ok,2]<-mn;vals[ok,3]<-mx;vals[ok,4]<-mx-mn
  stopifnot(all(a[ok,2]==a[ok,5]),all(a[ok,2]==a[ok,8]))
  vals[ok,5]<-as.numeric(a[ok,2]>0);vals[ok,6]<-rowSums(a[ok,c(3,6,9),drop=FALSE]>0)
 }
 land<-is.finite(mask) & mask==1;counts<-counts+c(sum(land & ok),sum(land & ok & vals[,5]>0,na.rm=TRUE),sum(land & ok & vals[,6]>0,na.rm=TRUE),sum(land & ok & vals[,6]==3,na.rm=TRUE),sum(land & ngcm<3))
 writeValues(out,vals,row,nr);gc(verbose=FALSE)
}
writeStop(out);readStop(inp);readStop(africa)
stopifnot(file.rename(file.path(dest,'ensemble.partial.tif'),file.path(dest,'ensemble.tif')))
v<-rast(file.path(dest,'ensemble.tif'))
# Plot-scale cell summaries only after native-cell model transformations and
# native-cell ensemble calculation; no average-climate prediction shortcut.
display<-c(aggregate(v[[1:4]],12,'mean',na.rm=TRUE),aggregate(v[[5:7]],12,'max',na.rm=TRUE));names(display)<-names(v)
mask<-aggregate(africa,12,'max',na.rm=TRUE);display<-terra::mask(display,mask,maskvalues=0)
writeRaster(display,file.path(dest,'display_6arcmin.tif'),overwrite=TRUE,wopt=list(gdal='COMPRESS=DEFLATE'))
write_json(list(status='PASS',scenario=s,period=p,counts=as.list(counts),
 count_unit='native grid cells, not equal-area land percentages',
 ensemble='Median, min and max of three independently transformed GCM offsets, after seed means, at each native cell',
 display='Means of native-cell offsets within 12x12 blocks; maximum native-cell flag counts within blocks. Only Africa continent displayed; no cultivation mask.'),
 file.path(dest,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines('COMPLETE',file.path(dest,'COMPLETE'));unlink(tmp,recursive=TRUE)
