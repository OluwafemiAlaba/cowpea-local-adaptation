Sys.setenv(GDAL_CACHEMAX='128')
# Location-level exploratory offsets; run from this directory after climate crops.
# Default GF extrapolation is retained and flagged, not interpreted as validated.
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(gradientForest))
source('projection_helpers.R');source('landscape_helpers.R')
cache<-readRDS('landscape_inputs/curves_and_ranges.rds')
args <- commandArgs(TRUE)
task <- as.integer(args[1]); stopifnot(task %in% 0:8)
gcm <- if(length(args)>1L) args[2] else 'ACCESS-CM2'
stopifnot(gcm %in% c('ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-HR'))
scenarios <- expand.grid(period=c('2021-2040','2041-2060','2081-2100'),
                         scenario=c('ssp126','ssp585','ssp245'),stringsAsFactors=FALSE)
scenario <- scenarios$scenario[task+1L]; period <- scenarios$period[task+1L]
dest <- file.path('projection_results',gcm,paste(scenario,period,sep='_'))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
fold <- read.csv('comparison_input/folds.csv',check.names=FALSE)
env <- read.csv('comparison_input/environment_six.csv',check.names=FALSE)
baseline <- read.csv('baseline_audit/historical_raster_extraction.csv',check.names=FALSE)
stopifnot(identical(fold$Accession,env$Accession),identical(env$Accession,baseline$Accession))
tolerances <- c(prec_08=0,prec_01=0,bio_15=.00051,bio_3=.000051,tmin_10=.000001,bio_19=0)
for(v in gf_predictors) stopifnot(all(abs(env[[v]]-baseline[[v]]) <= tolerances[v]))
# One row per collection coordinate. No density weighting of map summaries.
ix <- which(!duplicated(fold$site_id))
sites <- fold[ix,c('site_id','Longitude','Latitude')]
sites$n_plants <- as.integer(table(fold$site_id)[sites$site_id])
current <- check_environment(baseline[ix,])
future <- current
read_crop <- function(variable) {
  f <- file.path('projection_climate',paste0('wc2.1_30s_',variable,'_',gcm,'_',scenario,'_',period,'_Africa_selected.tif'))
  stopifnot(file.exists(sub('\\.tif$','.complete.json',f)))
  rast(f)
}
b <- read_crop('bioc'); p <- read_crop('prec'); t <- read_crop('tmin')
stopifnot(compareGeom(b,p,stopOnError=FALSE),compareGeom(b,t,stopOnError=FALSE))
addresses<-read.csv('landscape_inputs/canonical_climate_cells.csv')
addresses<-addresses[match(sites$site_id,addresses$site_id),]
stopifnot(identical(sites$site_id,addresses$site_id))
points <- project(vect(addresses,geom=c('climate_cell_longitude','climate_cell_latitude'),crs='EPSG:4326'),crs(b))
bv <- extract(b,points,ID=FALSE); pv <- extract(p,points,ID=FALSE); tv <- extract(t,points,ID=FALSE)
south <- sites$Latitude<0
future$prec_08 <- ifelse(south,pv[,2],pv[,4])
future$prec_01 <- ifelse(south,pv[,3],pv[,1])
future$bio_15 <- bv[,2]; future$bio_3 <- bv[,1]
future$tmin_10 <- ifelse(south,tv[,1],tv[,2]); future$bio_19 <- bv[,3]
ok <- complete.cases(future) & apply(as.matrix(future),1,function(z) all(is.finite(z)))
write.csv(cbind(sites,complete_future=ok,future),file.path(dest,'future_extractions_all_sites.csv'),row.names=FALSE)
if(!all(ok)) stop('Missing future climates: inspect extractions; no sites silently dropped.')
novel <- climate_novelty(current,future)
write.csv(cbind(sites,novel),file.path(dest,'future_climate_novelty.csv'),row.names=FALSE)
rows <- list()
for(seed in c(460,1460,2460)) for(set in c('partial_RDA','partial_RDA_eGWAS','all_three')) {
  ci<-cache$curves[[paste(set,seed,sep='__')]]
  offset<-sqrt(rowSums((transform_curves(ci,as.matrix(current))-transform_curves(ci,as.matrix(future)))^2))
  rows[[length(rows)+1L]] <- cbind(sites,data.frame(scenario=scenario,period=period,
    climate_model=gcm,set=set,seed=seed,offset=offset),novel)
}
result <- do.call(rbind,rows)
write.csv(result,file.path(dest,'offsets_all_sets_seeds.csv'),row.names=FALSE)
comparisons <- list()
for(seed in unique(result$seed)) for(pair in combn(unique(result$set),2,simplify=FALSE)) {
  a <- result[result$seed==seed & result$set==pair[1],]
  z <- result[result$seed==seed & result$set==pair[2],]
  z <- z[match(a$site_id,z$site_id),]
  comparisons[[length(comparisons)+1L]] <- data.frame(seed=seed,set1=pair[1],set2=pair[2],
    n=nrow(a),spearman=cor(a$offset,z$offset,method='spearman'))
}
write.csv(do.call(rbind,comparisons),file.path(dest,'between_set_offset_rank_agreement.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines('COMPLETE: site offsets only; continuous maps and manuscript update pending.',file.path(dest,'COMPLETE'))
