# Run from project root. Review saved site projections without refitting models.
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(jsonlite))
base <- 'manuscript_reanalysis/GF_validation_design_review'
out <- file.path(base,'projection_review')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(base,'projection_helpers.R'))
# Native macOS Quartz PDF avoids missing optional X11/cairo libraries.
pdf_device<-if(Sys.info()[['sysname']]=='Darwin')function(filename,width,height,...)grDevices::quartz(file=filename,type='pdf',width=width,height=height,...) else grDevices::cairo_pdf
models <- c('ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-HR')
ssps <- c('ssp126','ssp245','ssp585')
periods <- c('2021-2040','2041-2060','2081-2100')
sets <- c('partial_RDA','partial_RDA_eGWAS','all_three')
seeds <- c(460,1460,2460)
fold <- read.csv(file.path(base,'comparison_input/folds.csv'))
hist <- read.csv(file.path(base,'baseline_audit/historical_raster_extraction.csv'))
stopifnot(identical(fold$Accession,hist$Accession))
sites <- fold[!duplicated(fold$site_id),c('site_id','Longitude','Latitude')]
training <- hist[!duplicated(fold$site_id),gf_predictors]
stopifnot(nrow(sites)==348,nrow(fold)==460)
all <- list(); climates <- list(); checks <- list()
for(g in models) for(s in ssps) for(p in periods){
 d <- file.path(base,'projection_results',g,paste(s,p,sep='_'))
 stopifnot(file.exists(file.path(d,'COMPLETE')))
 x <- read.csv(file.path(d,'offsets_all_sets_seeds.csv'))
 f <- read.csv(file.path(d,'future_extractions_all_sites.csv'))
 n <- read.csv(file.path(d,'future_climate_novelty.csv'))
 stopifnot(nrow(x)==348*9,nrow(f)==348,!anyDuplicated(f$site_id),
   all(f$complete_future),all(is.finite(as.matrix(f[,gf_predictors]))),
   all(is.finite(x$offset)),all(x$offset>=0),all(x$scenario==s),
   all(x$period==p),all(x$climate_model==g),setequal(x$set,sets),
   setequal(x$seed,seeds),!anyDuplicated(x[,c('site_id','set','seed')]))
 for(st in sets)for(sd in seeds)stopifnot(setequal(x$site_id[x$set==st & x$seed==sd],sites$site_id))
 si <- sites[match(f$site_id,sites$site_id),]
 stopifnot(isTRUE(all.equal(f$Longitude,si$Longitude)),isTRUE(all.equal(f$Latitude,si$Latitude)))
 ni <- climate_novelty(training,f)
 n <- n[match(f$site_id,n$site_id),]
 stopifnot(identical(ni$outside_training_range,n$outside_training_range),
  all(ni$n_predictors_outside==n$n_predictors_outside),
  max(abs(ni$nearest_standardized_climate_distance-n$nearest_standardized_climate_distance))<1e-10)
 for(v in c('outside_training_range','n_predictors_outside','nearest_standardized_climate_distance')){
  expected <- n[[v]][match(x$site_id,n$site_id)]
  stopifnot(all(abs(as.numeric(x[[v]])-as.numeric(expected))<1e-10))
 }
 counts <- table(fold$site_id)
 stopifnot(all(x$n_plants==as.integer(counts[x$site_id])))
 f$climate_model <- g;f$scenario <- s;f$period <- p
 f <- cbind(f,ni)
 all[[length(all)+1L]] <- x; climates[[length(climates)+1L]] <- f
 checks[[length(checks)+1L]] <- data.frame(climate_model=g,scenario=s,period=p,
   locations=nrow(f),offset_rows=nrow(x),missing=0,nonfinite_offsets=0,
   outside_range=sum(ni$outside_training_range),percent_outside=100*mean(ni$outside_training_range))
}
x <- do.call(rbind,all);f <- do.call(rbind,climates);qc <- do.call(rbind,checks)
write.csv(qc,file.path(out,'combination_QC.csv'),row.names=FALSE)
# Collapse seeds within each fitted SNP set; never average raw climates.
id <- c('site_id','Longitude','Latitude','climate_model','scenario','period','set')
avg <- aggregate(x$offset,x[id],mean);names(avg)[ncol(avg)] <- 'offset'
write.csv(avg,file.path(out,'site_offsets_seed_mean.csv'),row.names=FALSE)
primary <- avg[avg$set=='partial_RDA',]
id2 <- c('site_id','Longitude','Latitude','scenario','period')
ens <- aggregate(primary$offset,primary[id2],median);names(ens)[ncol(ens)] <- 'median_offset'
mn <- aggregate(primary$offset,primary[id2],min);mx <- aggregate(primary$offset,primary[id2],max)
stopifnot(identical(ens[id2],mn[id2]),identical(ens[id2],mx[id2]))
ens$gcm_min <- mn$x;ens$gcm_max <- mx$x;ens$gcm_range <- mx$x-mn$x
prnov <- aggregate(f$outside_training_range,f[id2],sum);names(prnov)[ncol(prnov)] <- 'n_GCMs_outside'
ens <- merge(ens,prnov,by=id2,sort=FALSE)
write.csv(ens,file.path(out,'primary_three_GCM_site_summary.csv'),row.names=FALSE)
# Paired location rankings: sensitivity analyses do not establish predictive accuracy.
rho <- function(a,b)cor(a,b,method='spearman')
seedcorr <- list();setcorr <- list();gcmcorr <- list()
for(g in models)for(s in ssps)for(p in periods){
 z <- x[x$climate_model==g & x$scenario==s & x$period==p,]
 for(st in sets) for(pair in combn(seeds,2,simplify=FALSE)){
  a<-z[z$set==st & z$seed==pair[1],];b<-z[z$set==st & z$seed==pair[2],]
  seedcorr[[length(seedcorr)+1L]]<-data.frame(climate_model=g,scenario=s,period=p,set=st,seed1=pair[1],seed2=pair[2],rho=rho(a$offset,b$offset[match(a$site_id,b$site_id)]))
 }
 z <- avg[avg$climate_model==g & avg$scenario==s & avg$period==p,]
 for(pair in combn(sets,2,simplify=FALSE)){
  a<-z[z$set==pair[1],];b<-z[z$set==pair[2],]
  setcorr[[length(setcorr)+1L]]<-data.frame(climate_model=g,scenario=s,period=p,set1=pair[1],set2=pair[2],rho=rho(a$offset,b$offset[match(a$site_id,b$site_id)]))
 }
}
for(s in ssps)for(p in periods)for(pair in combn(models,2,simplify=FALSE)){
 a<-primary[primary$scenario==s & primary$period==p & primary$climate_model==pair[1],]
 b<-primary[primary$scenario==s & primary$period==p & primary$climate_model==pair[2],]
 gcmcorr[[length(gcmcorr)+1L]]<-data.frame(scenario=s,period=p,gcm1=pair[1],gcm2=pair[2],rho=rho(a$offset,b$offset[match(a$site_id,b$site_id)]))
}
seedcorr<-do.call(rbind,seedcorr);setcorr<-do.call(rbind,setcorr);gcmcorr<-do.call(rbind,gcmcorr)
write.csv(seedcorr,file.path(out,'seed_rank_agreement.csv'),row.names=FALSE)
write.csv(setcorr,file.path(out,'SNP_set_rank_agreement.csv'),row.names=FALSE)
write.csv(gcmcorr,file.path(out,'GCM_rank_agreement_primary.csv'),row.names=FALSE)
summaries <- do.call(rbind,lapply(split(primary,interaction(primary$climate_model,primary$scenario,primary$period,drop=TRUE)),function(z){
 data.frame(z[1,c('climate_model','scenario','period')],n=nrow(z),minimum=min(z$offset),q25=unname(quantile(z$offset,.25)),median=median(z$offset),q75=unname(quantile(z$offset,.75)),maximum=max(z$offset))
}))
write.csv(summaries,file.path(out,'primary_offset_distributions.csv'),row.names=FALSE)
outside_by_var<-do.call(rbind,lapply(gf_predictors,function(v){
 z<-f;z$outside<-z[[v]]<min(training[[v]]) | z[[v]]>max(training[[v]])
 a<-aggregate(z$outside,z[c('climate_model','scenario','period')],sum);names(a)[4]<-'n_outside';a$predictor<-v;a
}))
write.csv(outside_by_var,file.path(out,'extrapolation_by_predictor.csv'),row.names=FALSE)
span<-function(v)list(min=min(v),median=median(v),max=max(v))
report<-list(status='27 combinations passed saved-output consistency and recomputed climate-novelty checks',
 plants=460,locations=348,combinations=nrow(qc),offset_rows=nrow(x),SNP_sets=sets,seeds=seeds,
 seed_rank_agreement=span(seedcorr$rho),SNP_set_rank_agreement=span(setcorr$rho),GCM_rank_agreement_primary=span(gcmcorr$rho),
 outside_training_range_percent=span(qc$percent_outside),
 seasonal_alignment='Southern monthly prec_01=July, prec_08=February, tmin_10=April; bioclimatic scalars unshifted. Verified against extraction code and manifest; not a new independent raster extraction.',
 limitations=c('Site projections only; no spatial interpolation or unsampled-landscape prediction.',
 'Three GCM spread is a sensitivity range, not a probabilistic confidence interval.',
 'Marginal range checks do not detect all novel combinations of predictors.',
 'Geographic companion prediction remained below baseline overall; stable offsets do not validate fitness effects.',
 'Different fitted SNP sets have different GF scales; use ranks for between-set agreement.'))
write_json(report,file.path(out,'validation_summary.json'),pretty=TRUE,auto_unbox=TRUE)
# Plotting: exact coordinates, common scale across all primary panels.
borders <- read.csv('coordinate_investigation/working/gps_reference/map_polygons_africa.csv')
slabs<-c(ssp126='SSP1-2.6',ssp245='SSP2-4.5',ssp585='SSP5-8.5')
format_facets<-function(z){z$scenario<-factor(z$scenario,levels=ssps,labels=slabs);z$period<-factor(z$period,levels=periods);z}
ens<-format_facets(ens);primary<-format_facets(primary);qcplot<-format_facets(qc)
th<-theme_classic(base_size=12)+theme(legend.position='bottom',strip.background=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=10))
mapbase<-function(d)ggplot()+geom_polygon(data=borders,aes(longitude,latitude,group=ring),fill='grey96',colour='grey65',linewidth=.15)+coord_fixed(xlim=c(-20,55),ylim=c(-36,38),expand=FALSE)+facet_grid(period~scenario)+labs(x=NULL,y=NULL)+th+theme(axis.text=element_text(size=8),panel.spacing=grid::unit(.5,'lines'))
saveplot<-function(p,name,w=11,h=12){ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=220,bg='white');ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,bg='white',device=pdf_device)}
p<-mapbase(ens)+guides(colour=guide_colourbar(title.position='top',barwidth=grid::unit(5,'cm')))+geom_point(data=ens,aes(Longitude,Latitude,colour=median_offset),size=1.35)+scale_colour_viridis_c(option='C',begin=.08,end=.95,name='Genomic offset\n(modelled distance)',limits=c(0,max(primary$offset)))+
 labs(title='Exploratory genomic offset at 348 collection locations',subtitle='Primary partial-RDA SNP set; median across three climate models after averaging three GF seeds',caption='Exact collection coordinates; overlaps can conceal sites. One common colour scale across panels.\nPoints represent modelled genomic mismatch, not demonstrated fitness loss. No values are inferred between sites.')
saveplot(p,'Offset_primary_site_maps')
p<-mapbase(ens)+geom_point(data=ens,aes(Longitude,Latitude,colour=factor(n_GCMs_outside,levels=0:3)),size=1.35)+scale_colour_manual(values=c('#20708C','#E3B24A','#D66C30','#882255'),drop=FALSE,name='Climate models with ≥1\npredictor outside training range')+labs(title='Where future climates exceed sampled predictor ranges',subtitle='Range comparison uses the 348 current collection-site climates',caption='A value of zero does not rule out novel multivariate combinations. This is an extrapolation flag, not an exclusion rule.\nThe same climate flags apply to all SNP sets; no sites were silently removed.')
saveplot(p,'Offset_climate_extrapolation_maps')
p<-mapbase(ens)+geom_point(data=ens,aes(Longitude,Latitude,colour=gcm_range),size=1.35)+scale_colour_viridis_c(option='B',name='Maximum − minimum\nacross three GCMs')+labs(title='Climate-model spread in projected genomic offset',subtitle='Primary SNP set; each climate model averaged over three GF seeds',caption='The range reflects these three climate models, not a confidence interval or total projection uncertainty.')
saveplot(p,'Offset_GCM_spread_maps')
for(g in models){
 d<-primary[primary$climate_model==g,]
 p<-mapbase(d)+geom_point(data=d,aes(Longitude,Latitude,colour=offset),size=1.35)+scale_colour_viridis_c(option='C',begin=.08,end=.95,name='Genomic offset',limits=c(0,max(primary$offset)))+labs(title=paste('Exploratory site offsets:',g),subtitle='Primary SNP set; mean across three GF seeds',caption='Identical colour limits across all individual-GCM maps and the three-GCM median map. Exact coordinates; no interpolation.')
 saveplot(p,paste0('Offset_sites_',g))
}
cols<-c('ACCESS-CM2'='#276B8C','MRI-ESM2-0'='#BA6B25','MPI-ESM1-2-HR'='#655898')
p<-ggplot(primary,aes(period,offset,fill=climate_model))+geom_boxplot(outlier.size=.5,linewidth=.35)+facet_wrap(~scenario,nrow=1)+scale_fill_manual(values=cols)+labs(title='Projected offset distributions across collection locations',subtitle='Primary partial-RDA SNP set; seeds averaged within each climate model',x=NULL,y='Genomic offset (modelled distance)',fill=NULL,caption='Each location receives equal weight. Boxes show median and interquartile range; whiskers extend to 1.5 × IQR.\nLocations and scenarios are paired; these distributions are descriptive and are not independent replicates.')+th+theme(axis.text.x=element_text(angle=30,hjust=1))
saveplot(p,'Offset_primary_distributions',13,5.5)
a<-format_facets(gcmcorr);a$pair<-paste(a$gcm1,a$gcm2,sep=' / ')
p<-ggplot(a,aes(period,pair,fill=rho))+geom_tile()+geom_text(aes(label=sprintf('%.2f',rho)),size=3.5)+facet_wrap(~scenario,nrow=1)+scale_fill_viridis_c(limits=c(-1,1),name='Spearman ρ')+labs(title='Agreement between climate-model rankings of collection locations',x=NULL,y=NULL,caption='Primary SNP set; ranks use seed-averaged offsets at the same 348 locations. Agreement is not predictive validation.')+th+theme(axis.text.x=element_text(angle=30,hjust=1))
saveplot(p,'Offset_GCM_rank_agreement',14,4.8)
se<-seedcorr;se$set<-factor(se$set,levels=sets,labels=c('Partial RDA','Partial RDA + eGWAS','All three'))
a<-ggplot(se,aes(set,rho))+geom_boxplot(fill='#8CB5C8',outlier.size=.6)+labs(title='a  Random-seed agreement',x=NULL,y='Spearman correlation of site offsets')+th+theme(axis.text.x=element_text(angle=20,hjust=1))
sc<-setcorr;sc$pair<-factor(paste(sc$set1,sc$set2,sep=' / '),levels=c('partial_RDA / partial_RDA_eGWAS','partial_RDA / all_three','partial_RDA_eGWAS / all_three'),labels=c('Partial vs partial + eGWAS','Partial vs all three','Partial + eGWAS vs all three'))
b<-ggplot(sc,aes(pair,rho))+geom_boxplot(fill='#B9AED0',outlier.size=.6)+labs(title='b  Candidate-set agreement',x=NULL,y='Spearman correlation of site offsets')+th+theme(axis.text.x=element_text(angle=20,hjust=1))
saveplot(a+b+plot_annotation(title='Sensitivity of future-offset location rankings',caption='Each comparison pairs the same 348 locations. Different vertical scales expose the much smaller seed variation.\nSeed comparisons: 27 climates × 3 seed pairs per SNP set. SNP-set comparisons: 27 climates per pair, using seed means.\nThese are descriptive sensitivity checks, not independent replicates or validation against fitness.'),'Offset_seed_and_SNP_set_agreement',12,6)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE),'\n')
