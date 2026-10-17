#!/usr/bin/env Rscript
# Adapted from archived gwas_sum.r: corrected array panel, all 92 scans,
# both prespecified thresholds, complete category intersections, no model refit.
library(ggplot2);library(patchwork);library(jsonlite)
r <- 'manuscript_reanalysis/egwas_verified_460';o<-file.path(r,'gwas_summary')
meta<-read.delim(file.path(r,'post_gwas/parsed/cowpea460_phenotype_categories.txt'))
stopifnot(nrow(meta)==92,!anyDuplicated(meta$ID))
categories<-unique(meta$CATEGORY_1)
prettycat<-function(x) {
 x<-gsub('bio1_bio11_Temp_related_Bioclim','Temperature bioclim',x,fixed=TRUE)
 x<-gsub('bio12_bio19_Prec_related_Bioclim','Precipitation bioclim',x,fixed=TRUE)
 x<-gsub('_',' ',x,fixed=TRUE)
 vapply(x,function(z)paste(strwrap(z,width=52),collapse='\n'),character(1))
}
base<-theme_classic(base_size=12)+theme(axis.text=element_text(colour='black'),plot.title=element_text(face='bold',size=13))
savecsv<-function(x,n)write.csv(x,file.path(o,n),row.names=FALSE)
manifest<-list()
for(level in c('per_trait','complete_scan')) {
 h<-read.delim(file.path(r,'post_gwas/parsed',level,'cowpea460_gwas_hitlst.txt'))
 stopifnot(!anyDuplicated(paste(h$ID,h$rs)),all(h$ID %in% meta$ID))
 threshold<-if(level=='per_trait')0.05/41654 else 0.05/(41654*92)
 stopifnot(all(h$p_lrt<threshold))
 full<-read.delim(file.path(r,'post_gwas/parsed/per_trait/cowpea460_gwas_hitlst.txt'))
 stopifnot(setequal(paste(h$ID,h$rs),paste(full$ID[full$p_lrt<threshold],full$rs[full$p_lrt<threshold])))
 bytrait<-meta[,c('ID','CATEGORY_1')];bytrait$n_SNPs<-as.integer(table(factor(h$ID,levels=meta$ID)))
 bysnp<-as.data.frame(table(h$rs));names(bysnp)<-c('SNP','n_traits')
 hc<-merge(h,meta[,c('ID','CATEGORY_1')],by='ID');pairs<-unique(hc[,c('rs','CATEGORY_1')])
 bycat<-data.frame(category=categories,n_SNPs=as.integer(table(factor(pairs$CATEGORY_1,levels=categories))))
 nc<-as.data.frame(table(pairs$rs));names(nc)<-c('SNP','n_categories')
 memberships<-aggregate(CATEGORY_1~rs,pairs,function(x)paste(sort(unique(x)),collapse=' + '))
 intersections<-as.data.frame(table(memberships$CATEGORY_1));names(intersections)<-c('category_combination','n_SNPs');intersections<-intersections[order(-intersections$n_SNPs),]
 stopifnot(sum(bytrait$n_SNPs)==nrow(h),sum(bysnp$n_traits)==nrow(h),sum(intersections$n_SNPs)==nrow(bysnp))
 for(n in c('bytrait','bysnp','bycat','nc','memberships','intersections'))savecsv(get(n),paste0(level,'_',n,'.csv'))
 savecsv(h,paste0(level,'significant_associations.csv'))
 p1<-ggplot(bytrait,aes(n_SNPs))+geom_histogram(binwidth=1,boundary=-.5,fill='#0072B2',colour='white')+labs(title='a  Significant SNPs per climate scan',x='Significant SNPs',y='Climate scans (all 92)')+base
 p2<-ggplot(bysnp,aes(n_traits))+geom_histogram(binwidth=1,boundary=.5,fill='#D55E00',colour='white')+labs(title='b  Trait associations per significant SNP',x='Climate traits',y='Unique SNPs')+base
 p3<-ggplot(bycat,aes(reorder(category,n_SNPs),n_SNPs))+geom_col(fill='#0072B2')+scale_x_discrete(labels=prettycat)+coord_flip()+labs(title='c  Unique SNPs by climate category',x=NULL,y='Unique significant SNPs')+base
 p4<-ggplot(nc,aes(n_categories))+geom_histogram(binwidth=1,boundary=.5,fill='#D55E00',colour='white')+labs(title='d  Category breadth per SNP',x='Climate categories',y='Unique SNPs')+base
 title<-if(level=='per_trait')'Per-trait Bonferroni summary'else'Across-92-trait Bonferroni summary'
 main<-(p1+p2)/(p3+p4)+plot_annotation(title=title,subtitle=paste(nrow(h),'associations;',nrow(bysnp),'unique SNPs;',sum(bytrait$n_SNPs>0),'traits with hits'))
 # Exact mutually exclusive intersections; no top-N truncation.
 overlap<-ggplot(intersections,aes(reorder(category_combination,n_SNPs),n_SNPs))+geom_col(fill='#0072B2')+geom_text(aes(label=n_SNPs),hjust=-.15,size=3.5)+scale_x_discrete(labels=prettycat)+coord_flip()+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='Exact climate-category intersections',subtitle='Each significant SNP appears in exactly one combination',x=NULL,y='Unique SNPs')+base
 traitplot<-ggplot(bytrait,aes(factor(ID,levels=rev(meta$ID)),n_SNPs))+geom_col(fill='#0072B2')+coord_flip()+labs(title='Significant SNPs in all 92 climate scans',x=NULL,y='Significant SNPs')+base+theme(axis.text.y=element_text(size=6))
 snpplot<-ggplot(bysnp,aes(reorder(SNP,n_traits),n_traits))+geom_col(fill='#D55E00')+coord_flip()+labs(title='Trait associations for each significant SNP',x=NULL,y='Climate traits')+base+theme(axis.text.y=element_text(size=6))
 pdf(file.path(o,paste0(level,'_summary.pdf')),width=12,height=9,useDingbats=FALSE);print(main);print(overlap);print(traitplot);print(snpplot);dev.off()
 ggsave(file.path(o,paste0(level,'_summary.png')),main,width=12,height=9,dpi=180,bg='white')
 manifest[[level]]<-list(n_associations=nrow(h),n_unique_SNPs=nrow(bysnp),n_traits_with_hits=sum(bytrait$n_SNPs>0),n_traits_without_hits=sum(bytrait$n_SNPs==0),n_exact_intersections=nrow(intersections),threshold=threshold)
}
pve<-read.delim(file.path(r,'post_gwas/pve.txt'),header=FALSE,col.names=c('ID','PVE'))
stopifnot(nrow(pve)==92,!anyDuplicated(pve$ID),setequal(pve$ID,meta$ID),all(is.finite(pve$PVE)),all(pve$PVE>=0&pve$PVE<=1))
pve<-merge(meta,pve,by='ID');savecsv(pve,'GEMMA_PVE_by_trait.csv')
p1<-ggplot(pve,aes(PVE))+geom_histogram(bins=15,fill='#0072B2',colour='white')+labs(title='a  Distribution across 92 scans',x='GEMMA PVE',y='Climate traits')+base
p2<-ggplot(pve,aes(CATEGORY_1,PVE))+geom_boxplot(outlier.shape=NA,fill='#D55E00',alpha=.35)+geom_point(position=position_jitter(width=.12,height=0,seed=460),size=2)+scale_x_discrete(labels=prettycat)+coord_flip()+labs(title='b  Estimates by climate category',x=NULL,y='GEMMA PVE')+base
p<-p1+p2+plot_annotation(title='Model-based variance explained for collection-site climate',subtitle='Descriptive kinship-model estimates; not evidence of causal climatic adaptation')
ggsave(file.path(o,'GEMMA_PVE_summary.pdf'),p,width=12,height=6)
ggsave(file.path(o,'GEMMA_PVE_summary.png'),p,width=12,height=6,dpi=180,bg='white')
manifest$validation<-'PASS: 92 scans, unique trait-SNP pairs, threshold membership and all count totals reconciled';manifest$refit_models<-FALSE
write_json(manifest,file.path(o,'summary_validation.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
capture.output(sessionInfo(),file=file.path(o,'sessionInfo.txt'))
