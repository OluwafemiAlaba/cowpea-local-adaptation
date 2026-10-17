library(ggplot2)
root <- 'manuscript_reanalysis/go_window20kb/topGO/combined'
d <- read.csv(file.path(root,'results/without_simple_ALL_BP_topGO.csv'))
# Descriptive, explicitly theme-restricted display: no P-value filter.
pattern <- paste(c('stress','response to (abiotic|heat|cold|water|salt|osmotic|light|radiation|auxin|hormone|brassinosteroid|abscisic|ethylene)',
'photoperiod','flower','seed','germinat','circadian','vitamin E','tocopherol','photosynthesis','antioxid',
'hormone','abscisic','auxin','ethylene','brassinosteroid','jasmon','salicy','gibberell'),collapse='|')
d$theme_eligible <- grepl(pattern,d$Term,ignore.case=TRUE) & d$Significant>0
r <- d[d$theme_eligible,]; r <- r[order(-r$Significant,r$Term,r$GO_ID),]
target <- ceiling(nrow(r)*.2); cutoff <- r$Significant[target]
r$displayed <- r$Significant>=cutoff
write.csv(r,file.path(root,'environmental_process_counts_all.csv'),row.names=FALSE)
z <- r[r$displayed,]
write.csv(z,file.path(root,'environmental_process_counts_top20.csv'),row.names=FALSE)
z$Term <- factor(z$Term,levels=rev(z$Term))
p <- ggplot(z,aes(Significant,Term))+geom_col(fill='#267F8E',width=.65)+
 geom_text(aes(label=Significant),hjust=-.35,size=5)+
 scale_x_continuous(expand=expansion(mult=c(0,.13)),breaks=seq(0,40,10))+
 labs(x='Gene count',y=NULL)+theme_classic(base_size=16)+
 theme(axis.text=element_text(color='black'),axis.ticks.y=element_blank(),plot.margin=margin(15,25,15,15))
ggsave(file.path(root,'without_simple_environmental_counts_top20.png'),p,width=10,height=5,dpi=220)
caption <- sprintf('Descriptive counts of environmentally relevant processes among the pooled pRDA and eGWAS candidate genes, excluding simple RDA. Genes were mapped using +/-20 kb windows and counted once per term. The top 20%% of %d represented, theme-matched Biological Process terms were ranked by gene count, retaining ties: %d terms with at least %d genes. No significance filter was applied. Broad ontology headings were excluded. Genes can occur in several terms; these counts do not establish enrichment. The underlying pool contains 2,694 unique genes, of which 986 have usable Biological Process annotations.',nrow(r),nrow(z),cutoff)
writeLines(c(caption,'',paste('Reproducible thematic matching rule:',pattern)),file.path(root,'environmental_process_counts_caption.txt'))
cat(caption,'\n')
