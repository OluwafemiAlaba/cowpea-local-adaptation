library(ggplot2)
root <- 'manuscript_reanalysis/go_window20kb/topGO/combined'
d <- read.csv(file.path(root,'results/without_simple_ALL_BP_topGO.csv'))
# Explicit, reviewable selection based on function, not P values or count quantiles.
groups <- list(
 'Environmental response and protection'=c('GO:0006950','GO:0033554','GO:0006974','GO:0006281','GO:0006979','GO:0009628','GO:0009408','GO:0009415','GO:0000302','GO:0009314','GO:0009416','GO:0015979','GO:0009765','GO:0009636','GO:0009725','GO:0009733','GO:0009741','GO:0009742','GO:0009755','GO:0010189','GO:0042360','GO:0010167'),
 'Development and seasonal timing'=c('GO:0009791','GO:0048580','GO:0048573','GO:0009648','GO:2000028','GO:0009908','GO:0009909','GO:0048367','GO:0048831','GO:0090567','GO:0048366','GO:0048827','GO:0099402','GO:0010228','GO:0000902'))
Z<-list()
for(i in seq_along(groups)) {
 ids<-groups[[i]];stopifnot(all(ids%in%d$GO_ID))
 z<-d[match(ids,d$GO_ID),]; stopifnot(all(z$Significant>0))
 z$category<-names(groups)[i];z<-z[order(-z$Significant,z$Term),];Z[[i]]<-z
 z$Term<-factor(z$Term,levels=rev(z$Term))
 p<-ggplot(z,aes(Significant,Term))+geom_col(fill=c('#267F8E','#9B723E')[i],width=.68)+
 scale_x_continuous(expand=expansion(mult=c(0,.14)),breaks=if(i==1)seq(0,40,10)else 0:5)+
 scale_y_discrete(labels=function(x)vapply(x,function(t)paste(strwrap(t,49),collapse='\n'),''))+
 labs(title=paste0(LETTERS[i],'  ',names(groups)[i]),x='Gene count',y=NULL)+
 theme_bw(base_size=14)+theme(plot.title=element_text(face='bold',size=15),
 axis.text=element_text(color='black'),axis.ticks.y=element_blank(),panel.grid.major=element_line(colour='#E8E8E8',linewidth=.4),panel.grid.minor=element_blank(),panel.border=element_rect(colour='#777777',fill=NA,linewidth=.5),plot.margin=margin(15,25,15,15))
 ggsave(file.path(root,paste0('selected_processes_',LETTERS[i],'.png')),p,width=11,height=.32*nrow(z)+1.5,dpi=200)
}
write.csv(do.call(rbind,Z),file.path(root,'selected_environmental_development_processes.csv'),row.names=FALSE)
writeLines(c('Descriptive functional annotations of pooled pRDA and eGWAS candidates, excluding simple RDA. Genes were mapped with the approved +/-20 kb window rule and deduplicated. Selected represented environmental-response, protection, development and seasonal-timing terms are sorted by unique candidate-gene count within each panel. No P-value or top-percentage filter was applied. Selection is thematic, not exhaustive; all tested terms remain in the complete results. One gene may contribute to multiple, related terms, so bars must not be summed. Counts do not establish enrichment or adaptation. DNA damage and repair terms describe general cellular protection and are not specific evidence of climate responses. The pool contains 2,694 genes, with 986 BP-annotated genes. Seed development, germination and ABA-response terms were not added without matching annotations.','',
'Post-embryonic development: 5 genes. Regulation of post-embryonic development: 4 genes, a subset of these five.'),file.path(root,'selected_processes_caption.txt'))
