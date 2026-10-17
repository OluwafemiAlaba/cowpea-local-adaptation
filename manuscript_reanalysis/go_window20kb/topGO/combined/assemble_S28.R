library(ggplot2);library(patchwork)
root<-'manuscript_reanalysis/go_window20kb/topGO/combined'
a<-read.csv(file.path(root,'display_selection_audit.csv'));a<-a[a$analysis=='without_simple' & a$display,];a<-a[order(a$elim_P),];a$Term<-factor(a$Term,levels=rev(a$Term))
# Retain percentage-valued coordinates; define scaling explicitly in caption.
wrap<-function(x)vapply(x,function(t)paste(strwrap(t,36),collapse='\n'),'')
th<-theme_bw(base_size=17)+theme(axis.text=element_text(color='black',size=15),panel.grid.minor=element_blank(),panel.grid.major=element_line(color='#E8E8E8'),plot.title=element_text(size=17,face='bold'),plot.margin=margin(12,12,12,12))
p1<-ggplot(a,aes(gene_ratio_percent,Term))+geom_point(aes(size=Significant,color=-log10(elim_P)))+
 scale_size_continuous(name='Count',range=c(3,9),breaks=c(2,5,10,15))+
 scale_colour_gradientn(name=expression(-log[10](P)),colours=c('#65D226','#C8B600','#F17A13','#DC2020'))+
 scale_y_discrete(labels=wrap)+scale_x_continuous(breaks=c(0,.5,1,1.5),limits=c(0,1.7))+
 labs(title='A  Pathway & process enrichment',x='Gene ratio',y=NULL)+th+
 theme(legend.position='bottom',legend.box='vertical',legend.title=element_text(size=13),legend.text=element_text(size=12))
d<-read.csv(file.path(root,'selected_environmental_development_processes.csv'))
bar<-function(cat,title,col){z<-d[d$category==cat,];z<-z[order(-z$Significant,z$Term),];z$Term<-factor(z$Term,levels=rev(z$Term));ggplot(z,aes(Significant,Term))+geom_col(fill=col,width=.68)+scale_y_discrete(labels=wrap)+scale_x_continuous(expand=expansion(mult=c(0,.05)))+labs(title=title,x='Gene count',y=NULL)+th}
p2<-bar('Environmental response and protection','B  Environmental response','#267F8E')
p3<-bar('Development and seasonal timing','C  Development','#9B723E')
p<-p1+p2+p3+plot_layout(widths=c(.85,1,1))
ggsave(file.path(root,'Figure_S28_combined.png'),p,width=23,height=15,dpi=220)
