site_list=c("union_site_hyper_back_to_normal","union_site_hyper_hyper","union_site_hyper_hypo","union_site_hyper_stable","union_site_hypo_back_to_normal","union_site_hypo_hyper","union_site_hypo_hypo","union_site_hypo_stable" )
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChAMP)
library(plyr)
library(ggbreak) 
library(reshape2)
genefunc=function(i){
  data=read.table(paste0(i,'/',i,'.csv'))
  three=length(which(data$sum==3))
  two=length(which(data$sum==2))
  one=length(which(data$sum==1))
  return(c(three,two,one))
}

cgdata=data.frame(matrix(nrow=3))
for (i in site_list){
  assign(i,genefunc(i)) 
  eval(parse(text = paste0('cgdata$',i,'=get(i)')))
}
cgdata=cgdata[,-1]
plot_data=melt(cgdata)
plot_data$count=rep(c('2','2','1'),8)
library(ggplot2)
color=c('2'='#40D8DD','1'='#FF9595')
ggplot(plot_data,aes(variable,weight=value,fill=count))+ theme_bw()+scale_fill_manual(values = color)+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=15))+ labs(x="site",y="site_count")+ theme(legend.position = "bottom" ,legend.box = "horizontal")+scale_y_break(c(750, 3000)) 
ggsave('0328union_site.pdf',width = 9,height =15)


genefunc=function(i){
  data=read.table(paste0(i,'/',i,'.csv'))
  cgpoint=rownames(data) 
  gene=unique(na.omit(probe.features[cgpoint,]$gene))
  return(length(gene))
}
genedata=data.frame(matrix(nrow=1))
for (i in site_list){
  assign(i,genefunc(i)) 
  eval(parse(text = paste0('genedata$',i,'=get(i)')))
}
genedata=genedata[,-1]
plot_data=melt(genedata)
ggplot(plot_data,aes(variable,weight=value,fill='#135453'))+ theme_bw()+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=15))+ labs(x="site",y="gene_count")+ theme(legend.position = "bottom" ,legend.box = "horizontal") 
ggsave('0328union_gene.pdf',width = 9,height = 7)

