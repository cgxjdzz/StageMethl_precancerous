mean=apply(sum_norm,2,mean)
#delet=which(pd$label=='Primary Tumor' & mean<=0.35)
#mean=mean[-delet]
data1=data.frame(mean=mean,type=pd$label)

library(ggplot2)
plot.boxplot <- function(data,x,y,type,filename,title="boxplot"){
  a <- ggplot(data=data, aes(x =x, y =y ,color=type,group=type)) +
    geom_jitter(alpha = 0.3,size=3) +
    geom_boxplot(alpha = .5,size=1)+
    labs(x="sample",y="value",fill= "type")+
    ggtitle(title)+
    theme_bw() + 
    theme(panel.border = element_blank())+
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(panel.grid =element_blank())+  
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    scale_x_discrete(limit=c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum","Primary Tumor") )+
    scale_color_manual(limits=c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum","Primary Tumor"), values=c("#85B22E","#5F80B4","#E29827"))
  ggsave(paste0(filename, ".pdf"),plot=a,width=8,height=8)
}
plot.boxplot( data1,data1$type,data1$mean,data1$type,paste0("gastric cancer", "_boxplot"),title='gastric cancer')
