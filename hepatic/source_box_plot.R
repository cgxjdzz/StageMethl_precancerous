library(ggnewscale)
library(ggplot2)
require(GEOquery)
require(Biobase)
mean=apply(sum_norm,2,mean)
pd$branch='UCSC-TCGA'
load('eset.Rdata')
pD.all <- pData(eset[[1]])
pd$branch[which(pd$name %in% pD.all$geo_accession)]='GSE157973'
data2=data.frame(type=pd$label,mean=mean,batch=pd$branch)
data2$cancer='liver'
liver_box=data2
save(liver_box,file = 'liver_box.Rdata')
nana<-ggplot()+geom_jitter(alpha = 0.3,size=3,aes(x =data2$type, y =data2$mean ,color=data2$batch))+ #分组的点图
  scale_x_discrete(limit=c("Normal","Cirrhosis","Tumor") )+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('liver')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =data2$type, y =data2$mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 45))+
  labs(x="sample",y="value",fill= "type")

ggsave('liver_source_box_plot.pdf')



nana
