library(ggnewscale)
library(ggplot2)
load('1216_sum_norm.Rdata')
load('1109_GSE31835.Rdata')
pd$batch='GSE44661'
pd$batch[which(colnames(sum_norm) %in% colnames(beta_GSE31835))]='GSE31835'
mean=apply(sum_norm,2,mean)
data2=data.frame(type=pd$label,mean=mean,batch=pd$batch)

nana<-ggplot()+geom_jitter(alpha = 0.5,size=3,aes(x =data2$type, y =data2$mean ,color=data2$batch))+ #分组的点图
  scale_x_discrete(limit=c("skin punch biopsy","melanocytes","primary melanoma") )+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('melanoma')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =data2$type, y =data2$mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 45))+
  labs(x="sample",y="value",fill= "type")
nana
ggsave('0128melanoma_source_box_plot.pdf')





library(ggplot2)
library(magrittr)
library(WGCNA)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(GO.db)
library(dplyr)
library(LOLA)
library(GenomicRanges)
library(ChAMP)
pd$label=gsub(' ','_',pd$label)
sum_DMP=champ.DMP(sum_norm,pheno = pd$label)

sum_DMP$melanocytes_to_skin_punch_biopsy=sum_DMP$melanocytes_to_skin_punch_biopsy[which(abs(sum_DMP$melanocytes_to_skin_punch_biopsy$logFC)>0.1),]
sum_DMP$primary_melanoma_to_skin_punch_biopsy=sum_DMP$primary_melanoma_to_skin_punch_biopsy[which(abs(sum_DMP$primary_melanoma_to_skin_punch_biopsy$logFC)>0.1),]
melanocytes_to_skin_punch_biopsy_top50=sum_DMP$melanocytes_to_skin_punch_biopsy[order(sum_DMP$melanocytes_to_skin_punch_biopsy$P.Value),][c(1:50),]
primary_melanoma_to_skin_punch_biopsy_top50=sum_DMP$primary_melanoma_to_skin_punch_biopsy[order(sum_DMP$primary_melanoma_to_skin_punch_biopsy$P.Value),][c(1:50),]
melanocytes_to_skin_punch_biopsy=sum_norm[rownames(melanocytes_to_skin_punch_biopsy_top50),]
primary_melanoma_to_skin_punch_biopsy=sum_norm[rownames(primary_melanoma_to_skin_punch_biopsy_top50),]


pdf('melanoma_melanocytes_to_skin_punch_biopsy_heatmap_top50.pdf')
nana=HeatmapAnnotation(type=pd$label, col = list(type=c("skin_punch_biopsy" =  "#845EC2", "melanocytes" = "#4B4453",'primary_melanoma'='#B0A8B9')),name='type')
Heatmap(melanocytes_to_skin_punch_biopsy,show_row_names = F,top_annotation = nana,col=viridis::viridis(30),show_row_dend = F,name = 'beta',show_column_names = F)
dev.off()

pdf('melanoma_primary_melanoma_to_skin_punch_biopsy_heatmap_top50.pdf')
nana=HeatmapAnnotation(type=pd$label, col = list(type=c("skin_punch_biopsy" =  "#845EC2", "melanocytes" = "#4B4453",'primary_melanoma'='#B0A8B9')),name='type')
Heatmap(melanocytes_to_skin_punch_biopsy,show_row_names = F,top_annotation = nana,col=viridis::viridis(30),show_row_dend = F,name = 'beta',show_column_names = F)
dev.off()

