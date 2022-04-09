library(ggnewscale)
library(ggplot2)
mean=apply(sum_norm,2,mean)
data2=data.frame(type=sum_pd$label,mean=mean,batch=sum_pd$branch)
data2$cancer='colonic'
colonic_box=data2
save(colonic_box,file = 'colonic_box.Rdata')
nana<-ggplot()+geom_jitter(alpha = 0.3,size=3,aes(x =data2$type, y =data2$mean ,color=data2$batch))+ #分组的点图
  scale_x_discrete(limit=c("Normal","LGA","HGA",'Tumor') )+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('colonic')+#图标题
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
ggsave('0128colonic_source_box_plot.pdf')

source('zffdatamining_gene.R')





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


load('10.12myDMP.Rdata')
load('10.12myNorm.Rdata')
sum_DMP=myDMP
sum_combat=myNorm
sum_DMP$LGA_to_Normal=sum_DMP$LGA_to_Normal[which(abs(sum_DMP$LGA_to_Normal$logFC)>0.1),]
sum_DMP$Normal_to_HGA=sum_DMP$Normal_to_HGA[which(abs(sum_DMP$Normal_to_HGA$logFC)>0.1),]
sum_DMP$Normal_to_Tumor=sum_DMP$Normal_to_Tumor[which(abs(sum_DMP$Normal_to_Tumor$logFC)>0.1),]

LGA_to_Normal_top50=sum_DMP$LGA_to_Normal[order(sum_DMP$LGA_to_Normal$P.Value),][c(1:50),]
Normal_to_HGA_top50=sum_DMP$Normal_to_HGAr[order(sum_DMP$Normal_to_HGA$P.Value),][c(1:50),]
Normal_to_Tumor_top50=sum_DMP$Normal_to_Tumor[order(sum_DMP$Normal_to_Tumor$P.Value),][c(1:50),]
LGA_to_Normal=sum_combat[rownames(LGA_to_Normal_top50),]
Normal_to_HGA=sum_combat[rownames(Normal_to_HGA_top50),]
Normal_to_Tumor=sum_combat[rownames(Normal_to_Tumor_top50),]

label=pd$label
order=c(which(label=='Normal'),which(label=='LGA'),which(label=='HGA'),which(label=='Tumor'))
LGA_to_Normal=LGA_to_Normal[,order]

label=label[order]
  


pdf('4_colonic_LGA_to_Normal_heatmap_top50.pdf')
nana=HeatmapAnnotation(type=label, col = list(type=c("Normal" =  "#63BAAA", "LGA" = "#FCEAFF",'HGA'='#8F78A9','Tumor'='#5E0FDF')),name='type')
Heatmap(LGA_to_Normal,show_row_names = F,top_annotation = nana,col=viridis::viridis(30),name = 'beta',show_column_names = F,cluster_columns = F,column_dend_reorder = F,row_dend_reorder = F)
dev.off()

pdf('colonic_Normal_to_HGA_heatmap_top50.pdf')
nana=HeatmapAnnotation(type=pd$label, col = list(type=c("normal" =  "#845EC2", "cin3" = "#4B4453",'cancer'='#B0A8B9')),name='type')
Heatmap(Normal_to_HGA,show_row_names = F,top_annotation = nana,col=viridis::viridis(30),show_row_dend = F,name = 'beta',show_column_names = F)
dev.off()

pdf('colonic_Normal_to_Tumor_heatmap_top50.pdf')
nana=HeatmapAnnotation(type=pd$label, col = list(type=c("normal" =  "#845EC2", "cin3" = "#4B4453",'cancer'='#B0A8B9')),name='type')
Heatmap(Normal_to_Tumor,show_row_names = F,top_annotation = nana,col=viridis::viridis(30),show_row_dend = F,name = 'beta',show_column_names = F)
dev.off()

normal_norm=sum_norm[,which(label=='Normal')]
normal_mean=apply(normal_norm,1,mean)
sum_norm=myNorm
LGA_to_Normal_sum_norm=sum_norm[which(normal_mean>0.15),]
LGA_to_Normal_sum_DMP=sum_DMP$LGA_to_Normal[which(rownames(sum_DMP$LGA_to_Normal) %in% rownames(LGA_to_Normal_sum_norm)),]
ge.heatmap(LGA_to_Normal_sum_DMP,sum_norm =LGA_to_Normal_sum_norm,order=1,
                 filename = 'D:/甲基化/癌症发展/结肠癌/139404+癌68060/colonic_LGA_to_Normal_heatmap_top100_delete.pdf',label = pd$label,color = c("Normal" =  "#63BAAA", "LGA" = "#FCEAFF",'HGA'='#8F78A9','Tumor'='#5E0FDF'))
