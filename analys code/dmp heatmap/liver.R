rm(list = ls())
load('D:/甲基化/癌症发展/肝癌/0227_sum_norm.Rdata')
load('D:/甲基化/癌症发展/肝癌/sum_DMP.Rdata')

library(ChAMP)
library(plyr)

normal_v_pre=sum_DMP$Cirrhosis_to_Normal
normal_v_cancer=sum_DMP$Tumor_to_Normal

normal_v_cancer_return_to_normal=normal_v_cancer[which(abs(normal_v_cancer$logFC)<(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hyper=normal_v_pre[which(normal_v_pre$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hypo=normal_v_pre[which(normal_v_pre$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hyper=normal_v_cancer[which(normal_v_cancer$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hypo=normal_v_cancer[which(normal_v_cancer$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
save(normal_v_pre_hyper,normal_v_pre_hypo,normal_v_cancer_hyper,normal_v_cancer_hypo,file = 'liver_DMP.Rdata')


pre_to_cancer_stable=sum_DMP$Cirrhosis_to_Tumor[which(abs(sum_DMP$Cirrhosis_to_Tumor$logFC)<0.1),]
liver_hyper_stable=intersect(rownames(normal_v_pre_hyper),rownames(pre_to_cancer_stable))
liver_hypo_stable=intersect(rownames(normal_v_pre_hypo),rownames(pre_to_cancer_stable))
save(liver_hyper_stable,liver_hypo_stable,file = 'liver_hy_and_stable.Rdata')

liver_hyper_hyper=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hyper))
liver_hyper_hypo=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hypo))
liver_hypo_hyper=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hyper))
liver_hypo_hypo=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hypo))
liver_hyper_hyper=liver_hyper_hyper[-which(liver_hyper_hyper %in% liver_hyper_stable)]
liver_hypo_hypo=liver_hypo_hypo[-which(liver_hypo_hypo %in% liver_hypo_stable)]
save(liver_hyper_hyper,liver_hypo_hyper,liver_hyper_hypo,liver_hypo_hypo,file='liver_hyper_and_hypo.Rdata')


liver_hypo_return_to_normal=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_return_to_normal))
liver_hyper_return_to_normal=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_return_to_normal))
save(liver_hypo_return_to_normal,liver_hyper_return_to_normal,file = 'liver_hy_and_return_to_normal.Rdata')


####一下话前50图heatmap

liver_hyper_hypo=liver_hyper_hypo[order(normal_v_pre_hyper[liver_hyper_hypo,]$adj.P.Val)][c(1:50)]
liver_hypo_hyper=liver_hypo_hyper[order(normal_v_pre_hypo[liver_hypo_hyper,]$adj.P.Val)][c(1:50)]
liver_hyper_hyper=liver_hyper_hyper[order(normal_v_pre_hyper[liver_hyper_hyper,]$adj.P.Val)][c(1:50)]
liver_hypo_hypo=liver_hypo_hypo[order(normal_v_pre_hypo[liver_hypo_hypo,]$adj.P.Val)][c(1:50)]
liver_hyper_return_to_normal=liver_hyper_return_to_normal[order(normal_v_pre_hyper[liver_hyper_return_to_normal,]$adj.P.Val)][c(1:50)]
liver_hypo_return_to_normal=liver_hypo_return_to_normal[order(normal_v_pre_hypo[liver_hypo_return_to_normal,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_hyper',length(liver_hyper_hyper)),rep('hypo_hypo',length(liver_hypo_hypo)),rep('hypo_hyper',length(liver_hypo_hyper)),rep('hyper_hypo',length(liver_hyper_hypo)))

all_DMP=c(liver_hyper_hyper,liver_hypo_hypo,liver_hypo_hyper,liver_hyper_hypo)
draw_data=sum_norm[all_DMP,]
label=pd$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Cirrhosis'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]
library(ComplexHeatmap)
library(circlize)




color=c("Normal" =  "#BAB2BA", "Cirrhosis" = "#AA9E92",'Tumor'='#4D3D51')
sitecolor=c('hyper_hyper'='#CA4528','hypo_hypo'='#4D79A6','hypo_hyper'='#EFBA51','hyper_hypo'='#9EC6BE')
pdf('0212liver_heatmap.pdf')
ee=HeatmapAnnotation(type=site,col = list(type=sitecolor),name='site')
nana=HeatmapAnnotation(type=label, col = list(type=color),name='type')
p=Heatmap(draw_data,show_row_names = F,top_annotation = nana,
          col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
          row_order=c(1:length(site)),
          name = 'beta',show_column_names = F,cluster_columns = F,column_dend_reorder = F,
          row_dend_reorder = F,cluster_rows = F,left_annotation = 
            rowAnnotation(site = site, col = list(site=sitecolor)),row_split = (site))
p
dev.off()

site=c(rep('hyper_return_to_normal',length(liver_hyper_return_to_normal)),rep('hypo_return_to_normal',length(liver_hypo_return_to_normal)))
sitecolor=c('hyper_return_to_normal'='#CA4528','hypo_return_to_normal'='#4D79A6')
all_DMP=c(liver_hyper_return_to_normal,liver_hypo_return_to_normal)
draw_data=sum_norm[all_DMP,]
label=pd$label
###列按照组排序
draw_data=draw_data[,order]
label=label[order]
pdf('0212liver_return_to_normal_heatmap.pdf')
ee=HeatmapAnnotation(type=site,col = list(type=sitecolor),name='site')
nana=HeatmapAnnotation(type=label, col = list(type=color),name='type')
p=Heatmap(draw_data,show_row_names = F,top_annotation = nana,
          col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
          row_order=c(1:length(site)),
          name = 'beta',show_column_names = F,cluster_columns = F,column_dend_reorder = F,
          row_dend_reorder = F,cluster_rows = F,left_annotation = 
            rowAnnotation(site = site, col = list(site=sitecolor)),row_split = (site))
p
dev.off()






liver_hyper_stable=liver_hyper_stable[order(normal_v_pre_hyper[liver_hyper_stable,]$adj.P.Val)][c(1:50)]
liver_hypo_stable=liver_hypo_stable[order(normal_v_pre_hypo[liver_hypo_stable,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_stable',length(liver_hyper_stable)),rep('hypo_stable',length(liver_hypo_stable)))
sitecolor=c('hyper_stable'='#CA4528','hypo_stable'='#4D79A6')
all_DMP=c(liver_hyper_stable,liver_hypo_stable)
draw_data=sum_norm[all_DMP,]
label=pd$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Cirrhosis'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]
pdf('0212liver_stable_heatmap.pdf')
ee=HeatmapAnnotation(type=site,col = list(type=sitecolor),name='site')
nana=HeatmapAnnotation(type=label, col = list(type=color),name='type')
p=Heatmap(draw_data,show_row_names = F,top_annotation = nana,
          col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
          row_order=c(1:length(site)),
          name = 'beta',show_column_names = F,cluster_columns = F,column_dend_reorder = F,
          row_dend_reorder = F,cluster_rows = F,left_annotation = 
            rowAnnotation(site = site, col = list(site=sitecolor)),row_split = (site))
p
dev.off()

