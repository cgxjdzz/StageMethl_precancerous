
rm(list = ls())
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
pd$label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,pd$label)
library(plyr)

normal_v_pre=sum_DMP$CAF_to_Normal
normal_v_cancer=sum_DMP$Tumor_to_Normal

normal_v_cancer_back_to_normal=normal_v_cancer[which(abs(normal_v_cancer$logFC)<(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hyper=normal_v_pre[which(normal_v_pre$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hypo=normal_v_pre[which(normal_v_pre$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hyper=normal_v_cancer[which(normal_v_cancer$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hypo=normal_v_cancer[which(normal_v_cancer$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
save(normal_v_pre_hyper,normal_v_pre_hypo,normal_v_cancer_hyper,normal_v_cancer_hypo,file = 'prostatic_DMP.Rdata')


pre_to_cancer_stable=sum_DMP$CAF_to_Tumor[which(abs(sum_DMP$CAF_to_Tumor$logFC)<0.1),]
prostatic_hyper_stable=intersect(rownames(normal_v_pre_hyper),rownames(pre_to_cancer_stable))
prostatic_hypo_stable=intersect(rownames(normal_v_pre_hypo),rownames(pre_to_cancer_stable))
save(prostatic_hyper_stable,prostatic_hypo_stable,file = 'prostatic_hy_and_stable.Rdata')


prostatic_hyper_hyper=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hyper))
prostatic_hyper_hypo=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hypo))
prostatic_hypo_hyper=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hyper))
prostatic_hypo_hypo=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hypo))
prostatic_hyper_hyper=prostatic_hyper_hyper[-which(prostatic_hyper_hyper %in% prostatic_hyper_stable)]
prostatic_hypo_hypo=prostatic_hypo_hypo[-which(prostatic_hypo_hypo %in% prostatic_hypo_stable)]
save(prostatic_hyper_hyper,prostatic_hypo_hyper,prostatic_hyper_hypo,prostatic_hypo_hypo,file='prostatic_hyper_and_hypo.Rdata')


prostatic_hypo_back_to_normal=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_back_to_normal))
prostatic_hyper_back_to_normal=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_back_to_normal))
save(prostatic_hypo_back_to_normal,prostatic_hyper_back_to_normal,file = 'prostatic_hy_and_back_to_normal.Rdata')


####一下话前50图heatmap

prostatic_hyper_hypo=prostatic_hyper_hypo[order(normal_v_pre_hyper[prostatic_hyper_hypo,]$adj.P.Val)][c(1:50)]
prostatic_hypo_hyper=prostatic_hypo_hyper[order(normal_v_pre_hypo[prostatic_hypo_hyper,]$adj.P.Val)][c(1:50)]
prostatic_hyper_hyper=prostatic_hyper_hyper[order(normal_v_pre_hyper[prostatic_hyper_hyper,]$adj.P.Val)][c(1:50)]
prostatic_hypo_hypo=prostatic_hypo_hypo[order(normal_v_pre_hypo[prostatic_hypo_hypo,]$adj.P.Val)][c(1:50)]
prostatic_hyper_back_to_normal=prostatic_hyper_back_to_normal[order(normal_v_pre_hyper[prostatic_hyper_back_to_normal,]$adj.P.Val)][c(1:50)]
prostatic_hypo_back_to_normal=prostatic_hypo_back_to_normal[order(normal_v_pre_hypo[prostatic_hypo_back_to_normal,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_hyper',length(prostatic_hyper_hyper)),rep('hypo_hypo',length(prostatic_hypo_hypo)),rep('hypo_hyper',length(prostatic_hypo_hyper)),rep('hyper_hypo',length(prostatic_hyper_hypo)))

all_DMP=c(prostatic_hyper_hyper,prostatic_hypo_hypo,prostatic_hypo_hyper,prostatic_hyper_hypo)
draw_data=sum_norm[all_DMP,]

###列按照组排序
order=c(which(label=='Normal'),which(label=='CAF'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]

####行按照肿瘤FC值排序
#cancer_FC_up=all_DMP[which(cancer_FC<0)]
#cancer_FC_down=all_DMP[which(cancer_FC>0)]

#按beta值排序
#mean_up=normal_mean[cancer_FC_up]
#mean_down=normal_mean[cancer_FC_down]


#cancer_FC_up_preFC=pre_FC[which(cancer_FC<0)]
#cancer_FC_up=cancer_FC_up[order(cancer_FC_up_preFC)]
#cancer_FC_down_preFC=pre_FC[which(cancer_FC>0)]
#cancer_FC_down=cancer_FC_down[order(cancer_FC_down_preFC)]


###fold change
library(ComplexHeatmap)
library(circlize)




color=c("Normal" =  "#BAB2BA", "CAF" = "#AA9E92",'Tumor'='#4D3D51')
sitecolor=c('hyper_hyper'='#CA4528','hypo_hypo'='#4D79A6','hypo_hyper'='#EFBA51','hyper_hypo'='#9EC6BE')
pdf('0212prostatic_heatmap.pdf')
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

site=c(rep('hyper_back_to_normal',length(prostatic_hyper_back_to_normal)),rep('hypo_back_to_normal',length(prostatic_hypo_back_to_normal)))
sitecolor=c('hyper_back_to_normal'='#CA4528','hypo_back_to_normal'='#4D79A6')
all_DMP=c(prostatic_hyper_back_to_normal,prostatic_hypo_back_to_normal)
draw_data=sum_norm[all_DMP,]

###列按照组排序
draw_data=draw_data[,order]
label=label[order]
pdf('0212prostatic_back_to_normal_heatmap.pdf')
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






prostatic_hyper_stable=prostatic_hyper_stable[order(normal_v_pre_hyper[prostatic_hyper_stable,]$adj.P.Val)][c(1:50)]
prostatic_hypo_stable=prostatic_hypo_stable[order(normal_v_pre_hypo[prostatic_hypo_stable,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_stable',length(prostatic_hyper_stable)),rep('hypo_stable',length(prostatic_hypo_stable)))
sitecolor=c('hyper_stable'='#CA4528','hypo_stable'='#4D79A6')
all_DMP=c(prostatic_hyper_stable,prostatic_hypo_stable)
draw_data=sum_norm[all_DMP,]
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
###列按照组排序
order=c(which(label=='Normal'),which(label=='CAF'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]
pdf('0212prostatic_stable_heatmap.pdf')
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

