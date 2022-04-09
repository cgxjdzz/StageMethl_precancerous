rm(list = ls())
load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')

library(ChAMP)
library(plyr)

pd0204=data.frame(name=colnames(sum_norm),label=pd$`sum_impute$pd`)
pd0204$label=mapvalues(pd0204$label, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = sum_norm,pheno = pd0204$label)

normal_v_pre=sum_DMP$Normal_to_Precancer
normal_v_cancer=sum_DMP$Tumor_to_Normal

normal_v_cancer_back_to_normal=normal_v_cancer[which(abs(normal_v_cancer$logFC)<(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hyper=normal_v_pre[which(normal_v_pre$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hypo=normal_v_pre[which(normal_v_pre$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hyper=normal_v_cancer[which(normal_v_cancer$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hypo=normal_v_cancer[which(normal_v_cancer$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
save(normal_v_pre_hyper,normal_v_pre_hypo,normal_v_cancer_hyper,normal_v_cancer_hypo,file = 'colonic_DMP.Rdata')

pre_to_cancer_stable=sum_DMP$Normal_to_Precancer[which(abs(sum_DMP$Normal_to_Precancer$logFC)<0.1),]
colonic_hyper_stable=intersect(rownames(normal_v_pre_hyper),rownames(pre_to_cancer_stable))
colonic_hypo_stable=intersect(rownames(normal_v_pre_hypo),rownames(pre_to_cancer_stable))
save(colonic_hyper_stable,colonic_hypo_stable,file = 'colonic_hy_and_stable.Rdata')


colonic_hyper_hyper=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hyper))
colonic_hyper_hypo=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hypo))
colonic_hypo_hyper=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hyper))
colonic_hypo_hypo=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hypo))
#colonic_hyper_hyper=colonic_hyper_hyper[-which(colonic_hyper_hyper %in% colonic_hyper_stable)]
#colonic_hypo_hypo=colonic_hypo_hypo[-which(colonic_hypo_hypo %in% colonic_hypo_stable)]
save(colonic_hyper_hyper,colonic_hypo_hyper,colonic_hyper_hypo,colonic_hypo_hypo,file='colonic_hyper_and_hypo.Rdata')


colonic_hypo_back_to_normal=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_back_to_normal))
colonic_hyper_back_to_normal=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_back_to_normal))
save(colonic_hypo_back_to_normal,colonic_hyper_back_to_normal,file = 'colonic_hy_and_back_to_normal.Rdata')

label=pd0204$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Precancer'),which(label=='Tumor'))
sum_norm=sum_norm[,order]
pd0204=pd0204[order,]
####一下话前50图heatmap

colonic_hyper_hypo=colonic_hyper_hypo[order(normal_v_pre_hyper[colonic_hyper_hypo,]$adj.P.Val)]
colonic_hypo_hyper=colonic_hypo_hyper[order(normal_v_pre_hypo[colonic_hypo_hyper,]$adj.P.Val)]
#colonic_hyper_hyper=colonic_hyper_hyper[order(normal_v_pre_hyper[colonic_hyper_hyper,]$adj.P.Val)][c(1:20)]
colonic_hypo_hypo=colonic_hypo_hypo[order(normal_v_pre_hypo[colonic_hypo_hypo,]$adj.P.Val)][c(1:20)]
colonic_hyper_back_to_normal=colonic_hyper_back_to_normal[order(normal_v_pre_hyper[colonic_hyper_back_to_normal,]$adj.P.Val)]
colonic_hypo_back_to_normal=colonic_hypo_back_to_normal[order(normal_v_pre_hypo[colonic_hypo_back_to_normal,]$adj.P.Val)][c(1:20)]
# colonic_hypo_stable=colonic_hypo_stable[order(normal_v_pre_hypo[colonic_hypo_stable,]$adj.P.Val)][c(1:20)]
# colonic_hyper_stable=colonic_hyper_stable[order(normal_v_pre_hypo[colonic_hyper_stable,]$adj.P.Val)][c(1:20)]
identical(colnames(sum_norm),pd0204$name)

# site1=c(rep('hyper_hyper',length(colonic_hyper_hyper)),rep('hypo_hypo',length(colonic_hypo_hypo)),rep('hypo_hyper',length(colonic_hypo_hyper)),rep('hyper_hypo',length(colonic_hyper_hypo)))
# 
# all_DMP=c(colonic_hyper_hyper,colonic_hypo_hypo,colonic_hypo_hyper,colonic_hyper_hypo)
# draw_data1=sum_norm[all_DMP,]

library(ComplexHeatmap)
library(circlize)


# 
# 
color=c("Normal" =  "#BAB2BA", "Precancer" = "#AA9E92",'Tumor'='#4D3D51')
# sitecolor1=c('hyper_hyper'='#CA4528','hypo_hypo'='#4D79A6','hypo_hyper'='#EFBA51','hyper_hypo'='#9EC6BE')
# site2=c(rep('hyper_back_to_normal',length(colonic_hyper_back_to_normal)),rep('hypo_back_to_normal',length(colonic_hypo_back_to_normal)))
# sitecolor2=c('hyper_back_to_normal'='#CA4528','hypo_back_to_normal'='#4D79A6')
# all_DMP=c(colonic_hyper_back_to_normal,colonic_hypo_back_to_normal)
# draw_data2=sum_norm[all_DMP,]
# colonic_hyper_stable=colonic_hyper_stable[order(normal_v_pre_hyper[colonic_hyper_stable,]$adj.P.Val)]
# colonic_hypo_stable=colonic_hypo_stable[order(normal_v_pre_hypo[colonic_hypo_stable,]$adj.P.Val)]
# 
# site3=c(rep('hyper_stable',length(colonic_hyper_stable)),rep('hypo_stable',length(colonic_hypo_stable)))
# sitecolor3=c('hyper_stable'='#CA4528','hypo_stable'='#4D79A6')
# all_DMP=c(colonic_hyper_stable,colonic_hypo_stable)
# draw_data3=sum_norm[all_DMP,]
all_DMP=c(colonic_hyper_hyper,colonic_hyper_stable,colonic_hyper_back_to_normal,colonic_hyper_hypo,colonic_hypo_hypo,colonic_hypo_stable,colonic_hypo_back_to_normal,colonic_hypo_hyper)
draw_data=sum_norm[all_DMP,]
site=c(rep('hyper_hyper',length(colonic_hyper_hyper)),rep('hyper_stable',length(colonic_hyper_stable)),rep('hyper_back_to_normal',length(colonic_hyper_back_to_normal)),rep('hyper_hypo',length(colonic_hyper_hypo)),rep('hypo_hypo',length(colonic_hypo_hypo)),rep('hypo_stable',length(colonic_hypo_stable)),rep('hypo_back_to_normal',length(colonic_hypo_back_to_normal)),rep('hypo_hyper',length(colonic_hypo_hyper)))
sitecolor=c('hyper_hyper'='#E39A44','hypo_hypo'='#47864E','hypo_hyper'='#3290CF','hyper_hypo'='#C04C36',
            'hyper_back_to_normal'='#E8D2B7','hypo_back_to_normal'='#C1E5DF','hyper_stable'='#D3C856','hypo_stable'='#4CAC6E')
ee=HeatmapAnnotation(type=site,col = list(type=sitecolor),name='site',show_annotation_name=F,which = 'row',annotation_label=NULL, width = unit(1, "cm"),annotation_name_gp=gpar(fontsize = 10))
nana=HeatmapAnnotation(type=pd0204$label, col = list(type=color),name='type',show_annotation_name=F)
pdf('0320 8site_colonic_heatmap.pdf')
p=Heatmap(draw_data,show_row_names = F,top_annotation = nana,
          col=colorRamp2(c(0, 0.5, 1), c("#89C762", "white", "#F6AE60")),
          
          name = 'colonic',show_column_names = F,cluster_columns = F,column_dend_reorder = F,
          row_dend_reorder = F,cluster_rows = F)

ee+p
dev.off()

