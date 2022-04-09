rm(list = ls())
load('D:/甲基化/癌症发展/colonic/139404+68060/10.12myNorm.Rdata')
load('D:/甲基化/癌症发展/colonic/139404+68060/10.12impute.Rdata')
library(ChAMP)
library(plyr)
sum_norm=myNorm
pd0204=impute$pd
pd0204$label=mapvalues(impute$pd$label, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = myNorm,pheno = pd0204$label)

normal_v_pre=sum_DMP$Precancer_to_Normal
normal_v_cancer=sum_DMP$Normal_to_Tumor

normal_v_cancer_back_to_normal=normal_v_cancer[which(abs(normal_v_cancer$logFC)<(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hyper=normal_v_pre[which(normal_v_pre$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_pre_hypo=normal_v_pre[which(normal_v_pre$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hyper=normal_v_cancer[which(normal_v_cancer$logFC>(0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
normal_v_cancer_hypo=normal_v_cancer[which(normal_v_cancer$logFC<(-0.1)),c('logFC','adj.P.Val','gene','feature','cgi')]
save(normal_v_pre_hyper,normal_v_pre_hypo,normal_v_cancer_hyper,normal_v_cancer_hypo,file = 'colonic_DMP.Rdata')

pre_to_cancer_stable=sum_DMP$Precancer_to_Tumor[which(abs(sum_DMP$Precancer_to_Tumor$logFC)<0.1),]
colonic_hyper_stable=intersect(rownames(normal_v_pre_hyper),rownames(pre_to_cancer_stable))
colonic_hypo_stable=intersect(rownames(normal_v_pre_hypo),rownames(pre_to_cancer_stable))
save(colonic_hyper_stable,colonic_hypo_stable,file = 'colonic_hy_and_stable.Rdata')


colonic_hyper_hyper=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hyper))
colonic_hyper_hypo=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_hypo))
colonic_hypo_hyper=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hyper))
colonic_hypo_hypo=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_hypo))
colonic_hyper_hyper=colonic_hyper_hyper[-which(colonic_hyper_hyper %in% colonic_hyper_stable)]
colonic_hypo_hypo=colonic_hypo_hypo[-which(colonic_hypo_hypo %in% colonic_hypo_stable)]
save(colonic_hyper_hyper,colonic_hypo_hyper,colonic_hyper_hypo,colonic_hypo_hypo,file='colonic_hyper_and_hypo.Rdata')


colonic_hypo_back_to_normal=intersect(rownames(normal_v_pre_hypo),rownames(normal_v_cancer_back_to_normal))
colonic_hyper_back_to_normal=intersect(rownames(normal_v_pre_hyper),rownames(normal_v_cancer_back_to_normal))
save(colonic_hypo_back_to_normal,colonic_hyper_back_to_normal,file = 'colonic_hy_and_back_to_normal.Rdata')


####一下话前50图heatmap

#colonic_hyper_hypo=colonic_hyper_hypo[order(normal_v_pre_hyper[colonic_hyper_hypo,]$adj.P.Val)][c(1:50)]
#colonic_hypo_hyper=colonic_hypo_hyper[order(normal_v_pre_hypo[colonic_hypo_hyper,]$adj.P.Val)][c(1:50)]
colonic_hyper_hyper=colonic_hyper_hyper[order(normal_v_cancer_hyper[colonic_hyper_hyper,]$adj.P.Val)][c(1:50)]
colonic_hypo_hypo=colonic_hypo_hypo[order(normal_v_cancer_hypo[colonic_hypo_hypo,]$adj.P.Val)][c(1:50)]
colonic_hyper_back_to_normal=colonic_hyper_back_to_normal[order(normal_v_pre_hyper[colonic_hyper_back_to_normal,]$adj.P.Val)][c(1:50)]
colonic_hypo_back_to_normal=colonic_hypo_back_to_normal[order(normal_v_pre_hypo[colonic_hypo_back_to_normal,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_hyper',length(colonic_hyper_hyper)),rep('hypo_hypo',length(colonic_hypo_hypo)),rep('hypo_hyper',length(colonic_hypo_hyper)),rep('hyper_hypo',length(colonic_hyper_hypo)))

all_DMP=c(colonic_hyper_hyper,colonic_hypo_hypo,colonic_hypo_hyper,colonic_hyper_hypo)
draw_data=sum_norm[all_DMP,]
label=pd0204$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Precancer'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]



library(ComplexHeatmap)
library(circlize)




color=c("Normal" =  "#BAB2BA", "Precancer" = "#AA9E92",'Tumor'='#4D3D51')
sitecolor=c('hyper_hyper'='#CA4528','hypo_hypo'='#4D79A6','hypo_hyper'='#EFBA51','hyper_hypo'='#9EC6BE')
pdf('0212colonic_heatmap.pdf')
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

site=c(rep('hyper_back_to_normal',length(colonic_hyper_back_to_normal)),rep('hypo_back_to_normal',length(colonic_hypo_back_to_normal)))
sitecolor=c('hyper_back_to_normal'='#CA4528','hypo_back_to_normal'='#4D79A6')
all_DMP=c(colonic_hyper_back_to_normal,colonic_hypo_back_to_normal)
draw_data=sum_norm[all_DMP,]
label=pd0204$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Precancer'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]
pdf('0212colonic_back_to_normal_heatmap.pdf')
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


colonic_hyper_stable=colonic_hyper_stable[order(normal_v_pre_hyper[colonic_hyper_stable,]$adj.P.Val)][c(1:50)]
colonic_hypo_stable=colonic_hypo_stable[order(normal_v_pre_hypo[colonic_hypo_stable,]$adj.P.Val)][c(1:50)]



site=c(rep('hyper_stable',length(colonic_hyper_stable)),rep('hypo_stable',length(colonic_hypo_stable)))
sitecolor=c('hyper_stable'='#CA4528','hypo_stable'='#4D79A6')
all_DMP=c(colonic_hyper_stable,colonic_hypo_stable)
draw_data=sum_norm[all_DMP,]
label=pd0204$label
###列按照组排序
order=c(which(label=='Normal'),which(label=='Precancer'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]
pdf('0212colonic_stable_heatmap.pdf')
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

