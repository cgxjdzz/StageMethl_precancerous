library(ChAMP)
library(plyr)
label=mapvalues(pd$label, c("Solid Tissue Normal",'intestinal metaplasia biopsy from gastric antrum', "Primary Tumor"),c("Normal",'precancer', "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)

pre_DMP=sum_DMP$Normal_to_precancer[which(abs(sum_DMP$Normal_to_precancer$logFC)>0.1),]
cancer_DMP=sum_DMP$Normal_to_Tumor[which(abs(sum_DMP$Normal_to_Tumor$logFC)>0.1),]
all_DMP=intersect(rownames(pre_DMP),rownames(cancer_DMP))
draw_data=sum_norm[all_DMP,]

all_DMP=all_DMP[order(pre_DMP[all_DMP,]$P.Value)][c(1:150)]
###列按照组排序
order=c(which(label=='Normal'),which(label=='precancer'),which(label=='Tumor'))
draw_data=draw_data[,order]
label=label[order]

###行按照FC值排序
pre_FC=pre_DMP[all_DMP,]$logFC
cancer_FC=cancer_DMP[all_DMP,]$logFC

pre_FC_up=all_DMP[which(pre_FC>0)]
pre_FC_down=all_DMP[which(pre_FC<0)]
pre_FC_up_cancerFC=cancer_FC[which(pre_FC>0)]
pre_FC_up=pre_FC_up[order(pre_FC_up_cancerFC)]
pre_FC_down_cancerFC=cancer_FC[which(pre_FC<0)]
pre_FC_down=pre_FC_down[order(pre_FC_down_cancerFC)]
####行按照肿瘤FC值排序
cancer_FC_up=all_DMP[which(cancer_FC>0)]
cancer_FC_down=all_DMP[which(cancer_FC<0)]
cancer_FC_up_preFC=pre_FC[which(cancer_FC>0)]
cancer_FC_up=cancer_FC_up[order(cancer_FC_up_preFC)]
cancer_FC_down_preFC=pre_FC[which(cancer_FC<0)]
cancer_FC_down=cancer_FC_down[order(cancer_FC_down_preFC)]



###fold change
logFC_color=colorRamp2(c(-0.5,0, 0.5), c("#00C6CF",'white', "#AAAF0D"))
FC_pre=pre_DMP[c(cancer_FC_down,cancer_FC_up),]$logFC
FC_cancer=cancer_DMP[c(cancer_FC_down,cancer_FC_up),]$logFC
logFCdf=data.frame(pre=FC_pre,cancer=FC_cancer)



order2=order(pre_FC,cancer_FC)
draw_data=draw_data[c(cancer_FC_down,cancer_FC_up),]

color=c("Normal" =  "#63BAAA", "precancer" = "#FCEAFF",'Tumor'='#5E0FDF')
library(ComplexHeatmap)
library(circlize)
pdf('0204gastric_heatmap.pdf')
nana=HeatmapAnnotation(type=label, col = list(type=color),name='type')
p=Heatmap(draw_data,show_row_names = F,top_annotation = nana,col=colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),name = 'beta',show_column_names = F,cluster_columns = F,column_dend_reorder = F,row_dend_reorder = F,cluster_rows = F,row_split = c(rep(c("decrease"), length(cancer_FC_down)),rep(c('increase'),length(cancer_FC_up))))+
  rowAnnotation(df = logFCdf, col = list(pre=logFC_color,cancer=logFC_color), 
                width = unit(2, "cm"))
p
dev.off()
save(all_DMP,cancer_FC_down,cancer_FC_up,file='0204heatmapDMP.Rdata')