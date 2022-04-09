
load('D:/甲基化/癌症发展/0331改数据/黑色素瘤/0331_sum_DMP.Rdata')

library(ChAMP)

load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/melanoma_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/melanoma_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/melanoma_hyper_and_hypo.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)

DMP=sum_DMP$Precancer_to_Normal
# DMP=DMP[which(DMP$feature %in% c("3'UTR","5'UTR")),]
melanoma_hyper_hyper_gene=unique(DMP[melanoma_hyper_hyper,]$gene)
melanoma_hyper_hypo_gene=unique(DMP[melanoma_hyper_hypo,]$gene)
melanoma_hypo_hyper_gene=unique(DMP[melanoma_hypo_hyper,]$gene)
melanoma_hypo_hypo_gene=unique(DMP[melanoma_hypo_hypo,]$gene)
melanoma_hyper_stable_gene=unique(DMP[melanoma_hyper_stable,]$gene)
melanoma_hypo_stable_gene=unique(DMP[melanoma_hypo_stable,]$gene)
melanoma_hypo_back_to_normal_gene=unique(DMP[melanoma_hypo_back_to_normal,]$gene)
melanoma_hyper_back_to_normal_gene=unique(DMP[melanoma_hyper_back_to_normal,]$gene)
melanoma_hypo_back_to_normal_gene_bitr=bitr(melanoma_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hyper_back_to_normal_gene_bitr=bitr(melanoma_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hyper_hyper_gene_bitr=bitr(melanoma_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hyper_hypo_gene_bitr=bitr(melanoma_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hypo_hyper_gene_bitr=bitr(melanoma_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hypo_hypo_gene_bitr=bitr(melanoma_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hyper_stable_gene_bitr=bitr(melanoma_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hypo_stable_gene_bitr=bitr(melanoma_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(melanoma_hyper_hyper_gene_bitr,melanoma_hyper_hypo_gene_bitr,melanoma_hypo_hyper_gene_bitr,melanoma_hypo_hypo_gene_bitr,melanoma_hyper_stable_gene_bitr,melanoma_hypo_stable_gene_bitr,file = 'melanoma.Rdata')

