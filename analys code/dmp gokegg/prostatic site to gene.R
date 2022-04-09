
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/prostatic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/prostatic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/prostatic_hyper_and_hypo.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)

DMP=sum_DMP$CAF_to_Normal
DMP=DMP[which(DMP$feature %in% c('TSS200','TSS1500')),]
prostatic_hyper_hyper_gene=unique(DMP[prostatic_hyper_hyper,]$gene)
prostatic_hyper_hypo_gene=unique(DMP[prostatic_hyper_hypo,]$gene)
prostatic_hypo_hyper_gene=unique(DMP[prostatic_hypo_hyper,]$gene)
prostatic_hypo_hypo_gene=unique(DMP[prostatic_hypo_hypo,]$gene)
prostatic_hyper_stable_gene=unique(DMP[prostatic_hyper_stable,]$gene)
prostatic_hypo_stable_gene=unique(DMP[prostatic_hypo_stable,]$gene)
prostatic_hypo_back_to_normal_gene=unique(DMP[prostatic_hypo_back_to_normal,]$gene)
prostatic_hyper_back_to_normal_gene=unique(DMP[prostatic_hyper_back_to_normal,]$gene)
prostatic_hypo_back_to_normal_gene_bitr=bitr(prostatic_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hyper_back_to_normal_gene_bitr=bitr(prostatic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hyper_hyper_gene_bitr=bitr(prostatic_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hyper_hypo_gene_bitr=bitr(prostatic_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hypo_hyper_gene_bitr=bitr(prostatic_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hypo_hypo_gene_bitr=bitr(prostatic_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hyper_stable_gene_bitr=bitr(prostatic_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hypo_stable_gene_bitr=bitr(prostatic_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(prostatic_hyper_hyper_gene_bitr,prostatic_hyper_hypo_gene_bitr,prostatic_hypo_hyper_gene_bitr,prostatic_hypo_hypo_gene_bitr,prostatic_hyper_stable_gene_bitr,prostatic_hypo_stable_gene_bitr,file = 'prostatic.Rdata')
