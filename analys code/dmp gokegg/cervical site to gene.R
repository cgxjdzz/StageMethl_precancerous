
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/cervical_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/cervical_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/cervical_hyper_and_hypo.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)

DMP=sum_DMP$normal_to_cin3
DMP=DMP[which(DMP$feature %in% c('TSS200','TSS1500')),]
cervical_hyper_hyper_gene=unique(DMP[cervical_hyper_hyper,]$gene)
cervical_hyper_hypo_gene=unique(DMP[cervical_hyper_hypo,]$gene)
cervical_hypo_hyper_gene=unique(DMP[cervical_hypo_hyper,]$gene)
cervical_hypo_hypo_gene=unique(DMP[cervical_hypo_hypo,]$gene)
cervical_hyper_stable_gene=unique(DMP[cervical_hyper_stable,]$gene)
cervical_hypo_stable_gene=unique(DMP[cervical_hypo_stable,]$gene)
cervical_hypo_back_to_normal_gene=unique(DMP[cervical_hypo_back_to_normal,]$gene)
cervical_hyper_back_to_normal_gene=unique(DMP[cervical_hyper_back_to_normal,]$gene)
cervical_hypo_back_to_normal_gene_bitr=bitr(cervical_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hyper_back_to_normal_gene_bitr=bitr(cervical_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hyper_hyper_gene_bitr=bitr(cervical_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hyper_hypo_gene_bitr=bitr(cervical_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hypo_hyper_gene_bitr=bitr(cervical_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hypo_hypo_gene_bitr=bitr(cervical_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hyper_stable_gene_bitr=bitr(cervical_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hypo_stable_gene_bitr=bitr(cervical_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(cervical_hyper_hyper_gene_bitr,cervical_hyper_hypo_gene_bitr,cervical_hypo_hyper_gene_bitr,cervical_hypo_hypo_gene_bitr,cervical_hyper_stable_gene_bitr,cervical_hypo_stable_gene_bitr,file = 'cervical.Rdata')
