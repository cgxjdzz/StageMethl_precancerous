
load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/colonic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/colonic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/colonic_hyper_and_hypo.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
myNorm=sum_norm
pd=pd$`sum_impute$pd`
pd=mapvalues(pd, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = myNorm,pheno = pd)
DMP=sum_DMP$Normal_to_Precancer
DMP=DMP[which(DMP$feature %in% c('TSS200','TSS1500')),]
colonic_hyper_hyper_gene=unique(DMP[colonic_hyper_hyper,]$gene)
colonic_hyper_hypo_gene=unique(DMP[colonic_hyper_hypo,]$gene)
colonic_hypo_hyper_gene=unique(DMP[colonic_hypo_hyper,]$gene)
colonic_hypo_hypo_gene=unique(DMP[colonic_hypo_hypo,]$gene)
colonic_hyper_stable_gene=unique(DMP[colonic_hyper_stable,]$gene)
colonic_hypo_stable_gene=unique(DMP[colonic_hypo_stable,]$gene)
colonic_hypo_back_to_normal_gene=unique(DMP[colonic_hypo_back_to_normal,]$gene)
colonic_hyper_back_to_normal_gene=unique(DMP[colonic_hyper_back_to_normal,]$gene)
colonic_hypo_back_to_normal_gene_bitr=bitr(colonic_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hyper_back_to_normal_gene_bitr=bitr(colonic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hyper_hyper_gene_bitr=bitr(colonic_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hyper_hypo_gene_bitr=bitr(colonic_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hypo_hyper_gene_bitr=bitr(colonic_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hypo_hypo_gene_bitr=bitr(colonic_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hyper_stable_gene_bitr=bitr(colonic_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hypo_stable_gene_bitr=bitr(colonic_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(colonic_hyper_hyper_gene_bitr,colonic_hyper_hypo_gene_bitr,colonic_hypo_hyper_gene_bitr,colonic_hypo_hypo_gene_bitr,colonic_hyper_stable_gene_bitr,colonic_hypo_stable_gene_bitr,file = 'colonic.Rdata')
# 
# colonic_hyper_hyper_go=enrichGO(colonic_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# 
# colonic_hyper_hyper_kegg=enrichKEGG(colonic_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
# 
# 
# colonic_hyper_hypo_go=enrichGO(colonic_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# colonic_hyper_hypo_kegg=enrichKEGG(colonic_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
# 
# 
# colonic_hypo_hyper_go=enrichGO(colonic_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# colonic_hypo_hyper_kegg=enrichKEGG(colonic_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
# 
# 
# colonic_hypo_hypo_go=enrichGO(colonic_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# colonic_hypo_hypo_kegg=enrichKEGG(colonic_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
# 
# 
# colonic_hypo_stable_go=enrichGO(colonic_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# colonic_hypo_stable_kegg=enrichKEGG(colonic_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
# 
# 
# 
# colonic_hyper_stable_go=enrichGO(colonic_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
#                qvalueCutoff = 0.2,keyType = 'ENTREZID')
# colonic_hyper_stable_kegg=enrichKEGG(colonic_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
#                    qvalueCutoff = 0.2, use_internal_data = FALSE)
