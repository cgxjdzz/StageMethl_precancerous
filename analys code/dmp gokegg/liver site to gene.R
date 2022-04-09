
load('D:/甲基化/癌症发展/肝癌/sum_DMP.Rdata')

load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/liver_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/liver_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/liver_hyper_and_hypo.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)

DMP=sum_DMP$Cirrhosis_to_Normal
DMP=DMP[which(DMP$feature %in% c('Body','1stExon')),]
liver_hyper_hyper_gene=unique(DMP[liver_hyper_hyper,]$gene)
liver_hyper_hypo_gene=unique(DMP[liver_hyper_hypo,]$gene)
liver_hypo_hyper_gene=unique(DMP[liver_hypo_hyper,]$gene)
liver_hypo_hypo_gene=unique(DMP[liver_hypo_hypo,]$gene)
liver_hyper_stable_gene=unique(DMP[liver_hyper_stable,]$gene)
liver_hypo_stable_gene=unique(DMP[liver_hypo_stable,]$gene)
liver_hypo_back_to_normal_gene=unique(DMP[liver_hypo_back_to_normal,]$gene)
liver_hyper_back_to_normal_gene=unique(DMP[liver_hyper_back_to_normal,]$gene)
liver_hypo_back_to_normal_gene_bitr=bitr(liver_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hyper_back_to_normal_gene_bitr=bitr(liver_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hyper_hyper_gene_bitr=bitr(liver_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hyper_hypo_gene_bitr=bitr(liver_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hypo_hyper_gene_bitr=bitr(liver_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hypo_hypo_gene_bitr=bitr(liver_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hyper_stable_gene_bitr=bitr(liver_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
liver_hypo_stable_gene_bitr=bitr(liver_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(liver_hyper_hyper_gene_bitr,liver_hyper_hypo_gene_bitr,liver_hypo_hyper_gene_bitr,liver_hypo_hypo_gene_bitr,liver_hyper_stable_gene_bitr,liver_hypo_stable_gene_bitr,file = 'liver.Rdata')


load('liver.Rdata')
liver_hyper_hyper_go=enrichGO(liver_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

liver_hyper_hyper_kegg=enrichKEGG(liver_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


liver_hyper_hypo_go=enrichGO(liver_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
liver_hyper_hypo_kegg=enrichKEGG(liver_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


liver_hypo_hyper_go=enrichGO(liver_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
liver_hypo_hyper_kegg=enrichKEGG(liver_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


liver_hypo_hypo_go=enrichGO(liver_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
liver_hypo_hypo_kegg=enrichKEGG(liver_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


liver_hypo_stable_go=enrichGO(liver_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
liver_hypo_stable_kegg=enrichKEGG(liver_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



liver_hyper_stable_go=enrichGO(liver_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
liver_hyper_stable_kegg=enrichKEGG(liver_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(liver_hyper_hyper_go,liver_hyper_hyper_kegg,liver_hyper_hypo_go,liver_hyper_hypo_kegg,liver_hypo_hyper_go,liver_hypo_hyper_kegg,liver_hypo_hypo_go,liver_hypo_hypo_kegg,liver_hypo_stable_go,liver_hypo_stable_kegg,liver_hyper_stable_go,liver_hyper_stable_kegg,file = 'liver_gokegg.Rdata')