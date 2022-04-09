load('colonic_hy_and_back_to_normal.Rdata')
load('cervical_hy_and_back_to_normal.Rdata')
load('melanoma_hy_and_back_to_normal.Rdata')
load('prostatic_hy_and_back_to_normal.Rdata')
load('gastric_hy_and_back_to_normal.Rdata')
library(ChAMP)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
data(probe.features.epic)
cervical_hyper_back_to_normal_gene=unique(probe.features[cervical_hyper_back_to_normal,]$gene)
cervical_hypo_back_to_normal_gene=unique(probe.features[cervical_hypo_back_to_normal,]$gene)

colonic_hyper_back_to_normal_gene=unique(probe.features[colonic_hyper_back_to_normal,]$gene)
colonic_hypo_back_to_normal_gene=unique(probe.features[colonic_hypo_back_to_normal,]$gene)

melanoma_hyper_back_to_normal_gene=unique(probe.features[melanoma_hyper_back_to_normal,]$gene)
melanoma_hypo_back_to_normal_gene=unique(probe.features[melanoma_hypo_back_to_normal,]$gene)

prostatic_hyper_back_to_normal_gene=unique(probe.features[prostatic_hyper_back_to_normal,]$gene)
prostatic_hypo_back_to_normal_gene=unique(probe.features[prostatic_hypo_back_to_normal,]$gene)

gastric_hyper_back_to_normal_gene=unique(probe.features[gastric_hyper_back_to_normal,]$gene)
gastric_hypo_back_to_normal_gene=unique(probe.features[gastric_hypo_back_to_normal,]$gene)


cervical_hyper_back_to_normal_gene_bitr=bitr(cervical_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
cervical_hyper_back_to_normal_gene_bitr=bitr(cervical_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID

cervical_hyper_back_to_normal_go=enrichGO(cervical_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hyper_back_to_normal_kegg=enrichKEGG(cervical_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



cervical_hyper_back_to_normal_go=enrichGO(cervical_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hyper_back_to_normal_kegg=enrichKEGG(cervical_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

colonic_hyper_back_to_normal_gene_bitr=bitr(colonic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
colonic_hyper_back_to_normal_gene_bitr=bitr(colonic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID

colonic_hyper_back_to_normal_go=enrichGO(colonic_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hyper_back_to_normal_kegg=enrichKEGG(colonic_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)



colonic_hyper_back_to_normal_go=enrichGO(colonic_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hyper_back_to_normal_kegg=enrichKEGG(colonic_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)

melanoma_hyper_back_to_normal_gene_bitr=bitr(melanoma_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
melanoma_hyper_back_to_normal_gene_bitr=bitr(melanoma_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID

melanoma_hyper_back_to_normal_go=enrichGO(melanoma_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hyper_back_to_normal_kegg=enrichKEGG(melanoma_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)



melanoma_hyper_back_to_normal_go=enrichGO(melanoma_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hyper_back_to_normal_kegg=enrichKEGG(melanoma_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)


prostatic_hyper_back_to_normal_gene_bitr=bitr(prostatic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
prostatic_hyper_back_to_normal_gene_bitr=bitr(prostatic_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID

prostatic_hyper_back_to_normal_go=enrichGO(prostatic_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hyper_back_to_normal_kegg=enrichKEGG(prostatic_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)



prostatic_hyper_back_to_normal_go=enrichGO(prostatic_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hyper_back_to_normal_kegg=enrichKEGG(prostatic_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)

gastric_hyper_back_to_normal_gene_bitr=bitr(gastric_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hyper_back_to_normal_gene_bitr=bitr(gastric_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID

gastric_hyper_back_to_normal_go=enrichGO(gastric_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hyper_back_to_normal_kegg=enrichKEGG(gastric_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)



gastric_hyper_back_to_normal_go=enrichGO(gastric_hyper_back_to_normal_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                      qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hyper_back_to_normal_kegg=enrichKEGG(gastric_hyper_back_to_normal_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                          qvalueCutoff = 0.2, use_internal_data = FALSE)