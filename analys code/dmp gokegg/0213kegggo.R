
library(clusterProfiler)
load('cervical.Rdata')
cervical_hyper_hyper_go=enrichGO(cervical_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

cervical_hyper_hyper_kegg=enrichKEGG(cervical_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


cervical_hyper_hypo_go=enrichGO(cervical_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hyper_hypo_kegg=enrichKEGG(cervical_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


cervical_hypo_hyper_go=enrichGO(cervical_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hypo_hyper_kegg=enrichKEGG(cervical_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


cervical_hypo_hypo_go=enrichGO(cervical_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hypo_hypo_kegg=enrichKEGG(cervical_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


cervical_hypo_stable_go=enrichGO(cervical_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hypo_stable_kegg=enrichKEGG(cervical_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



cervical_hyper_stable_go=enrichGO(cervical_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
cervical_hyper_stable_kegg=enrichKEGG(cervical_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(cervical_hyper_hyper_go,cervical_hyper_hyper_kegg,cervical_hyper_hypo_go,cervical_hyper_hypo_kegg,cervical_hypo_hyper_go,cervical_hypo_hyper_kegg,cervical_hypo_hypo_go,cervical_hypo_hypo_kegg,cervical_hypo_stable_go,cervical_hypo_stable_kegg,cervical_hyper_stable_go,cervical_hyper_stable_kegg,file = 'cervical_gokegg.Rdata')


load('melanoma.Rdata')
melanoma_hyper_hyper_go=enrichGO(melanoma_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

melanoma_hyper_hyper_kegg=enrichKEGG(melanoma_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


melanoma_hyper_hypo_go=enrichGO(melanoma_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hyper_hypo_kegg=enrichKEGG(melanoma_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


melanoma_hypo_hyper_go=enrichGO(melanoma_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hypo_hyper_kegg=enrichKEGG(melanoma_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


melanoma_hypo_hypo_go=enrichGO(melanoma_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hypo_hypo_kegg=enrichKEGG(melanoma_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


melanoma_hypo_stable_go=enrichGO(melanoma_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hypo_stable_kegg=enrichKEGG(melanoma_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



melanoma_hyper_stable_go=enrichGO(melanoma_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
melanoma_hyper_stable_kegg=enrichKEGG(melanoma_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(melanoma_hyper_hyper_go,melanoma_hyper_hyper_kegg,melanoma_hyper_hypo_go,melanoma_hyper_hypo_kegg,melanoma_hypo_hyper_go,melanoma_hypo_hyper_kegg,melanoma_hypo_hypo_go,melanoma_hypo_hypo_kegg,melanoma_hypo_stable_go,melanoma_hypo_stable_kegg,melanoma_hyper_stable_go,melanoma_hyper_stable_kegg,file = 'melanoma_gokegg.Rdata')

load('gastric.Rdata')
gastric_hyper_hyper_go=enrichGO(gastric_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

gastric_hyper_hyper_kegg=enrichKEGG(gastric_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


gastric_hyper_hypo_go=enrichGO(gastric_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hyper_hypo_kegg=enrichKEGG(gastric_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


gastric_hypo_hyper_go=enrichGO(gastric_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hypo_hyper_kegg=enrichKEGG(gastric_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


gastric_hypo_hypo_go=enrichGO(gastric_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hypo_hypo_kegg=enrichKEGG(gastric_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


gastric_hypo_stable_go=enrichGO(gastric_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hypo_stable_kegg=enrichKEGG(gastric_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



gastric_hyper_stable_go=enrichGO(gastric_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
gastric_hyper_stable_kegg=enrichKEGG(gastric_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(gastric_hyper_hyper_go,gastric_hyper_hyper_kegg,gastric_hyper_hypo_go,gastric_hyper_hypo_kegg,gastric_hypo_hyper_go,gastric_hypo_hyper_kegg,gastric_hypo_hypo_go,gastric_hypo_hypo_kegg,gastric_hypo_stable_go,gastric_hypo_stable_kegg,gastric_hyper_stable_go,gastric_hyper_stable_kegg,file = 'gastric_gokegg.Rdata')


load('prostatic.Rdata')
prostatic_hyper_hyper_go=enrichGO(prostatic_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

prostatic_hyper_hyper_kegg=enrichKEGG(prostatic_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


prostatic_hyper_hypo_go=enrichGO(prostatic_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hyper_hypo_kegg=enrichKEGG(prostatic_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


prostatic_hypo_hyper_go=enrichGO(prostatic_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hypo_hyper_kegg=enrichKEGG(prostatic_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


prostatic_hypo_hypo_go=enrichGO(prostatic_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hypo_hypo_kegg=enrichKEGG(prostatic_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


prostatic_hypo_stable_go=enrichGO(prostatic_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hypo_stable_kegg=enrichKEGG(prostatic_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



prostatic_hyper_stable_go=enrichGO(prostatic_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
prostatic_hyper_stable_kegg=enrichKEGG(prostatic_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(prostatic_hyper_hyper_go,prostatic_hyper_hyper_kegg,prostatic_hyper_hypo_go,prostatic_hyper_hypo_kegg,prostatic_hypo_hyper_go,prostatic_hypo_hyper_kegg,prostatic_hypo_hypo_go,prostatic_hypo_hypo_kegg,prostatic_hypo_stable_go,prostatic_hypo_stable_kegg,prostatic_hyper_stable_go,prostatic_hyper_stable_kegg,file = 'prostatic_gokegg.Rdata')

load('colonic.Rdata')
colonic_hyper_hyper_go=enrichGO(colonic_hyper_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')

colonic_hyper_hyper_kegg=enrichKEGG(colonic_hyper_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)


colonic_hyper_hypo_go=enrichGO(colonic_hyper_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hyper_hypo_kegg=enrichKEGG(colonic_hyper_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


colonic_hypo_hyper_go=enrichGO(colonic_hypo_hyper_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hypo_hyper_kegg=enrichKEGG(colonic_hypo_hyper_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)


colonic_hypo_hypo_go=enrichGO(colonic_hypo_hypo_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hypo_hypo_kegg=enrichKEGG(colonic_hypo_hypo_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2, use_internal_data = FALSE)


colonic_hypo_stable_go=enrichGO(colonic_hypo_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hypo_stable_kegg=enrichKEGG(colonic_hypo_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)



colonic_hyper_stable_go=enrichGO(colonic_hyper_stable_gene_bitr, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
colonic_hyper_stable_kegg=enrichKEGG(colonic_hyper_stable_gene_bitr, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                     qvalueCutoff = 0.2, use_internal_data = FALSE)

save(colonic_hyper_hyper_go,colonic_hyper_hyper_kegg,colonic_hyper_hypo_go,colonic_hyper_hypo_kegg,colonic_hypo_hyper_go,colonic_hypo_hyper_kegg,colonic_hypo_hypo_go,colonic_hypo_hypo_kegg,colonic_hypo_stable_go,colonic_hypo_stable_kegg,colonic_hyper_stable_go,colonic_hyper_stable_kegg,file = 'colonic_gokegg.Rdata')