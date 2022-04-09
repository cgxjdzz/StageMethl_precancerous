
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum" ,"Primary Tumor"),c("Normal",'Precancer',"Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/gastric_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/gastric_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0212新新大纲/dmp的热图/gastric_hyper_and_hypo.Rdata')
                             

DMP=sum_DMP$Normal_to_Precancer
DMP=DMP[which(DMP$feature %in% c('TSS200','TSS1500')),]
gastric_hyper_hyper_gene=unique(DMP[gastric_hyper_hyper,]$gene)
gastric_hyper_hypo_gene=unique(DMP[gastric_hyper_hypo,]$gene)
gastric_hypo_hyper_gene=unique(DMP[gastric_hypo_hyper,]$gene)
gastric_hypo_hypo_gene=unique(DMP[gastric_hypo_hypo,]$gene)
gastric_hyper_stable_gene=unique(DMP[gastric_hyper_stable,]$gene)
gastric_hypo_stable_gene=unique(DMP[gastric_hypo_stable,]$gene)
gastric_hypo_back_to_normal_gene=unique(DMP[gastric_hypo_back_to_normal,]$gene)
gastric_hyper_back_to_normal_gene=unique(DMP[gastric_hyper_back_to_normal,]$gene)

gastric_hypo_back_to_normal_gene_bitr=bitr(gastric_hypo_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hyper_back_to_normal_gene_bitr=bitr(gastric_hyper_back_to_normal_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hyper_hyper_gene_bitr=bitr(gastric_hyper_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hyper_hypo_gene_bitr=bitr(gastric_hyper_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hypo_hyper_gene_bitr=bitr(gastric_hypo_hyper_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hypo_hypo_gene_bitr=bitr(gastric_hypo_hypo_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hyper_stable_gene_bitr=bitr(gastric_hyper_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
gastric_hypo_stable_gene_bitr=bitr(gastric_hypo_stable_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)$ENTREZID
save(gastric_hyper_hyper_gene_bitr,gastric_hyper_hypo_gene_bitr,gastric_hypo_hyper_gene_bitr,gastric_hypo_hypo_gene_bitr,gastric_hyper_stable_gene_bitr,gastric_hypo_stable_gene_bitr,file = 'gastric.Rdata')
