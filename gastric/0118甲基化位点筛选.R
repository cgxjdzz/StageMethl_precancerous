Solid_Tissue_Normal_norm=sum_norm[,which(pd$label=='Solid Tissue Normal')]
intestinal_metaplasia_biopsy_from_gastric_antrum_norm=sum_norm[,which(pd$label=='intestinal metaplasia biopsy from gastric antrum')]
cancer_norm=sum_norm[,which(pd$label=='Primary Tumor')]
source('zffdatamining_gene.R')
#筛选从Solid_Tissue_Normal到intestinal_metaplasia_biopsy_from_gastric_antrum的位点
Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down=c()
Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up=c()

for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  if(is.na(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal[cgname,]$logFC)){
    next
  }
  
  temp=ge.groupdifference(data=as.numeric(Solid_Tissue_Normal_norm[i,]),cut = 0.9)
  
  if((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal[cgname,]$logFC)<temp[1]){
    Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down=c(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down,cgname)
  }
  if((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal[cgname,]$logFC)>temp[2]){
    Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up=c(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up,cgname)
  }
}
save(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down,file='0121Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down.Rdata')
save(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up,file='0121Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up.Rdata')

#筛选从intestinal_metaplasia_biopsy_from_gastric_antrum到Primary_Tumor的位点
intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down=c()
intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up=c()

for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  if(is.na(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor[cgname,]$logFC)){
    next
  }
  temp=ge.groupdifference(data=as.numeric(intestinal_metaplasia_biopsy_from_gastric_antrum_norm[i,]),cut = 0.9)
  if((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor[cgname,]$logFC)<temp[1]){
    intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up=c(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up,cgname)
  }
  if((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor[cgname,]$logFC)>temp[2]){
    intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down=c(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down,cgname)
  }
}
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down,file='0121intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down.Rdata')
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up,file='0121intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up.Rdata')


gastric_up_up=intersect((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up),(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up))
save(gastric_up_up,file = '0121DMPgastric_up_up.Rdata')

gastric_down_down=intersect((intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down),(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down))
save(gastric_down_down,file ='0121DMPgastric_down_down.Rdata' )





library(stringr)
library(org.Hs.eg.db)
library(minfi)
library(clusterProfiler)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)

load('DMP gastric_down_down_gene.Rdata')
load('DMP gastric_up_up_gene.Rdata')
###
gastric_up_up_gene_bitr <- bitr(unique(gastric_up_up_gene$Var1), fromType = "SYMBOL",
                                  toType = c( "ENTREZID"),
                                  OrgDb = org.Hs.eg.db)

gastric_up_up_go <- enrichGO(gastric_up_up_gene_bitr$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.01,
                               qvalueCutoff = 0.2,keyType = 'ENTREZID')
png(file='gastric_up_up_go_barplot.png', width = 1000,height = 600)
barplot(gastric_up_up_go,showCategory=20,drop=T)

dev.off()

dotplot(gastric_up_up_go,showCategory=50)

gastric_up_up_kegg=enrichKEGG(gastric_up_up_gene_bitr$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                qvalueCutoff = 0.2, use_internal_data = FALSE)
save(gastric_up_up_go,gastric_up_up_kegg,file='gastric_up_up_go_kegg.Rdata')
#####################

gastric_down_down_gene_bitr <- bitr(unique(gastric_down_down_gene$Var1), fromType = "SYMBOL",
                                      toType = c( "ENTREZID"),
                                      OrgDb = org.Hs.eg.db)

gastric_down_down_go <- enrichGO(gastric_down_down_gene_bitr$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.2,keyType = 'ENTREZID')
png(file='gastric_down_down_go_barplot.png', width = 1000,height = 600)
barplot(gastric_down_down_go,showCategory=20,drop=T)
dev.off()

dotplot(go,showCategory=50)

gastric_down_down_kegg=enrichKEGG(gastric_down_down_gene_bitr$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                    qvalueCutoff = 0.2, use_internal_data = FALSE)
save(gastric_down_down_go,gastric_down_down_kegg,file='gastric_down_down_go_kegg.Rdata')