
load('D:/甲基化/癌症发展/肝癌/sum_DMP.Rdata')
liver_DMP=sum_DMP$Cirrhosis_to_Normal
load('D:/甲基化/癌症发展/0331改数据/黑色素瘤/0331_sum_DMP.Rdata')
melanoma_DMP=sum_DMP$Precancer_to_Normal
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
prostatic_DMP=sum_DMP$CAF_to_Normal
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum" ,"Primary Tumor"),c("Normal",'Precancer',"Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
gastric_DMP=sum_DMP$Normal_to_Precancer
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
cervical_DMP=sum_DMP$normal_to_cin3
load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')
myNorm=sum_norm
pd=pd$`sum_impute$pd`
pd=mapvalues(pd, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = myNorm,pheno = pd)
colonic_DMP=sum_DMP$Normal_to_Precancer
###################################
site_list=c("union_site_hyper_back_to_normal","union_site_hyper_hyper","union_site_hyper_hypo","union_site_hyper_stable","union_site_hypo_back_to_normal","union_site_hypo_hyper","union_site_hypo_hypo","union_site_hypo_stable" )
library(ChAMP)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
data(probe.features.epic)
for (i in site_list){
  data=read.table(paste0(i,'/',i,'.csv'))
  cgpoint=rownames(data)
  gene=unique(na.omit(probe.features[cgpoint,]$gene))
  go=enrichGO(gene, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,keyType = "SYMBOL",
                                   qvalueCutoff = 0.2)
  geneid=bitr(gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb =org.Hs.eg.db )
  kegg=enrichKEGG(geneid$ENTREZID)
  
  go_result=go@result
  go_result=go_result[order(go_result$p.adjust),]
  write.csv(go@result,file= paste0(i,'/go_result.csv'))
  go_result$GeneRatio=factor(go_result$GeneRatio,levels=c(unique(go_result$GeneRatio)[order(as.numeric(as.data.frame(strsplit(unique(go_result$GeneRatio),'/'))[1,]))]))
  go_result=go_result[c(1:30),]
  
  p <- ggplot(go_result, mapping = aes(x=(GeneRatio), y=Description, size=Count, color=p.adjust))
  p <- p + geom_point()
  p <- p + xlab("Cancer category") + ylab("Terms")
  p <- p + scale_colour_gradient(high="#990000", low="#FF9999")
  p <- p + scale_size_continuous(range = c(5,7))#点大小缩放
  p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                              axis.title = element_text(size = 20))
  p
  ggsave(paste0(i,'/go_result.pdf'),width =9, height =6)
  
  kegg_result=kegg@result
  kegg_result=kegg_result[order(kegg_result$p.adjust),]
  write.csv(kegg@result,file= paste0(i,'/kegg_result.csv'))
  kegg_result$GeneRatio=factor(kegg_result$GeneRatio,levels=c(unique(kegg_result$GeneRatio)[order(as.numeric(as.data.frame(strsplit(unique(kegg_result$GeneRatio),'/'))[1,]))]))
  kegg_result=kegg_result[c(1:30),]
  
  p <- ggplot(kegg_result, mapping = aes(x=(GeneRatio), y=Description, size=Count, color=p.adjust))
  p <- p + geom_point()
  p <- p + xlab("Cancer catekeggry") + ylab("Terms")
  p <- p + scale_colour_gradient(high="#990000", low="#FF9999")
  p <- p + scale_size_continuous(range = c(5,7))#点大小缩放
  p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                              axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                              axis.title = element_text(size = 20))
  p
  ggsave(paste0(i,'/kegg_result.pdf'),width =9, height =6)
  }