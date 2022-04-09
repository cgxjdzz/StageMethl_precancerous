




library(stringr)
library(org.Hs.eg.db)
library(minfi)
library(clusterProfiler)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)



###将pD和beta按组排序
ee=myCombat[,order(c(pD$label))]
pD=pD[,order(c(pD$label))]
identical(colnames(ee),pD$geo_accession)


###QC流程
library("FactoMineR")#鐢讳富鎴愬垎鍒嗘瀽鍥鹃渶瑕佸姞杞借繖涓や釜鍖?
library("factoextra")  
dat=t(ee)
group_list=pD$label

dat.pca <- PCA(dat , graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('COmbat_sort_all_PCA.png')


####heatmap
dat=ee
dat[1:4,1:4] 
cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=group_list)
rownames(ac)=colnames(n)  
pheatmap(n,show_colnames =F,show_rownames = F,cluster_rows = F,
         annotation_col=ac,filename = 'Combat_sort_heatmap_top1000_sd.png')
dev.off()



#####相关系数热图

exprSet=ee
pheatmap::pheatmap(cor(exprSet)) 
# 组内的样本的相似性应该是要高于组间的！
colD=data.frame(group_list=group_list)
rownames(colD)=colnames(exprSet)
pheatmap::pheatmap(cor(exprSet),
                   annotation_row =colD,
                   annotation_col = colD,
                   show_rownames = F,
                   show_colnames = F,
                   filename = 'Combat_sort_cor_all.png')

####champ的QC图

champ.QC(beta = myCombat,pheno = pD$label)
QC.GUI()

####各组平均值的密度图
library(dplyr)
library(ggplot2)
data=myCombat
group_list=impute$pd
group_name=unique(group_list$label)
sum=list()
k=1
for (i in c(group_name)){
  temp_pd=filter(group_list,label==i)
  temp_data=data[,which(colnames(data) %in% temp_pd$geo_accession)]
  nn=apply(as.array(temp_data), MARGIN=1, FUN = sum)
  sum$Tumor=nn
  k=k+1
  #temp_data=data %>% select(!starts_with("GSM"))????
}
sum=as.data.frame(sum)
ggplot(sum,aes(x=sum,fill=new,alpha = 1/10))+geom_density()
png('average_density.png')
plot(density(as.numeric(myCombat[,c(which(pd=='LGA'))])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,4.5),
     xlab="β",ylab="Density",main="Density Plot")+
  lines(density(as.numeric(myCombat[,c(which(pd=='Normal'))])),col="blue",lwd=2.5)+
  lines(density(as.numeric(myCombat[,c(which(pd=="HGA"))])),col="green",lwd=2.5)+
  lines(density(as.numeric(myCombat[,c(which(pd=="Tumor"))])),col="orange",lwd=2.5)
legend("top",                legend=c(unique(pd)),        #图例内容
      col=c("red","blue","green","orange"),                 #图例颜色
       lty=1,lwd=1)                                          #图例大小

ggsave(px,'average_density.png')
####平均甲基化图

aa=apply(as.array(sum_norm), MARGIN=2, FUN=sum)
data=list()
data$average_beta=aa/dim(sum_norm)[1]
data$sample=pd$label
data=as.data.frame(data)
ggplot(data=data,aes(x=sample,y=average_beta,colour=sample)) +geom_boxplot()
ggsave('average.png')

 

##DMP位点截取
myDMP=champ.DMP(beta=myCombat,pheno = impute$pd$label)
DMP.GUI(DMP = myDMP[[1]],beta=myCombat,pheno = impute$pd$label,cutgroupnumber = )








###go和kegg

cervical_up_up_gene_bitr <- bitr(unique(cervical_up_up_gene$Var1), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)

cervical_up_up_go <- enrichGO(cervical_up_up_gene_bitr$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2,keyType = 'ENTREZID')
png(file='cervical_up_up_go_barplot.png', width = 1000,height = 600)
barplot(cervical_up_up_go,showCategory=20,drop=T)

dev.off()

dotplot(cervical_up_up_go,showCategory=50)

cervical_up_up_kegg=enrichKEGG(cervical_up_up_gene_bitr$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
           pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
           qvalueCutoff = 0.2, use_internal_data = FALSE)
save(cervical_up_up_go,cervical_up_up_kegg,file='cervical_up_up_go_kegg.Rdata')


cervical_down_down_gene_bitr <- bitr(unique(cervical_down_down_gene$Var1), fromType = "SYMBOL",
                                     toType = c( "ENTREZID"),
                                     OrgDb = org.Hs.eg.db)

cervical_down_down_go <- enrichGO(cervical_down_down_gene_bitr$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2,keyType = 'ENTREZID')
png(file='cervical_down_down_go_barplot.png', width = 1000,height = 600)
barplot(cervical_down_down_go,showCategory=20,drop=T)
dev.off()

dotplot(go,showCategory=50)

cervical_down_down_kegg=enrichKEGG(cervical_down_down_gene_bitr$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                                   qvalueCutoff = 0.2, use_internal_data = FALSE)
save(cervical_down_down_go,cervical_down_down_kegg,file='cervical_down_down_go_kegg.Rdata')