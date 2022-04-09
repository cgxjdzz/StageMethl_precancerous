require(GEOquery)
require(Biobase)


load('eset.Rdata')
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
pd_GSE157973=pD.all[, c("geo_accession",  "diagnosis:ch1")][which(pD.all$`diagnosis:ch1`=='Cirrhosis'),]
beta_GSE157973=beta.m[,which(colnames(beta.m) %in% pd_GSE157973$geo_accession)]


require(readr)
temp_pd=read_tsv('TCGA-LIHC.GDC_phenotype.tsv')
temp=read_tsv('TCGA-LIHC.methylation450.tsv')
temp=as.data.frame(temp)
row.names(temp) <- temp[,1]
temp=temp[,-1]

pd_TCGA=temp_pd[,c("submitter_id.samples","sample_type.samples")]
beta_TCGA=temp[,which(colnames(temp) %in% pd_TCGA$submitter_id.samples)]
pd_TCGA=pd_TCGA[which(pd_TCGA$submitter_id.samples %in% colnames(temp)),]
beta_TCGA=beta_TCGA[,c(pd_TCGA$submitter_id.samples)]


identical(colnames(beta_TCGA),(pd_TCGA$submitter_id.samples))
save(beta_TCGA,pd_TCGA,file = '0227TCGA.Rdata')

beta_tm=beta_GSE157973[which(rownames(beta_GSE157973) %in% rownames(beta_TCGA)),]
TCGA_beta_tm=beta_TCGA[which(rownames(beta_TCGA) %in% rownames(beta_GSE157973)),]
TCGA_beta_tm=TCGA_beta_tm[which(TCGA_beta_tm[,1] !=''),]
beta_tm=beta_tm[which(rownames(beta_tm) %in% rownames(TCGA_beta_tm)),]
beta_tm=beta_tm[c(rownames(TCGA_beta_tm)),]
identical(rownames(TCGA_beta_tm),rownames(beta_tm))
sum_beta=cbind(beta_tm,TCGA_beta_tm)
pD_tm=pd_GSE157973
colnames(pD_tm)=c('name','label')
colnames(pd_TCGA)=c('name','label')
sum_pd=rbind(pD_tm,pd_TCGA)
identical(colnames(sum_beta),sum_pd$name)
save(sum_beta,sum_pd,file = '0227_sum_data.Rdata')






load('0227_sum_data.Rdata')

library(ChAMP)

sum_filter=champ.filter(beta=as.matrix(sum_beta),pd=as.matrix(sum_pd))
save(sum_filter,file = '0227_sum_filter.Rdata')

load('0227_sum_filter.Rdata')
sum_impute=champ.impute(beta=sum_filter$beta,pd=sum_filter$pd,method='Delete')
save(sum_impute,file = '0227_sum_impute.Rdata')

gl=sum_impute$pd[,2]
champ.QC(beta=sum_impute$beta,pheno =gl )
sum_norm=champ.norm(beta=sum_impute$beta,cores=5,method='PBC')

pd=as.data.frame(sum_impute$pd)
save(sum_norm,pd,file = '0227_sum_norm.Rdata')

pd$label[which(pd$label=='Primary Tumor' | pd$label=='Recurrent Tumor')]='Tumor'
pd$label[which(pd$label=='Solid Tissue Normal')]='Normal'

sum_DMP=champ.DMP(beta=sum_norm,pheno=pd$label)
save(sum_DMP,file='sum_DMP.Rdata')


exprSet=sum_norm
pheatmap::pheatmap(cor(exprSet)) 
# 组内的样本的相似性应该是要高于组间的！
group_list=pd$label
colD=data.frame(group_list=group_list)
rownames(colD)=colnames(exprSet)
pheatmap::pheatmap(cor(exprSet),
                   annotation_row =colD,
                   annotation_col = colD,
                   show_rownames = F,
                   show_colnames = F,
                   filename = '0227_sumnorm_sort_cor_all.png')


library(dplyr)
library(ggplot2)

pdlabel=pd$label
png('0227_sumnorm_average_density.png')
plot(density(as.numeric(sum_norm[,c(which(pdlabel=="Cirrhosis"))])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,10),
     xlab="β",ylab="Density",main="Density Plot")+
  lines(density(as.numeric(sum_norm[,c(which(pdlabel=="Tumor"))])),col="blue",lwd=2.5)+
  lines(density(as.numeric(sum_norm[,c(which(pdlabel=='Normal'))])),col="orange",lwd=2.5)
legend("top",                legend=c(unique(pdlabel)),        #图例内容
       col=c("red","blue","orange"),                 #图例颜色
       lty=1,lwd=1)                                          #图例大小
dev.off()


aa=apply(as.array(sum_norm), MARGIN=2, FUN=sum)
data=list()
data$average_beta=aa/dim(sum_norm)[1]
data$sample=factor(pd$label,levels=c(c("Normal","Cirrhosis","Tumor")))
data=as.data.frame(data)
ggplot(data=data,aes(x=sample,y=average_beta,colour=sample)) +geom_boxplot()
ggsave('0227_sumnorm_sort_average.png')


library("FactoMineR")#画主成分分析图需要加载这两个???
library("factoextra")  
dat=t(sum_norm)
group_list=pd$label

dat.pca <- PCA(dat , graph = FALSE)
save(dat.pca,file='dat_pca.Rdata')
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('0227_sumnorm_sort_all_PCA.png')

