
require(GEOquery)
require(Biobase)
eset <- getGEO("GSE103186",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
save(beta,pD.all,file = '0107_origin_data.Rdata')
pd_GSE103186=pD.all[, c("title", "geo_accession",  "source_name_ch1")][which(pD.all$source_name_ch1=='CAF'),]
beta_GSE103186=beta.m[,which(colnames(beta.m) %in% pd$geo_accession)]
identical(colnames(beta_GSE103186),pd_GSE103186$geo_accession)
save(beta_GSE103186,pd_GSE103186,file = '0107GSE103186.Rdata')

require(readr)
temp_pd=read_tsv('TCGA.STAD.sampleMap_STAD_clinicalMatrix.tsv')
temp=read_tsv('HumanMethylation450.tsv')
temp=as.data.frame(temp)
rownames(temp) <- temp[,1]
temp=temp[,-1]

#这里还有别的但只选了这个 
pd_GSE103186=pd_GSE103186[,c('geo_accession',"tissue:ch1")][which(pd_GSE103186$`tissue:ch1`=='intestinal metaplasia biopsy from gastric antrum'),]
pd_TCGA=temp_pd[,c("sampleID","sample_type")]
#pd_TCGA=pd_TCGA[-which(pd_TCGA$sample_type.samples=='Metastatic'),]
beta_TCGA=temp[,which(colnames(temp) %in% pd_TCGA$sampleID)]
pd_TCGA=pd_TCGA[which(pd_TCGA$sampleID %in% colnames(temp)),]
beta_TCGA=beta_TCGA[,c(pd_TCGA$sampleID)]


identical(colnames(beta_TCGA),(pd_TCGA$sampleID))
save(beta_TCGA,pd_TCGA,file = '0107TCGA.Rdata')





beta_tm=beta_GSE103186[which(rownames(beta_GSE103186) %in% rownames(beta_TCGA)),]
TCGA_beta_tm=beta_TCGA[which(rownames(beta_TCGA) %in% rownames(beta_GSE103186)),]
TCGA_beta_tm=TCGA_beta_tm[which(TCGA_beta_tm[,1] !=''),]
beta_tm=beta_tm[which(rownames(beta_tm) %in% rownames(TCGA_beta_tm)),]
beta_tm=beta_tm[c(rownames(TCGA_beta_tm)),]
identical(rownames(TCGA_beta_tm),rownames(beta_tm))
sum_beta=cbind(beta_tm,TCGA_beta_tm)

pD_tm=pd_GSE103186[,c("geo_accession","tissue:ch1" )]
colnames(pD_tm)=c('name','label')
colnames(pd_TCGA)=c('name','label')
sum_pd=rbind(pD_tm,pd_TCGA)

save(sum_beta,sum_pd,file = '0107_sum_data.Rdata')
identical(colnames(sum_beta),sum_pd$name)

#### 筛选一下位点让所有的癌症都是同样的位点
#sum_beta=sum_beta[which(rownames(sum_beta) %in% cgname_ZgSgQlx),]


library(ChAMP)
sum_filter$beta=as.list(sum_filter$beta)
sum_filter$pd=as.list(sum_filter$pd)
sum_filter=champ.filter(beta=as.matrix(sum_beta),pd=as.matrix(sum_pd))
save(sum_filter,file = '0107_sum_filter.Rdata')
gc()


Sys.setenv("VROOM_CONNECTION_SIZE"=131072*13)


sum_impute=champ.impute(beta=sum_filter$beta,pd=sum_filter$pd)
save(sum_impute,file = '0107_sum_impute.Rdata')

gl=sum_impute$pd[,2]
champ.QC(beta=sum_impute$beta,pheno =gl )
sum_norm=champ.norm(beta=sum_impute$beta,cores=5,method='PBC')

pd=as.data.frame(sum_impute$pd)
save(sum_norm,pd,file = '0107_sum_norm.Rdata')







###







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
                   filename = '0107_sumnorm_sort_cor_all.png')
library(dplyr)
library(ggplot2)

pd=pd$label
png('11.1_sumnorm_average_density.png')
plot(density(as.numeric(sum_norm[,c(which(pd=="CAF"))])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,10),
     xlab="β",ylab="Density",main="Density Plot")+
  lines(density(as.numeric(sum_norm[,c(which(pd=="Primary Tumor"))])),col="blue",lwd=2.5)+
  lines(density(as.numeric(sum_norm[,c(which(pd=='Solid Tissue Normal'))])),col="orange",lwd=2.5)
legend("top",                legend=c(unique(pd)),        #图例内容
       col=c("red","blue","orange"),                 #图例颜色
       lty=1,lwd=1)                                          #图例大小
dev.off()


#sum_norm=sum_norm[,order(factor(pd$label,levels=c(c("Solid Tissue Normal","CAF","Primary Tumor"))))]
#pd=pd[order(factor(pd$label,levels=c(c("Solid Tissue Normal","CAF","Primary Tumor")))),]
#identical(colnames(sum_norm),pd)

aa=apply(as.array(sum_norm), MARGIN=2, FUN=sum)
data=list()
data$average_beta=aa/dim(sum_norm)[1] 
data$sample=factor(pd$label,levels=c(c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum","Primary Tumor")))
data=as.data.frame(data)
ggplot(data=data,aes(x=sample,y=average_beta,colour=sample)) +geom_boxplot()
ggsave('0107_sumnorm_sort_average.png')

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
ggsave('0107_sumnorm_sort_all_PCA.png')

