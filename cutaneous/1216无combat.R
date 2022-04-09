require(GEOquery)
require(Biobase)
eset <- getGEO("GSE44661",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
save(beta,pD.all,file = 'GSE44661_orin.Rdata')
pd_GSE44661=pD.all[, c("geo_accession",  "tissue:ch1")][which(pD.all$`tissue:ch1` %in% c('melanocytes','primary melanoma')),]
beta_GSE44661=beta.m[,which(colnames(beta.m) %in% pd_GSE44661$geo_accession)]

save(pd_GSE44661,beta_GSE44661,file='1216_GSE44661.Rdata')


eset <- getGEO("GSE31835",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
save(beta,pD.all,file = 'GSE31835_orin.Rdata')
pd_GSE31835=pD.all[, c("geo_accession", "source_name_ch1")][c(1:10),]
beta_GSE31835=beta.m[,which(colnames(beta.m) %in% pd_GSE31835$geo_accession)]
save(pd_GSE31835,beta_GSE31835,file='1216_GSE31835.Rdata')

##################
beta_tm=beta_GSE44661[which(rownames(beta_GSE44661) %in% rownames(beta_GSE31835)),]
TCGA_beta_tm=beta_GSE31835[which(rownames(beta_GSE31835) %in% rownames(beta_GSE44661)),]
TCGA_beta_tm=TCGA_beta_tm[which(TCGA_beta_tm[,1] !=''),]
beta_tm=beta_tm[which(rownames(beta_tm) %in% rownames(TCGA_beta_tm)),]
beta_tm=beta_tm[c(rownames(TCGA_beta_tm)),]
identical(rownames(TCGA_beta_tm),rownames(beta_tm))
sum_beta=cbind(beta_tm,TCGA_beta_tm)

pD_tm=pd_GSE44661
colnames(pD_tm)=c('name','label')
colnames(pd_GSE31835)=c('name','label')
sum_pd=rbind(pD_tm,pd_GSE31835)
save(sum_beta,sum_pd,file = '1216_sum_data_beforebatch.Rdata')


library(ChAMP)
sum_filter=champ.filter(beta=as.matrix(sum_beta),pd=as.matrix(sum_pd))
save(sum_filter,file = '1216_sum_filter.Rdata')
gc()

sum_impute=champ.impute(beta=sum_filter$beta,pd=sum_filter$pd)
save(sum_impute,file = '1216_sum_impute.Rdata')

pd=as.data.frame(sum_impute$pd)
save(sum_norm,pd,file = '1216_sum_norm.Rdata')

gl=sum_impute$pd[,2]
champ.QC(beta=sum_impute$beta,pheno =gl )
sum_norm=champ.norm(beta=sum_impute$beta,cores=5,method='PBC')


#### 筛选一下位点让所有的癌症都是同样的位点
sum_norm=sum_norm[which(rownames(sum_norm) %in% cgname_ZgSgQlx),]

library(dplyr)
library(ggplot2)
aa=apply(as.array(sum_norm), MARGIN=2, FUN=sum)
data=list()
data$average_beta=aa/dim(sum_norm)[1]
data$sample=factor(pd$label,levels=c(c("skin punch biopsy","melanocytes","primary melanoma")))
data=as.data.frame(data)
ggplot(data=data,aes(x=sample,y=average_beta,colour=sample)) +geom_boxplot()
ggsave('1216_sumnorm_sort_average.png')