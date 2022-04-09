
require(GEOquery)
require(Biobase)
eset <- getGEO("GSE99553",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
save(beta,pD.all,file = '0107——GSE99553_origin_data.Rdata')
pd_GSE99553=pD.all[, c("title", "geo_accession",  "cancer status:ch1")][which(pD.all$`cancer status:ch1`=='Control'),]
beta_GSE99553=beta.m[,which(colnames(beta.m) %in% pd_GSE99553$geo_accession)]
identical(colnames(beta_GSE99553),pd_GSE99553$geo_accession)
save(beta_GSE99553,pd_GSE99553,file = '0107GSE99553.Rdata')

beta_tm=beta_GSE99553[which(rownames(beta_GSE99553) %in% rownames(sum_beta)),]
sum_beta_tm=sum_beta[which(rownames(sum_beta) %in% rownames(beta_GSE99553)),]
sum_beta_tm=sum_beta_tm[which(sum_beta_tm[,1] !=''),]
beta_tm=beta_tm[which(rownames(beta_tm) %in% rownames(sum_beta_tm)),]
beta_tm=beta_tm[c(rownames(sum_beta_tm)),]
identical(rownames(sum_beta_tm),rownames(beta_tm))
sum_beta=cbind(beta_tm,sum_beta_tm)

pd_GSE99553$`cancer status:ch1`='Solid Tissue Normal'
pD_tm=pd_GSE99553[,c("geo_accession","cancer status:ch1" )]
colnames(pD_tm)=c('name','label')
sum_pd=rbind(pD_tm,sum_pd)
