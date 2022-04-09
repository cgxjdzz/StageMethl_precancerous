#去掉41384再来一次


inter=intersect(rownames(GSE68060_beta),rownames(GSE139404_beta))
GSE46306_beta_tm=GSE68060_beta[which(rownames(GSE68060_beta) %in% inter),]
TCGA_beta_tm=GSE139404_beta[which(rownames(GSE139404_beta) %in% inter),]
sum_beta=cbind(TCGA_beta_tm,GSE46306_beta_tm)
colnames(GSE68060_pd)=c('name','label','branch')
colnames(GSE139404_pd)=c('name','label','branch')
sum_pd=rbind(GSE139404_pd,GSE68060_pd)

save(sum_beta,sum_pd,file = '0316_sum_data.Rdata')

table(sum_pd$label)
sum_pd=sum_pd[-which(sum_pd$label %in% c('Metastatic','cervical carcinoma in situ')),]
sum_pd$label[which(sum_pd$label %in% c('Invasive cervical carcinoma','Primary Tumor'))] = 'cancer'
sum_pd$label[which(sum_pd$label %in% c('Solid Tissue Normal'))] = 'normal'
sum_beta=sum_beta[,c(sum_pd$name)]
table(sum_pd$label)


library(limma)
batch = sum_pd$branch
sum_beta=removeBatchEffect(sum_beta, batch)

save(sum_beta,sum_pd,file = '0316_sum_data_afterBatch.Rdata')

library(ChAMP)

sum_filter=champ.filter(beta=as.matrix(sum_beta),pd=as.matrix(sum_pd$branch))
save(sum_filter,file = '0316_sum_filter.Rdata')
gc()

sum_impute=champ.impute(beta=sum_filter$beta,pd=sum_filter$pd)
save(sum_impute,file = '0316_sum_impute.Rdata')

gl=sum_impute$pd[,2]
champ.QC(beta=sum_impute$beta,pheno =gl )
sum_norm=champ.norm(beta=sum_impute$beta,cores=5,method='PBC')

pd=as.data.frame(sum_impute$pd)
save(sum_norm,pd,file = '0316_sum_norm.Rdata')

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
ggsave('0316_sumnorm_sort_all_PCA.png')

aa=apply(as.array(sum_norm), MARGIN=2, FUN=sum)
data=list()
data$average_beta=aa/dim(sum_norm)[1]
data$sample=factor(pd$label,levels=c(c("normal",'cin3','cancer')))
data=as.data.frame(data)
ggplot(data=data,aes(x=sample,y=average_beta,colour=sample)) +geom_boxplot()
ggsave('0316_sumnorm_sort_average.png')
