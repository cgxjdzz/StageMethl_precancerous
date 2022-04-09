hsa_kegg <- clusterProfiler::download_KEGG("hsa")
hsa_kegg=hsa_kegg$KEGGPATHID2EXTID
pathway='hsa04151'
pathway_geneid=hsa_kegg[which(hsa_kegg$from==pathway),]
library(clusterProfiler)
library(org.Hs.eg.db)
data(probe.features.epic)
library(ChAMP)
library(readr)
pathway_gene=bitr(pathway_geneid$to,fromType = 'ENTREZID',toType ='SYMBOL',OrgDb =org.Hs.eg.db)
pathway_gene=pathway_gene$SYMBOL
probe.features=probe.features[which(probe.features$gene %in% pathway_gene),]
pathway_probe=rownames(probe.features)
save(pathway_probe,file='pathway_probe.Rdata')
load('pathway_probe.Rdata')

cervical_data=read_csv('0215cervical_survival_data.csv')
rownames(cervical_data)=cervical_data[,1]
cervical_data=cervical_data[,-1]
cervical_data=cervical_data[pathway_probe[which(pathway_probe %in% rownames(cervical_data))],]

save(cervical_data,file='cervical_pathwaydata.Rdata')
cervical_data[is.na(cervical_data)]=0
different_data=data.frame(matrix(ncol = dim(cervical_data)[2],nrow = 0))
for (i in c(1:dim(cervical_data)[1])){
  row=as.numeric(cervical_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(cervical_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(cervical_data)
save(sum_different,file = 'cervical_different.Rdata')
##########################################
colonic_data=read_csv('colonic_survival_data.csv')
rownames(colonic_data)=colonic_data[,1]
colonic_data=colonic_data[,-1]
colonic_data=colonic_data[pathway_probe[which(pathway_probe %in% rownames(colonic_data))],]

save(colonic_data,file='colonic_pathwaydata.Rdata')
colonic_data[is.na(colonic_data)]=0
different_data=data.frame(matrix(ncol = dim(colonic_data)[2],nrow = 0))
for (i in c(1:dim(colonic_data)[1])){
  row=as.numeric(colonic_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(colonic_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(colonic_data)
save(sum_different,file = 'colonic_different.Rdata')
##########################################
gastric_data=read_csv('gastric_survival_data.csv')
gastric_data=data.frame(gastric_data)
rownames(gastric_data)=gastric_data[,1]
gastric_data=gastric_data[,-1]
gastric_data=gastric_data[pathway_probe[which(pathway_probe %in% rownames(gastric_data))],]

save(gastric_data,file='gastric_pathwaydata.Rdata')
gastric_data[is.na(gastric_data)]=0
different_data=data.frame(matrix(ncol = dim(gastric_data)[2],nrow = 0))
for (i in c(1:dim(gastric_data)[1])){
  row=as.numeric(gastric_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(gastric_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(gastric_data)
save(sum_different,file = 'gastric_different.Rdata')
##########################################
liver_data=read_csv('liver_survival_data.csv')
liver_data=data.frame(liver_data)
rownames(liver_data)=liver_data[,1]
liver_data=liver_data[,-1]
liver_data=liver_data[pathway_probe[which(pathway_probe %in% rownames(liver_data))],]

save(liver_data,file='liver_pathwaydata.Rdata')
liver_data[is.na(liver_data)]=0
different_data=data.frame(matrix(ncol = dim(liver_data)[2],nrow = 0))
for (i in c(1:dim(liver_data)[1])){
  row=as.numeric(liver_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(liver_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(liver_data)
save(sum_different,file = 'liver_different.Rdata')
##########################################
melanoma_data=read_csv('melanoma_survival_data.csv')
melanoma_data=data.frame(melanoma_data)
rownames(melanoma_data)=melanoma_data[,1]
melanoma_data=melanoma_data[,-1]
melanoma_data=melanoma_data[pathway_probe[which(pathway_probe %in% rownames(melanoma_data))],]

save(melanoma_data,file='melanoma_pathwaydata.Rdata')
melanoma_data[is.na(melanoma_data)]=0
different_data=data.frame(matrix(ncol = dim(melanoma_data)[2],nrow = 0))
for (i in c(1:dim(melanoma_data)[1])){
  row=as.numeric(melanoma_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(melanoma_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(melanoma_data)
save(sum_different,file = 'melanoma_different.Rdata')
##########################################
prostatic_data=read_tsv('TCGA-PRAD.methylation450.tsv')
prostatic_data=data.frame(prostatic_data)
rownames(prostatic_data)=prostatic_data[,1]
prostatic_data=prostatic_data[,-1]
prostatic_data=prostatic_data[pathway_probe[which(pathway_probe %in% rownames(prostatic_data))],]

save(prostatic_data,file='prostatic_pathwaydata.Rdata')
prostatic_data[is.na(prostatic_data)]=0
different_data=data.frame(matrix(ncol = dim(prostatic_data)[2],nrow = 0))
for (i in c(1:dim(prostatic_data)[1])){
  row=as.numeric(prostatic_data[i,])
  sd=sd(row)
  mean=mean(row)
  out=which(row>(mean+3*sd)|row<(mean-3*sd))
  different_row_data=rep(0,dim(prostatic_data)[2])
  different_row_data[out]=1
  different_data=rbind(different_data,different_row_data)
}
sum_different=colSums(prostatic_data)
save(sum_different,file = 'prostatic_different.Rdata')


load('cervical_different.Rdata')
cervical_different=sum_different
load('gastric_different.Rdata')
gastric_different=sum_different
load('prostatic_different.Rdata')
prostatic_different=sum_different
load('colonic_different.Rdata')
colonic_different=sum_different
load('melanoma_different.Rdata')
melanoma_different=sum_different
load('liver_different.Rdata')
liver_different=sum_different

boxplot_different=list(cervical=cervical_different,gastric=gastric_different,prostatic=prostatic_different,colonic=colonic_different,melanoma=melanoma_different,liver=liver_different)
boxplot(boxplot_different)
save(cervical_different,gastric_different,prostatic_different,colonic_different,melanoma_different,liver_different,file = '0321different.Rdata')
