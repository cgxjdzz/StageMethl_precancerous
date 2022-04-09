library(readr)
survival=read_tsv('TCGA-STAD.survival.tsv')
temp=read_tsv('TCGA-STAD.methylation450.tsv')
temp=as.data.frame(temp)
rownames(temp)=temp[,1]
temp=temp[,-1]

rownames(survival)=survival$sample
survival=survival[colnames(temp),]
survival=data.frame(sample=survival$sample,OS.time=survival$OS.time)
survival=na.omit(survival)

write.csv(survival,file = 'gastric_survival_OSdata.csv')

temp=temp[,survival$sample]
write.csv(temp,file = 'gastric_survival_data.csv')


# nana=survival[which(survival$OS==1),c('_PATIENT','OS.time')]
# nana=unique(nana)
# rownames(nana)=nana$X_PATIENT
# 
# temp_sample=c((substr(colnames(temp),0,12)))
# temp_sample=temp_sample[-1]
# OS.data=nana[temp_sample,]
# rownames(OS.data)=colnames(temp)[-1]
# OS.data=na.omit(OS.data)
# 
# sample_list=rownames(OS.data)
# data=temp