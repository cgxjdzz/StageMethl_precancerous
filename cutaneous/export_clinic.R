library(readr)
survival=read.table('survival.TXT',fill=T)
colnames(survival)=survival[1,]
survival=survival[-1,]
temp=read_tsv('TCGA-SKCM.methylation450.tsv')
temp=as.data.frame(temp)
rownames(temp)=temp[,1]
temp=temp[,-1]

rownames(survival)=survival$sample
survival=survival[colnames(temp),]
survival=data.frame(sample=survival$sample,OS.time=survival$OS.time)
write.csv(survival,file = 'prostatic_survival_OSdata.csv')

write.csv(temp,file = 'prostatic_survival_data.csv')


nana=survival[which(survival$OS==1),c('_PATIENT','OS.time')]
nana=unique(nana)
colnames(nana)=c('X_PATIENT','OS.time')
rownames(nana)=nana$X_PATIENT

temp_sample=c((substr(colnames(temp),0,12)))
temp_sample=temp_sample[-1]
OS.data=nana[temp_sample,]
rownames(OS.data)=colnames(temp)[-1]
OS.data=na.omit(OS.data)

sample_list=rownames(OS.data)
data=temp[,sample_list]
rownames(data)=temp[,1]
write.csv(OS.data,file = 'melanoma_survival_OSdata.csv')
write.csv(data,file='melanoma_survival_data.csv')
