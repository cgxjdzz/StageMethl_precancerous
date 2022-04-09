
library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
load('D:/甲基化/癌症发展/colonic/139404+68060/sum-champ_myNorm.Rdata')
load('0319colonic_hyper_hypo.Rdata')


background=rownames(myNorm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
colonic_background_feature=table(background$feature)/dim(background)[1]
colonic_background_cgi=table(background$cgi)/dim(background)[1]
colonic_background_feature=melt(as.data.frame(colonic_background_feature))
rownames(colonic_background_feature)=colonic_background_feature[,1]
colonic_background_feature=colonic_background_feature[,-c(1,2)]
colonic_background_cgi=melt(as.data.frame(colonic_background_cgi))
rownames(colonic_background_cgi)=colonic_background_cgi[,1]
colonic_background_cgi=colonic_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[colonic_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
colonic_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
colonic_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
colonic_norm_pre_hyper_feature=melt(as.data.frame(colonic_norm_pre_hyper_feature))
colonic_norm_pre_hyper_feature$value=colonic_norm_pre_hyper_feature$value/colonic_background_feature
colonic_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',colonic_norm_pre_hyper_feature[,2])
colonic_norm_pre_hyper_cgi=melt(as.data.frame(colonic_norm_pre_hyper_cgi))
colonic_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',colonic_norm_pre_hyper_cgi[,2])
colonic_norm_pre_hyper_cgi$value=colonic_norm_pre_hyper_cgi$value/colonic_background_cgi

norm_pre_hypo=feature[colonic_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
colonic_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
colonic_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
colonic_norm_pre_hypo_feature=melt(as.data.frame(colonic_norm_pre_hypo_feature))
colonic_norm_pre_hypo_feature$value=colonic_norm_pre_hypo_feature$value/colonic_background_feature
colonic_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',colonic_norm_pre_hypo_feature[,2])
colonic_norm_pre_hypo_cgi=melt(as.data.frame(colonic_norm_pre_hypo_cgi))
colonic_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',colonic_norm_pre_hypo_cgi[,2])
colonic_norm_pre_hypo_cgi$value=colonic_norm_pre_hypo_cgi$value/colonic_background_cgi

norm_cancer_hyper=feature[colonic_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
colonic_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
colonic_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
colonic_norm_cancer_hyper_feature=melt(as.data.frame(colonic_norm_cancer_hyper_feature))
colonic_norm_cancer_hyper_feature$value=colonic_norm_cancer_hyper_feature$value/colonic_background_feature
colonic_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',colonic_norm_cancer_hyper_feature[,2])
colonic_norm_cancer_hyper_cgi=melt(as.data.frame(colonic_norm_cancer_hyper_cgi))
colonic_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',colonic_norm_cancer_hyper_cgi[,2])
colonic_norm_cancer_hyper_cgi$value=colonic_norm_cancer_hyper_cgi$value/colonic_background_cgi

norm_cancer_hypo=feature[colonic_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
colonic_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
colonic_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
colonic_norm_cancer_hypo_feature=melt(as.data.frame(colonic_norm_cancer_hypo_feature))
colonic_norm_cancer_hypo_feature$value=colonic_norm_cancer_hypo_feature$value/colonic_background_feature
colonic_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',colonic_norm_cancer_hypo_feature[,2])
colonic_norm_cancer_hypo_cgi=melt(as.data.frame(colonic_norm_cancer_hypo_cgi))
colonic_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',colonic_norm_cancer_hypo_cgi[,2])
colonic_norm_cancer_hypo_cgi$value=colonic_norm_cancer_hypo_cgi$value/colonic_background_cgi

pre_cancer_hyper=feature[colonic_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
colonic_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
colonic_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
colonic_pre_cancer_hyper_feature=melt(as.data.frame(colonic_pre_cancer_hyper_feature))
colonic_pre_cancer_hyper_feature$value=colonic_pre_cancer_hyper_feature$value/colonic_background_feature
colonic_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',colonic_pre_cancer_hyper_feature[,2])
colonic_pre_cancer_hyper_cgi=melt(as.data.frame(colonic_pre_cancer_hyper_cgi))
colonic_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',colonic_pre_cancer_hyper_cgi[,2])
colonic_pre_cancer_hyper_cgi$value=colonic_pre_cancer_hyper_cgi$value/colonic_background_cgi

pre_cancer_hypo=feature[colonic_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
colonic_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
colonic_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
colonic_pre_cancer_hypo_feature=melt(as.data.frame(colonic_pre_cancer_hypo_feature))
colonic_pre_cancer_hypo_feature$value=colonic_pre_cancer_hypo_feature$value/colonic_background_feature
colonic_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',colonic_pre_cancer_hypo_feature[,2])
colonic_pre_cancer_hypo_cgi=melt(as.data.frame(colonic_pre_cancer_hypo_cgi))
colonic_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',colonic_pre_cancer_hypo_cgi[,2])
colonic_pre_cancer_hypo_cgi$value=colonic_pre_cancer_hypo_cgi$value/colonic_background_cgi

colonic_feature=rbind(colonic_norm_pre_hyper_feature,colonic_norm_pre_hypo_feature,colonic_norm_cancer_hyper_feature,colonic_norm_cancer_hypo_feature,colonic_pre_cancer_hyper_feature,colonic_pre_cancer_hypo_feature)
colonic_feature$value[is.na(colonic_feature$value)]=0
colonic_feature_duijiplot=colonic_feature
colonic_feature_duijiplot$cancer='colonic'
colonic_feature=dcast(colonic_feature,Var1~variable,value.var = 'value')
rownames(colonic_feature)=colonic_feature$Var1
colonic_feature=colonic_feature[,-1]
pdf('colonic_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(colonic_feature),main = "colonic",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(colonic_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


colonic_cgi=rbind(colonic_norm_pre_hyper_cgi,colonic_norm_pre_hypo_cgi,colonic_norm_cancer_hyper_cgi,colonic_norm_cancer_hypo_cgi,colonic_pre_cancer_hyper_cgi,colonic_pre_cancer_hypo_cgi)
colonic_cgi$value[is.na(colonic_cgi$value)]=0
colonic_cgi_duijiplot=colonic_cgi
colonic_cgi_duijiplot$cancer='colonic'
colonic_cgi=dcast(colonic_cgi,Var1~variable,value.var = 'value')
rownames(colonic_cgi)=colonic_cgi$Var1
colonic_cgi=colonic_cgi[,-1]
pdf('colonic_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(colonic_cgi),main = "colonic",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(colonic_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()



library(ggplot2)
ggplot(colonic_cgi)
