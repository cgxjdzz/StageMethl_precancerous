
library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
load('0319liver_hyper_hypo.Rdata')
load('D:/甲基化/癌症发展/肝癌/0227_sum_norm.Rdata')

background=rownames(sum_norm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
liver_background_feature=table(background$feature)/dim(background)[1]
liver_background_cgi=table(background$cgi)/dim(background)[1]
liver_background_feature=melt(as.data.frame(liver_background_feature))
rownames(liver_background_feature)=liver_background_feature[,1]
liver_background_feature=liver_background_feature[,-c(1,2)]
liver_background_cgi=melt(as.data.frame(liver_background_cgi))
rownames(liver_background_cgi)=liver_background_cgi[,1]
liver_background_cgi=liver_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[liver_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
liver_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
liver_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
liver_norm_pre_hyper_feature=melt(as.data.frame(liver_norm_pre_hyper_feature))
liver_norm_pre_hyper_feature$value=liver_norm_pre_hyper_feature$value/liver_background_feature
liver_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',liver_norm_pre_hyper_feature[,2])
liver_norm_pre_hyper_cgi=melt(as.data.frame(liver_norm_pre_hyper_cgi))
liver_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',liver_norm_pre_hyper_cgi[,2])
liver_norm_pre_hyper_cgi$value=liver_norm_pre_hyper_cgi$value/liver_background_cgi

norm_pre_hypo=feature[liver_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
liver_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
liver_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
liver_norm_pre_hypo_feature=melt(as.data.frame(liver_norm_pre_hypo_feature))
liver_norm_pre_hypo_feature$value=liver_norm_pre_hypo_feature$value/liver_background_feature
liver_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',liver_norm_pre_hypo_feature[,2])
liver_norm_pre_hypo_cgi=melt(as.data.frame(liver_norm_pre_hypo_cgi))
liver_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',liver_norm_pre_hypo_cgi[,2])
liver_norm_pre_hypo_cgi$value=liver_norm_pre_hypo_cgi$value/liver_background_cgi

norm_cancer_hyper=feature[liver_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
liver_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
liver_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
liver_norm_cancer_hyper_feature=melt(as.data.frame(liver_norm_cancer_hyper_feature))
liver_norm_cancer_hyper_feature$value=liver_norm_cancer_hyper_feature$value/liver_background_feature
liver_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',liver_norm_cancer_hyper_feature[,2])
liver_norm_cancer_hyper_cgi=melt(as.data.frame(liver_norm_cancer_hyper_cgi))
liver_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',liver_norm_cancer_hyper_cgi[,2])
liver_norm_cancer_hyper_cgi$value=liver_norm_cancer_hyper_cgi$value/liver_background_cgi

norm_cancer_hypo=feature[liver_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
liver_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
liver_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
liver_norm_cancer_hypo_feature=melt(as.data.frame(liver_norm_cancer_hypo_feature))
liver_norm_cancer_hypo_feature$value=liver_norm_cancer_hypo_feature$value/liver_background_feature
liver_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',liver_norm_cancer_hypo_feature[,2])
liver_norm_cancer_hypo_cgi=melt(as.data.frame(liver_norm_cancer_hypo_cgi))
liver_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',liver_norm_cancer_hypo_cgi[,2])
liver_norm_cancer_hypo_cgi$value=liver_norm_cancer_hypo_cgi$value/liver_background_cgi

pre_cancer_hyper=feature[liver_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
liver_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
liver_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
liver_pre_cancer_hyper_feature=melt(as.data.frame(liver_pre_cancer_hyper_feature))
liver_pre_cancer_hyper_feature$value=liver_pre_cancer_hyper_feature$value/liver_background_feature
liver_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',liver_pre_cancer_hyper_feature[,2])
liver_pre_cancer_hyper_cgi=melt(as.data.frame(liver_pre_cancer_hyper_cgi))
liver_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',liver_pre_cancer_hyper_cgi[,2])
liver_pre_cancer_hyper_cgi$value=liver_pre_cancer_hyper_cgi$value/liver_background_cgi

pre_cancer_hypo=feature[liver_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
liver_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
liver_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
liver_pre_cancer_hypo_feature=melt(as.data.frame(liver_pre_cancer_hypo_feature))
liver_pre_cancer_hypo_feature$value=liver_pre_cancer_hypo_feature$value/liver_background_feature
liver_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',liver_pre_cancer_hypo_feature[,2])
liver_pre_cancer_hypo_cgi=melt(as.data.frame(liver_pre_cancer_hypo_cgi))
liver_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',liver_pre_cancer_hypo_cgi[,2])
liver_pre_cancer_hypo_cgi$value=liver_pre_cancer_hypo_cgi$value/liver_background_cgi


liver_feature=rbind(liver_norm_pre_hyper_feature,liver_norm_pre_hypo_feature,liver_norm_cancer_hyper_feature,liver_norm_cancer_hypo_feature,liver_pre_cancer_hyper_feature,liver_pre_cancer_hypo_feature)
liver_feature$value[is.na(liver_feature$value)]=0
liver_feature_duijiplot=liver_feature
liver_feature_duijiplot$cancer='liver'
liver_feature=dcast(liver_feature,Var1~variable,value.var = 'value')
rownames(liver_feature)=liver_feature$Var1
liver_feature=liver_feature[,-1]
pdf('liver_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(liver_feature),main = "liver",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(liver_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


liver_cgi=rbind(liver_norm_pre_hyper_cgi,liver_norm_pre_hypo_cgi,liver_norm_cancer_hyper_cgi,liver_norm_cancer_hypo_cgi,liver_pre_cancer_hyper_cgi,liver_pre_cancer_hypo_cgi)
liver_cgi$value[is.na(liver_cgi$value)]=0
liver_cgi_duijiplot=liver_cgi
liver_cgi_duijiplot$cancer='liver'
liver_cgi=dcast(liver_cgi,Var1~variable,value.var = 'value')
rownames(liver_cgi)=liver_cgi$Var1
liver_cgi=liver_cgi[,-1]
pdf('liver_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(liver_cgi),main = "liver",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(liver_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()




