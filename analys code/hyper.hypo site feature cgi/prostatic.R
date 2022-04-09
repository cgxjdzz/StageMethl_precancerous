
library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
load('0319prostatic_hyper_hypo.Rdata')
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')

background=rownames(sum_norm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
prostatic_background_feature=table(background$feature)/dim(background)[1]
prostatic_background_cgi=table(background$cgi)/dim(background)[1]
prostatic_background_feature=melt(as.data.frame(prostatic_background_feature))
rownames(prostatic_background_feature)=prostatic_background_feature[,1]
prostatic_background_feature=prostatic_background_feature[,-c(1,2)]
prostatic_background_cgi=melt(as.data.frame(prostatic_background_cgi))
rownames(prostatic_background_cgi)=prostatic_background_cgi[,1]
prostatic_background_cgi=prostatic_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[prostatic_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
prostatic_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
prostatic_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
prostatic_norm_pre_hyper_feature=melt(as.data.frame(prostatic_norm_pre_hyper_feature))
prostatic_norm_pre_hyper_feature$value=prostatic_norm_pre_hyper_feature$value/prostatic_background_feature
prostatic_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',prostatic_norm_pre_hyper_feature[,2])
prostatic_norm_pre_hyper_cgi=melt(as.data.frame(prostatic_norm_pre_hyper_cgi))
prostatic_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',prostatic_norm_pre_hyper_cgi[,2])
prostatic_norm_pre_hyper_cgi$value=prostatic_norm_pre_hyper_cgi$value/prostatic_background_cgi

norm_pre_hypo=feature[prostatic_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
prostatic_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
prostatic_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
prostatic_norm_pre_hypo_feature=melt(as.data.frame(prostatic_norm_pre_hypo_feature))
prostatic_norm_pre_hypo_feature$value=prostatic_norm_pre_hypo_feature$value/prostatic_background_feature
prostatic_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',prostatic_norm_pre_hypo_feature[,2])
prostatic_norm_pre_hypo_cgi=melt(as.data.frame(prostatic_norm_pre_hypo_cgi))
prostatic_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',prostatic_norm_pre_hypo_cgi[,2])
prostatic_norm_pre_hypo_cgi$value=prostatic_norm_pre_hypo_cgi$value/prostatic_background_cgi

norm_cancer_hyper=feature[prostatic_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
prostatic_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
prostatic_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
prostatic_norm_cancer_hyper_feature=melt(as.data.frame(prostatic_norm_cancer_hyper_feature))
prostatic_norm_cancer_hyper_feature$value=prostatic_norm_cancer_hyper_feature$value/prostatic_background_feature
prostatic_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',prostatic_norm_cancer_hyper_feature[,2])
prostatic_norm_cancer_hyper_cgi=melt(as.data.frame(prostatic_norm_cancer_hyper_cgi))
prostatic_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',prostatic_norm_cancer_hyper_cgi[,2])
prostatic_norm_cancer_hyper_cgi$value=prostatic_norm_cancer_hyper_cgi$value/prostatic_background_cgi

norm_cancer_hypo=feature[prostatic_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
prostatic_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
prostatic_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
prostatic_norm_cancer_hypo_feature=melt(as.data.frame(prostatic_norm_cancer_hypo_feature))
prostatic_norm_cancer_hypo_feature$value=prostatic_norm_cancer_hypo_feature$value/prostatic_background_feature
prostatic_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',prostatic_norm_cancer_hypo_feature[,2])
prostatic_norm_cancer_hypo_cgi=melt(as.data.frame(prostatic_norm_cancer_hypo_cgi))
prostatic_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',prostatic_norm_cancer_hypo_cgi[,2])
prostatic_norm_cancer_hypo_cgi$value=prostatic_norm_cancer_hypo_cgi$value/prostatic_background_cgi

pre_cancer_hyper=feature[prostatic_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
prostatic_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
prostatic_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
prostatic_pre_cancer_hyper_feature=melt(as.data.frame(prostatic_pre_cancer_hyper_feature))
prostatic_pre_cancer_hyper_feature$value=prostatic_pre_cancer_hyper_feature$value/prostatic_background_feature
prostatic_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',prostatic_pre_cancer_hyper_feature[,2])
prostatic_pre_cancer_hyper_cgi=melt(as.data.frame(prostatic_pre_cancer_hyper_cgi))
prostatic_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',prostatic_pre_cancer_hyper_cgi[,2])
prostatic_pre_cancer_hyper_cgi$value=prostatic_pre_cancer_hyper_cgi$value/prostatic_background_cgi

pre_cancer_hypo=feature[prostatic_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
prostatic_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
prostatic_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
prostatic_pre_cancer_hypo_feature=melt(as.data.frame(prostatic_pre_cancer_hypo_feature))
prostatic_pre_cancer_hypo_feature$value=prostatic_pre_cancer_hypo_feature$value/prostatic_background_feature
prostatic_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',prostatic_pre_cancer_hypo_feature[,2])
prostatic_pre_cancer_hypo_cgi=melt(as.data.frame(prostatic_pre_cancer_hypo_cgi))
prostatic_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',prostatic_pre_cancer_hypo_cgi[,2])
prostatic_pre_cancer_hypo_cgi$value=prostatic_pre_cancer_hypo_cgi$value/prostatic_background_cgi


prostatic_feature=rbind(prostatic_norm_pre_hyper_feature,prostatic_norm_pre_hypo_feature,prostatic_norm_cancer_hyper_feature,prostatic_norm_cancer_hypo_feature,prostatic_pre_cancer_hyper_feature,prostatic_pre_cancer_hypo_feature)
prostatic_feature$value[is.na(prostatic_feature$value)]=0
prostatic_feature_duijiplot=prostatic_feature
prostatic_feature_duijiplot$cancer='prostatic'
prostatic_feature=dcast(prostatic_feature,Var1~variable,value.var = 'value')
rownames(prostatic_feature)=prostatic_feature$Var1
prostatic_feature=prostatic_feature[,-1]
pdf('prostatic_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(prostatic_feature),main = "prostatic",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(prostatic_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


prostatic_cgi=rbind(prostatic_norm_pre_hyper_cgi,prostatic_norm_pre_hypo_cgi,prostatic_norm_cancer_hyper_cgi,prostatic_norm_cancer_hypo_cgi,prostatic_pre_cancer_hyper_cgi,prostatic_pre_cancer_hypo_cgi)
prostatic_cgi$value[is.na(prostatic_cgi$value)]=0
prostatic_cgi_duijiplot=prostatic_cgi
prostatic_cgi_duijiplot$cancer='prostatic'
prostatic_cgi=dcast(prostatic_cgi,Var1~variable,value.var = 'value')
rownames(prostatic_cgi)=prostatic_cgi$Var1
prostatic_cgi=prostatic_cgi[,-1]
pdf('prostatic_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(prostatic_cgi),main = "prostatic",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(prostatic_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()





