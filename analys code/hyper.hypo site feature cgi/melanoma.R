
library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
load('0319melanoma_hyper_hypo.Rdata')
load('D:/甲基化/癌症发展/黑色素瘤/1216_sum_norm.Rdata')

background=rownames(sum_norm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
melanoma_background_feature=table(background$feature)/dim(background)[1]
melanoma_background_cgi=table(background$cgi)/dim(background)[1]
melanoma_background_feature=melt(as.data.frame(melanoma_background_feature))
rownames(melanoma_background_feature)=melanoma_background_feature[,1]
melanoma_background_feature=melanoma_background_feature[,-c(1,2)]
melanoma_background_cgi=melt(as.data.frame(melanoma_background_cgi))
rownames(melanoma_background_cgi)=melanoma_background_cgi[,1]
melanoma_background_cgi=melanoma_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[melanoma_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
melanoma_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
melanoma_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
melanoma_norm_pre_hyper_feature=melt(as.data.frame(melanoma_norm_pre_hyper_feature))
melanoma_norm_pre_hyper_feature$value=melanoma_norm_pre_hyper_feature$value/melanoma_background_feature
melanoma_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',melanoma_norm_pre_hyper_feature[,2])
melanoma_norm_pre_hyper_cgi=melt(as.data.frame(melanoma_norm_pre_hyper_cgi))
melanoma_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',melanoma_norm_pre_hyper_cgi[,2])
melanoma_norm_pre_hyper_cgi$value=melanoma_norm_pre_hyper_cgi$value/melanoma_background_cgi

norm_pre_hypo=feature[melanoma_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
melanoma_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
melanoma_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
melanoma_norm_pre_hypo_feature=melt(as.data.frame(melanoma_norm_pre_hypo_feature))
melanoma_norm_pre_hypo_feature$value=melanoma_norm_pre_hypo_feature$value/melanoma_background_feature
melanoma_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',melanoma_norm_pre_hypo_feature[,2])
melanoma_norm_pre_hypo_cgi=melt(as.data.frame(melanoma_norm_pre_hypo_cgi))
melanoma_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',melanoma_norm_pre_hypo_cgi[,2])
melanoma_norm_pre_hypo_cgi$value=melanoma_norm_pre_hypo_cgi$value/melanoma_background_cgi

norm_cancer_hyper=feature[melanoma_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
melanoma_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
melanoma_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
melanoma_norm_cancer_hyper_feature=melt(as.data.frame(melanoma_norm_cancer_hyper_feature))
melanoma_norm_cancer_hyper_feature$value=melanoma_norm_cancer_hyper_feature$value/melanoma_background_feature
melanoma_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',melanoma_norm_cancer_hyper_feature[,2])
melanoma_norm_cancer_hyper_cgi=melt(as.data.frame(melanoma_norm_cancer_hyper_cgi))
melanoma_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',melanoma_norm_cancer_hyper_cgi[,2])
melanoma_norm_cancer_hyper_cgi$value=melanoma_norm_cancer_hyper_cgi$value/melanoma_background_cgi

norm_cancer_hypo=feature[melanoma_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
melanoma_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
melanoma_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
melanoma_norm_cancer_hypo_feature=melt(as.data.frame(melanoma_norm_cancer_hypo_feature))
melanoma_norm_cancer_hypo_feature$value=melanoma_norm_cancer_hypo_feature$value/melanoma_background_feature
melanoma_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',melanoma_norm_cancer_hypo_feature[,2])
melanoma_norm_cancer_hypo_cgi=melt(as.data.frame(melanoma_norm_cancer_hypo_cgi))
melanoma_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',melanoma_norm_cancer_hypo_cgi[,2])
melanoma_norm_cancer_hypo_cgi$value=melanoma_norm_cancer_hypo_cgi$value/melanoma_background_cgi

pre_cancer_hyper=feature[melanoma_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
melanoma_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
melanoma_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
melanoma_pre_cancer_hyper_feature=melt(as.data.frame(melanoma_pre_cancer_hyper_feature))
melanoma_pre_cancer_hyper_feature$value=melanoma_pre_cancer_hyper_feature$value/melanoma_background_feature
melanoma_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',melanoma_pre_cancer_hyper_feature[,2])
melanoma_pre_cancer_hyper_cgi=melt(as.data.frame(melanoma_pre_cancer_hyper_cgi))
melanoma_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',melanoma_pre_cancer_hyper_cgi[,2])
melanoma_pre_cancer_hyper_cgi$value=melanoma_pre_cancer_hyper_cgi$value/melanoma_background_cgi

pre_cancer_hypo=feature[melanoma_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
melanoma_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
melanoma_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
melanoma_pre_cancer_hypo_feature=melt(as.data.frame(melanoma_pre_cancer_hypo_feature))
melanoma_pre_cancer_hypo_feature$value=melanoma_pre_cancer_hypo_feature$value/melanoma_background_feature
melanoma_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',melanoma_pre_cancer_hypo_feature[,2])
melanoma_pre_cancer_hypo_cgi=melt(as.data.frame(melanoma_pre_cancer_hypo_cgi))
melanoma_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',melanoma_pre_cancer_hypo_cgi[,2])
melanoma_pre_cancer_hypo_cgi$value=melanoma_pre_cancer_hypo_cgi$value/melanoma_background_cgi

melanoma_feature=rbind(melanoma_norm_pre_hyper_feature,melanoma_norm_pre_hypo_feature,melanoma_norm_cancer_hyper_feature,melanoma_norm_cancer_hypo_feature,melanoma_pre_cancer_hyper_feature,melanoma_pre_cancer_hypo_feature)
melanoma_feature$value[is.na(melanoma_feature$value)]=0
melanoma_feature_duijiplot=melanoma_feature
melanoma_feature_duijiplot$cancer='melanoma'
melanoma_feature=dcast(melanoma_feature,Var1~variable,value.var = 'value')
rownames(melanoma_feature)=melanoma_feature$Var1
melanoma_feature=melanoma_feature[,-1]
pdf('melanoma_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(melanoma_feature),main = "melanoma",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(melanoma_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


melanoma_cgi=rbind(melanoma_norm_pre_hyper_cgi,melanoma_norm_pre_hypo_cgi,melanoma_norm_cancer_hyper_cgi,melanoma_norm_cancer_hypo_cgi,melanoma_pre_cancer_hyper_cgi,melanoma_pre_cancer_hypo_cgi)
melanoma_cgi$value[is.na(melanoma_cgi$value)]=0
melanoma_cgi_duijiplot=melanoma_cgi
melanoma_cgi_duijiplot$cancer='melanoma'
melanoma_cgi=dcast(melanoma_cgi,Var1~variable,value.var = 'value')
rownames(melanoma_cgi)=melanoma_cgi$Var1
melanoma_cgi=melanoma_cgi[,-1]
pdf('melanoma_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(melanoma_cgi),main = "melanoma",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(melanoma_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()



