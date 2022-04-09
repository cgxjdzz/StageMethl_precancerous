
library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
load('0319gastric_hyper_hypo.Rdata')
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')

background=rownames(sum_norm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
gastric_background_feature=table(background$feature)/dim(background)[1]
gastric_background_cgi=table(background$cgi)/dim(background)[1]
gastric_background_feature=melt(as.data.frame(gastric_background_feature))
rownames(gastric_background_feature)=gastric_background_feature[,1]
gastric_background_feature=gastric_background_feature[,-c(1,2)]
gastric_background_cgi=melt(as.data.frame(gastric_background_cgi))
rownames(gastric_background_cgi)=gastric_background_cgi[,1]
gastric_background_cgi=gastric_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[gastric_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
gastric_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
gastric_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
gastric_norm_pre_hyper_feature=melt(as.data.frame(gastric_norm_pre_hyper_feature))
gastric_norm_pre_hyper_feature$value=gastric_norm_pre_hyper_feature$value/gastric_background_feature
gastric_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',gastric_norm_pre_hyper_feature[,2])
gastric_norm_pre_hyper_cgi=melt(as.data.frame(gastric_norm_pre_hyper_cgi))
gastric_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',gastric_norm_pre_hyper_cgi[,2])
gastric_norm_pre_hyper_cgi$value=gastric_norm_pre_hyper_cgi$value/gastric_background_cgi

norm_pre_hypo=feature[gastric_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
gastric_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
gastric_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
gastric_norm_pre_hypo_feature=melt(as.data.frame(gastric_norm_pre_hypo_feature))
gastric_norm_pre_hypo_feature$value=gastric_norm_pre_hypo_feature$value/gastric_background_feature
gastric_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',gastric_norm_pre_hypo_feature[,2])
gastric_norm_pre_hypo_cgi=melt(as.data.frame(gastric_norm_pre_hypo_cgi))
gastric_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',gastric_norm_pre_hypo_cgi[,2])
gastric_norm_pre_hypo_cgi$value=gastric_norm_pre_hypo_cgi$value/gastric_background_cgi

norm_cancer_hyper=feature[gastric_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
gastric_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
gastric_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
gastric_norm_cancer_hyper_feature=melt(as.data.frame(gastric_norm_cancer_hyper_feature))
gastric_norm_cancer_hyper_feature$value=gastric_norm_cancer_hyper_feature$value/gastric_background_feature
gastric_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',gastric_norm_cancer_hyper_feature[,2])
gastric_norm_cancer_hyper_cgi=melt(as.data.frame(gastric_norm_cancer_hyper_cgi))
gastric_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',gastric_norm_cancer_hyper_cgi[,2])
gastric_norm_cancer_hyper_cgi$value=gastric_norm_cancer_hyper_cgi$value/gastric_background_cgi

norm_cancer_hypo=feature[gastric_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
gastric_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
gastric_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
gastric_norm_cancer_hypo_feature=melt(as.data.frame(gastric_norm_cancer_hypo_feature))
gastric_norm_cancer_hypo_feature$value=gastric_norm_cancer_hypo_feature$value/gastric_background_feature
gastric_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',gastric_norm_cancer_hypo_feature[,2])
gastric_norm_cancer_hypo_cgi=melt(as.data.frame(gastric_norm_cancer_hypo_cgi))
gastric_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',gastric_norm_cancer_hypo_cgi[,2])
gastric_norm_cancer_hypo_cgi$value=gastric_norm_cancer_hypo_cgi$value/gastric_background_cgi

pre_cancer_hyper=feature[gastric_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
gastric_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
gastric_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
gastric_pre_cancer_hyper_feature=melt(as.data.frame(gastric_pre_cancer_hyper_feature))
gastric_pre_cancer_hyper_feature$value=gastric_pre_cancer_hyper_feature$value/gastric_background_feature
gastric_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',gastric_pre_cancer_hyper_feature[,2])
gastric_pre_cancer_hyper_cgi=melt(as.data.frame(gastric_pre_cancer_hyper_cgi))
gastric_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',gastric_pre_cancer_hyper_cgi[,2])
gastric_pre_cancer_hyper_cgi$value=gastric_pre_cancer_hyper_cgi$value/gastric_background_cgi

pre_cancer_hypo=feature[gastric_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
gastric_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
gastric_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
gastric_pre_cancer_hypo_feature=melt(as.data.frame(gastric_pre_cancer_hypo_feature))
gastric_pre_cancer_hypo_feature$value=gastric_pre_cancer_hypo_feature$value/gastric_background_feature
gastric_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',gastric_pre_cancer_hypo_feature[,2])
gastric_pre_cancer_hypo_cgi=melt(as.data.frame(gastric_pre_cancer_hypo_cgi))
gastric_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',gastric_pre_cancer_hypo_cgi[,2])
gastric_pre_cancer_hypo_cgi$value=gastric_pre_cancer_hypo_cgi$value/gastric_background_cgi


gastric_feature=rbind(gastric_norm_pre_hyper_feature,gastric_norm_pre_hypo_feature,gastric_norm_cancer_hyper_feature,gastric_norm_cancer_hypo_feature,gastric_pre_cancer_hyper_feature,gastric_pre_cancer_hypo_feature)
gastric_feature$value[is.na(gastric_feature$value)]=0
gastric_feature_duijiplot=gastric_feature
gastric_feature_duijiplot$cancer='gastric'
gastric_feature=dcast(gastric_feature,Var1~variable,value.var = 'value')
rownames(gastric_feature)=gastric_feature$Var1
gastric_feature=gastric_feature[,-1]
pdf('gastric_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(gastric_feature),main = "gastric",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(gastric_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


gastric_cgi=rbind(gastric_norm_pre_hyper_cgi,gastric_norm_pre_hypo_cgi,gastric_norm_cancer_hyper_cgi,gastric_norm_cancer_hypo_cgi,gastric_pre_cancer_hyper_cgi,gastric_pre_cancer_hypo_cgi)
gastric_cgi$value[is.na(gastric_cgi$value)]=0
gastric_cgi_duijiplot=gastric_cgi
gastric_cgi_duijiplot$cancer='gastric'
gastric_cgi=dcast(gastric_cgi,Var1~variable,value.var = 'value')
rownames(gastric_cgi)=gastric_cgi$Var1
gastric_cgi=gastric_cgi[,-1]
pdf('gastric_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(gastric_cgi),main = "gastric",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(gastric_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()




