library(ChAMP)
library(reshape)
data(probe.features.epic)#850K
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
load('D:/甲基化/癌症发展/宫颈癌/0128sum_combat.Rdata')
load('0319cervical_hyper_hypo.Rdata')

myNorm=sum_combat
background=rownames(myNorm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
cervical_background_feature=table(background$feature)/dim(background)[1]
cervical_background_cgi=table(background$cgi)/dim(background)[1]
cervical_background_feature=melt(as.data.frame(cervical_background_feature))
rownames(cervical_background_feature)=cervical_background_feature[,1]
cervical_background_feature=cervical_background_feature[,-c(1,2)]
cervical_background_cgi=melt(as.data.frame(cervical_background_cgi))
rownames(cervical_background_cgi)=cervical_background_cgi[,1]
cervical_background_cgi=cervical_background_cgi[,-c(1,2)]


norm_pre_hyper=feature[cervical_norm_pre_hyper,]
norm_pre_hyper=na.omit(data.frame(feature=norm_pre_hyper$feature,cgi=norm_pre_hyper$cgi))
cervical_norm_pre_hyper_feature=table(norm_pre_hyper$feature)/dim(norm_pre_hyper)[1]
cervical_norm_pre_hyper_cgi=table(norm_pre_hyper$cgi)/dim(norm_pre_hyper)[1]
cervical_norm_pre_hyper_feature=melt(as.data.frame(cervical_norm_pre_hyper_feature))
cervical_norm_pre_hyper_feature$value=cervical_norm_pre_hyper_feature$value/cervical_background_feature
cervical_norm_pre_hyper_feature[,2]=gsub('Freq','norm_pre_hyper',cervical_norm_pre_hyper_feature[,2])
cervical_norm_pre_hyper_cgi=melt(as.data.frame(cervical_norm_pre_hyper_cgi))
cervical_norm_pre_hyper_cgi[,2]=gsub('Freq','norm_pre_hyper',cervical_norm_pre_hyper_cgi[,2])
cervical_norm_pre_hyper_cgi$value=cervical_norm_pre_hyper_cgi$value/cervical_background_cgi

norm_pre_hypo=feature[cervical_norm_pre_hypo,]
norm_pre_hypo=na.omit(data.frame(feature=norm_pre_hypo$feature,cgi=norm_pre_hypo$cgi))
cervical_norm_pre_hypo_feature=table(norm_pre_hypo$feature)/dim(norm_pre_hypo)[1]
cervical_norm_pre_hypo_cgi=table(norm_pre_hypo$cgi)/dim(norm_pre_hypo)[1]
cervical_norm_pre_hypo_feature=melt(as.data.frame(cervical_norm_pre_hypo_feature))
cervical_norm_pre_hypo_feature$value=cervical_norm_pre_hypo_feature$value/cervical_background_feature
cervical_norm_pre_hypo_feature[,2]=gsub('Freq','norm_pre_hypo',cervical_norm_pre_hypo_feature[,2])
cervical_norm_pre_hypo_cgi=melt(as.data.frame(cervical_norm_pre_hypo_cgi))
cervical_norm_pre_hypo_cgi[,2]=gsub('Freq','norm_pre_hypo',cervical_norm_pre_hypo_cgi[,2])
cervical_norm_pre_hypo_cgi$value=cervical_norm_pre_hypo_cgi$value/cervical_background_cgi

norm_cancer_hyper=feature[cervical_norm_cancer_hyper,]
norm_cancer_hyper=na.omit(data.frame(feature=norm_cancer_hyper$feature,cgi=norm_cancer_hyper$cgi))
cervical_norm_cancer_hyper_feature=table(norm_cancer_hyper$feature)/dim(norm_cancer_hyper)[1]
cervical_norm_cancer_hyper_cgi=table(norm_cancer_hyper$cgi)/dim(norm_cancer_hyper)[1]
cervical_norm_cancer_hyper_feature=melt(as.data.frame(cervical_norm_cancer_hyper_feature))
cervical_norm_cancer_hyper_feature$value=cervical_norm_cancer_hyper_feature$value/cervical_background_feature
cervical_norm_cancer_hyper_feature[,2]=gsub('Freq','norm_cancer_hyper',cervical_norm_cancer_hyper_feature[,2])
cervical_norm_cancer_hyper_cgi=melt(as.data.frame(cervical_norm_cancer_hyper_cgi))
cervical_norm_cancer_hyper_cgi[,2]=gsub('Freq','norm_cancer_hyper',cervical_norm_cancer_hyper_cgi[,2])
cervical_norm_cancer_hyper_cgi$value=cervical_norm_cancer_hyper_cgi$value/cervical_background_cgi

norm_cancer_hypo=feature[cervical_norm_cancer_hypo,]
norm_cancer_hypo=na.omit(data.frame(feature=norm_cancer_hypo$feature,cgi=norm_cancer_hypo$cgi))
cervical_norm_cancer_hypo_feature=table(norm_cancer_hypo$feature)/dim(norm_cancer_hypo)[1]
cervical_norm_cancer_hypo_cgi=table(norm_cancer_hypo$cgi)/dim(norm_cancer_hypo)[1]
cervical_norm_cancer_hypo_feature=melt(as.data.frame(cervical_norm_cancer_hypo_feature))
cervical_norm_cancer_hypo_feature$value=cervical_norm_cancer_hypo_feature$value/cervical_background_feature
cervical_norm_cancer_hypo_feature[,2]=gsub('Freq','norm_cancer_hypo',cervical_norm_cancer_hypo_feature[,2])
cervical_norm_cancer_hypo_cgi=melt(as.data.frame(cervical_norm_cancer_hypo_cgi))
cervical_norm_cancer_hypo_cgi[,2]=gsub('Freq','norm_cancer_hypo',cervical_norm_cancer_hypo_cgi[,2])
cervical_norm_cancer_hypo_cgi$value=cervical_norm_cancer_hypo_cgi$value/cervical_background_cgi

pre_cancer_hyper=feature[cervical_pre_cancer_hyper,]
pre_cancer_hyper=na.omit(data.frame(feature=pre_cancer_hyper$feature,cgi=pre_cancer_hyper$cgi))
cervical_pre_cancer_hyper_feature=table(pre_cancer_hyper$feature)/dim(pre_cancer_hyper)[1]
cervical_pre_cancer_hyper_cgi=table(pre_cancer_hyper$cgi)/dim(pre_cancer_hyper)[1]
cervical_pre_cancer_hyper_feature=melt(as.data.frame(cervical_pre_cancer_hyper_feature))
cervical_pre_cancer_hyper_feature$value=cervical_pre_cancer_hyper_feature$value/cervical_background_feature
cervical_pre_cancer_hyper_feature[,2]=gsub('Freq','pre_cancer_hyper',cervical_pre_cancer_hyper_feature[,2])
cervical_pre_cancer_hyper_cgi=melt(as.data.frame(cervical_pre_cancer_hyper_cgi))
cervical_pre_cancer_hyper_cgi[,2]=gsub('Freq','pre_cancer_hyper',cervical_pre_cancer_hyper_cgi[,2])
cervical_pre_cancer_hyper_cgi$value=cervical_pre_cancer_hyper_cgi$value/cervical_background_cgi

pre_cancer_hypo=feature[cervical_pre_cancer_hypo,]
pre_cancer_hypo=na.omit(data.frame(feature=pre_cancer_hypo$feature,cgi=pre_cancer_hypo$cgi))
cervical_pre_cancer_hypo_feature=table(pre_cancer_hypo$feature)/dim(pre_cancer_hypo)[1]
cervical_pre_cancer_hypo_cgi=table(pre_cancer_hypo$cgi)/dim(pre_cancer_hypo)[1]
cervical_pre_cancer_hypo_feature=melt(as.data.frame(cervical_pre_cancer_hypo_feature))
cervical_pre_cancer_hypo_feature$value=cervical_pre_cancer_hypo_feature$value/cervical_background_feature
cervical_pre_cancer_hypo_feature[,2]=gsub('Freq','pre_cancer_hypo',cervical_pre_cancer_hypo_feature[,2])
cervical_pre_cancer_hypo_cgi=melt(as.data.frame(cervical_pre_cancer_hypo_cgi))
cervical_pre_cancer_hypo_cgi[,2]=gsub('Freq','pre_cancer_hypo',cervical_pre_cancer_hypo_cgi[,2])
cervical_pre_cancer_hypo_cgi$value=cervical_pre_cancer_hypo_cgi$value/cervical_background_cgi


cervical_feature=rbind(cervical_norm_pre_hyper_feature,cervical_norm_pre_hypo_feature,cervical_norm_cancer_hyper_feature,cervical_norm_cancer_hypo_feature,cervical_pre_cancer_hyper_feature,cervical_pre_cancer_hypo_feature)
cervical_feature$value[is.na(cervical_feature$value)]=0
cervical_feature_duijiplot=cervical_feature
cervical_feature_duijiplot$cancer='cervical'
cervical_feature=dcast(cervical_feature,Var1~variable,value.var = 'value')
rownames(cervical_feature)=cervical_feature$Var1
cervical_feature=cervical_feature[,-1]
pdf('cervical_feature.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(cervical_feature),main = "cervical",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(cervical_feature),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()


cervical_cgi=rbind(cervical_norm_pre_hyper_cgi,cervical_norm_pre_hypo_cgi,cervical_norm_cancer_hyper_cgi,cervical_norm_cancer_hypo_cgi,cervical_pre_cancer_hyper_cgi,cervical_pre_cancer_hypo_cgi)
cervical_cgi$value[is.na(cervical_cgi$value)]=0
cervical_cgi_duijiplot=cervical_cgi
cervical_cgi_duijiplot$cancer='cervical'
cervical_cgi=dcast(cervical_cgi,Var1~variable,value.var = 'value')
rownames(cervical_cgi)=cervical_cgi$Var1
cervical_cgi=cervical_cgi[,-1]
pdf('cervical_cgi.pdf', width=8.51, height=5.49)
par(mar=c(8,4,1,2))
barplot(as.matrix(cervical_cgi),main = "cervical",
        xlab = "", ylab = "Score",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(cervical_cgi),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7,args.legend=list(y=500))
dev.off()




