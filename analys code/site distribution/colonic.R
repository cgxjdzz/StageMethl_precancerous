library(ChAMP)
library(reshape2)
data(probe.features.epic)#850K
#data(probe.features)#450K

load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hyper_and_hypo.Rdata')
colnames(pd)='label'
library(ChAMP)
library(plyr)
pd0204=pd
pd0204$label=mapvalues(pd$label, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = sum_norm,pheno = pd0204$label)
hypo_hyper=colonic_hypo_hyper
hyper_hypo=colonic_hyper_hypo
hyper_hyper=colonic_hyper_hyper
hypo_hypo=colonic_hypo_hypo
hypo_stable=colonic_hypo_stable
hyper_stable=colonic_hyper_stable
hyper_back_to_normal=colonic_hyper_back_to_normal
hypo_back_to_normal=colonic_hypo_back_to_normal

normal_pre=rownames(sum_DMP$Normal_to_Precancer)

normal_cancer=rownames(sum_DMP$Tumor_to_Normal)


background=rownames(sum_norm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

normal_pre=feature[normal_pre,]
normal_pre=na.omit(data.frame(feature=normal_pre$feature,cgi=normal_pre$cgi))
colonic_normal_pre_feature=table(normal_pre$feature)/dim(normal_pre)[1]
colonic_normal_pre_cgi=table(normal_pre$cgi)/dim(normal_pre)[1]
colonic_normal_pre_feature=melt(as.data.frame(colonic_normal_pre_feature))
colonic_normal_pre_feature[,2]=gsub('Freq','colonic',colonic_normal_pre_feature[,2])
colonic_normal_pre_cgi=melt(as.data.frame(colonic_normal_pre_cgi))
colonic_normal_pre_cgi[,2]=gsub('Freq','colonic',colonic_normal_pre_cgi[,2])

normal_cancer=feature[normal_cancer,]
normal_cancer=na.omit(data.frame(feature=normal_cancer$feature,cgi=normal_cancer$cgi))
colonic_normal_cancer_feature=table(normal_cancer$feature)/dim(normal_cancer)[1]
colonic_normal_cancer_cgi=table(normal_cancer$cgi)/dim(normal_cancer)[1]
colonic_normal_cancer_feature=melt(as.data.frame(colonic_normal_cancer_feature))
colonic_normal_cancer_feature[,2]=gsub('Freq','colonic',colonic_normal_cancer_feature[,2])
colonic_normal_cancer_cgi=melt(as.data.frame(colonic_normal_cancer_cgi))
colonic_normal_cancer_cgi[,2]=gsub('Freq','colonic',colonic_normal_cancer_cgi[,2])


hypo_hyper=feature[hypo_hyper,]
hypo_hyper=na.omit(data.frame(feature=hypo_hyper$feature,cgi=hypo_hyper$cgi))
colonic_hypo_hyper_feature=table(hypo_hyper$feature)/dim(hypo_hyper)[1]
colonic_hypo_hyper_cgi=table(hypo_hyper$cgi)/dim(hypo_hyper)[1]
colonic_hypo_hyper_feature=melt(as.data.frame(colonic_hypo_hyper_feature))
colonic_hypo_hyper_feature[,2]=gsub('Freq','colonic',colonic_hypo_hyper_feature[,2])
colonic_hypo_hyper_cgi=melt(as.data.frame(colonic_hypo_hyper_cgi))
colonic_hypo_hyper_cgi[,2]=gsub('Freq','colonic',colonic_hypo_hyper_cgi[,2])

hyper_hypo=feature[hyper_hypo,]
hyper_hypo=na.omit(data.frame(feature=hyper_hypo$feature,cgi=hyper_hypo$cgi))
colonic_hyper_hypo_feature=table(hyper_hypo$feature)/dim(hyper_hypo)[1]
colonic_hyper_hypo_cgi=table(hyper_hypo$cgi)/dim(hyper_hypo)[1]
colonic_hyper_hypo_feature=melt(as.data.frame(colonic_hyper_hypo_feature))
colonic_hyper_hypo_feature[,2]=gsub('Freq','colonic',colonic_hyper_hypo_feature[,2])
colonic_hyper_hypo_cgi=melt(as.data.frame(colonic_hyper_hypo_cgi))
colonic_hyper_hypo_cgi[,2]=gsub('Freq','colonic',colonic_hyper_hypo_cgi[,2])

hyper_hyper=feature[hyper_hyper,]
hyper_hyper=na.omit(data.frame(feature=hyper_hyper$feature,cgi=hyper_hyper$cgi))
colonic_hyper_hyper_feature=table(hyper_hyper$feature)/dim(hyper_hyper)[1]
colonic_hyper_hyper_cgi=table(hyper_hyper$cgi)/dim(hyper_hyper)[1]
colonic_hyper_hyper_feature=melt(as.data.frame(colonic_hyper_hyper_feature))
colonic_hyper_hyper_feature[,2]=gsub('Freq','colonic',colonic_hyper_hyper_feature[,2])
colonic_hyper_hyper_cgi=melt(as.data.frame(colonic_hyper_hyper_cgi))
colonic_hyper_hyper_cgi[,2]=gsub('Freq','colonic',colonic_hyper_hyper_cgi[,2])

hypo_hypo=feature[hypo_hypo,]
hypo_hypo=na.omit(data.frame(feature=hypo_hypo$feature,cgi=hypo_hypo$cgi))
colonic_hypo_hypo_feature=table(hypo_hypo$feature)/dim(hypo_hypo)[1]
colonic_hypo_hypo_cgi=table(hypo_hypo$cgi)/dim(hypo_hypo)[1]
colonic_hypo_hypo_feature=melt(as.data.frame(colonic_hypo_hypo_feature))
colonic_hypo_hypo_feature[,2]=gsub('Freq','colonic',colonic_hypo_hypo_feature[,2])
colonic_hypo_hypo_cgi=melt(as.data.frame(colonic_hypo_hypo_cgi))
colonic_hypo_hypo_cgi[,2]=gsub('Freq','colonic',colonic_hypo_hypo_cgi[,2])


background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
colonic_background_feature=table(background$feature)/dim(background)[1]
colonic_background_cgi=table(background$cgi)/dim(background)[1]
colonic_background_feature=melt(as.data.frame(colonic_background_feature))
colonic_background_feature[,2]=gsub('Freq','colonic',colonic_background_feature[,2])
colonic_background_cgi=melt(as.data.frame(colonic_background_cgi))
colonic_background_cgi[,2]=gsub('Freq','colonic',colonic_background_cgi[,2])


hyper_stable=feature[hyper_stable,]
hyper_stable=na.omit(data.frame(feature=hyper_stable$feature,cgi=hyper_stable$cgi))
colonic_hyper_stable_feature=table(hyper_stable$feature)/dim(hyper_stable)[1]
colonic_hyper_stable_cgi=table(hyper_stable$cgi)/dim(hyper_stable)[1]
colonic_hyper_stable_feature=melt(as.data.frame(colonic_hyper_stable_feature))
colonic_hyper_stable_feature[,2]=gsub('Freq','colonic',colonic_hyper_stable_feature[,2])
colonic_hyper_stable_cgi=melt(as.data.frame(colonic_hyper_stable_cgi))
colonic_hyper_stable_cgi[,2]=gsub('Freq','colonic',colonic_hyper_stable_cgi[,2])

hypo_stable=feature[hypo_stable,]
hypo_stable=na.omit(data.frame(feature=hypo_stable$feature,cgi=hypo_stable$cgi))
colonic_hypo_stable_feature=table(hypo_stable$feature)/dim(hypo_stable)[1]
colonic_hypo_stable_cgi=table(hypo_stable$cgi)/dim(hypo_stable)[1]
colonic_hypo_stable_feature=melt(as.data.frame(colonic_hypo_stable_feature))
colonic_hypo_stable_feature[,2]=gsub('Freq','colonic',colonic_hypo_stable_feature[,2])
colonic_hypo_stable_cgi=melt(as.data.frame(colonic_hypo_stable_cgi))
colonic_hypo_stable_cgi[,2]=gsub('Freq','colonic',colonic_hypo_stable_cgi[,2])

hypo_back_to_normal=feature[hypo_back_to_normal,]
hypo_back_to_normal=na.omit(data.frame(feature=hypo_back_to_normal$feature,cgi=hypo_back_to_normal$cgi))
colonic_hypo_back_to_normal_feature=table(hypo_back_to_normal$feature)/dim(hypo_back_to_normal)[1]
colonic_hypo_back_to_normal_cgi=table(hypo_back_to_normal$cgi)/dim(hypo_back_to_normal)[1]
colonic_hypo_back_to_normal_feature=melt(as.data.frame(colonic_hypo_back_to_normal_feature))
colonic_hypo_back_to_normal_feature[,2]=gsub('Freq','colonic',colonic_hypo_back_to_normal_feature[,2])
colonic_hypo_back_to_normal_cgi=melt(as.data.frame(colonic_hypo_back_to_normal_cgi))
colonic_hypo_back_to_normal_cgi[,2]=gsub('Freq','colonic',colonic_hypo_back_to_normal_cgi[,2])

hyper_back_to_normal=feature[hyper_back_to_normal,]
hyper_back_to_normal=na.omit(data.frame(feature=hyper_back_to_normal$feature,cgi=hyper_back_to_normal$cgi))
colonic_hyper_back_to_normal_feature=table(hyper_back_to_normal$feature)/dim(hyper_back_to_normal)[1]
colonic_hyper_back_to_normal_cgi=table(hyper_back_to_normal$cgi)/dim(hyper_back_to_normal)[1]
colonic_hyper_back_to_normal_feature=melt(as.data.frame(colonic_hyper_back_to_normal_feature))
colonic_hyper_back_to_normal_feature[,2]=gsub('Freq','colonic',colonic_hyper_back_to_normal_feature[,2])
colonic_hyper_back_to_normal_cgi=melt(as.data.frame(colonic_hyper_back_to_normal_cgi))
colonic_hyper_back_to_normal_cgi[,2]=gsub('Freq','colonic',colonic_hyper_back_to_normal_cgi[,2])






save(colonic_background_cgi,colonic_normal_cancer_cgi,colonic_hypo_hypo_cgi,colonic_hypo_hyper_cgi,colonic_normal_pre_cgi,colonic_hyper_hyper_cgi,colonic_hyper_hypo_cgi,colonic_hypo_stable_cgi,colonic_hyper_stable_cgi,colonic_hyper_back_to_normal_cgi,colonic_hypo_back_to_normal_cgi,file='colonic_cgi_rate.Rdata')
save(colonic_background_feature,colonic_normal_cancer_feature,colonic_hypo_hypo_feature,colonic_hypo_hyper_feature,colonic_normal_pre_feature,colonic_hyper_hyper_feature,colonic_hyper_hypo_feature,colonic_hypo_stable_feature,colonic_hyper_stable_feature,colonic_hyper_back_to_normal_feature,colonic_hypo_back_to_normal_feature,file='colonic_feature_rate.Rdata')

colonic_background_cgi$variable='background'
colonic_normal_pre_cgi$variable='normal_v_pre'
colonic_normal_cancer_cgi$variable='pre_v_cancer'
colonic_hypo_hypo_cgi$variable='hypo_hypo'
colonic_hyper_hyper_cgi$variable='hyper_hyper'
colonic_hypo_hyper_cgi$variable='hypo_hyper'
colonic_hyper_hypo_cgi$variable='hyper_hypo'
colonic_hypo_stable_cgi$variable='hypo_stable'
colonic_hyper_stable_cgi$variable='hyper_stable'
colonic_hypo_back_to_normal_cgi$variable='hypo_RtoN'
colonic_hyper_back_to_normal_cgi$variable='hyper_RtoN'
colonic_barplotdata=rbind(colonic_background_cgi,colonic_normal_pre_cgi,colonic_normal_cancer_cgi,colonic_hypo_hypo_cgi,colonic_hyper_hyper_cgi,colonic_hypo_hyper_cgi,colonic_hyper_hypo_cgi,colonic_hypo_stable_cgi,colonic_hyper_stable_cgi,colonic_hypo_back_to_normal_cgi,colonic_hyper_back_to_normal_cgi)
colonic_cgi_barplotdata=dcast(colonic_barplotdata,Var1~variable,value.var = 'value')
rownames(colonic_cgi_barplotdata)=colonic_cgi_barplotdata[,1]
colonic_cgi_barplotdata=as.matrix(colonic_cgi_barplotdata[,-1])
pdf('colonic_cgi.pdf', width=8.51, height=5.49)
barplot(colonic_cgi_barplotdata,main = "colonic",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(colonic_cgi_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7)
dev.off()


colonic_background_feature$variable='background'
colonic_normal_pre_feature$variable='normal_v_pre'
colonic_normal_cancer_feature$variable='pre_v_cancer'
colonic_hypo_hypo_feature$variable='hypo_hypo'
colonic_hyper_hyper_feature$variable='hyper_hyper'
colonic_hypo_hyper_feature$variable='hypo_hyper'
colonic_hyper_hypo_feature$variable='hyper_hypo'
colonic_hypo_stable_feature$variable='hypo_stable'
colonic_hyper_stable_feature$variable='hyper_stable'
colonic_hypo_back_to_normal_feature$variable='hypo_RtoN'
colonic_hyper_back_to_normal_feature$variable='hyper_RtoN'
colonic_barplotdata=rbind(colonic_background_feature,colonic_normal_pre_feature,colonic_normal_cancer_feature,colonic_hypo_hypo_feature,colonic_hyper_hyper_feature,colonic_hypo_hyper_feature,colonic_hyper_hypo_feature,colonic_hypo_stable_feature,colonic_hyper_stable_feature,colonic_hypo_back_to_normal_feature,colonic_hyper_back_to_normal_feature)
colonic_feature_barplotdata=dcast(colonic_barplotdata,Var1~variable,value.var = 'value')
rownames(colonic_feature_barplotdata)=colonic_feature_barplotdata[,1]
colonic_feature_barplotdata=as.matrix(colonic_feature_barplotdata[,-1])
pdf('colonic_feature.pdf', width=8.51, height=5.49)
barplot(colonic_feature_barplotdata,main = "colonic",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(colonic_feature_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7, args.legend = list(y=0.5,cex=0.8))
#locator(1) 
dev.off()
