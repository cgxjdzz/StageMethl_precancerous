library(ChAMP)
library(plyr)
data(probe.features.epic)#850K
#data(probe.features)#450K

load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hyper_and_hypo.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
hypo_hyper=prostatic_hypo_hyper
hyper_hypo=prostatic_hyper_hypo
hyper_hyper=prostatic_hyper_hyper
hypo_hypo=prostatic_hypo_hypo
hypo_stable=prostatic_hypo_stable
hyper_stable=prostatic_hyper_stable
hyper_back_to_normal=prostatic_hyper_back_to_normal
hypo_back_to_normal=prostatic_hypo_back_to_normal


myNorm=sum_norm
normal_pre=rownames(sum_DMP$CAF_to_Normal)
normal_cancer=rownames(sum_DMP$Tumor_to_Normal)

myNorm=sum_norm
background=rownames(myNorm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

normal_pre=feature[normal_pre,]
normal_pre=na.omit(data.frame(feature=normal_pre$feature,cgi=normal_pre$cgi))
prostatic_normal_pre_feature=table(normal_pre$feature)/dim(normal_pre)[1]
prostatic_normal_pre_cgi=table(normal_pre$cgi)/dim(normal_pre)[1]
prostatic_normal_pre_feature=melt(as.data.frame(prostatic_normal_pre_feature))
prostatic_normal_pre_feature[,2]=gsub('Freq','prostatic',prostatic_normal_pre_feature[,2])
prostatic_normal_pre_cgi=melt(as.data.frame(prostatic_normal_pre_cgi))
prostatic_normal_pre_cgi[,2]=gsub('Freq','prostatic',prostatic_normal_pre_cgi[,2])

normal_cancer=feature[normal_cancer,]
normal_cancer=na.omit(data.frame(feature=normal_cancer$feature,cgi=normal_cancer$cgi))
prostatic_normal_cancer_feature=table(normal_cancer$feature)/dim(normal_cancer)[1]
prostatic_normal_cancer_cgi=table(normal_cancer$cgi)/dim(normal_cancer)[1]
prostatic_normal_cancer_feature=melt(as.data.frame(prostatic_normal_cancer_feature))
prostatic_normal_cancer_feature[,2]=gsub('Freq','prostatic',prostatic_normal_cancer_feature[,2])
prostatic_normal_cancer_cgi=melt(as.data.frame(prostatic_normal_cancer_cgi))
prostatic_normal_cancer_cgi[,2]=gsub('Freq','prostatic',prostatic_normal_cancer_cgi[,2])


hypo_hyper=feature[hypo_hyper,]
hypo_hyper=na.omit(data.frame(feature=hypo_hyper$feature,cgi=hypo_hyper$cgi))
prostatic_hypo_hyper_feature=table(hypo_hyper$feature)/dim(hypo_hyper)[1]
prostatic_hypo_hyper_cgi=table(hypo_hyper$cgi)/dim(hypo_hyper)[1]
prostatic_hypo_hyper_feature=melt(as.data.frame(prostatic_hypo_hyper_feature))
prostatic_hypo_hyper_feature[,2]=gsub('Freq','prostatic',prostatic_hypo_hyper_feature[,2])
prostatic_hypo_hyper_cgi=melt(as.data.frame(prostatic_hypo_hyper_cgi))
prostatic_hypo_hyper_cgi[,2]=gsub('Freq','prostatic',prostatic_hypo_hyper_cgi[,2])

hyper_hypo=feature[hyper_hypo,]
hyper_hypo=na.omit(data.frame(feature=hyper_hypo$feature,cgi=hyper_hypo$cgi))
prostatic_hyper_hypo_feature=table(hyper_hypo$feature)/dim(hyper_hypo)[1]
prostatic_hyper_hypo_cgi=table(hyper_hypo$cgi)/dim(hyper_hypo)[1]
prostatic_hyper_hypo_feature=melt(as.data.frame(prostatic_hyper_hypo_feature))
prostatic_hyper_hypo_feature[,2]=gsub('Freq','prostatic',prostatic_hyper_hypo_feature[,2])
prostatic_hyper_hypo_cgi=melt(as.data.frame(prostatic_hyper_hypo_cgi))
prostatic_hyper_hypo_cgi[,2]=gsub('Freq','prostatic',prostatic_hyper_hypo_cgi[,2])

hyper_hyper=feature[hyper_hyper,]
hyper_hyper=na.omit(data.frame(feature=hyper_hyper$feature,cgi=hyper_hyper$cgi))
prostatic_hyper_hyper_feature=table(hyper_hyper$feature)/dim(hyper_hyper)[1]
prostatic_hyper_hyper_cgi=table(hyper_hyper$cgi)/dim(hyper_hyper)[1]
prostatic_hyper_hyper_feature=melt(as.data.frame(prostatic_hyper_hyper_feature))
prostatic_hyper_hyper_feature[,2]=gsub('Freq','prostatic',prostatic_hyper_hyper_feature[,2])
prostatic_hyper_hyper_cgi=melt(as.data.frame(prostatic_hyper_hyper_cgi))
prostatic_hyper_hyper_cgi[,2]=gsub('Freq','prostatic',prostatic_hyper_hyper_cgi[,2])

hypo_hypo=feature[hypo_hypo,]
hypo_hypo=na.omit(data.frame(feature=hypo_hypo$feature,cgi=hypo_hypo$cgi))
prostatic_hypo_hypo_feature=table(hypo_hypo$feature)/dim(hypo_hypo)[1]
prostatic_hypo_hypo_cgi=table(hypo_hypo$cgi)/dim(hypo_hypo)[1]
prostatic_hypo_hypo_feature=melt(as.data.frame(prostatic_hypo_hypo_feature))
prostatic_hypo_hypo_feature[,2]=gsub('Freq','prostatic',prostatic_hypo_hypo_feature[,2])
prostatic_hypo_hypo_cgi=melt(as.data.frame(prostatic_hypo_hypo_cgi))
prostatic_hypo_hypo_cgi[,2]=gsub('Freq','prostatic',prostatic_hypo_hypo_cgi[,2])


background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
prostatic_background_feature=table(background$feature)/dim(background)[1]
prostatic_background_cgi=table(background$cgi)/dim(background)[1]
prostatic_background_feature=melt(as.data.frame(prostatic_background_feature))
prostatic_background_feature[,2]=gsub('Freq','prostatic',prostatic_background_feature[,2])
prostatic_background_cgi=melt(as.data.frame(prostatic_background_cgi))
prostatic_background_cgi[,2]=gsub('Freq','prostatic',prostatic_background_cgi[,2])



hyper_stable=feature[hyper_stable,]
hyper_stable=na.omit(data.frame(feature=hyper_stable$feature,cgi=hyper_stable$cgi))
prostatic_hyper_stable_feature=table(hyper_stable$feature)/dim(hyper_stable)[1]
prostatic_hyper_stable_cgi=table(hyper_stable$cgi)/dim(hyper_stable)[1]
prostatic_hyper_stable_feature=melt(as.data.frame(prostatic_hyper_stable_feature))
prostatic_hyper_stable_feature[,2]=gsub('Freq','prostatic',prostatic_hyper_stable_feature[,2])
prostatic_hyper_stable_cgi=melt(as.data.frame(prostatic_hyper_stable_cgi))
prostatic_hyper_stable_cgi[,2]=gsub('Freq','prostatic',prostatic_hyper_stable_cgi[,2])

hypo_stable=feature[hypo_stable,]
hypo_stable=na.omit(data.frame(feature=hypo_stable$feature,cgi=hypo_stable$cgi))
prostatic_hypo_stable_feature=table(hypo_stable$feature)/dim(hypo_stable)[1]
prostatic_hypo_stable_cgi=table(hypo_stable$cgi)/dim(hypo_stable)[1]
prostatic_hypo_stable_feature=melt(as.data.frame(prostatic_hypo_stable_feature))
prostatic_hypo_stable_feature[,2]=gsub('Freq','prostatic',prostatic_hypo_stable_feature[,2])
prostatic_hypo_stable_cgi=melt(as.data.frame(prostatic_hypo_stable_cgi))
prostatic_hypo_stable_cgi[,2]=gsub('Freq','prostatic',prostatic_hypo_stable_cgi[,2])

hypo_back_to_normal=feature[hypo_back_to_normal,]
hypo_back_to_normal=na.omit(data.frame(feature=hypo_back_to_normal$feature,cgi=hypo_back_to_normal$cgi))
prostatic_hypo_back_to_normal_feature=table(hypo_back_to_normal$feature)/dim(hypo_back_to_normal)[1]
prostatic_hypo_back_to_normal_cgi=table(hypo_back_to_normal$cgi)/dim(hypo_back_to_normal)[1]
prostatic_hypo_back_to_normal_feature=melt(as.data.frame(prostatic_hypo_back_to_normal_feature))
prostatic_hypo_back_to_normal_feature[,2]=gsub('Freq','prostatic',prostatic_hypo_back_to_normal_feature[,2])
prostatic_hypo_back_to_normal_cgi=melt(as.data.frame(prostatic_hypo_back_to_normal_cgi))
prostatic_hypo_back_to_normal_cgi[,2]=gsub('Freq','prostatic',prostatic_hypo_back_to_normal_cgi[,2])

hyper_back_to_normal=feature[hyper_back_to_normal,]
hyper_back_to_normal=na.omit(data.frame(feature=hyper_back_to_normal$feature,cgi=hyper_back_to_normal$cgi))
prostatic_hyper_back_to_normal_feature=table(hyper_back_to_normal$feature)/dim(hyper_back_to_normal)[1]
prostatic_hyper_back_to_normal_cgi=table(hyper_back_to_normal$cgi)/dim(hyper_back_to_normal)[1]
prostatic_hyper_back_to_normal_feature=melt(as.data.frame(prostatic_hyper_back_to_normal_feature))
prostatic_hyper_back_to_normal_feature[,2]=gsub('Freq','prostatic',prostatic_hyper_back_to_normal_feature[,2])
prostatic_hyper_back_to_normal_cgi=melt(as.data.frame(prostatic_hyper_back_to_normal_cgi))
prostatic_hyper_back_to_normal_cgi[,2]=gsub('Freq','prostatic',prostatic_hyper_back_to_normal_cgi[,2])




save(prostatic_background_cgi,prostatic_normal_cancer_cgi,prostatic_hypo_hypo_cgi,prostatic_hypo_hyper_cgi,prostatic_normal_pre_cgi,prostatic_hyper_hyper_cgi,prostatic_hyper_hypo_cgi,prostatic_hypo_stable_cgi,prostatic_hyper_stable_cgi,prostatic_hyper_back_to_normal_cgi,prostatic_hypo_back_to_normal_cgi,file='prostatic_cgi_rate.Rdata')
save(prostatic_background_feature,prostatic_normal_cancer_feature,prostatic_hypo_hypo_feature,prostatic_hypo_hyper_feature,prostatic_normal_pre_feature,prostatic_hyper_hyper_feature,prostatic_hyper_hypo_feature,prostatic_hypo_stable_feature,prostatic_hyper_stable_feature,prostatic_hyper_back_to_normal_feature,prostatic_hypo_back_to_normal_feature,file='prostatic_feature_rate.Rdata')

prostatic_background_cgi$variable='background'
prostatic_normal_pre_cgi$variable='normal_v_pre'
prostatic_normal_cancer_cgi$variable='pre_v_cancer'
prostatic_hypo_hypo_cgi$variable='hypo_hypo'
prostatic_hyper_hyper_cgi$variable='hyper_hyper'
prostatic_hypo_hyper_cgi$variable='hypo_hyper'
prostatic_hyper_hypo_cgi$variable='hyper_hypo'
prostatic_hypo_stable_cgi$variable='hypo_stable'
prostatic_hyper_stable_cgi$variable='hyper_stable'
prostatic_hypo_back_to_normal_cgi$variable='hypo_RtoN'
prostatic_hyper_back_to_normal_cgi$variable='hyper_RtoN'
prostatic_barplotdata=rbind(prostatic_background_cgi,prostatic_normal_pre_cgi,prostatic_normal_cancer_cgi,prostatic_hypo_hypo_cgi,prostatic_hyper_hyper_cgi,prostatic_hypo_hyper_cgi,prostatic_hyper_hypo_cgi,prostatic_hypo_stable_cgi,prostatic_hyper_stable_cgi,prostatic_hypo_back_to_normal_cgi,prostatic_hyper_back_to_normal_cgi)
prostatic_cgi_barplotdata=dcast(prostatic_barplotdata,Var1~variable,value.var = 'value')
rownames(prostatic_cgi_barplotdata)=prostatic_cgi_barplotdata[,1]
prostatic_cgi_barplotdata=as.matrix(prostatic_cgi_barplotdata[,-1])
pdf('prostatic_cgi.pdf', width=8.51, height=5.49)
barplot(prostatic_cgi_barplotdata,main = "prostatic",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(prostatic_cgi_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7)
dev.off()


prostatic_background_feature$variable='background'
prostatic_normal_pre_feature$variable='normal_v_pre'
prostatic_normal_cancer_feature$variable='pre_v_cancer'
prostatic_hypo_hypo_feature$variable='hypo_hypo'
prostatic_hyper_hyper_feature$variable='hyper_hyper'
prostatic_hypo_hyper_feature$variable='hypo_hyper'
prostatic_hyper_hypo_feature$variable='hyper_hypo'
prostatic_hypo_stable_feature$variable='hypo_stable'
prostatic_hyper_stable_feature$variable='hyper_stable'
prostatic_hypo_back_to_normal_feature$variable='hypo_RtoN'
prostatic_hyper_back_to_normal_feature$variable='hyper_RtoN'
prostatic_barplotdata=rbind(prostatic_background_feature,prostatic_normal_pre_feature,prostatic_normal_cancer_feature,prostatic_hypo_hypo_feature,prostatic_hyper_hyper_feature,prostatic_hypo_hyper_feature,prostatic_hyper_hypo_feature,prostatic_hypo_stable_feature,prostatic_hyper_stable_feature,prostatic_hypo_back_to_normal_feature,prostatic_hyper_back_to_normal_feature)
prostatic_feature_barplotdata=dcast(prostatic_barplotdata,Var1~variable,value.var = 'value')
rownames(prostatic_feature_barplotdata)=prostatic_feature_barplotdata[,1]
prostatic_feature_barplotdata=as.matrix(prostatic_feature_barplotdata[,-1])
pdf('prostatic_feature.pdf', width=8.51, height=5.49)
barplot(prostatic_feature_barplotdata,main = "prostatic",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(prostatic_feature_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7, args.legend = list(y=0.5,cex=0.8))
#locator(1) 
dev.off()