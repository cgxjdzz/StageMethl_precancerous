library(ChAMP)
library(reshape)
data(probe.features.epic)#850K
#data(probe.features)#450K
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
load('D:/甲基化/癌症发展/宫颈癌/0128sum_combat.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hyper_and_hypo.Rdata')
hypo_hyper=cervical_hypo_hyper
hyper_hypo=cervical_hyper_hypo
hyper_hyper=cervical_hyper_hyper
hypo_hypo=cervical_hypo_hypo
hypo_stable=cervical_hypo_stable
hyper_stable=cervical_hyper_stable
hyper_back_to_normal=cervical_hyper_back_to_normal
hypo_back_to_normal=cervical_hypo_back_to_normal

normal_pre=rownames(sum_DMP$normal_to_cin3)
normal_cancer=rownames(sum_DMP$normal_to_cancer)

myNorm=sum_combat
background=rownames(myNorm)
feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

normal_pre=feature[normal_pre,]
normal_pre=na.omit(data.frame(feature=normal_pre$feature,cgi=normal_pre$cgi))
cervical_normal_pre_feature=table(normal_pre$feature)/dim(normal_pre)[1]
cervical_normal_pre_cgi=table(normal_pre$cgi)/dim(normal_pre)[1]
cervical_normal_pre_feature=melt(as.data.frame(cervical_normal_pre_feature))
cervical_normal_pre_feature[,2]=gsub('Freq','cervical',cervical_normal_pre_feature[,2])
cervical_normal_pre_cgi=melt(as.data.frame(cervical_normal_pre_cgi))
cervical_normal_pre_cgi[,2]=gsub('Freq','cervical',cervical_normal_pre_cgi[,2])

normal_cancer=feature[normal_cancer,]
normal_cancer=na.omit(data.frame(feature=normal_cancer$feature,cgi=normal_cancer$cgi))
cervical_normal_cancer_feature=table(normal_cancer$feature)/dim(normal_cancer)[1]
cervical_normal_cancer_cgi=table(normal_cancer$cgi)/dim(normal_cancer)[1]
cervical_normal_cancer_feature=melt(as.data.frame(cervical_normal_cancer_feature))
cervical_normal_cancer_feature[,2]=gsub('Freq','cervical',cervical_normal_cancer_feature[,2])
cervical_normal_cancer_cgi=melt(as.data.frame(cervical_normal_cancer_cgi))
cervical_normal_cancer_cgi[,2]=gsub('Freq','cervical',cervical_normal_cancer_cgi[,2])


hypo_hyper=feature[hypo_hyper,]
hypo_hyper=na.omit(data.frame(feature=hypo_hyper$feature,cgi=hypo_hyper$cgi))
cervical_hypo_hyper_feature=table(hypo_hyper$feature)/dim(hypo_hyper)[1]
cervical_hypo_hyper_cgi=table(hypo_hyper$cgi)/dim(hypo_hyper)[1]
cervical_hypo_hyper_feature=melt(as.data.frame(cervical_hypo_hyper_feature))
cervical_hypo_hyper_feature[,2]=gsub('Freq','cervical',cervical_hypo_hyper_feature[,2])
cervical_hypo_hyper_cgi=melt(as.data.frame(cervical_hypo_hyper_cgi))
cervical_hypo_hyper_cgi[,2]=gsub('Freq','cervical',cervical_hypo_hyper_cgi[,2])

hyper_hypo=feature[hyper_hypo,]
hyper_hypo=na.omit(data.frame(feature=hyper_hypo$feature,cgi=hyper_hypo$cgi))
cervical_hyper_hypo_feature=table(hyper_hypo$feature)/dim(hyper_hypo)[1]
cervical_hyper_hypo_cgi=table(hyper_hypo$cgi)/dim(hyper_hypo)[1]
cervical_hyper_hypo_feature=melt(as.data.frame(cervical_hyper_hypo_feature))
cervical_hyper_hypo_feature[,2]=gsub('Freq','cervical',cervical_hyper_hypo_feature[,2])
cervical_hyper_hypo_cgi=melt(as.data.frame(cervical_hyper_hypo_cgi))
cervical_hyper_hypo_cgi[,2]=gsub('Freq','cervical',cervical_hyper_hypo_cgi[,2])

hyper_hyper=feature[hyper_hyper,]
hyper_hyper=na.omit(data.frame(feature=hyper_hyper$feature,cgi=hyper_hyper$cgi))
cervical_hyper_hyper_feature=table(hyper_hyper$feature)/dim(hyper_hyper)[1]
cervical_hyper_hyper_cgi=table(hyper_hyper$cgi)/dim(hyper_hyper)[1]
cervical_hyper_hyper_feature=melt(as.data.frame(cervical_hyper_hyper_feature))
cervical_hyper_hyper_feature[,2]=gsub('Freq','cervical',cervical_hyper_hyper_feature[,2])
cervical_hyper_hyper_cgi=melt(as.data.frame(cervical_hyper_hyper_cgi))
cervical_hyper_hyper_cgi[,2]=gsub('Freq','cervical',cervical_hyper_hyper_cgi[,2])

hypo_hypo=feature[hypo_hypo,]
hypo_hypo=na.omit(data.frame(feature=hypo_hypo$feature,cgi=hypo_hypo$cgi))
cervical_hypo_hypo_feature=table(hypo_hypo$feature)/dim(hypo_hypo)[1]
cervical_hypo_hypo_cgi=table(hypo_hypo$cgi)/dim(hypo_hypo)[1]
cervical_hypo_hypo_feature=melt(as.data.frame(cervical_hypo_hypo_feature))
cervical_hypo_hypo_feature[,2]=gsub('Freq','cervical',cervical_hypo_hypo_feature[,2])
cervical_hypo_hypo_cgi=melt(as.data.frame(cervical_hypo_hypo_cgi))
cervical_hypo_hypo_cgi[,2]=gsub('Freq','cervical',cervical_hypo_hypo_cgi[,2])


background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
cervical_background_feature=table(background$feature)/dim(background)[1]
cervical_background_cgi=table(background$cgi)/dim(background)[1]
cervical_background_feature=melt(as.data.frame(cervical_background_feature))
cervical_background_feature[,2]=gsub('Freq','cervical',cervical_background_feature[,2])
cervical_background_cgi=melt(as.data.frame(cervical_background_cgi))
cervical_background_cgi[,2]=gsub('Freq','cervical',cervical_background_cgi[,2])

hyper_stable=feature[hyper_stable,]
hyper_stable=na.omit(data.frame(feature=hyper_stable$feature,cgi=hyper_stable$cgi))
cervical_hyper_stable_feature=table(hyper_stable$feature)/dim(hyper_stable)[1]
cervical_hyper_stable_cgi=table(hyper_stable$cgi)/dim(hyper_stable)[1]
cervical_hyper_stable_feature=melt(as.data.frame(cervical_hyper_stable_feature))
cervical_hyper_stable_feature[,2]=gsub('Freq','cervical',cervical_hyper_stable_feature[,2])
cervical_hyper_stable_cgi=melt(as.data.frame(cervical_hyper_stable_cgi))
cervical_hyper_stable_cgi[,2]=gsub('Freq','cervical',cervical_hyper_stable_cgi[,2])

hypo_stable=feature[hypo_stable,]
hypo_stable=na.omit(data.frame(feature=hypo_stable$feature,cgi=hypo_stable$cgi))
cervical_hypo_stable_feature=table(hypo_stable$feature)/dim(hypo_stable)[1]
cervical_hypo_stable_cgi=table(hypo_stable$cgi)/dim(hypo_stable)[1]
cervical_hypo_stable_feature=melt(as.data.frame(cervical_hypo_stable_feature))
cervical_hypo_stable_feature[,2]=gsub('Freq','cervical',cervical_hypo_stable_feature[,2])
cervical_hypo_stable_cgi=melt(as.data.frame(cervical_hypo_stable_cgi))
cervical_hypo_stable_cgi[,2]=gsub('Freq','cervical',cervical_hypo_stable_cgi[,2])

hypo_back_to_normal=feature[hypo_back_to_normal,]
hypo_back_to_normal=na.omit(data.frame(feature=hypo_back_to_normal$feature,cgi=hypo_back_to_normal$cgi))
cervical_hypo_back_to_normal_feature=table(hypo_back_to_normal$feature)/dim(hypo_back_to_normal)[1]
cervical_hypo_back_to_normal_cgi=table(hypo_back_to_normal$cgi)/dim(hypo_back_to_normal)[1]
cervical_hypo_back_to_normal_feature=melt(as.data.frame(cervical_hypo_back_to_normal_feature))
cervical_hypo_back_to_normal_feature[,2]=gsub('Freq','cervical',cervical_hypo_back_to_normal_feature[,2])
cervical_hypo_back_to_normal_cgi=melt(as.data.frame(cervical_hypo_back_to_normal_cgi))
cervical_hypo_back_to_normal_cgi[,2]=gsub('Freq','cervical',cervical_hypo_back_to_normal_cgi[,2])

hyper_back_to_normal=feature[hyper_back_to_normal,]
hyper_back_to_normal=na.omit(data.frame(feature=hyper_back_to_normal$feature,cgi=hyper_back_to_normal$cgi))
cervical_hyper_back_to_normal_feature=table(hyper_back_to_normal$feature)/dim(hyper_back_to_normal)[1]
cervical_hyper_back_to_normal_cgi=table(hyper_back_to_normal$cgi)/dim(hyper_back_to_normal)[1]
cervical_hyper_back_to_normal_feature=melt(as.data.frame(cervical_hyper_back_to_normal_feature))
cervical_hyper_back_to_normal_feature[,2]=gsub('Freq','cervical',cervical_hyper_back_to_normal_feature[,2])
cervical_hyper_back_to_normal_cgi=melt(as.data.frame(cervical_hyper_back_to_normal_cgi))
cervical_hyper_back_to_normal_cgi[,2]=gsub('Freq','cervical',cervical_hyper_back_to_normal_cgi[,2])



save(cervical_background_cgi,cervical_normal_cancer_cgi,cervical_hypo_hypo_cgi,cervical_hypo_hyper_cgi,cervical_normal_pre_cgi,cervical_hyper_hyper_cgi,cervical_hyper_hypo_cgi,cervical_hypo_stable_cgi,cervical_hyper_stable_cgi,cervical_hyper_back_to_normal_cgi,cervical_hypo_back_to_normal_cgi,file='cervical_cgi_rate.Rdata')
save(cervical_background_feature,cervical_normal_cancer_feature,cervical_hypo_hypo_feature,cervical_hypo_hyper_feature,cervical_normal_pre_feature,cervical_hyper_hyper_feature,cervical_hyper_hypo_feature,cervical_hypo_stable_feature,cervical_hyper_stable_feature,cervical_hyper_back_to_normal_feature,cervical_hypo_back_to_normal_feature,file='cervical_feature_rate.Rdata')

library(reshape2)

cervical_background_cgi$variable='background'
cervical_normal_pre_cgi$variable='normal_v_pre'
cervical_normal_cancer_cgi$variable='pre_v_cancer'
cervical_hypo_hypo_cgi$variable='hypo_hypo'
cervical_hyper_hyper_cgi$variable='hyper_hyper'
cervical_hypo_hyper_cgi$variable='hypo_hyper'
cervical_hyper_hypo_cgi$variable='hyper_hypo'
cervical_hypo_stable_cgi$variable='hypo_stable'
cervical_hyper_stable_cgi$variable='hyper_stable'
cervical_hypo_back_to_normal_cgi$variable='hypo_RtoN'
cervical_hyper_back_to_normal_cgi$variable='hyper_RtoN'
cervical_barplotdata=rbind(cervical_background_cgi,cervical_normal_pre_cgi,cervical_normal_cancer_cgi,cervical_hypo_hypo_cgi,cervical_hyper_hyper_cgi,cervical_hypo_hyper_cgi,cervical_hyper_hypo_cgi,cervical_hypo_stable_cgi,cervical_hyper_stable_cgi,cervical_hypo_back_to_normal_cgi,cervical_hyper_back_to_normal_cgi)
cervical_cgi_barplotdata=dcast(cervical_barplotdata,Var1~variable,value.var = 'value')
rownames(cervical_cgi_barplotdata)=cervical_cgi_barplotdata[,1]
cervical_cgi_barplotdata=as.matrix(cervical_cgi_barplotdata[,-1])
pdf('Cervical_cgi.pdf', width=8.51, height=5.49)
barplot(cervical_cgi_barplotdata,main = "Cervical",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(cervical_cgi_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7)
dev.off()


cervical_background_feature$variable='background'
cervical_normal_pre_feature$variable='normal_v_pre'
cervical_normal_cancer_feature$variable='pre_v_cancer'
cervical_hypo_hypo_feature$variable='hypo_hypo'
cervical_hyper_hyper_feature$variable='hyper_hyper'
cervical_hypo_hyper_feature$variable='hypo_hyper'
cervical_hyper_hypo_feature$variable='hyper_hypo'
cervical_hypo_stable_feature$variable='hypo_stable'
cervical_hyper_stable_feature$variable='hyper_stable'
cervical_hypo_back_to_normal_feature$variable='hypo_RtoN'
cervical_hyper_back_to_normal_feature$variable='hyper_RtoN'
cervical_barplotdata=rbind(cervical_background_feature,cervical_normal_pre_feature,cervical_normal_cancer_feature,cervical_hypo_hypo_feature,cervical_hyper_hyper_feature,cervical_hypo_hyper_feature,cervical_hyper_hypo_feature,cervical_hypo_stable_feature,cervical_hyper_stable_feature,cervical_hypo_back_to_normal_feature,cervical_hyper_back_to_normal_feature)
cervical_feature_barplotdata=dcast(cervical_barplotdata,Var1~variable,value.var = 'value')
rownames(cervical_feature_barplotdata)=cervical_feature_barplotdata[,1]
cervical_feature_barplotdata=as.matrix(cervical_feature_barplotdata[,-1])
pdf('Cervical_feature.pdf', width=8.51, height=5.49)
barplot(cervical_feature_barplotdata,main = "Cervical",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#BAB2BA','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(cervical_feature_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7, args.legend = list(y=0.5,cex=0.8))
#locator(1) 
dev.off()

