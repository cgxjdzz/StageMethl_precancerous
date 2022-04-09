library(ChAMP)
data(probe.features.epic)#850K
#data(probe.features)#450K

load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hyper_and_hypo.Rdata')
load("D:/甲基化/癌症发展/胃癌/sum_dmp.Rdata")
hypo_hyper=gastric_hypo_hyper
hyper_hypo=gastric_hyper_hypo
hyper_hyper=gastric_hyper_hyper
hypo_hypo=gastric_hypo_hypo
hypo_stable=gastric_hypo_stable
hyper_stable=gastric_hyper_stable
hyper_back_to_normal=gastric_hyper_back_to_normal
hypo_back_to_normal=gastric_hypo_back_to_normal

myNorm=sum_norm
normal_pre=rownames(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal)
normal_cancer=rownames(sum_DMP$Primary_Tumor_to_Solid_Tissue_Normal)

myNorm=sum_norm
background=rownames(myNorm)

feature=data.frame(feature=probe.features$feature,cgi=probe.features$cgi)
rownames(feature)=rownames(probe.features)

normal_pre=feature[normal_pre,]
normal_pre=na.omit(data.frame(feature=normal_pre$feature,cgi=normal_pre$cgi))
gastric_normal_pre_feature=table(normal_pre$feature)/dim(normal_pre)[1]
gastric_normal_pre_cgi=table(normal_pre$cgi)/dim(normal_pre)[1]
gastric_normal_pre_feature=melt(as.data.frame(gastric_normal_pre_feature))
gastric_normal_pre_feature[,2]=gsub('Freq','gastric',gastric_normal_pre_feature[,2])
gastric_normal_pre_cgi=melt(as.data.frame(gastric_normal_pre_cgi))
gastric_normal_pre_cgi[,2]=gsub('Freq','gastric',gastric_normal_pre_cgi[,2])

normal_cancer=feature[normal_cancer,]
normal_cancer=na.omit(data.frame(feature=normal_cancer$feature,cgi=normal_cancer$cgi))
gastric_normal_cancer_feature=table(normal_cancer$feature)/dim(normal_cancer)[1]
gastric_normal_cancer_cgi=table(normal_cancer$cgi)/dim(normal_cancer)[1]
gastric_normal_cancer_feature=melt(as.data.frame(gastric_normal_cancer_feature))
gastric_normal_cancer_feature[,2]=gsub('Freq','gastric',gastric_normal_cancer_feature[,2])
gastric_normal_cancer_cgi=melt(as.data.frame(gastric_normal_cancer_cgi))
gastric_normal_cancer_cgi[,2]=gsub('Freq','gastric',gastric_normal_cancer_cgi[,2])


hypo_hyper=feature[hypo_hyper,]
hypo_hyper=na.omit(data.frame(feature=hypo_hyper$feature,cgi=hypo_hyper$cgi))
gastric_hypo_hyper_feature=table(hypo_hyper$feature)/dim(hypo_hyper)[1]
gastric_hypo_hyper_cgi=table(hypo_hyper$cgi)/dim(hypo_hyper)[1]
gastric_hypo_hyper_feature=melt(as.data.frame(gastric_hypo_hyper_feature))
gastric_hypo_hyper_feature[,2]=gsub('Freq','gastric',gastric_hypo_hyper_feature[,2])
gastric_hypo_hyper_cgi=melt(as.data.frame(gastric_hypo_hyper_cgi))
gastric_hypo_hyper_cgi[,2]=gsub('Freq','gastric',gastric_hypo_hyper_cgi[,2])

hyper_hypo=feature[hyper_hypo,]
hyper_hypo=na.omit(data.frame(feature=hyper_hypo$feature,cgi=hyper_hypo$cgi))
gastric_hyper_hypo_feature=table(hyper_hypo$feature)/dim(hyper_hypo)[1]
gastric_hyper_hypo_cgi=table(hyper_hypo$cgi)/dim(hyper_hypo)[1]
gastric_hyper_hypo_feature=melt(as.data.frame(gastric_hyper_hypo_feature))
gastric_hyper_hypo_feature[,2]=gsub('Freq','gastric',gastric_hyper_hypo_feature[,2])
gastric_hyper_hypo_cgi=melt(as.data.frame(gastric_hyper_hypo_cgi))
gastric_hyper_hypo_cgi[,2]=gsub('Freq','gastric',gastric_hyper_hypo_cgi[,2])

hyper_hyper=feature[hyper_hyper,]
hyper_hyper=na.omit(data.frame(feature=hyper_hyper$feature,cgi=hyper_hyper$cgi))
gastric_hyper_hyper_feature=table(hyper_hyper$feature)/dim(hyper_hyper)[1]
gastric_hyper_hyper_cgi=table(hyper_hyper$cgi)/dim(hyper_hyper)[1]
gastric_hyper_hyper_feature=melt(as.data.frame(gastric_hyper_hyper_feature))
gastric_hyper_hyper_feature[,2]=gsub('Freq','gastric',gastric_hyper_hyper_feature[,2])
gastric_hyper_hyper_cgi=melt(as.data.frame(gastric_hyper_hyper_cgi))
gastric_hyper_hyper_cgi[,2]=gsub('Freq','gastric',gastric_hyper_hyper_cgi[,2])

hypo_hypo=feature[hypo_hypo,]
hypo_hypo=na.omit(data.frame(feature=hypo_hypo$feature,cgi=hypo_hypo$cgi))
gastric_hypo_hypo_feature=table(hypo_hypo$feature)/dim(hypo_hypo)[1]
gastric_hypo_hypo_cgi=table(hypo_hypo$cgi)/dim(hypo_hypo)[1]
gastric_hypo_hypo_feature=melt(as.data.frame(gastric_hypo_hypo_feature))
gastric_hypo_hypo_feature[,2]=gsub('Freq','gastric',gastric_hypo_hypo_feature[,2])
gastric_hypo_hypo_cgi=melt(as.data.frame(gastric_hypo_hypo_cgi))
gastric_hypo_hypo_cgi[,2]=gsub('Freq','gastric',gastric_hypo_hypo_cgi[,2])


background=feature[background,]
background=na.omit(data.frame(feature=background$feature,cgi=background$cgi))
gastric_background_feature=table(background$feature)/dim(background)[1]
gastric_background_cgi=table(background$cgi)/dim(background)[1]
gastric_background_feature=melt(as.data.frame(gastric_background_feature))
gastric_background_feature[,2]=gsub('Freq','gastric',gastric_background_feature[,2])
gastric_background_cgi=melt(as.data.frame(gastric_background_cgi))
gastric_background_cgi[,2]=gsub('Freq','gastric',gastric_background_cgi[,2])



hyper_stable=feature[hyper_stable,]
hyper_stable=na.omit(data.frame(feature=hyper_stable$feature,cgi=hyper_stable$cgi))
gastric_hyper_stable_feature=table(hyper_stable$feature)/dim(hyper_stable)[1]
gastric_hyper_stable_cgi=table(hyper_stable$cgi)/dim(hyper_stable)[1]
gastric_hyper_stable_feature=melt(as.data.frame(gastric_hyper_stable_feature))
gastric_hyper_stable_feature[,2]=gsub('Freq','gastric',gastric_hyper_stable_feature[,2])
gastric_hyper_stable_cgi=melt(as.data.frame(gastric_hyper_stable_cgi))
gastric_hyper_stable_cgi[,2]=gsub('Freq','gastric',gastric_hyper_stable_cgi[,2])

hypo_stable=feature[hypo_stable,]
hypo_stable=na.omit(data.frame(feature=hypo_stable$feature,cgi=hypo_stable$cgi))
gastric_hypo_stable_feature=table(hypo_stable$feature)/dim(hypo_stable)[1]
gastric_hypo_stable_cgi=table(hypo_stable$cgi)/dim(hypo_stable)[1]
gastric_hypo_stable_feature=melt(as.data.frame(gastric_hypo_stable_feature))
gastric_hypo_stable_feature[,2]=gsub('Freq','gastric',gastric_hypo_stable_feature[,2])
gastric_hypo_stable_cgi=melt(as.data.frame(gastric_hypo_stable_cgi))
gastric_hypo_stable_cgi[,2]=gsub('Freq','gastric',gastric_hypo_stable_cgi[,2])

hypo_back_to_normal=feature[hypo_back_to_normal,]
hypo_back_to_normal=na.omit(data.frame(feature=hypo_back_to_normal$feature,cgi=hypo_back_to_normal$cgi))
gastric_hypo_back_to_normal_feature=table(hypo_back_to_normal$feature)/dim(hypo_back_to_normal)[1]
gastric_hypo_back_to_normal_cgi=table(hypo_back_to_normal$cgi)/dim(hypo_back_to_normal)[1]
gastric_hypo_back_to_normal_feature=melt(as.data.frame(gastric_hypo_back_to_normal_feature))
gastric_hypo_back_to_normal_feature[,2]=gsub('Freq','gastric',gastric_hypo_back_to_normal_feature[,2])
gastric_hypo_back_to_normal_cgi=melt(as.data.frame(gastric_hypo_back_to_normal_cgi))
gastric_hypo_back_to_normal_cgi[,2]=gsub('Freq','gastric',gastric_hypo_back_to_normal_cgi[,2])

hyper_back_to_normal=feature[hyper_back_to_normal,]
hyper_back_to_normal=na.omit(data.frame(feature=hyper_back_to_normal$feature,cgi=hyper_back_to_normal$cgi))
gastric_hyper_back_to_normal_feature=table(hyper_back_to_normal$feature)/dim(hyper_back_to_normal)[1]
gastric_hyper_back_to_normal_cgi=table(hyper_back_to_normal$cgi)/dim(hyper_back_to_normal)[1]
gastric_hyper_back_to_normal_feature=melt(as.data.frame(gastric_hyper_back_to_normal_feature))
gastric_hyper_back_to_normal_feature[,2]=gsub('Freq','gastric',gastric_hyper_back_to_normal_feature[,2])
gastric_hyper_back_to_normal_cgi=melt(as.data.frame(gastric_hyper_back_to_normal_cgi))
gastric_hyper_back_to_normal_cgi[,2]=gsub('Freq','gastric',gastric_hyper_back_to_normal_cgi[,2])




save(gastric_background_cgi,gastric_normal_cancer_cgi,gastric_hypo_hypo_cgi,gastric_hypo_hyper_cgi,gastric_normal_pre_cgi,gastric_hyper_hyper_cgi,gastric_hyper_hypo_cgi,gastric_hypo_stable_cgi,gastric_hyper_stable_cgi,gastric_hyper_back_to_normal_cgi,gastric_hypo_back_to_normal_cgi,file='gastric_cgi_rate.Rdata')
save(gastric_background_feature,gastric_normal_cancer_feature,gastric_hypo_hypo_feature,gastric_hypo_hyper_feature,gastric_normal_pre_feature,gastric_hyper_hyper_feature,gastric_hyper_hypo_feature,gastric_hypo_stable_feature,gastric_hyper_stable_feature,gastric_hyper_back_to_normal_feature,gastric_hypo_back_to_normal_feature,file='gastric_feature_rate.Rdata')

gastric_background_cgi$variable='background'
gastric_normal_pre_cgi$variable='normal_v_pre'
gastric_normal_cancer_cgi$variable='pre_v_cancer'
gastric_hypo_hypo_cgi$variable='hypo_hypo'
gastric_hyper_hyper_cgi$variable='hyper_hyper'
gastric_hypo_hyper_cgi$variable='hypo_hyper'
gastric_hyper_hypo_cgi$variable='hyper_hypo'
gastric_hypo_stable_cgi$variable='hypo_stable'
gastric_hyper_stable_cgi$variable='hyper_stable'
gastric_hypo_back_to_normal_cgi$variable='hypo_RtoN'
gastric_hyper_back_to_normal_cgi$variable='hyper_RtoN'
gastric_barplotdata=rbind(gastric_background_cgi,gastric_normal_pre_cgi,gastric_normal_cancer_cgi,gastric_hypo_hypo_cgi,gastric_hyper_hyper_cgi,gastric_hypo_hyper_cgi,gastric_hyper_hypo_cgi,gastric_hypo_stable_cgi,gastric_hyper_stable_cgi,gastric_hypo_back_to_normal_cgi,gastric_hyper_back_to_normal_cgi)
gastric_cgi_barplotdata=dcast(gastric_barplotdata,Var1~variable,value.var = 'value')
rownames(gastric_cgi_barplotdata)=gastric_cgi_barplotdata[,1]
gastric_cgi_barplotdata=as.matrix(gastric_cgi_barplotdata[,-1])
pdf('gastric_cgi.pdf', width=8.51, height=5.49)
barplot(gastric_cgi_barplotdata,main = "gastric",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#1699EF'),
        legend.text = rownames(gastric_cgi_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7)
dev.off()


gastric_background_feature$variable='background'
gastric_normal_pre_feature$variable='normal_v_pre'
gastric_normal_cancer_feature$variable='pre_v_cancer'
gastric_hypo_hypo_feature$variable='hypo_hypo'
gastric_hyper_hyper_feature$variable='hyper_hyper'
gastric_hypo_hyper_feature$variable='hypo_hyper'
gastric_hyper_hypo_feature$variable='hyper_hypo'
gastric_hypo_stable_feature$variable='hypo_stable'
gastric_hyper_stable_feature$variable='hyper_stable'
gastric_hypo_back_to_normal_feature$variable='hypo_RtoN'
gastric_hyper_back_to_normal_feature$variable='hyper_RtoN'
gastric_barplotdata=rbind(gastric_background_feature,gastric_normal_pre_feature,gastric_normal_cancer_feature,gastric_hypo_hypo_feature,gastric_hyper_hyper_feature,gastric_hypo_hyper_feature,gastric_hyper_hypo_feature,gastric_hypo_stable_feature,gastric_hyper_stable_feature,gastric_hypo_back_to_normal_feature,gastric_hyper_back_to_normal_feature)
gastric_feature_barplotdata=dcast(gastric_barplotdata,Var1~variable,value.var = 'value')
rownames(gastric_feature_barplotdata)=gastric_feature_barplotdata[,1]
gastric_feature_barplotdata=as.matrix(gastric_feature_barplotdata[,-1])
pdf('gastric_feature.pdf', width=8.51, height=5.49)
barplot(gastric_feature_barplotdata,main = "gastric",
        xlab = "", ylab = "Frequency",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#BAB2BA','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(gastric_feature_barplotdata),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 0.7, args.legend = list(x=92,y=0.52,cex=0.8))
#locator(1) 
dev.off()
  