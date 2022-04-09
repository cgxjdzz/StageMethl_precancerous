load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hyper_and_hypo.Rdata')
library(ChAMP)
library(ggplot2)
library(ggthemes)
cancerlist=c('cervical','colonic','liver','prostatic','gastric','melanoma')
plot_data=data.frame(cervical=c(length(cervical_hyper_hyper),length(cervical_hyper_hypo),length(cervical_hyper_back_to_normal),length(cervical_hyper_stable),length(cervical_hypo_hypo),length(cervical_hypo_hyper),length(cervical_hypo_back_to_normal),length(cervical_hypo_stable)),
                     colonic=c(length(colonic_hyper_hyper),length(colonic_hyper_hypo),length(colonic_hyper_back_to_normal),length(colonic_hyper_stable),length(colonic_hypo_hypo),length(colonic_hypo_hyper),length(colonic_hypo_back_to_normal),length(colonic_hypo_stable)),
                     gastric=c(length(gastric_hyper_hyper),length(gastric_hyper_hypo),length(gastric_hyper_back_to_normal),length(gastric_hyper_stable),length(gastric_hypo_hypo),length(gastric_hypo_hyper),length(gastric_hypo_back_to_normal),length(gastric_hypo_stable)),
                     prostatic=c(length(prostatic_hyper_hyper),length(prostatic_hyper_hypo),length(prostatic_hyper_back_to_normal),length(prostatic_hyper_stable),length(prostatic_hypo_hypo),length(prostatic_hypo_hyper),length(prostatic_hypo_back_to_normal),length(prostatic_hypo_stable)),
                     melanoma=c(length(melanoma_hyper_hyper),length(melanoma_hyper_hypo),length(melanoma_hyper_back_to_normal),length(melanoma_hyper_stable),length(melanoma_hypo_hypo),length(melanoma_hypo_hyper),length(melanoma_hypo_back_to_normal),length(melanoma_hypo_stable)),
                     liver=c(length(liver_hyper_hyper),length(liver_hyper_hypo),length(liver_hyper_back_to_normal),length(liver_hyper_stable),length(liver_hypo_hypo),length(liver_hypo_hyper),length(liver_hypo_back_to_normal),length(liver_hypo_stable)))
rownames(plot_data)=c('hyper_hyper','hyper_hypo','hyper_RtoN','hyper_stable','hypo_hypo','hypo_hyper','hypo_RtoN','hypo_stable')

pdf('site_count.pdf', width=8.51, height=5.49)
barplot(as.matrix(plot_data),main = "Counts",
        xlab = "", ylab = "Counts",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#BAB2BA','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(plot_data),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 1.2, args.legend = list(x=55,y=35000,cex=1.2))
#locator(1) 
dev.off()

plot_data_2=data.frame(cervical=plot_data$cervical/sum(plot_data$cervical),
                       colonic=plot_data$colonic/sum(plot_data$colonic),
                       gastric=plot_data$gastric/sum(plot_data$gastric),
                       prostatic=plot_data$prostatic/sum(plot_data$prostatic),
                       melanoma=plot_data$melanoma/sum(plot_data$melanoma),
                       liver=plot_data$liver/sum(plot_data$liver)
                       )
rownames(plot_data_2)=rownames(plot_data)
pdf('site_rate.pdf', width=8.51, height=5.49)
barplot(as.matrix(plot_data_2),main = "Rate",
        xlab = "", ylab = "Rate",
        col = c('#4D3D51','#CA4528','#EFBA51','#4D79A6','#BAB2BA','#1699EF','#9EC6BE','#AA9E92'),
        legend.text = rownames(plot_data),beside=TRUE,las=2,cex.axis = 0.8,cex.names = 1.2, args.legend = list(x=55,y=0.6,cex=1.2))
#locator(1) 
dev.off()
