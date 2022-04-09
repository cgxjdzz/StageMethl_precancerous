library(qqman)           #加载qqman包
library(RColorBrewer)   
library(ChAMP)
library(ggplot2)
library(dplyr)
library(CMplot)
###################
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
cervical_DMP=sum_DMP$cancer_to_cin3
load('D:/甲基化/癌症发展/0331改数据//colonic/139404+68060/0331_sum_norm.Rdata')
library(ChAMP)
library(plyr)
pd0204=data.frame(label=pd)
pd0204$sum_impute.pd=mapvalues(pd0204$sum_impute.pd, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = sum_norm,pheno = pd0204$sum_impute.pd)
colonic_DMP=sum_DMP$Tumor_to_Precancer
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum" ,"Primary Tumor"),c("Normal",'Precancer',"Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
gastric_DMP=sum_DMP$Precancer_to_Tumor
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
prostatic_DMP=sum_DMP$CAF_to_Tumor
load('D:/甲基化/癌症发展/0331改数据/黑色素瘤/0331_sum_DMP.Rdata')
library(ChAMP)
melanoma_DMP=sum_DMP$Tumor_to_Precancer
load('D:/甲基化/癌症发展/肝癌/sum_DMP.Rdata')
liver_DMP=sum_DMP$Cirrhosis_to_Tumor
########################
data("probe.features.epic")
epic=probe.features
data("probe.features")
probe=rbind(epic[,c("gene",'CHR',"MAPINFO")],probe.features[,c("gene",'CHR',"MAPINFO")])
probe=unique(probe)
cancerlist=c('cervical','colonic','gastric','melanoma','liver','prostatic')
for (i in cancerlist){
  load(paste0('0319',i,'_hyper_hypo.Rdata'))
  i_hyper=get(paste0(i,'_pre_cancer_hyper'))
  probe_hyper=probe[i_hyper,]
  p_hyper=get(paste0(i,'_DMP'))[i_hyper,'adj.P.Val']
  probe_hyper$p=p_hyper+1e-320
  probe_hyper$CHR=as.numeric(probe_hyper$CHR)
  colnames(probe_hyper)=c('SNP','Chromosome','Position','trait')
  pdf(paste0(i,'_hyper_manhattan.pdf'),height = 3.3)
  CMplot(probe_hyper,plot.type = "m",
         threshold = c(0.01,0.05)/nrow(probe_hypo),
         threshold.col=c('grey','black'),
         threshold.lty = c(1,2),threshold.lwd = c(1,1), amplify = T,
         signal.cex = c(0.5,0.5), signal.pch = c(20,20),signal.col = c("red","orange"),file.output = F,cex.axis=0.3,cex=0.1)
  dev.off()
  i_hypo=get(paste0(i,'_pre_cancer_hypo'))
  probe_hypo=probe[i_hypo,]
  p_hypo=get(paste0(i,'_DMP'))[i_hypo,'adj.P.Val']
  probe_hypo$p=p_hypo+1e-320
  probe_hypo$CHR=as.numeric(probe_hypo$CHR)
  colnames(probe_hypo)=c('SNP','Chromosome','Position','trait')
  pdf(paste0(i,'_hypo_manhattan.pdf'),height = 3.3)
  CMplot(probe_hypo,plot.type = "m",
         threshold = c(0.01,0.05)/nrow(probe_hypo),
         threshold.col=c('grey','black'),
         threshold.lty = c(1,2),threshold.lwd = c(1,1), amplify = T,
         signal.cex =  c(0.5,0.5), signal.pch = c(20,20),signal.col = c("red","orange"),file.output =F ,cex.axis=0.3,cex=0.3)
  dev.off()
}
# 
# axisdf = probe_hyper %>% group_by(CHR) %>% summarize(center=( max(MAPINFO) + min(MAPINFO) ) / 2 )
# ggplot(probe_hyper, aes(x=MAPINFO, y=-log10(p))) +
#   
#   geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
#   scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
#   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
#   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
#   theme_bw() +
#   theme( 
#     legend.position="none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
