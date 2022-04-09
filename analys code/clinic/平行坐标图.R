library(GGally)
load('D:/甲基化/癌症发展/宫颈癌/0128sum_combat.Rdata')
cervical_beta=sum_combat
cervical_label=pd$label
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
prostatic_beta=sum_norm
prostatic_label=pd$label
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
gastric_beta=sum_norm
gastric_label=pd$label
load('D:/甲基化/癌症发展/肝癌/0227_sum_norm.Rdata')
liver_beta=sum_norm
liver_label=pd$label
load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')
colonic_beta=sum_norm
colonic_label=pd$`sum_impute$pd`
cancerlist=c('prostatic','cervical','liver','gastric','colonic')
countlist=c('100000','200000','300000','400000')
sitelist=c('hypo_hypo','hyper_hyper','hypo_hyper','hyper_hypo','hypo_stable','hyper_stable','hypo_back_to_normal','hyper_back_to_normal')
library(reshape2)

for (i in cancerlist){
  temp_beta=t(get(paste0(i,'_beta')))
  temp_label=get(paste0(i,'_label'))
  ee=aggregate(temp_beta,list(temp_label),mean)
  rownames(ee)=ee[,1]
  ee=ee[,-1]
  nn=data.frame(t(ee))
  colnames(nn)=rownames(ee)
  nn$site='ee'
  for (j in sitelist){
    valuesite=get(paste0(i,'_',j))
    if(rownames(nn) %in% valuesite){
      nn[which(rownames(nn) %in% valuesite),]$site=j
    }
  
  }
  nn=nn[-which(nn$site==ee),]
}