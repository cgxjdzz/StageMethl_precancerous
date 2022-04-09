library(ChAMP)
library(plyr)
load('D:/甲基化/癌症发展/肝癌/sum_DMP.Rdata')
liver_DMP=sum_DMP$Cirrhosis_to_Normal
load('D:/甲基化/癌症发展/0331改数据/黑色素瘤/0331_sum_DMP.Rdata')
melanoma_DMP=sum_DMP$Precancer_to_Normal
load('D:/甲基化/癌症发展/前列腺癌/1217_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal", "Primary Tumor"),c("Normal", "Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
prostatic_DMP=sum_DMP$CAF_to_Normal
load('D:/甲基化/癌症发展/胃癌/0107_sum_norm.Rdata')
label=mapvalues(pd$label, c("Solid Tissue Normal","intestinal metaplasia biopsy from gastric antrum" ,"Primary Tumor"),c("Normal",'Precancer',"Tumor"))
sum_DMP=champ.DMP(sum_norm,label)
gastric_DMP=sum_DMP$Normal_to_Precancer
load("D:/甲基化/癌症发展/宫颈癌/0128combat_sum_DMP.Rdata")
cervical_DMP=sum_DMP$normal_to_cin3
load('D:/甲基化/癌症发展/0331改数据/colonic/139404+68060/0331_sum_norm.Rdata')
myNorm=sum_norm
pd=pd$`sum_impute$pd`
pd=mapvalues(pd, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = myNorm,pheno = pd)
colonic_DMP=sum_DMP$Normal_to_Precancer


######################################################



load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/prostatic_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/melanoma_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/cervical_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/liver_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/gastric_hyper_and_hypo.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_back_to_normal.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hy_and_stable.Rdata')
load('D:/甲基化/癌症发展/0317更改位点定义/dmp的热图/colonic_hyper_and_hypo.Rdata')
cancerlist=c('prostatic','cervical','liver','gastric','colonic')
sitelist=c('hypo_hypo','hyper_hyper','hypo_hyper','hyper_hypo','hypo_stable','hyper_stable','hypo_back_to_normal','hyper_back_to_normal')
for( i in sitelist){
  dir.create(i)
  for(j in cancerlist){
    name=paste0(j,'_',i)
    write.table(get(name),file = paste0(i,'/',name,'.csv'),row.names = F,col.names = F,sep=' ',eol = "\n")
  }
}


for( i in sitelist){
  dir.create(paste0('union_site_',i))
  eval(parse(text=paste0(i,'=c()')))
  for(j in cancerlist){
    name=paste0(j,'_',i)
    if (length(get(name))>1000){
      eval(parse(text=paste0(name,'=',name,'[order(',j,'_DMP[',j,'_',i,',]$adj.P.Val)]')))
      eval(parse(text=paste0(i,'=c(',i,',get(name)[c(1:1000)])')))
    }else{
      eval(parse(text=paste0(i,'=c(',i,',get(name))')))
    }
      
    

  }
  eval(parse(text=paste0(i,'=unique(',i,')')))
  data=data.frame(row.names =get( i),site=get(i))
  for (j in cancerlist){
    eval(parse(text=paste0('data$',j,'=0')))
    name=paste0(j,'_',i)
    eval(parse(text=paste0('data$',j,'[which(rownames(data) %in% get(name))]=1')))
  }
  data=data[,-1]
  data$sum=rowSums(data)
  table(data$sum)
  data=data[order(-data$sum),]
  print(i)
  table(data$sum)
  print(max(data$sum))
  write.table(data,file = paste0('union_site/union_site_',i,'/','union_site_',i,'.csv'),row.names = T,col.names = T,sep=' ',eol = "\n")
}
