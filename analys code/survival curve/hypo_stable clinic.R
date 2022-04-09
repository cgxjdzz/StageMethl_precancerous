library(readr)
library(survival)
library(survminer)
library(ggplot2)
cancerlist=c('prostatic','cervical','liver','gastric','colonic')

for(i in cancerlist){
  ee=read.table(paste0('hypo_stable/',i,'_coef.txt'),fill=T)
  cg_data=ee[which(ee[,3]<0.05),]
  cgpoint=ee[which(ee[,3]<0.05),1]
  assign(paste0(i,'_cgpoint'),cgpoint)
  
  
  memory.limit(16000)
  clinic=data.frame(read_tsv(paste0('D:/甲基化/癌症发展/0317更改位点定义/0402 clinic/',i,'/clinical.tsv')))
  clinic_dead=clinic[which(clinic$vital_status=='Dead'),]
  clinic_alive=clinic[-which(clinic$vital_status=='Dead'),]
  survival_dead=data.frame(ID=clinic_dead$case_submitter_id,status=clinic_dead$vital_status,time=as.numeric(clinic_dead$days_to_death))
  survival_alive=data.frame(ID=clinic_alive$case_submitter_id,status=clinic_alive$vital_status,time=as.numeric(clinic_alive$days_to_last_follow_up))
  survival_data=unique(rbind(survival_dead,survival_alive))
  #1是生0是死
  survival_data$status[-which(survival_data$status=='Dead')]=1
  survival_data$status[which(survival_data$status=='Dead')]=0
  dir.create(i)
  write.csv(survival_data,file = paste0(i,'/survival_data.csv'))
  beta_data=data.frame(read_tsv(paste0('D:/甲基化/癌症发展/0317更改位点定义/0402 clinic/',i,'/beta.tsv')))
  rownames(beta_data)=beta_data[,1]
  beta_data=beta_data[,-1]
  beta_data=beta_data[-which(is.na(beta_data[,1])),]
  beta_data=na.omit(beta_data)
  id=gsub('.','-',substr(colnames(beta_data),1,12), fixed = TRUE)
  colnames(beta_data)=id
  rownames(survival_data)=survival_data$ID
  id=intersect(id ,survival_data$ID)
  beta_data=beta_data[,id]
  survival_data=survival_data[id,]
  beta_data=beta_data[cgpoint,]
  
  different_data=data.frame(matrix(ncol = dim(beta_data)[2],nrow = 0))
  predictor <- numeric(dim(beta_data)[2])
  
  for (k in 1:length(cgpoint)) {
    predictor <- predictor + cg_data[k,2] * as.numeric(beta_data[k,])
  }
  score_median <- median(predictor, na.rm = TRUE)

  label=rep(0,dim(survival_data)[1])
  label[which(predictor>score_median)]=1
  label[which(predictor<score_median)]=-1
  data=data.frame(survival=as.numeric(survival_data$time),status=as.numeric(survival_data$status),label=label)
  attach(data)
  fit <- survfit(Surv(survival,status)~label,data=data)
  ggsurvplot(fit, data = data,pval = TRUE, conf.int = TRUE,palette = "hue")
  ggsave(paste0('predict_divide_hypo_stable/',i,'.pdf'))
}

