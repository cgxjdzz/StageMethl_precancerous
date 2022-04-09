cancerlist=c('prostatic','cervical','liver','gastric','colonic')
sitelist=c('hypo_hypo','hyper_hyper','hypo_hyper','hyper_hypo','hypo_stable','hyper_stable','hypo_back_to_normal','hyper_back_to_normal')
library(readr)
library(survival)
library(getopt)
library(readr)
library(survival)
library(pROC)


do_coxph <- function(df, name) {
  if (name=='status'|name=='time') {
    return()
  }
  # 
  # if (sum(is.na(df[[name]])) / nrow(df) > opt$valid) {
  #   return()
  # }
  # 
  #----------------------------#
  # coxph
  #----------------------------#
  library(survival)
  formula_string = paste('Surv(time,status) ~ ', name, sep = '')
  rex.cox <- coxph(as.formula(formula_string), df)
  summary <- summary(rex.cox)
  coef <- summary$coefficients[1]
  p_value <- summary$coefficients[5]
  #k-m
  predictor <- coef * df[[name]]
  score_median <- median(predictor, na.rm = TRUE)
  group <- predictor < score_median
  diff <- survdiff(Surv(df$time, df$status) ~ group, rho = 0)
  kmp <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = FALSE)
  ####ROC
  df_tmp <- data.frame(outcome, predictor)
  df_roc <-
    df_tmp[which(!is.na(df_tmp$outcome) &
                   !is.na(df_tmp$predictor)), ]
  rocauc <-
    pROC::roc(df_roc$outcome,
              df_roc$predictor,
              levels = c(0, 1),
              na.rm = TRUE)$auc
  # sum=rbind(sum,data.frame(        name,
  #                 coef,
  #                 p_value,
  #                 rocauc,
  #                 score_median,
  #                 kmp))
  sprintf(
    "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
    name,
    coef,
    p_value,
    rocauc,
    score_median,
    kmp)
  
}







for (i in cancerlist){
  memory.limit(16000)
  clinic=data.frame(read_tsv(paste0(i,'/clinical.tsv')))
  clinic_dead=clinic[which(clinic$vital_status=='Dead'),]
  clinic_alive=clinic[-which(clinic$vital_status=='Dead'),]
  survival_dead=data.frame(ID=clinic_dead$case_submitter_id,status=clinic_dead$vital_status,time=as.numeric(clinic_dead$days_to_death))
  survival_alive=data.frame(ID=clinic_alive$case_submitter_id,status=clinic_alive$vital_status,time=as.numeric(clinic_alive$days_to_last_follow_up))
  survival_data=unique(rbind(survival_dead,survival_alive))
  #1是生0是死
  survival_data$status[-which(survival_data$status=='Dead')]=1
  survival_data$status[which(survival_data$status=='Dead')]=0
  write.csv(survival_data,file=paste0(i,'/0402survival.csv'))
  beta_data=data.frame(read_tsv(paste0(i,'/beta.tsv')))
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
  cox_data=data.frame(t(beta_data),row.names = colnames(beta_data))
  ###是不是有可能转置时出了问题
  colnames(cox_data)=rownames(beta_data)
  cox_data=cox_data[,c(300001:dim(cox_data)[2])]
  cox_data$time=survival_data$time
  cox_data$status=as.numeric(survival_data$status)
  outcome <-
    ifelse(survival_data$time < 365 * 3, ifelse(survival_data$status == 1, 0, NA), 1)
  ########################
  # cor_abs=abs(cor(cox_data,survival_data$time))
  # name=rownames(cor_abs)[which(cor_abs>0.1)]
  # cox_data=cox_data[,name]
  ###########################
  
  library("foreach")
  library("doParallel")
  cl<-makeCluster(10)
  registerDoParallel(cl)
  ee=foreach(name = colnames(cox_data)) %dopar% {
    (do_coxph(cox_data, name))
  }
  stopCluster(cl)
  write.table(unlist(ee),paste0(i,'/',i,'_400000_coef.txt'),quote =F,row.names = F,
              col.names = F,eol = "")
  
  
  

  # summary=summary(cox_result)
  # # save(summary)
  # coef <- summary$coefficients[1]
  # p_value <- summary$coefficients[5]
  }