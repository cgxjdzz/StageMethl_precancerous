suppressPackageStartupMessages({
  library(getopt)
  library(readr)
  library(survival)
  library(pROC)
})
outcome <-
  ifelse(survival_data$time < 365 * 3, ifelse(df$status == 1, 0, NA), 1)

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


library("foreach")
library("doParallel")
cl<-makeCluster(11)
registerDoParallel(cl)
ee=foreach(name = colnames(cox_data_test)) %dopar% {
  (do_coxph(cox_data_test, name))
}
stopCluster(cl)
write.table(unlist(ee),'ee.txt',quote =F,row.names = F,
            col.names = F,eol = "")
