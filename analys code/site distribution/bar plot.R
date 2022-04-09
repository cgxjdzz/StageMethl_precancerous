cervical=read.csv('cervical_model_score.csv')
colonic=read.csv('colonic_model_score.csv')
gastric=read.csv('gastric_model_score.csv')
melanoma=read.csv('melanoma_model_score.csv')
prostatic=read.csv('prostatic_model_score.csv')

data=rbind(cervical,colonic,gastric,melanoma,prostatic)
split=strsplit(data$X,'_')
split=as.data.frame(split)
data$cancer=c(rep('cervical',8),rep('colonic',8),rep('gastric',8),rep('melanoma',8),rep('prostatic',8))
data$X=c(rep(c('hypo_return_to_normal','hyper_return_to_normal','hyper_stable','hypo_stable','hyper_hyper','hypo_hyper','hyper_hypo','hypo_hypo'),5))
data$X0=-data$X0
library(reshape2)
data2=dcast(data,X~cancer,value.var = 'X0')
data2=as.matrix(data2)
rownames(data2)=data2[,1]
data2=data2[,-1]

am <- factor(am)


# Change factor levels

# other_table <- xtabs(~cyl + am , data = mtcars) # Equivalent

data2_num <- apply(data2,1, as.numeric)
rownames(data2_num)=colnames(data2)
data2_num=t(data2_num)

barplot(data2_num,main = "Grouped barchart",
              xlab = "", ylab = "Score",
              col = c("darkgrey", "darkblue", "red"),
           legend.text = rownames(data2_num),beside=TRUE)
        
        