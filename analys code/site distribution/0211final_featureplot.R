load('cervical_feature_rate.Rdata')
load('colonic_feature_rate.Rdata')
load('gastric_feature_rate.Rdata')
load('melanoma_feature_rate.Rdata')
load('prostatic_feature_rate.Rdata')
load('liver_feature_rate.Rdata')
hyper_back_to_normal_temp=data.frame(Var1=c("3'UTR","5'UTR"),variable=c('colonic','colonic'),value=c(0,0))
colonic_hyper_back_to_normal_feature=rbind(colonic_hyper_back_to_normal_feature,hyper_back_to_normal_temp)
rownames(colonic_hyper_back_to_normal_feature)=colonic_hyper_back_to_normal_feature$Var1
colonic_hyper_back_to_normal_feature=colonic_hyper_back_to_normal_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hypo_back_to_normal_temp=data.frame(Var1=c('1stExon','TSS200'),variable=c('colonic','colonic'),value=c(0,0))
colonic_hypo_back_to_normal_feature=rbind(colonic_hypo_back_to_normal_feature,hypo_back_to_normal_temp)
rownames(colonic_hypo_back_to_normal_feature)=colonic_hypo_back_to_normal_feature$Var1
colonic_hypo_back_to_normal_feature=colonic_hypo_back_to_normal_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hyper_hyper_temp=data.frame(Var1=c('1stExon' ,"3'UTR",'TSS200'),variable=c('colonic','colonic','colonic'),value=c(0,0,0))
colonic_hyper_hyper_feature=rbind(colonic_hyper_hyper_feature,hyper_hyper_temp)
rownames(colonic_hyper_hyper_feature)=colonic_hyper_hyper_feature$Var1
colonic_hyper_hyper_feature=colonic_hyper_hyper_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hypo_hyper_temp=data.frame(Var1=c("3'UTR",'TSS1500'),variable=c('colonic','colonic'),value=c(0,0))
colonic_hypo_hyper_feature=rbind(colonic_hypo_hyper_feature,hypo_hyper_temp)
rownames(colonic_hypo_hyper_feature)=colonic_hypo_hyper_feature$Var1
colonic_hypo_hyper_feature=colonic_hypo_hyper_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hyper_hypo_temp=data.frame(Var1=c("3'UTR",'IGR','1stExon','TSS200'),variable=c('colonic','colonic','colonic','colonic'),value=c(0,0,0,0))
colonic_hyper_hypo_feature=rbind(colonic_hyper_hypo_feature,hyper_hypo_temp)
rownames(colonic_hyper_hypo_feature)=colonic_hyper_hypo_feature$Var1
colonic_hyper_hypo_feature=colonic_hyper_hypo_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hypo_hypo_temp=data.frame(Var1=c("3'UTR"  ,  'IGR'),variable=rep('melanoma',2),value=c(0,0))
melanoma_hypo_hypo_feature=rbind(melanoma_hypo_hypo_feature,hypo_hypo_temp)
rownames(melanoma_hypo_hypo_feature)=melanoma_hypo_hypo_feature$Var1
melanoma_hypo_hypo_feature=melanoma_hypo_hypo_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]



hyper_stable_temp=data.frame(Var1=c("3'UTR", "5'UTR"  , 'Body' ,  'IGR'  ),variable=c('melanoma','melanoma','melanoma','melanoma'),value=c(0,0,0,0))
melanoma_hyper_stable_feature=rbind(melanoma_hyper_stable_feature,hyper_stable_temp)
rownames(melanoma_hyper_stable_feature)=melanoma_hyper_stable_feature$Var1
melanoma_hyper_stable_feature=melanoma_hyper_stable_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hypo_stable_temp=data.frame(Var1=c('1stExon' ,"3'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),variable=rep('melanoma',6),value=c(0,0,0,0,0,0))
melanoma_hypo_stable_feature=rbind(melanoma_hypo_stable_feature,hypo_stable_temp)
rownames(melanoma_hypo_stable_feature)=melanoma_hypo_stable_feature$Var1
melanoma_hypo_stable_feature=melanoma_hypo_stable_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

hyper_back_to_normal_temp=data.frame(Var1=c("3'UTR"),variable=c('melanoma'),value=c(0))
melanoma_hyper_back_to_normal_feature=rbind(melanoma_hyper_back_to_normal_feature,hyper_back_to_normal_temp)
rownames(melanoma_hyper_back_to_normal_feature)=melanoma_hyper_back_to_normal_feature$Var1
melanoma_hyper_back_to_normal_feature=melanoma_hyper_back_to_normal_feature[c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),]

zero=data.frame(Var1=c('1stExon' ,"3'UTR", "5'UTR"  , 'Body' ,  'IGR'  , 'TSS1500' ,'TSS200'),variable=rep('',7),value=rep(7,7))


colonic_hyper_stable_feature=zero
colonic_hypo_stable_feature=zero

library(plyr)

background_feature=rbind(cervical_background_feature,colonic_background_feature,gastric_background_feature,prostatic_background_feature,melanoma_background_feature,liver_background_feature)

normal_pre_feature=rbind(cervical_normal_pre_feature,colonic_normal_pre_feature,gastric_normal_pre_feature,prostatic_normal_pre_feature,melanoma_normal_pre_feature,liver_normal_pre_feature)
normal_pre_feature$value[is.na(normal_pre_feature$value)]=0
normal_pre_feature$value=normal_pre_feature$value/background_feature$value
colnames(normal_pre_feature)=c('feature','Group','Frequency')
normal_pre_feature$site=rep('normal_pre',dim(normal_pre_feature)[1])

normal_cancer_feature=rbind(cervical_normal_cancer_feature,colonic_normal_cancer_feature,gastric_normal_cancer_feature,prostatic_normal_cancer_feature,melanoma_normal_cancer_feature,liver_normal_cancer_feature)
normal_cancer_feature$value[is.na(normal_cancer_feature$value)]=0
normal_cancer_feature$value=normal_cancer_feature$value/background_feature$value
colnames(normal_cancer_feature)=c('feature','Group','Frequency')
normal_cancer_feature$site=rep('normal_cancer',dim(normal_cancer_feature)[1])

hypo_hypo_feature=rbind(cervical_hypo_hypo_feature,colonic_hypo_hypo_feature,gastric_hypo_hypo_feature,prostatic_hypo_hypo_feature,melanoma_hypo_hypo_feature,liver_hypo_hypo_feature)
hypo_hypo_feature$value[is.na(hypo_hypo_feature$value)]=0
hypo_hypo_feature$value=hypo_hypo_feature$value/background_feature$value
colnames(hypo_hypo_feature)=c('feature','Group','Frequency')
hypo_hypo_feature$site=rep('hypo_hypo',dim(hypo_hypo_feature)[1])

hyper_hyper_feature=rbind(cervical_hyper_hyper_feature,colonic_hyper_hyper_feature,gastric_hyper_hyper_feature,prostatic_hyper_hyper_feature,melanoma_hyper_hyper_feature,liver_hyper_hyper_feature)
hyper_hyper_feature$value[is.na(hyper_hyper_feature$value)]=0
hyper_hyper_feature$value=hyper_hyper_feature$value/background_feature$value
colnames(hyper_hyper_feature)=c('feature','Group','Frequency')
hyper_hyper_feature$site=rep('hyper_hyper',dim(hyper_hyper_feature)[1])

hyper_hypo_feature=rbind(cervical_hyper_hypo_feature,colonic_hyper_hypo_feature,gastric_hyper_hypo_feature,prostatic_hyper_hypo_feature,melanoma_hyper_hypo_feature,liver_hyper_hypo_feature)
hyper_hypo_feature$value[is.na(hyper_hypo_feature$value)]=0
hyper_hypo_feature$value=hyper_hypo_feature$value/background_feature$value
colnames(hyper_hypo_feature)=c('feature','Group','Frequency')
hyper_hypo_feature$site=rep('hyper_hypo',dim(hyper_hypo_feature)[1])

hypo_hyper_feature=rbind(cervical_hypo_hyper_feature,colonic_hypo_hyper_feature,gastric_hypo_hyper_feature,prostatic_hypo_hyper_feature,melanoma_hypo_hyper_feature,liver_hypo_hyper_feature)
hypo_hyper_feature$value[is.na(hypo_hyper_feature$value)]=0
hypo_hyper_feature$value=hypo_hyper_feature$value/background_feature$value
colnames(hypo_hyper_feature)=c('feature','Group','Frequency')
hypo_hyper_feature$site=rep('hypo_hyper',dim(hypo_hyper_feature)[1])

hyper_stable_feature=rbind(cervical_hyper_stable_feature,colonic_hyper_stable_feature,gastric_hyper_stable_feature,prostatic_hyper_stable_feature,melanoma_hyper_stable_feature,liver_hyper_stable_feature)
hyper_stable_feature$value[is.na(hyper_stable_feature$value)]=0
hyper_stable_feature$value=hyper_stable_feature$value/background_feature$value
colnames(hyper_stable_feature)=c('feature','Group','Frequency')
hyper_stable_feature$site=rep('hyper_stable',dim(hyper_stable_feature)[1])

hypo_stable_feature=rbind(cervical_hypo_stable_feature,colonic_hypo_stable_feature,gastric_hypo_stable_feature,prostatic_hypo_stable_feature,melanoma_hypo_stable_feature,liver_hypo_stable_feature)
hypo_stable_feature$value[is.na(hypo_stable_feature$value)]=0
hypo_stable_feature$value=hypo_stable_feature$value/background_feature$value
colnames(hypo_stable_feature)=c('feature','Group','Frequency')
hypo_stable_feature$site=rep('hypo_stable',dim(hypo_stable_feature)[1])

hypo_back_to_normal_feature=rbind(cervical_hypo_back_to_normal_feature,colonic_hypo_back_to_normal_feature,gastric_hypo_back_to_normal_feature,prostatic_hypo_back_to_normal_feature,melanoma_hypo_back_to_normal_feature,liver_hypo_back_to_normal_feature)
hypo_back_to_normal_feature$value[is.na(hypo_back_to_normal_feature$value)]=0
hypo_back_to_normal_feature$value=hypo_back_to_normal_feature$value/background_feature$value
colnames(hypo_back_to_normal_feature)=c('feature','Group','Frequency')
hypo_back_to_normal_feature$site=rep('hypo_back_to_normal',dim(hypo_back_to_normal_feature)[1])

hyper_back_to_normal_feature=rbind(cervical_hyper_back_to_normal_feature,colonic_hyper_back_to_normal_feature,gastric_hyper_back_to_normal_feature,prostatic_hyper_back_to_normal_feature,melanoma_hyper_back_to_normal_feature,liver_hyper_back_to_normal_feature)
hyper_back_to_normal_feature$value[is.na(hyper_back_to_normal_feature$value)]=0
hyper_back_to_normal_feature$value=hyper_back_to_normal_feature$value/background_feature$value
colnames(hyper_back_to_normal_feature)=c('feature','Group','Frequency')
hyper_back_to_normal_feature$site=rep('hyper_back_to_normal',dim(hyper_back_to_normal_feature)[1])


feature_data=rbind(normal_pre_feature,normal_cancer_feature,hypo_hypo_feature,hyper_hyper_feature,hyper_hypo_feature,hypo_hyper_feature,hypo_stable_feature,hyper_stable_feature,hypo_back_to_normal_feature,hyper_back_to_normal_feature)
feature_data$site=mapvalues(feature_data$site,c('normal_pre','normal_cancer'),c('normal_v_pre','normal_v_cancer'))

feature_data$site=mapvalues(feature_data$site,c('hypo_back_to_normal','hyper_back_to_normal'),c('hypo_RtoN','hyper_RtoN'))
feature_data$site=factor(feature_data$site,levels = c('normal_v_pre','normal_v_cancer','hyper_hyper', 'hyper_stable','hyper_RtoN','hyper_hypo', 'hypo_hypo','hypo_stable','hypo_RtoN', 'hypo_hyper'))
feature_data$Group=mapvalues(feature_data$Group,c('melanoma','liver'),c('cutaneous','hepatic'))
feature_data$Group=factor(feature_data$Group,levels=c('cervical','colonic','gastric','prostatic','cutaneous','hepatic'))
library(ggplot2)
library(ggthemes)
feature_data=na.omit(feature_data)
color=c('1stExon'='#4D3D51',"5'UTR"='#CA4528',"3'UTR"='#EFBA51','Body'='#4D79A6','ExonBnd'='#BAB2BA','IGR'='4D79A6','TSS1500'='#9EC6BE','TSS200'='#AA9E92')
ggplot(feature_data,aes(site,weight=Frequency,fill=feature))+ theme_bw()+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=15))+scale_fill_manual(values = color)+facet_wrap(~Group)+ labs(x="site",y="score")+ scale_y_continuous(breaks = seq(0, 13, len = 14))+theme(legend.position = "bottom" ,legend.box = "horizontal") 
ggsave('0320feature_duiji.pdf',width = 9,height = 7)
