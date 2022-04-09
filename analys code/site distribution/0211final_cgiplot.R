load('cervical_cgi_rate.Rdata')
load('colonic_cgi_rate.Rdata')
load('gastric_cgi_rate.Rdata')
load('melanoma_cgi_rate.Rdata')
load('prostatic_cgi_rate.Rdata')
load('liver_cgi_rate.Rdata')
library(plyr)
background_cgi=rbind(cervical_background_cgi,colonic_background_cgi,gastric_background_cgi,prostatic_background_cgi,melanoma_background_cgi,liver_background_cgi)

normal_pre_cgi=rbind(cervical_normal_pre_cgi,colonic_normal_pre_cgi,gastric_normal_pre_cgi,prostatic_normal_pre_cgi,melanoma_normal_pre_cgi,liver_normal_pre_cgi)
normal_pre_cgi$value[is.na(normal_pre_cgi$value)]=0
normal_pre_cgi$value=normal_pre_cgi$value/background_cgi$value
colnames(normal_pre_cgi)=c('cgi','Group','Frequency')
normal_pre_cgi$site=rep('normal_pre',dim(normal_pre_cgi)[1])

normal_cancer_cgi=rbind(cervical_normal_cancer_cgi,colonic_normal_cancer_cgi,gastric_normal_cancer_cgi,prostatic_normal_cancer_cgi,melanoma_normal_cancer_cgi,liver_normal_cancer_cgi)
normal_cancer_cgi$value[is.na(normal_cancer_cgi$value)]=0
normal_cancer_cgi$value=normal_cancer_cgi$value/background_cgi$value
colnames(normal_cancer_cgi)=c('cgi','Group','Frequency')
normal_cancer_cgi$site=rep('normal_cancer',dim(normal_cancer_cgi)[1])

hypo_hypo_cgi=rbind(cervical_hypo_hypo_cgi,colonic_hypo_hypo_cgi,gastric_hypo_hypo_cgi,prostatic_hypo_hypo_cgi,melanoma_hypo_hypo_cgi,liver_hypo_hypo_cgi)
hypo_hypo_cgi$value[is.na(hypo_hypo_cgi$value)]=0
hypo_hypo_cgi$value=hypo_hypo_cgi$value/background_cgi$value
colnames(hypo_hypo_cgi)=c('cgi','Group','Frequency')
hypo_hypo_cgi$site=rep('hypo_hypo',dim(hypo_hypo_cgi)[1])

hyper_hyper_cgi=rbind(cervical_hyper_hyper_cgi,colonic_hyper_hyper_cgi,gastric_hyper_hyper_cgi,prostatic_hyper_hyper_cgi,melanoma_hyper_hyper_cgi,liver_hyper_hyper_cgi)
hyper_hyper_cgi$value[is.na(hyper_hyper_cgi$value)]=0
hyper_hyper_cgi$value=hyper_hyper_cgi$value/background_cgi$value
colnames(hyper_hyper_cgi)=c('cgi','Group','Frequency')
hyper_hyper_cgi$site=rep('hyper_hyper',dim(hyper_hyper_cgi)[1])

hyper_hypo_cgi=rbind(cervical_hyper_hypo_cgi,colonic_hyper_hypo_cgi,gastric_hyper_hypo_cgi,prostatic_hyper_hypo_cgi,melanoma_hyper_hypo_cgi,liver_hyper_hypo_cgi)
hyper_hypo_cgi$value[is.na(hyper_hypo_cgi$value)]=0
hyper_hypo_cgi$value=hyper_hypo_cgi$value/background_cgi$value
colnames(hyper_hypo_cgi)=c('cgi','Group','Frequency')
hyper_hypo_cgi$site=rep('hyper_hypo',dim(hyper_hypo_cgi)[1])

hypo_hyper_cgi=rbind(cervical_hypo_hyper_cgi,colonic_hypo_hyper_cgi,gastric_hypo_hyper_cgi,prostatic_hypo_hyper_cgi,melanoma_hypo_hyper_cgi,liver_hypo_hyper_cgi)
hypo_hyper_cgi$value[is.na(hypo_hyper_cgi$value)]=0
hypo_hyper_cgi$value=hypo_hyper_cgi$value/background_cgi$value
colnames(hypo_hyper_cgi)=c('cgi','Group','Frequency')
hypo_hyper_cgi$site=rep('hypo_hyper',dim(hypo_hyper_cgi)[1])

hyper_stable_cgi=rbind(cervical_hyper_stable_cgi,colonic_hyper_stable_cgi,gastric_hyper_stable_cgi,prostatic_hyper_stable_cgi,melanoma_hyper_stable_cgi,liver_hyper_stable_cgi)
hyper_stable_cgi$value[is.na(hyper_stable_cgi$value)]=0
hyper_stable_cgi$value=hyper_stable_cgi$value/background_cgi$value
colnames(hyper_stable_cgi)=c('cgi','Group','Frequency')
hyper_stable_cgi$site=rep('hyper_stable',dim(hyper_stable_cgi)[1])

hypo_stable_cgi=rbind(cervical_hypo_stable_cgi,colonic_hypo_stable_cgi,gastric_hypo_stable_cgi,prostatic_hypo_stable_cgi,melanoma_hypo_stable_cgi,liver_hypo_stable_cgi)
hypo_stable_cgi$value[is.na(hypo_stable_cgi$value)]=0
hypo_stable_cgi$value=hypo_stable_cgi$value/background_cgi$value
colnames(hypo_stable_cgi)=c('cgi','Group','Frequency')
hypo_stable_cgi$site=rep('hypo_stable',dim(hypo_stable_cgi)[1])

hypo_back_to_normal_cgi=rbind(cervical_hypo_back_to_normal_cgi,colonic_hypo_back_to_normal_cgi,gastric_hypo_back_to_normal_cgi,prostatic_hypo_back_to_normal_cgi,melanoma_hypo_back_to_normal_cgi,liver_hypo_back_to_normal_cgi)
hypo_back_to_normal_cgi$value[is.na(hypo_back_to_normal_cgi$value)]=0
hypo_back_to_normal_cgi$value=hypo_back_to_normal_cgi$value/background_cgi$value
colnames(hypo_back_to_normal_cgi)=c('cgi','Group','Frequency')
hypo_back_to_normal_cgi$site=rep('hypo_back_to_normal',dim(hypo_back_to_normal_cgi)[1])

hyper_back_to_normal_cgi=rbind(cervical_hyper_back_to_normal_cgi,colonic_hyper_back_to_normal_cgi,gastric_hyper_back_to_normal_cgi,prostatic_hyper_back_to_normal_cgi,melanoma_hyper_back_to_normal_cgi,liver_hyper_back_to_normal_cgi)
hyper_back_to_normal_cgi$value[is.na(hyper_back_to_normal_cgi$value)]=0
hyper_back_to_normal_cgi$value=hyper_back_to_normal_cgi$value/background_cgi$value
colnames(hyper_back_to_normal_cgi)=c('cgi','Group','Frequency')
hyper_back_to_normal_cgi$site=rep('hyper_back_to_normal',dim(hyper_back_to_normal_cgi)[1])


cgi_data=rbind(normal_pre_cgi,normal_cancer_cgi,hypo_hypo_cgi,hyper_hyper_cgi,hyper_hypo_cgi,hypo_hyper_cgi,hypo_stable_cgi,hyper_stable_cgi,hypo_back_to_normal_cgi,hyper_back_to_normal_cgi)
cgi_data$site=mapvalues(cgi_data$site,c('normal_pre','normal_cancer'),c('normal_v_pre','normal_v_cancer'))

cgi_data$site=mapvalues(cgi_data$site,c('hypo_back_to_normal','hyper_back_to_normal'),c('hypo_RtoN','hyper_RtoN'))
cgi_data$site=factor(cgi_data$site,levels = c('normal_v_pre','normal_v_cancer','hyper_hyper', 'hyper_stable','hyper_RtoN','hyper_hypo', 'hypo_hypo','hypo_stable','hypo_RtoN', 'hypo_hyper'))
cgi_data$Group=mapvalues(cgi_data$Group,c('melanoma','liver'),c('cutaneous','hepatic'))
cgi_data$Group=factor(cgi_data$Group,levels=c('cervical','colonic','gastric','prostatic','cutaneous','hepatic'))
library(ggplot2)
library(ggthemes)
color=c('island'='#4D3D51','opensea'='#CA4528','shelf'='#EFBA51','shore'='#1699EF')
ggplot(cgi_data,aes(site,weight=Frequency,fill=cgi))+ theme_bw()+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle =90, hjust = 1,size=15),strip.text=element_text(size=15))+scale_fill_manual(values = color)+facet_wrap(~Group)+ labs(x="site",y="score")+ scale_y_continuous(breaks = seq(0, 13, len = 14))+theme(legend.position = "bottom" ,legend.box = "horizontal")  
ggsave('0320cgi_duiji.pdf',width = 9,height = 7)

