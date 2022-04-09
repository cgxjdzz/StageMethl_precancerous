load('cervical_cgi_rate.Rdata')
load('colonic_cgi_rate.Rdata')
load('gastric_cgi_rate.Rdata')
load('melanoma_cgi_rate.Rdata')
load('prostatic_cgi_rate.Rdata')

background_cgi=rbind(cervical_background_cgi,colonic_background_cgi,gastric_background_cgi,melanoma_background_cgi,prostatic_background_cgi)
colnames(background_cgi)=c('cgi','Group','Frequency')


normal_pre_cgi=rbind(cervical_normal_pre_cgi,colonic_normal_pre_cgi,gastric_normal_pre_cgi,melanoma_normal_pre_cgi,prostatic_normal_pre_cgi)
colnames(normal_pre_cgi)=c('cgi','Group','Frequency')

pre_cancer_cgi=rbind(cervical_pre_cancer_cgi,colonic_pre_cancer_cgi,gastric_pre_cancer_cgi,melanoma_pre_cancer_cgi,prostatic_pre_cancer_cgi)
colnames(pre_cancer_cgi)=c('cgi','Group','Frequency')

down_down_cgi=rbind(cervical_down_down_cgi,colonic_down_down_cgi,gastric_down_down_cgi,melanoma_down_down_cgi,prostatic_down_down_cgi)
colnames(down_down_cgi)=c('cgi','Group','Frequency')

up_up_cgi=rbind(cervical_up_up_cgi,colonic_up_up_cgi,gastric_up_up_cgi,melanoma_up_up_cgi,prostatic_up_up_cgi)
colnames(up_up_cgi)=c('cgi','Group','Frequency')

up_down_cgi=rbind(cervical_up_down_cgi,colonic_up_down_cgi,gastric_up_down_cgi,melanoma_up_down_cgi,prostatic_up_down_cgi)
colnames(up_down_cgi)=c('cgi','Group','Frequency')

down_up_cgi=rbind(cervical_down_up_cgi,colonic_down_up_cgi,gastric_down_up_cgi,melanoma_down_up_cgi,prostatic_down_up_cgi)
colnames(down_up_cgi)=c('cgi','Group','Frequency')

library(ggplot2)

p1=ggplot(background_cgi, aes(x = Group, y = Frequency, fill = cgi)) +
  geom_bar(position = position_stack(), stat ="identity", width = .7) +
  coord_flip()
p1
