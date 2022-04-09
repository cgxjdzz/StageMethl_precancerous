load('cervical_rate.Rdata')
load('colonic_rate.Rdata')
load('gastric_rate.Rdata')
load('melanoma_rate.Rdata')
load('prostatic_rate.Rdata')
cgidata=rbind(cervical_cgi,colonic_cgi,melanoma_cgi,prostatic_cgi,gastric_cgi)

featuredata=rbind(cervical_feature,colonic_feature,melanoma_feature,prostatic_feature,gastric_feature)


p1=ggplot(cgidata, aes(x = group, y = counts, fill = cgi)) +
  geom_bar(position = position_stack(), stat ="identity", width = .7) +
  coord_flip()
ggsave(p1,filename = 'cgi_count.pdf')

p2=ggplot(featuredata, aes(x = group, y = counts, fill = feature)) +
  geom_bar(position = position_stack(), stat ="identity", width = .7) +
  coord_flip()
ggsave(p2,filename = 'feature_count.pdf')




library(ggplot2)
