cgi_duiji=rbind(colonic_cgi_duijiplot,cervical_cgi_duijiplot,gastric_cgi_duijiplot,prostatic_cgi_duijiplot,melanoma_cgi_duijiplot,liver_cgi_duijiplot)
feature_duiji=rbind(colonic_feature_duijiplot,cervical_feature_duijiplot,gastric_feature_duijiplot,prostatic_feature_duijiplot,melanoma_feature_duijiplot,liver_feature_duijiplot)
cgi_duiji$cancer=factor(cgi_duiji$cancer,levels=c('cervical','colonic','gastric','prostatic','melanoma','liver'))
color=c('island'='#4D3D51','opensea'='#CA4528','shelf'='#EFBA51','shore'='#1699EF')
ggplot(cgi_duiji,aes(variable,weight=value,fill=Var1))+ theme_bw()+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=20))+scale_fill_manual(values = color)+facet_wrap(~cancer)+ labs(x="site",y="score")+ scale_y_continuous(breaks = seq(0, 8, len = 9))+theme(legend.position = "bottom" ,legend.box = "horizontal") 
ggsave('cgi_duiji.pdf',height = 7,width = 9)

feature_duiji$cancer=factor(feature_duiji$cancer,levels=c('cervical','colonic','gastric','prostatic','melanoma','liver'))
color=c('1stExon'='#4D3D51',"5'UTR"='#CA4528',"3'UTR"='#EFBA51','Body'='#4D79A6','ExonBnd'='#BAB2BA','IGR'='4D79A6','TSS1500'='#9EC6BE','TSS200'='#AA9E92')
ggplot(feature_duiji,aes(variable,weight=value,fill=Var1))+ theme_bw()+ theme(panel.grid=element_blank())+geom_bar(position="stack")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=20))+scale_fill_manual(values = color)+facet_wrap(~cancer)+ labs(x="site",y="score")+ scale_y_continuous(breaks = seq(0, 10, len = 11))+theme(legend.position = "bottom" ,legend.box = "horizontal")
ggsave('feature_duiji.pdf',height = 7,width = 9)

