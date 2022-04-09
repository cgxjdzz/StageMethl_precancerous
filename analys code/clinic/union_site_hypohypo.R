cancerlist=c('prostatic','cervical','liver','gastric','colonic')
countlist=c('100000','200000','300000','400000')
sitelist=c('hypo_hypo','hyper_hyper','hypo_hyper','hyper_hypo','hypo_stable','hyper_stable','hypo_back_to_normal','hyper_back_to_normal')
library(ggplot2)
plot_data_value=data.frame(matrix(ncol=4,nrow = 0))
colnames(plot_data_value)=c('cancer','value','site','rate_type')
plot_data_all=data.frame(matrix(ncol=4,nrow = 0))
colnames(plot_data_all)=c('cancer','value','site','rate_type')
for (k in sitelist){
  ee=read.table(paste0('union_site_',k,'/union_site_',k,'.csv'))
  cgpoint=rownames(ee)[which(ee$sum>1)]
  for (i in cancerlist){
    for (j in countlist ) {
      if (j=='100000'){
        data=read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T)
      }else{
        data=rbind(data,read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T))
      }
      
    }
    data=data[which(data[,1] %in% cgpoint),]
    plot_data_value=rbind(plot_data_value,data.frame(cancer=i,value=length(which(data[,3]<0.05)),site=k,rate_type='value'))
    plot_data_all=rbind(plot_data_all,data.frame(cancer=i,value=nrow(data),site=k,rate_type='all'))
  }
}
for (i in cancerlist){
  plot_data_all[which(plot_data_all$cancer==i),'value']=plot_data_all[which(plot_data_all$cancer==i),'value']/sum(plot_data_all[which(plot_data_all$cancer==i),'value'])
}
for (i in cancerlist){
  plot_data_value[which(plot_data_value$cancer==i),'value']=plot_data_value[which(plot_data_value$cancer==i),'value']/sum(plot_data_value[which(plot_data_value$cancer==i),'value'])
}
plot_data=rbind(plot_data_value,plot_data_all)
# color=c('1'='#B51C58','0'='#0094FD')
library(plyr)
plot_data$rate_type=factor(plot_data$rate_type,levels = c('value','all'))
plot_data$site=mapvalues(plot_data$site,c('hyper_back_to_normal','hypo_back_to_normal'),c('hyper_RtoN','hypo_RtoN'))
plot_data$site=factor(plot_data$site,levels = c('hyper_hyper', 'hyper_stable','hyper_RtoN','hyper_hypo', 'hypo_hypo','hypo_stable','hypo_RtoN', 'hypo_hyper'))
sitecolor=c('hyper_hyper'='#E39A44','hypo_hypo'='#47864E','hypo_hyper'='#3290CF','hyper_hypo'='#C04C36',
            'hyper_RtoN'='#E8D2B7','hypo_RtoN'='#C1E5DF','hyper_stable'='#D3C856','hypo_stable'='#4CAC6E')
ggplot(data.frame(plot_data),aes(x=rate_type,weight=value,fill=site))+geom_bar(position="stack")+facet_wrap(~cancer,ncol=5)+scale_fill_manual(values = sitecolor)+ labs(x="DMP",y="rate")+
  theme_bw()+ theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=18))+
  theme(legend.position = "bottom" ,legend.box = "horizontal")+   coord_cartesian(ylim = c(0.5, 1))
ggsave('valuable_site_rate_barplot.pdf')
