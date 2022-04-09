cancerlist=c('prostatic','cervical','liver','gastric','colonic')
countlist=c('100000','200000','300000','400000')
sitelist=c('hypo_hypo','hyper_hyper','hypo_hyper','hyper_hypo','hypo_stable','hyper_stable','hypo_back_to_normal','hyper_back_to_normal')
library(ggplot2)
plot_data=data.frame(matrix(ncol=4,nrow = 0))
colnames(plot_data)=c('cancer','value','site','rate_type')
for (i in cancerlist){
  for (j in countlist ) {
    if (j=='100000'){
      data=read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T)
    }else{
      data=rbind(data,read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T))
    }
  }
  data=rbind(data,data.frame(cancer=i,value=(length(which(data[,3]<0.05))/nrow(data)),site='all',rate_type='1'))
  data=rbind(data,data.frame(cancer=i,value=1-(length(which(data[,3]<0.05))/nrow(data)),site='all',rate_type='0'))

}



for (i in cancerlist){
  data=read.table(paste0('hypo_hypo/',i,'_coef.txt'))
  print(paste0('hypo_hypo_',i,' ',(length(which(data[,3]<0.05))),'/',nrow(data),'  ',(length(which(data[,3]<0.05))/nrow(data))))
}
##每个位点
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
    assign(paste0(i,'_',k),data[which(data[,3]<0.05),1])
    plot_data=rbind(plot_data,data.frame(cancer=i,value=(length(which(data[,3]<0.05))/nrow(data)),site=k,rate_type='1'))
    plot_data=rbind(plot_data,data.frame(cancer=i,value=1-(length(which(data[,3]<0.05))/nrow(data)),site=k,rate_type='0'))
  }
}
color=c('1'='#B51C58','0'='#0094FD')
library(plyr)
plot_data$site=mapvalues(plot_data$site,c('hyper_back_to_normal','hypo_back_to_normal'),c('hyper_RtoN','hypo_RtoN'))
plot_data$site=factor(plot_data$site,levels = c('hyper_hyper', 'hyper_stable','hyper_RtoN','hyper_hypo', 'hypo_hypo','hypo_stable','hypo_RtoN', 'hypo_hyper'))
plot_data$rate_type=factor(plot_data$rate_type,levels = c('1','0'))
ggplot(data.frame(plot_data),aes(x=cancer,weight=value,fill=rate_type))+geom_bar(position="stack")+scale_fill_manual(values = color)+facet_wrap(~site,ncol=4)+ labs(x="DMP",y="rate")+
  theme_bw()+ theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15),strip.text=element_text(size=18))+
  theme(legend.position = "bottom" ,legend.box = "horizontal")+   coord_cartesian(ylim = c(0.5, 1))
ggsave('valuable_site_barplot.pdf')


hyper_hypo_sitelist=c('norm_cancer_hyper','norm_cancer_hypo','norm_pre_hyper','norm_pre_hypo')
#癌前hyper hypo
for (k in hyper_hypo_sitelist){
  
  for (i in cancerlist){
    cgpoint=get(paste0(i,'_',k))
    for (j in countlist ) {
      if (j=='100000'){
        data=read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T)
      }else{
        data=rbind(data,read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T))
      }
    }
    data=data[which(data[,1] %in% cgpoint),]
    print(paste0(k,'_',i,' ',(length(which(data[,3]<0.05))),'/',nrow(data),'  ',(length(which(data[,3]<0.05))/nrow(data))))
  }
  print('\n')
}
