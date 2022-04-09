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
    plot_data_value=rbind(plot_data_value,data.frame(cancer=i,value=(length(which(data[,3]<0.05))),site=k,rate_type='value'))
    plot_data_all=rbind(plot_data_all,data.frame(cancer=i,value=nrow(data),site=k,rate_type='all'))
  }
}

for (i in cancerlist) {
  for (j in countlist ) {
    if (j=='100000'){
      data=read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T)
    }else{
      data=rbind(data,read.table(paste0(i,'/',i,'_',j,'_coef.txt'),fill=T))
    }
  }
  data$site='ee'
  for (k in sitelist) {
    ee=read.table(paste0('union_site_',k,'/union_site_',k,'.csv'))
    cgpoint=rownames(ee)[which(ee$sum>1)]
    data$site[which(data[,1] %in% cgpoint)]=k
  }
  data=data[-which(data$site=='ee'),]
  data=data[order(data[,3]),]
  print(i)
  print(table(data[c(1:50),'site']))
}
