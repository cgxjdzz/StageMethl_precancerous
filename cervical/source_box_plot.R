library(ggnewscale)
library(ggplot2)
mean=apply(sum_combat,2,mean)
data2=data.frame(type=pd$label,mean=mean,batch=pd$branch)
data2$cancer='cervical'
cervical_box=data2
save(cervical_box,file = 'cervical_box.Rdata')
nana<-ggplot()+geom_jitter(alpha = 0.3,size=3,aes(x =data2$type, y =data2$mean ,color=data2$batch))+ #分组的点图
  scale_x_discrete(limit=c("normal","cin3","cancer") )+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('cervical carcinoma')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =data2$type, y =data2$mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 45))+
  labs(x="sample",y="value",fill= "type")

ggsave('0128cervical_source_box_plot.pdf')



nana

