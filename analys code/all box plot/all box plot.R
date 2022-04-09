load('cervical_box.Rdata')
load('colonic_box.Rdata')
load('gastric_box.Rdata')
load('liver_box.Rdata')
load('melanoma_box.Rdata')
load('prostatic_box.Rdata')
plot_data=rbind(cervical_box,colonic_box,gastric_box,liver_box,melanoma_box,prostatic_box)
library(ggplot2)
library(ggnewscale)
p1<-ggplot(cervical_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("normal",'cin3','cancer') ,labels=c('normal\n(n=23)','precancer\n(n=17)','cancer\n(n=316)'))+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('cervical')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
p1
ggsave('cervical source box plot.pdf')

p2<-ggplot(colonic_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("Normal",'LGA','HGA','Tumor') ,labels=c('normal\n(n=45)','LGA\n(n=18)','HGA\n(n=22)','cancer\n(n=147)'))+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('colonic')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
ggsave('colonic source box plot.pdf')

p3<-ggplot(gastric_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("Solid Tissue Normal",'intestinal metaplasia','Primary Tumor') ,labels=c('normal\n(n=44)','precancer\n(n=76)','cancer\n(n=396)'))+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('gastric')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
ggsave('gastric source box plot.pdf')

p4<-ggplot(prostatic_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("Solid Tissue Normal",'CAF','Primary Tumor'),labels=c('normal\n(n=50)','precancer\n(n=18)','cancer\n(n=502)') )+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('prostatic')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
ggsave('prostatic source box plot.pdf')

p5<-ggplot(melanoma_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("Normal",'Precancer','Tumor') ,labels=c('normal\n(n=10)','precancer\n(n=3)','cancer\n(n=233)'))+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('cutaneous')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
ggsave('cutaneous source box plot.pdf')

p6<-ggplot(liver_box)+geom_jitter(alpha = 0.3,size=3,aes(x =type, y =mean ,color=batch))+ #分组的点图
  scale_x_discrete(limit=c("Normal",'Cirrhosis','Tumor') ,labels=c('normal\n(n=50)','precancer\n(n=130)','cancer\n(n=380)'))+#x轴的顺序，discrete离散
  
  scale_color_discrete('source')+#设置颜色及图例，离散
  ggtitle('hepatic')+#图标题
  new_scale('color')+#新scale设置
  geom_boxplot(alpha = .5,size=1,aes(x =type, y =mean))+#箱图
  theme_bw() + 
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(panel.grid =element_blank())+  
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 0,size=10))+
  labs(x="sample",y="value",fill= "type")
ggsave('hepatic source box plot.pdf')
  
p<- cowplot::plot_grid(p1, p2, p3, p4,p5,p6, nrow = 2, labels = LETTERS[1:4])

p= ggpubr::ggarrange(p1, p2, p3, p4,p5,p6, nrow = 2, ncol = 3,vjust = 1)#将p1-p4四幅图组合成一幅图，按照两行两列排列，标签分别为A、B、C、D，颜色为红色(通过font.label = list()修改)，无法通过label.color = 'red'或其他方式修改。
ggsave(filename = 'boxplot_all_n.pdf',plot=p,width = 7.31*3,height = 5.75*2,limitsize = FALSE)
p
library(grid)
grid.newpage()
# pushViewport函数提供了添加视图以及在树中的视图之间导航的方法。
pushViewport(viewport(layout = grid.layout(2,3)))
# viewport函数创建视图，描述图形设备上的矩形区域，并在这些区域中定义许多坐标系统。

vplayout <- function(x,y) +
  viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
print(p3, vp = vplayout(1, 3))
print(p4, vp = vplayout(2, 1))
print(p5, vp = vplayout(2, 2))
print(p6, vp = vplayout(2, 3))