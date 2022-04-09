LGA_2_Normal=myDMP$LGA_to_Normal
save(LGA_2_Normal,file = '10.18DMP_LGA_2_Normal.Rdata')
LGA_2_HGA=myDMP$LGA_to_HGA
save(LGA_2_HGA,file = '10.18DMP_LGA_2_HGA.Rdata')
LGA_2_Tumor=myDMP$LGA_to_Tumor
save(LGA_2_Tumor,file = '10.18DMP_LGA_2_Tumor.Rdata')
Normal_2_HGA=myDMP$Normal_to_HGA
save(Normal_2_HGA,file = '10.18DMPNormal_2_HGA.Rdata')
Normal_2_Tumor=myDMP$Normal_to_Tumor
save(Normal_2_Tumor,file = '10.18DMPNormal_2_Tumor.Rdata')
HGA_2_Tumor=myDMP$HGA_to_Tumor
save(HGA_2_Tumor,file = '10.18DMPHGA_2_Tumor.Rdata')

library(dplyr)
deg=LGA_2_HGA
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
LGA_2_HGA_down=filter(deg,g=='DOWN')
LGA_2_HGA_up=filter(deg,g=='UP')
save(LGA_2_HGA_down,file = '10.18DMP_LGA_2_HGA_down.Rdata')
save(LGA_2_HGA_up,file = '10.18DMP_LGA_2_HGA_up.Rdata')

deg=HGA_2_Tumor
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
HGA_2_Tumor_down=filter(deg,g=='DOWN')
HGA_2_Tumor_up=filter(deg,g=='UP')
save(HGA_2_Tumor_down,file = '10.18DMP_HGA_2_Tumor_down.Rdata')
save(HGA_2_Tumor_up,file = '10.18DMP_HGA_2_Tumor_up.Rdata')


up_up=intersect(rownames(HGA_2_Tumor_up),rownames(LGA_2_HGA_up))
save(up_up,file = '10.18DMPup_up.Rdata')

down_down=intersect(rownames(HGA_2_Tumor_down),rownames(LGA_2_HGA_down))
save(down_down,file ='10.18DMPdown_down.Rdata' )

colonic_up_down=intersect(rownames(LGA_2_HGA_up),rownames(HGA_2_Tumor_down))
colonic_up_down_gene=as.data.frame(table(droplevels((LGA_2_HGA_up[colonic_up_down,])$gene)))
colonic_down_up=intersect(rownames(LGA_2_HGA_down),rownames(HGA_2_Tumor_up))
colonic_down_up_gene=as.data.frame(table(droplevels((LGA_2_HGA_down[colonic_down_up,])$gene)))
save(colonic_up_down,colonic_up_down_gene,colonic_down_up,colonic_down_up_gene,file = '0207colonic_downandup.Rdata')




down_down_gene=as.data.frame(table(droplevels((HGA_2_Tumor_down[down_down,])$gene)))
up_up_gene=as.data.frame(table(droplevels((HGA_2_Tumor_up[up_up,])$gene)))
down_down_gene=down_down_gene$Var1
up_up_gene=up_up_gene$Var1
list=c('down_down','up_up','down_down_gene','up_up_gene')
library(readr)
for (i in 1:length(list)){
        nn=list[i]
        getnn=data.frame(get(nn))
        colnames(getnn)=nn
        write_csv(as.data.frame(nn),append = T,file =paste('0206_colonic',nn,'.csv'))
        write_csv(getnn,append = T,file =paste('0206_colonic',nn,'.csv'))
}

plot(density(as.numeric(HGA_2_Tumor$logFC[which(HGA_2_Tumor$logFC>0)])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,10),
     xlab="LogFC",ylab="Density",main="Density Plot")
plot(density(as.numeric(HGA_2_Tumor$logFC[which(HGA_2_Tumor$logFC<0)])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,20),
     xlab="LogFC",ylab="Density",main="Density Plot")

plot(density(as.numeric(LGA_2_HGA$logFC[which(LGA_2_HGA$logFC>0)])),
     col="red",lwd=2.5,bty="l" ,ylim=c(0,10),
     xlab="LogFC",ylab="Density",main="Density Plot")
