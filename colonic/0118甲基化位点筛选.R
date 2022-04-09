LGA_norm=sum_norm[,which(pd$label=='LGA')]
HGA_norm=sum_norm[,which(pd$label=='HGA')]
cancer_norm=sum_norm[,which(pd$label=='Tumor')]
source('zffdatamining_gene.R')
#筛选从LGA到HGA的位点
LGA_to_HGA_down=c()
LGA_to_HGA_up=c()

for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  temp=ge.groupdifference(data=as.numeric(LGA_norm[i,]),cut = 0.7)
  if(is.na(LGA_2_HGA[cgname,]$logFC)){
    next
  }
  if((LGA_2_HGA[cgname,]$logFC)<temp[1]){
    LGA_to_HGA_down=c(LGA_to_HGA_down,cgname)
  }
  if((LGA_2_HGA[cgname,]$logFC)>temp[2]){
    LGA_to_HGA_up=c(LGA_to_HGA_up,cgname)
  }
}
save(LGA_to_HGA_down,file='LGA_to_HGA_down.Rdata')
save(LGA_to_HGA_up,file='LGA_to_HGA_up.Rdata')

#筛选从HGA到Tumor的位点
HGA_to_Tumor_down=c()
HGA_to_Tumor_up=c()

for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  temp=ge.groupdifference(data=as.numeric(HGA_norm[i,]),cut = 0.7)
  if(is.na(HGA_2_Tumor[cgname,]$logFC)){
    next
  }
  if((HGA_2_Tumor[cgname,]$logFC)<temp[1]){
    HGA_to_Tumor_down=c(HGA_to_Tumor_down,cgname)
  }
  if((HGA_2_Tumor[cgname,]$logFC)>temp[2]){
    HGA_to_Tumor_up=c(HGA_to_Tumor_up,cgname)
  }
}
save(HGA_to_Tumor_down,file='HGA_to_Tumor_down.Rdata')
save(HGA_to_Tumor_up,file='HGA_to_Tumor_up.Rdata')


colonic_up_up=intersect(rownames(HGA_to_Tumor_up),rownames(LGA_to_HGA_up))
save(up_up,file = '10.18DMPcolonic_up_up.Rdata')

colonic_down_down=intersect(rownames(HGA_to_Tumor_down),rownames(LGA_to_HGA_down))
save(down_down,file ='10.18DMPcolonic_down_down.Rdata' )
