normal=sum_combat[,which(pd$label=='normal')]
source('zffdatamining_gene.R')

normal_to_cin3_down=c()
normal_to_cin3_up=c()
normal_norm=sum_combat[,which(pd$label=='normal')]
for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  temp=ge.groupdifference(data=as.numeric(cin3_norm[i,]),cut = 0.9)
  if(is.na(normal_to_cin3[cgname,]$logFC)){
    next
  }
  if((normal_to_cin3[cgname,]$logFC)<temp[1]){
    normal_to_cin3_down=c(normal_to_cin3_down,cgname)
  }
  if((normal_to_cin3[cgname,]$logFC)>temp[2]){
    normal_to_cin3_up=c(normal_to_cin3_up,cgname)
  }
}


cin3_to_cancer_down=c()
cin3_to_cancer_up=c()
cin3_norm=sum_combat[,which(pd$label=='cin3')]
cancer_norm=sum_combat[,which(pd$label=='cancer')]
for (i in 1:dim(sum_norm)[1]) {
  cgname=rownames(sum_norm)[i]
  temp=ge.groupdifference(data=as.numeric(cin3_norm[i,]),cut = 0.9)
  if(is.na(cancer_to_cin3[cgname,]$logFC)){
    next
  }
  if((cancer_to_cin3[cgname,]$logFC)<temp[1]){
    cin3_to_cancer_up=c(cin3_to_cancer_up,cgname)
  }
  if((cancer_to_cin3[cgname,]$logFC)>temp[2]){
    cin3_to_cancer_down=c(cin3_to_cancer_down,cgname)
  }
}

select_down_down=intersect(cin3_to_cancer_down,normal_to_cin3_down)
select_up_up=intersect(cin3_to_cancer_up,normal_to_cin3_up)
