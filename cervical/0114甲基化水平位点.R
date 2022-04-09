load('1118_sum_norm.Rdata')
library(ChAMP)
library(readr)
sum_DMP <- champ.DMP(beta = sum_norm,pheno=pd$label)
save(sum_DMP,file='sum_DMP.Rdata')

normal_to_cin3=sum_DMP$normal_to_cin3
save(normal_to_cin3,file='0128normal_to_cin3.Rdata')
cancer_to_cin3=sum_DMP$cancer_to_cin3
save(cancer_to_cin3,file='0128cancer_to_cin3.Rdata')


library(dplyr)
deg=normal_to_cin3
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
normal_to_cin3_down=filter(deg,g=='DOWN')
normal_to_cin3_up=filter(deg,g=='UP')
save(normal_to_cin3_down,file = '0128DMP normal_to_cin3_down.Rdata')
save(normal_to_cin3_up,file = '0128DMP normal_to_cin3_up.Rdata')

deg=cancer_to_cin3
#这里应该是找从cin3到cancer的低甲基化位点，但是deg结果是相反的，所以之后的up和down也就调换了一下
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
cin3_to_cancer_down=filter(deg,g=='UP')
cin3_to_cancer_up=filter(deg,g=='DOWN')
save(cin3_to_cancer_down,file = '0128DMP cancer_to_cin3_down.Rdata')
save(cin3_to_cancer_up,file = '0128DMP cancer_to_cin3_up.Rdata')

cervical_down_down=intersect(rownames(normal_to_cin3_down),rownames(cin3_to_cancer_down))
save(cervical_down_down,file ='0128DMP cervical_down_down.Rdata' )
cervical_down_down_gene=as.data.frame(table(droplevels((normal_to_cin3_down[cervical_down_down,])$gene)))
save(cervical_down_down_gene,file ='0128DMP cervical_down_down_gene.Rdata' )

cervical_up_up=intersect(rownames(normal_to_cin3_up),rownames(cin3_to_cancer_up))
save(cervical_up_up,file ='0128DMP cervical_up_up.Rdata' )
cervical_up_up_gene=as.data.frame(table(droplevels((normal_to_cin3_up[cervical_up_up,])$gene)))
save(cervical_up_up_gene,file ='0128DMP cervical_up_up_gene.Rdata' )

cervical_up_down=intersect(rownames(normal_to_cin3_up),rownames(cin3_to_cancer_down))
cervical_up_down_gene=as.data.frame(table(droplevels((normal_to_cin3_up[cervical_up_down,])$gene)))
cervical_down_up=intersect(rownames(normal_to_cin3_down),rownames(cin3_to_cancer_up))
cervical_down_up_gene=as.data.frame(table(droplevels((normal_to_cin3_down[cervical_down_up,])$gene)))
save(cervical_up_down,cervical_up_down_gene,cervical_down_up,cervical_down_up_gene,file = '0207cervical_downandup.Rdata')

cervical_down_down_gene=cervical_down_down_gene$Var1
cervical_up_up_gene=cervical_up_up_gene$Var1
  list=c('cervical_down_down','cervical_up_up','cervical_down_down_gene','cervical_up_up_gene')
for (i in 1:length(list)){
  nn=list[i]
  getnn=data.frame(get(nn))
  colnames(getnn)=nn
  write_csv(as.data.frame(nn),append = T,file =paste('0206_cervical',nn,'.csv'))
  write_csv(getnn,append = T,file =paste('0206_cervical',nn,'.csv'))
  }


