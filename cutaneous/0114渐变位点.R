load('1109_sum_norm.Rdata')
library(ChAMP)
pd$label=gsub(' ','_',pd$label)
sum_DMP <- champ.DMP(beta = sum_norm,pheno=pd$label)



melanocytes_to_skin_punch_biopsy=sum_DMP$melanocytes_to_skin_punch_biopsy
save(melanocytes_to_skin_punch_biopsy,file='melanocytes_to_skin_punch_biopsy.Rdata')
melanocytes_to_primary_melanoma=sum_DMP$melanocytes_to_primary_melanoma
save(melanocytes_to_primary_melanoma,file='melanocytes_to_primary_melanoma.Rdata')


library(dplyr)
deg=melanocytes_to_skin_punch_biopsy
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
skin_punch_biopsy_to_melanocytes_stable=filter(deg,g=='stable')
skin_punch_biopsy_to_melanocytes_down=filter(deg,g=='UP')
skin_punch_biopsy_to_melanocytes_up=filter(deg,g=='DOWN')
save(skin_punch_biopsy_to_melanocytes_down,file = 'DMP skin_punch_biopsy_to_melanocytes_down.Rdata')
save(skin_punch_biopsy_to_melanocytes_up,file = 'DMP skin_punch_biopsy_to_melanocytes_up.Rdata')

deg=melanocytes_to_primary_melanoma
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
melanocytes_to_primary_melanoma_down=filter(deg,g=='DOWN')
melanocytes_to_primary_melanoma_up=filter(deg,g=='UP')
save(melanocytes_to_primary_melanoma_down,file = 'DMP melanocytes_to_primary_melanoma_down.Rdata')
save(melanocytes_to_primary_melanoma_up,file = 'DMP melanocytes_to_primary_melanoma_up.Rdata')

melanoma_down_down=intersect(rownames(skin_punch_biopsy_to_melanocytes_down),rownames(melanocytes_to_primary_melanoma_down))
#根本没有啊，持续下降的位点
print(dim(melanoma_down_down))
save(melanoma_down_down,file ='DMP cervical_down_down.Rdata' )

melanoma_up_up=intersect(rownames(skin_punch_biopsy_to_melanocytes_up),rownames(melanocytes_to_primary_melanoma_up))
save(melanoma_up_up,file='DMP cervical_up_up.Rdata')#


low_low_cancer=intersect(rownames(skin_punch_biopsy_to_melanocytes_down),rownames(melanocytes_to_primary_melanoma_down))
high_high_cancer=intersect(rownames(skin_punch_biopsy_to_melanocytes_up),rownames(melanocytes_to_primary_melanoma_up))
save(low_low_cancer,high_high_cancer,file='0122癌前和癌.R')

down_down_gene=as.data.frame(table(droplevels((melanocytes_to_primary_melanoma_down[melanoma_down_down,])$gene)))
up_up_gene=as.data.frame(table(droplevels((melanocytes_to_primary_melanoma_up[melanoma_up_up,])$gene)))

melanoma_up_down=intersect(rownames(skin_punch_biopsy_to_melanocytes_up),rownames(melanocytes_to_primary_melanoma_down))
melanoma_up_down_gene=as.data.frame(table(droplevels((skin_punch_biopsy_to_melanocytes_up[melanoma_up_down,])$gene)))
melanoma_down_up=intersect(rownames(skin_punch_biopsy_to_melanocytes_down),rownames(melanocytes_to_primary_melanoma_up))
melanoma_down_up_gene=as.data.frame(table(droplevels((skin_punch_biopsy_to_melanocytes_down[melanoma_down_up,])$gene)))
save(melanoma_up_down,melanoma_up_down_gene,melanoma_down_up,melanoma_down_up_gene,file = '0207melanoma_downandup.Rdata')



down_down_gene=down_down_gene$Var1
up_up_gene=up_up_gene$Var1
list=c('melanoma_down_down','melanoma_up_up','down_down_gene','up_up_gene')
library(readr)
for (i in 1:length(list)){
  nn=list[i]
  getnn=data.frame(get(nn))
  colnames(getnn)=nn
  write_csv(as.data.frame(nn),append = T,file =paste('0206_melanoma',nn,'.csv'))
  write_csv(getnn,append = T,file =paste('0206_melanoma',nn,'.csv'))
}

stable_low_cancer=intersect(rownames(melanocytes_to_primary_melanoma_down),rownames(skin_punch_biopsy_to_melanocytes_stable))
stable_high_cancer=intersect(rownames(melanocytes_to_primary_melanoma_up),rownames(skin_punch_biopsy_to_melanocytes_stable))
save(stable_high_cancer,stable_low_cancer,file='0122癌突变_melanoma_0.2.Rdata')
