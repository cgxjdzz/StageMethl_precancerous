load('1217_sum_norm.Rdata')
library(ChAMP)
pd$label=gsub(' ','_',pd$label)
sum_DMP <- champ.DMP(beta = sum_norm,pheno=pd$label)
save(sum_DMP,file='sum_dmp.Rdata')


intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal=sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal,file='intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal.Rdata')
intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor=sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor,file='intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor.Rdata')


library(dplyr)
deg=intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down=filter(deg,g=='UP')
Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up=filter(deg,g=='DOWN')
save(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down,file = 'DMP Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down.Rdata')
save(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up,file = 'DMP Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up.Rdata')

deg=intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor
deg$g=ifelse(abs(deg$logFC) < 0.02,'stable',
             ifelse(deg$logFC > 0.02,'UP','DOWN'))
intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down=filter(deg,g=='DOWN')
intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up=filter(deg,g=='UP')
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down,file = 'DMP intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down.Rdata')
save(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up,file = 'DMP intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up.Rdata')

gastric_down_down=intersect(rownames(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down),rownames(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down))
save(gastric_down_down,file ='DMP gastric_down_down.Rdata' )
gastric_down_down_gene=as.data.frame(table(droplevels((Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down[gastric_down_down,])$gene)))
save(gastric_down_down_gene,file ='DMP gastric_down_down_gene.Rdata' )

gastric_up_up=intersect(rownames(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up),rownames(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up))
save(gastric_up_up,file ='DMP gastric_up_up.Rdata' )
gastric_up_up_gene=as.data.frame(table(droplevels((Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up[gastric_up_up,])$gene)))
save(gastric_up_up_gene,file ='DMP gastric_up_up_gene.Rdata')

gastric_up_down=intersect(rownames(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up),rownames(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_down))
gastric_up_down_gene=as.data.frame(table(droplevels((Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_up[gastric_up_down,])$gene)))
gastric_down_up=intersect(rownames(Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down),rownames(intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor_up))
gastric_down_up_gene=as.data.frame(table(droplevels((Solid_Tissue_Normal_to_intestinal_metaplasia_biopsy_from_gastric_antrum_down[gastric_down_up,])$gene)))
save(gastric_up_down,gastric_up_down_gene,gastric_down_up,gastric_down_up_gene,file = '0207gastric_downandup.Rdata')




gastric_down_down_gene=gastric_down_down_gene$Var1
gastric_up_up_gene=gastric_up_up_gene$Var1
list=c('gastric_down_down','gastric_up_up','gastric_down_down_gene','gastric_up_up_gene')
library(readr)
for (i in 1:length(list)){
  nn=list[i]
  getnn=data.frame(get(nn))
  colnames(getnn)=nn
  write_csv(as.data.frame(nn),append = T,file =paste('0206_gastric',nn,'.csv'))
  write_csv(getnn,append = T,file =paste('0206_gastric',nn,'.csv'))
}
