source('zffdatamining_gene.R')
library(ggplot2)
load('0107_sum_impute.Rdata')
pd=data.frame(sum_impute$pd)
beta=sum_impute$beta
load('sum_dmp.Rdata')


normal_dmp=intersect(rownames(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal),rownames(sum_DMP$Primary_Tumor_to_Solid_Tissue_Normal))
IMB_dmp=intersect(rownames(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal),rownames(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor))
cancer_dmp=intersect(rownames(sum_DMP$Primary_Tumor_to_Solid_Tissue_Normal),rownames(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Primary_Tumor))

beta_normal=beta[normal_dmp,which(pd$label=='Solid Tissue Normal')]
IMB_dmp=IMB_dmp[which(IMB_dmp %in% rownames(beta))]
beta_IMB=beta[IMB_dmp,which(pd$label=='intestinal metaplasia biopsy from gastric antrum')]
beta_cancer=beta[cancer_dmp,which(pd$label=='Primary Tumor')]

normal_scorelist=ge.disordercalculate(beta_normal)
normal_disorder_scorelist=normal_scorelist[which(normal_scorelist >= 4)]
normal_disorder_dmp=normal_dmp[which(normal_scorelist >= 4)]
normal_final=ge.disorderfinal(normal_disorder_scorelist,length(normal_dmp))
pdf("gastic_normal_scoreplot.pdf")
myplot <- ggplot(data.frame(y=normal_scorelist,x=c(1:length(normal_scorelist))),aes(x,y)) + geom_point(alpha = 0.3)
print(myplot)
dev.off()

IMB_scorelist=ge.disordercalculate(beta_IMB)
IMB_disorder_scorelist=IMB_scorelist[which(IMB_scorelist >= 4)]
IMB_disorder_dmp=IMB_dmp[which(IMB_scorelist >= 4)]
IMB_final=ge.disorderfinal(IMB_disorder_scorelist,length(IMB_dmp))
pdf("gastic_IMB_scoreplot.pdf")
myplot <- ggplot(data.frame(y=IMB_scorelist,x=c(1:length(IMB_scorelist))),aes(x,y)) + geom_point(alpha = 0.3)
print(myplot)
dev.off()

cancer_scorelist=ge.disordercalculate(beta_cancer)
cancer_disorder_scorelist=cancer_scorelist[which(cancer_scorelist >= 4)]
cancer_disorder_dmp=cancer_dmp[which(cancer_scorelist >= 4)]
cancer_final=ge.disorderfinal(cancer_disorder_scorelist,length(cancer_dmp))
pdf("gastic_cancer_scoreplot.pdf")
myplot <- ggplot(data.frame(y=cancer_scorelist,x=c(1:length(cancer_scorelist))),aes(x,y)) + geom_point(alpha = 0.3)
print(myplot)
dev.off()

save(normal_disorder_dmp,IMB_disorder_dmp,cancer_disorder_dmp,file='gastic_disorder_dmp.Rdata')
c(normal_final,
  IMB_final,
  cancer_final)

normal_disorder_gene=(droplevels(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal[normal_disorder_dmp,]$gene))
CAF_disorder_gene=(droplevels(sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal[IMB_disorder_dmp,]$gene))
cancer_disorder_gene=(droplevels(sum_DMP$Primary_Tumor_to_Solid_Tissue_Normal[cancer_disorder_dmp,]$gene))
gastric_disorder_gene=table(c(normal_disorder_gene,CAF_disorder_gene,cancer_disorder_gene))
save(gastric_disorder_gene,file='gastric_disorder_gene.Rdata')
