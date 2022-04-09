
library(clusterProfiler)
cancer_list=c('cervical','colonic','gastric','prostatic','melanoma','liver')
site_list=c('hyper_hyper','hyper_hypo','hyper_back_to_normal','hyper_stable','hypo_hypo','hypo_hyper','hypo_back_to_normal','hypo_stable')




# cervical_hyper_hypo_gene_bitr=c()
cervical_hyper_hyper_gene_bitr=c()
 cervical_hyper_stable_gene_bitr=c()
gastric_hyper_hyper_gene_bitr=c()
# colonic_hyper_back_to_normal_gene_bitr=c()
# gastric_hyper_back_to_normal_gene_bitr=c()
for (i in site_list){
  
  
  cp <- list(cervical=get(paste0('cervical_',i,'_gene_bitr')),colonic=get(paste0('colonic_',i,'_gene_bitr')),gastric=get(paste0('gastric_',i,'_gene_bitr')), prostatic=get(paste0('prostatic_',i,'_gene_bitr')),hepatic=get(paste0('liver_',i,'_gene_bitr')))
  go_compare <- compareCluster(cp, fun = "enrichGO", ont= 'BP', OrgDb="org.Hs.eg.db", pvalueCutoff=0.01, pAdjustMethod = "BH",qvalueCutoff = 0.01)
  go_result=go_compare@compareClusterResult
  pdf(paste0('0409_TSS_cut0.01_',i,'_go.pdf'), width = 10, height = 20)
  dotplot(go_compare, showCategory=10,font.size=6)
  dev.off()
  
  
  kegg_compare <- compareCluster(cp, fun = "enrichKEGG",pvalueCutoff=0.01,qvalueCutoff = 0.01)
  kegg_result=kegg_compare@compareClusterResult
  pdf(paste0('0409_TSS_cut0.01_',i,'_kegg.pdf'), width = 10, height = 10)
  dotplot(kegg_compare, showCategory=10,font.size=8)
  dev.off()
}
