colonic_up_up_gene=as.data.frame(table(droplevels((HGA_2_Tumor_up[up_up,])$gene)))
save(colonic_up_up_gene,file ='DMP colonic_up_up_gene.Rdata')
colonic_down_down_gene=as.data.frame(table(droplevels((HGA_2_Tumor_down[down_down,])$gene)))
save(colonic_down_down_gene,file ='DMP colonic_down_down_gene.Rdata')
