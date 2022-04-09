library(ggplot2)
cancerlist=c('prostatic','cervical','liver','gastric','colonic')
sitelist=c('hypo_hypo','hyper_hyper','hypo_stable','hyper_stable')
dir.create('0407venn')
for (i in sitelist){
  for (j in cancerlist){
    ee=read.table(paste0(i,'/',j,'_coef.txt'))
    assign(j,ee[which(ee[,3]<0.05),1])
  }
  library(VennDiagram)
  veen <-venn.diagram(list(prostatic=prostatic,cervical=cervical,liver=liver,gastric=gastric,colonic=colonic), fill=c("red","green","blue",'yellow','gray'), alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, filename=NULL)
  pdf(paste0('0407venn/',i,".pdf"))
  grid.draw(veen)
  dev.off()
}

