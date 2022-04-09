load('D:/甲基化/癌症发展/colonic/139404+68060/10.12myNorm.Rdata')
load('D:/甲基化/癌症发展/colonic/139404+68060/10.12impute.Rdata')
load('colonic_hy_and_back_to_normal.Rdata')
load('colonic_hy_and_stable.Rdata')
load('colonic_hyper_and_hypo.Rdata')
sum_norm=myNorm
pd0204=impute$pd
pd0204$label=mapvalues(impute$pd$label, c("LGA", "HGA"),c("Precancer", "Precancer"))
sum_DMP=champ.DMP(beta = myNorm,pheno = pd0204$label)
colonic_hyper_hyper_gene=c