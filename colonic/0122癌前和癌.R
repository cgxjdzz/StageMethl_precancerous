LGA_to_Normal=myDMP$Normal_to_HGA
Normal_to_HGA=myDMP$Normal_to_Tumor
library(dplyr)

deg=LGA_to_Normal
deg$g=ifelse(abs(deg$logFC) < 0.2,'stable',
             ifelse(deg$logFC > 0.2,'UP','DOWN'))
Normal_to_LGA_stable=filter(deg,g=='stable')
Normal_to_LGA_down=filter(deg,g=='UP')
Normal_to_LGA_up=filter(deg,g=='DOWN')

deg=Normal_to_HGA
deg$g=ifelse(abs(deg$logFC) < 0.2,'stable',
             ifelse(deg$logFC > 0.2,'UP','DOWN'))
Normal_to_HGA_down=filter(deg,g=='DOWN')
Normal_to_HGA_up=filter(deg,g=='UP')


low_low_cancer=intersect(rownames(Normal_to_HGA_down),rownames(Normal_to_LGA_down))
high_high_cancer=intersect(rownames(Normal_to_HGA_up),rownames(Normal_to_LGA_up))
save(low_low_cancer,high_high_cancer,file='0122癌前和癌_colonic_0.2.Rdata')

stable_low_cancer=intersect(rownames(Normal_to_HGA_down),rownames(Normal_to_LGA_stable))
stable_high_cancer=intersect(rownames(Normal_to_HGA_up),rownames(Normal_to_LGA_stable))
save(stable_high_cancer,stable_low_cancer,file='0122stable突变_colonic_0.2.Rdata')

