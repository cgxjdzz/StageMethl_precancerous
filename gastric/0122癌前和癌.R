
LGA_to_Normal=sum_DMP$Primary_Tumor_to_Solid_Tissue_Normal
Normal_to_HGA=sum_DMP$intestinal_metaplasia_biopsy_from_gastric_antrum_to_Solid_Tissue_Normal
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
Normal_to_HGA_down=filter(deg,g=='UP')
Normal_to_HGA_up=filter(deg,g=='DOWN')


low_low_cancer=intersect(rownames(Normal_to_HGA_down),rownames(Normal_to_LGA_down))
high_high_cancer=intersect(rownames(Normal_to_HGA_up),rownames(Normal_to_LGA_up))
save(low_low_cancer,high_high_cancer,file='0122癌前和癌_gastric_0.2.Rdata')

stable_low_cancer=intersect(rownames(Normal_to_HGA_down),rownames(Normal_to_LGA_stable))
stable_high_cancer=intersect(rownames(Normal_to_HGA_up),rownames(Normal_to_LGA_stable))
