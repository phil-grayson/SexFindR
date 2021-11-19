library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("~/SexFindR/Supplemental_Code/Mosquito/DifCover/")

##### --- plots --- #####

together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases spanned" = stop-base) %>% filter(grepl("NC_03510",scaf)) 

scaffold_lengths <- read_tsv("GCF_002204515.2_AaegL5.0_genomic.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("NC_03510",scaf))

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(scaf=="NC_035107.1", "Y-like", "Autosome")) %>% mutate(proportion = `bases spanned`/length)

x <- ggscatter(proportion, 
               x = "log2(Male coverage/Female coverage)", 
               y = "proportion", 
               color = "Chromosome", 
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
               ylab = "proportion of chromosome",
               size = 2) 
x

filtered_proportion <- proportion %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656 | `log2(Male coverage/Female coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

filtered_proportion %>% ggplot(aes(x=scaf,y=`total chromosome proportion with significantly different coverage`)) + geom_point(aes(color=Chromosome)) + theme(axis.text.x = element_text(angle = 90))

y <- ggdotchart(filtered_proportion, 
                x = "scaf", 
                y = "total chromosome proportion with significantly different coverage", 
                color = "Chromosome", 
                palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
                xlab = "Chromosome",
                #sorting = "descending", 
                add = "segments", 
                add.params = list(color = "lightgray", size = 1), 
                group = "Chromosome", 
                dot.size = 4 ) 
x+y

