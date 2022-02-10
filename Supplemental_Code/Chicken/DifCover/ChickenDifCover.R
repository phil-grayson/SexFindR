library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("~/SexFindR/Supplemental_Code/Chicken/DifCover/")

together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_1.222222222.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases spanned" = stop-base) %>% filter(grepl("NC",scaf)) 

scaffold_lengths <- read_tsv("chick_scaffold_lengths.txt", col_names = c("scaf","length")) %>% filter(grepl("NC",scaf))

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(scaf=="NC_006127.5", "Z", ifelse(scaf=="NC_006126.5", "W", "Autosome"))) %>% mutate(proportion = `bases spanned`/length)

proportion <- proportion %>% mutate(`log2(Male coverage/Female coverage)` = `log2(Female coverage/Male coverage)`*-1)

s = 35

x <- ggscatter(proportion, 
               x = "log2(Male coverage/Female coverage)",
               y = "proportion",
               color = "Chromosome", 
               ylab = "Proportion of chromosome (in region)",
               palette = c("#9BA4A9", "#FFC107", "#029880"), 
               size = 2) +
                geom_vline(xintercept =0,linetype="dotted") +
                font("legend.title", size = s) +
                font("legend.text", size = s) +
                font("xlab", size = s) +
                font("ylab", size = s) +
                font("xy.text", size = s)


filtered_proportion <- proportion %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656 | `log2(Female coverage/Male coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

y <- ggdotchart(filtered_proportion, 
           x = "scaf", 
           y = "total chromosome proportion with significantly different coverage", 
           color = "Chromosome", 
           palette = c("#9BA4A9", "#FFC107", "#029880"), 
           xlab = "Chromosome",
           ylab = "Proportion of chromosome (significant)", 
           #sorting = "descending", 
           add = "segments", 
           add.params = list(color = "lightgray", size = 1), 
           group = "Chromosome", 
           dot.size = 4 )  +
            font("legend.title", size = s) +
            font("legend.text", size = s) +
            font("xlab", size = s) +
            font("ylab", size = s) +
            font("xy.text", size = s) + theme(axis.text.x=element_blank()) 

x+y

ggsave("chicken_sex_chromosome.pdf", width = 23, height = 11.13349)
