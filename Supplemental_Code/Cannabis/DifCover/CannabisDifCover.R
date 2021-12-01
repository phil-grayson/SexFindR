library(tidyverse)
library(patchwork)
library(ggpubr)


##### --- Original "mother plus Y" --- #####
setwd("~/SexFindR/Supplemental_Code/Cannabis/DifCover/")


together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_2.468996284.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases spanned" = stop-base) 

scaffold_lengths <- read_tsv("/Users/phil/Google Drive/Postdoc/Other_Datasets/Cannabis/male_ref/male_ref_scaffold_lengths.txt", col_names = c("scaf","length"))

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(grepl("Y_",scaf), "Y-region", "Autosome")) %>% mutate(proportion = `bases spanned`/length)

x <- ggscatter(proportion, 
               x = "log2(Female coverage/Male coverage)", 
               y = "proportion", 
               color = "Chromosome", 
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
               ylab = "proportion of chromosome",
               size = 2)  + geom_vline(xintercept =0,linetype="dotted")

filtered_proportion <- proportion %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656 | `log2(Female coverage/Male coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

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
                dot.size = 4 ) + theme(axis.text.x=element_blank())

x+y
