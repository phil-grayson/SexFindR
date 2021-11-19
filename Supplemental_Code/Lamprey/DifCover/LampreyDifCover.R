library(tidyverse)
library(patchwork)
library(ggpubr)

setwd("~/SexFindR/Supplemental_Code/Lamprey/DifCover/")
scaffold_lengths <- read_tsv("vgp_scaffold_sizes.txt", col_names = c("scaf","length"))

##### --- F2_M6 --- #####

together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_1.080840086.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases spanned" = stop-base) 

scaffold_lengths <- read_tsv("/Users/phil/Desktop/Apr_May_Lamprey/vgp_scaffold_sizes.txt", col_names = c("scaf","length"))

proportion_full <- full_join(together%>% filter(grepl("NC_0",scaf)),scaffold_lengths) %>% mutate(proportion_full = `bases spanned`/length) 

x <- ggscatter(proportion_full, 
               x = "log2(Female coverage/Male coverage)", 
               y = "proportion_full", 
               color = "#9BA4A9", 
               #palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
               ylab = "proportion of chromosome",
               size = 2)  + geom_vline(xintercept =0,linetype="dotted")
x

filtered_proportion_full <- proportion_full %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656 | `log2(Female coverage/Male coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

y <- ggdotchart(filtered_proportion_full, 
                x = "scaf", 
                y = "total chromosome proportion with significantly different coverage", 
                color = "#9BA4A9", 
                #palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
                xlab = "Chromosome",
                #sorting = "descending", 
                add = "segments", 
                add.params = list(color = "lightgray", size = 1), 
                #group = "Chromosome", 
                dot.size = 4 ) + theme(axis.text.x=element_blank())
y

x+y

