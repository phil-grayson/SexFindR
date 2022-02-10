library(tidyverse)
library(patchwork)
library(ggpubr)


##### --- Original "mother plus Y" --- #####
setwd("~/SexFindR/Supplemental_Code/Cannabis/DifCover/")


together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_2.468996284.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases spanned" = stop-base) 

scaffold_lengths <- read_tsv("male_ref_scaffold_lengths.txt", col_names = c("scaf","length"))

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(grepl("Y_",scaf), "Y-region", "Autosome")) %>% mutate(proportion = `bases spanned`/length)

#want to flip axis for consistency between samples
proportion <- proportion %>% mutate(`log2(Male coverage/Female coverage)` = `log2(Female coverage/Male coverage)`*-1)

#font size (s)
s = 35

x <- ggscatter(proportion, 
               #x = "log2(Female coverage/Male coverage)", 
               x = "log2(Male coverage/Female coverage)",
               y = "proportion", 
               color = "Chromosome", 
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
               ylab = "Proportion of chromosome (in region)",
               size = 2) +
               geom_vline(xintercept =0,linetype="dotted") +
               font("legend.title", size = s) +
               font("legend.text", size = s) +
               font("xlab", size = s) +
               font("ylab", size = s) +
               font("xy.text", size = s)

filtered_proportion <- proportion %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656 | `log2(Female coverage/Male coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("Proportion of chromosome (significant)" = sum(`bases spanned`)/length)

y <- ggdotchart(filtered_proportion, 
                x = "scaf", 
                y = "Proportion of chromosome (significant)", 
                color = "Chromosome", 
                palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
                xlab = "Chromosome",
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

ggsave("cannabis_sex_chromosome.pdf", width = 23, height = 11.13349)
