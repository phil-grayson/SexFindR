library(tidyverse)
library(patchwork)
#install.packages("ggpubr")
library(ggpubr)

setwd("~/Desktop/chicken_difCov_figures/")

together <- read_tsv(file = "sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_1.222222222.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases spanned" = stop-base) %>% filter(grepl("NC",scaf)) 

scaffold_lengths <- read_tsv("/Users/phil/Google\ Drive/Postdoc/Other_Datasets/Chickens/chick_scaffold_lengths.txt", col_names = c("scaf","length")) %>% filter(grepl("NC",scaf))

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(scaf=="NC_006127.5", "Z", ifelse(scaf=="NC_006126.5", "W", "Autosome"))) %>% mutate(proportion = `bases spanned`/length)

proportion %>% ggplot(aes(x=`log2(Female coverage/Male coverage)`,y=proportion)) +geom_point(aes(color=Chromosome))

x <- ggscatter(proportion, 
           x = "log2(Female coverage/Male coverage)", 
           y = "proportion", 
           color = "Chromosome", 
           palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
           ylab = "proportion of chromosome",
           size = "bases spanned") 
           #add = "segments", 
           #add.params = list(color = "lightgray", size = 1), 
           #group = "Chromosome", 

ggsave("chicken_634vs635_proportion_sized.pdf")


x <- ggscatter(proportion, 
               x = "log2(Female coverage/Male coverage)", 
               y = "proportion", 
               color = "Chromosome", 
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"), 
               ylab = "proportion of chromosome",
               size = 2)


ggsave("chicken_634vs635_proportion.pdf")

filtered_proportion <- proportion %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656 | `log2(Female coverage/Male coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

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
ggsave("chicken_634vs635_total_proportion_outside_0p7369656.pdf")

x+y
ggsave("chicken_634vs635_together.pdf")
