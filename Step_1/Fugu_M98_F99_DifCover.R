library(tidyverse)
library(patchwork)
library(ggpubr)

together <- read_tsv(file = "~/SexFindR/Step_1/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases spanned" = stop-base) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")

scaffold_lengths <- read_tsv("~/SexFindR/Step_3/GCF_901000725.2_fTakRub1.2_genomic.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")

proportion <- full_join(together,scaffold_lengths) %>% mutate(Chromosome = ifelse(grepl("NC_042303.1",scaf), "XY-like", "Autosome")) %>% mutate(proportion = `bases spanned`/length)

s=35

x <- ggscatter(proportion,
               x = "log2(Male coverage/Female coverage)",
               y = "proportion",
               color = "Chromosome",
               palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"),
               ylab = "Proportion of chromosome (in region)",
               size = 2)  + geom_vline(xintercept =0,linetype="dotted")  +
                font("legend.title", size = s) +
                font("legend.text", size = s) +
                font("xlab", size = s) +
                font("ylab", size = s) +
                font("xy.text", size = s)


filtered_proportion <- proportion %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656 | `log2(Male coverage/Female coverage)` <= -0.7369656) %>% group_by(scaf) %>% mutate("total chromosome proportion with significantly different coverage" = sum(`bases spanned`)/length)

y <- ggdotchart(filtered_proportion,
                x = "scaf",
                y = "total chromosome proportion with significantly different coverage",
                color = "Chromosome",
                palette = c("#9BA4A9", "#34ADE8", "#8F3ED8"),
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
                font("xy.text", size = s) + 
                theme(axis.text.x=element_blank())

x+y


ggsave("fugu_sex_chromosome.pdf", width = 23, height = 11.13349)

