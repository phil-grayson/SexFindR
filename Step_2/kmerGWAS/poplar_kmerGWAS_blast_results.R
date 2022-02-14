library(tidyverse)
library(patchwork)
setwd("~/SexFindR/Step_2/kmerGWAS/") #update path based on SexFindR folder location

# mapped to V2.2 to referenece back to Geraldes et al. 2015
poplar_kmerV2 <- read_tsv("poplar_abyss_v2.outfmt6",col_names = F)

# X1 is unique for each k-mer, so we want to select the best hit for each one and plot those:
poplar_kmerV2_filter <- poplar_kmerV2 %>% group_by(X1) %>% top_n(n = 1, wt = -X11)

a<-poplar_kmerV2_filter %>% ungroup() %>% count(X2) %>% rename(`Abyss contig blastn hits`=n,Chromosome=X2) %>% ggplot(aes(x=Chromosome,y=`Abyss contig blastn hits`)) + geom_col() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("A. Male-specific k-mers from kmerGWAS in poplar assembled with Abyss and linked to v2.2 of the reference genome using blastn")

b<-poplar_kmerV2_filter%>% filter(X2 == "scaffold_19") %>% ggplot(aes(x=X9,y=X3)) + geom_point() + xlab('Position on chromosome (bp)') + ylab('Percent identity for blastn') + scale_x_continuous(breaks=seq(0,16000000,1000000)) + ggtitle("B. Distribution of male-specific k-mers from kmerGWAS across poplar v2.2 chromosome 19 from blastn")

a/b
ggsave("kmersGWAS_poplar_blast.pdf", width = 11.5, height = 7)
