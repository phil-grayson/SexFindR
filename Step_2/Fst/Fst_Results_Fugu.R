library(tidyverse)

##### --- Fugu --- #####
setwd("~/SexFindR/Step_2/Fst")
Fstugu <- read_tsv("biallelic_fst_noMissing.weir.fst") %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0)) %>% rename(scaf = CHROM, base=POS) %>% mutate(base=as.numeric(base)) %>% filter(WEIR_AND_COCKERHAM_FST >=0) %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0))
scaffold_lengths <- read_tsv("GCF_901000725.2_fTakRub1.2_genomic.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1") 

cutoff1 <- round(nrow(Fstugu)*0.05)
Fstugu_sort_cut1 <- head(Fstugu %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff1) 
min(Fstugu_sort_cut1$`WEIR_AND_COCKERHAM_FST`)

# lowest value there is 0.157874 (5% cutoff)

cutoff2 <- round(nrow(Fstugu)*0.01)
Fstugu_sort_cut2 <- head(Fstugu %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff2) 
min(Fstugu_sort_cut2$`WEIR_AND_COCKERHAM_FST`)
# lowest value there is 0.247786 (1% cutoff)

table(Fstugu_sort_cut2$scaf)

quantile(Fstugu$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)

mean(Fstugu$WEIR_AND_COCKERHAM_FST)

FstFugu <- left_join(scaffold_lengths,Fstugu) %>% filter(WEIR_AND_COCKERHAM_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(FstFugu$scaf))
FstFugu$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(FstFugu$scaf)){
  nbp[i] <- max(FstFugu[FstFugu$scaf == i,]$base)
  FstFugu[FstFugu$scaf == i,"BPcum"] <- FstFugu[FstFugu$scaf == i,"base"] + s
  s <- s + nbp[i]
}

axis.set <- FstFugu %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

FstFugu %>% ggplot(aes(x=BPcum,y=WEIR_AND_COCKERHAM_FST,color=as.factor(scaf))) +geom_point(alpha = 0.75) +scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + theme(axis.text.x = element_text(angle = 90)) + labs(x="CHROMOSOME",color="") + guides(color = guide_legend(ncol = 1, byrow = F)) + geom_hline(yintercept = 0.157874, linetype="dotted", color = "black", size=0.75) + geom_hline(yintercept = 0.247786, linetype="dotted", color = "black", size=1.5)  + geom_vline(xintercept = 308269038,linetype="dotted",color="red",size=1) + labs(title="Fst for fugu with 5% and 1% cutoff and the FuguSDR as dotted lines")
ggsave("fugu_Fst_manhattan_with95_99_noMissing.pdf")