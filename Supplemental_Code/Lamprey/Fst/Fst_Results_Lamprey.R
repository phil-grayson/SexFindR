library(tidyverse)

setwd("~/SexFindR/Supplemental_Code/Lamprey/")

##### --- Lamprey --- #####

lamprey_index <- read_tsv("DifCover/vgp_scaffold_sizes.txt",col_names = F) %>% rename(scaf=X1,length=X2)

lamprey_index_big <- lamprey_index %>% filter(grepl("NC_",scaf), length > 500000)

Fst <- read_delim("Fst/sex_biallelic_noMissing_Jan.weir.fst.zip",delim="\t",col_names = T) %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0))
mean(Fst$WEIR_AND_COCKERHAM_FST)

quantile(Fst$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)
#for 1% and 5% hlines in plot

#can use full genome 
FstLamprey_all <- left_join(lamprey_index,Fst%>%rename(scaf=CHROM,base=POS)) %>% filter(WEIR_AND_COCKERHAM_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(FstLamprey_all$scaf))
FstLamprey_all$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(FstLamprey_all$scaf)){
  nbp[i] <- max(FstLamprey_all[FstLamprey_all$scaf == i,]$base)
  FstLamprey_all[FstLamprey_all$scaf == i,"BPcum"] <- FstLamprey_all[FstLamprey_all$scaf == i,"base"] + s
  s <- s + nbp[i]
}

axis.set.all <- FstLamprey_all %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

#or just the large scaffolds
FstLamprey <- left_join(lamprey_index_big,Fst%>%rename(scaf=CHROM,base=POS)) %>% filter(WEIR_AND_COCKERHAM_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(FstLamprey$scaf))
FstLamprey$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(FstLamprey$scaf)){
  nbp[i] <- max(FstLamprey[FstLamprey$scaf == i,]$base)
  FstLamprey[FstLamprey$scaf == i,"BPcum"] <- FstLamprey[FstLamprey$scaf == i,"base"] + s
  s <- s + nbp[i]
}

axis.set <- FstLamprey %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

#all scaffolds
FstLamprey_all %>% ggplot(aes(x=BPcum,y=WEIR_AND_COCKERHAM_FST,color=as.factor(scaf))) +geom_point(alpha = 0.75) +scale_x_continuous(label = axis.set.all$scaf, breaks = axis.set.all$center) + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept = 0.008299406, linetype="dotted", color = "black", size=0.75) + geom_hline(yintercept = 0.017770800, linetype="dotted", color = "black", size=1.5) + labs(x="CHROMOSOME",color="") + theme(legend.position = "none")
#ggsave("lamprey_all_noMissing_Fst_manhattan_with95_99.pdf")

#big scaffolds
FstLamprey %>% ggplot(aes(x=BPcum,y=WEIR_AND_COCKERHAM_FST,color=as.factor(scaf))) +geom_point(alpha = 0.75) +scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept = 0.008299406, linetype="dotted", color = "black", size=0.75) + geom_hline(yintercept = 0.017770800, linetype="dotted", color = "black", size=1.5) + labs(x="CHROMOSOME",color="") + theme(legend.position = "none")
#ggsave("lamprey_all_noMissing_big_Fst_manhattan_with95_99.pdf")

#might want to put it into context since Fst is so low for lamprey (y axis 0-1 would be a closer comparison to fugu which goes above 0.8)
FstLamprey %>% ggplot(aes(x=BPcum,y=WEIR_AND_COCKERHAM_FST,color=as.factor(scaf))) +geom_point(alpha = 0.75) +scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept = 0.008299406, linetype="dotted", color = "black", size=0.75) + geom_hline(yintercept = 0.017770800, linetype="dotted", color = "black", size=1.5) + labs(x="CHROMOSOME",color="") +ylim(0,1) + theme(legend.position = "none")
#ggsave("lamprey_all_noMissing_0to1_big_Fst_manhattan_with95_99.pdf")

