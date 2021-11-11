library(tidyverse)
library(patchwork)
library(ggthemes)
# Focus on 22 fully assembled chromosomes (contain "NC_")

# SNP density with permutations, remove any that don't have p <= 0.05. sort the windows by p-value first and then by absolute value of the change
# then provide a rank order 
fugu_SNP <- read_tsv("~/SexFindR/Step_3/SNPdensity_perm_with_true_p_fugu.txt.zip") %>% select(scaf,base,mean_MvF_dif,Pvalue)
fugu_SNP_filter <- fugu_SNP %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1") %>% filter(Pvalue <= 0.05) %>% mutate(base=base+10000) # adding window size for consistency between analyses
fugu_SNP_filter_arrange <- fugu_SNP_filter %>% arrange(Pvalue,-abs(mean_MvF_dif)) %>% rownames_to_column() %>% rename(SNPdensity_rank=rowname) %>% mutate(SNPdensity_rank=as.numeric(SNPdensity_rank))

# GWAS. grab top 5% of sites and write them out for downstream python
fugu_gemma <- read_tsv("~/SexFindR/Step_3/test_gemma_out.assoc.txt.zip") %>% separate(rs, c("scaf","base"),sep=':')%>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")
cutoff <- round(nrow(fugu_gemma)*0.05)
fugu_gemma_sort_cut <- head(fugu_gemma %>% arrange(p_lrt), n=cutoff) 
#write_tsv(path="fugu_gemma_sort_cut.txt",fugu_gemma_sort_cut %>% select(scaf,base),col_names = F)
gemma_window_count_fugu <- read_tsv("~/SexFindR/Step_3/gemma_window_count.txt") %>% arrange(-count) %>% rownames_to_column() %>% rename(GWAScount_rank=rowname) %>% mutate(GWAScount_rank = as.numeric(GWAScount_rank))

# combine SNP density and GWAS 
SNP_gemma_fugu <- full_join(fugu_SNP_filter_arrange,gemma_window_count_fugu) %>% select(scaf,base,GWAScount_rank,SNPdensity_rank)
top_SNP_gemma_fugu <- SNP_gemma_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100)

# now we want Fst
Fstugu <- read_tsv("~/SexFindR/Step_3/biallelic_fst.weir.fst") %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0)) %>% rename(scaf = CHROM, base=POS) %>% mutate(base=as.numeric(base)) %>% filter(WEIR_AND_COCKERHAM_FST >=0) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")
cutoff <- round(nrow(Fstugu)*0.05)
Fstugu_sort_cut <- head(Fstugu %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff) 
#write_tsv(path="Fstugu_sort_cut.txt",Fstugu_sort_cut %>% select(scaf,base),col_names = F)

Fstugu_window_count <- read_tsv("~/SexFindR/Step_3/Fstugu_window_count.txt") %>% arrange(-count) %>% rownames_to_column() %>% rename(Fstcount_rank=rowname) %>% mutate(Fstcount_rank = as.numeric(Fstcount_rank))

# combine SNP_gemma with Fstugu_window_count
all_SNP_based_fugu <- full_join(SNP_gemma_fugu,Fstugu_window_count %>% select(-count)) %>% mutate(total_rank = GWAScount_rank+SNPdensity_rank+Fstcount_rank)
top_SNP_based_fugu <- all_SNP_based_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100,Fstcount_rank <=100) %>% mutate(total_rank = GWAScount_rank+SNPdensity_rank+Fstcount_rank)
#top_SNP_based_fugu %>% select(-total_rank) %>% write_tsv("Fugu_SexFindR_candidate_windows.txt")

# want to make a plot for the top region

##### --- Fugu only NC_042303.1 --- #####
scaffold_lengths <- read_tsv("~/SexFindR/Step_3/GCF_901000725.2_fTakRub1.2_genomic.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1") 
# 95% and 99% outliers
#quantile(Fstugu$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)
#test <- left_join(scaffold_lengths,Fstugu) %>% filter(WEIR_AND_COCKERHAM_FST>=0) %>% filter(scaf == "NC_042303.1")
#color="dodgerblue" + geom_vline(xintercept = 12708100,linetype="dotted",color="red",size=1) 
#a <- test %>% ggplot(aes(x=base, y=WEIR_AND_COCKERHAM_FST, color=WEIR_AND_COCKERHAM_FST)) + geom_point() + theme(axis.text.x = element_text(angle = 90)) + labs(x="",color="") + labs(y="Fst (Weir and Cockerham)",color="") + geom_vline(xintercept = 12710000,linetype="dotted",color="red",size=1) + labs(title = "A. Fst")+ geom_vline(xintercept = 9270000,linetype="dotted",color="coral4",size=0.5) + geom_vline(xintercept = 2990000,linetype="dotted",color="coral4",size=0.5) + geom_vline(xintercept = 9290000,linetype="dotted",color="coral4",size=0.5) + geom_vline(xintercept = 10600000,linetype="dotted",color="coral4",size=0.5) + scale_x_continuous(limits = c(0,max(test$base)*1.05), expand = c(0, 0)) +  scale_y_continuous(limits = c(0,max(test$WEIR_AND_COCKERHAM_FST)*1.05), expand = c(0, 0)) + theme_igray() + scale_color_gradient(low = "#e5e5e5",high = "#0072B2")
#a

x_length = scaffold_lengths %>% filter(scaf=="NC_042303.1")

fugu_fst_count <- left_join(scaffold_lengths,Fstugu_window_count) %>% filter(scaf == "NC_042303.1")

a <- fugu_fst_count %>% ggplot(aes(x=base, y=count, color=count)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = 90)) + labs(x="",color="") + labs(y="Fst (Weir and Cockerham) outliers per 10k window",color="") + geom_vline(xintercept = 12710000,linetype="dotted",color="red",size=1) + labs(title = "A. Fst")+ geom_vline(xintercept = 9270000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 2990000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 9290000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 10600000,linetype="dotted",color="black",size=0.5) +
  scale_x_continuous(n.breaks = 10,labels = scales::comma, limits = c(0,x_length$length), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,max(fugu_fst_count$count)*1.05), expand = c(0, 0)) + theme_igray() + scale_color_gradient(low = "#e5e5e5",high = "#0072B2") +theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title = element_text(size=30))


mid=0
b <- fugu_SNP%>% filter(scaf == "NC_042303.1") %>% ggplot(aes(x=base,y=mean_MvF_dif,color=mean_MvF_dif)) +geom_point(size=6)+ geom_vline(xintercept = 12710000,linetype="dotted",color="red",size=1)  + labs(x="",color="") + labs(y="10kb SNP Density (Male mean - Female mean)",color="")+ labs(title = "B. SNP Density")+ geom_vline(xintercept = 9270000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 2990000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 9290000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 10600000,linetype="dotted",color="black",size=0.5)+ theme_igray() + scale_x_continuous(n.breaks = 10,labels = scales::comma, limits = c(0,x_length$length), expand = c(0, 0)) + scale_color_gradient2(midpoint=mid, low = "#332288", mid = "#e5e5e5",high = "#332288") +theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title = element_text(size=30))

c <- gemma_window_count_fugu %>% filter(scaf == "NC_042303.1") %>% ggplot(aes(x=base,y=count,color=count)) +geom_point(size=6)+ geom_vline(xintercept = 12710000,linetype="dotted",color="red",size=1)  + labs(x="Position (bp) on fugu NC_042303.1",color="") + labs(y="GWAS outliers per 10kb window",color="") + labs(title = "C. Gemma GWAS") + geom_vline(xintercept = 9270000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 2990000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 9290000,linetype="dotted",color="black",size=0.5) + geom_vline(xintercept = 10600000,linetype="dotted",color="black",size=0.5)+ theme_igray() +scale_x_continuous(n.breaks = 10,labels = scales::comma, limits = c(0,x_length$length), expand = c(0, 0))+ scale_color_gradient(low = "#e5e5e5",high = "#882255") +theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),  
        axis.title.x = element_text(size = 30, vjust = -0.5, face = "plain"),
        axis.title.y = element_text(size = 25, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title = element_text(size=30))

a/b/c
ggsave("NC_042303_fugu_sexFindR_results_figure.pdf")
ggsave("NC_042303_fugu_sexFindR_results_figure.png")
