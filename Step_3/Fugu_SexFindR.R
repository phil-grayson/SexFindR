library(tidyverse)
setwd("~/Desktop/fugu/")
# Focus on 22 fully assembled chromosomes for the test

# SNP density with permutations, remove any that don't have p <= 0.05. sort the windows by p-value first and then by absolute value of the change
# then provide a rank order 
fugu_SNP <- read_tsv("~/Desktop/fugu/snpdensity/SNPdensity_perm_with_true_p_fugu.txt") %>% select(scaf,base,mean_MvF_dif,Pvalue)
fugu_SNP_filter <- fugu_SNP %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1") %>% filter(Pvalue <= 0.05) %>% mutate(base=base+10000) # adding window size for consistency between analyses
fugu_SNP_filter_arrange <- fugu_SNP_filter %>% arrange(Pvalue,-abs(mean_MvF_dif)) %>% rownames_to_column() %>% rename(SNPdensity_rank=rowname) %>% mutate(SNPdensity_rank=as.numeric(SNPdensity_rank))

# GWAS. grab top 5% of sites and write them out for downstream python
fugu_gemma <- read_tsv("~/Desktop/GEMMA/Fugu/test_gemma_out.assoc.txt") %>% separate(rs, c("scaf","base"),sep=':')%>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")
cutoff <- round(nrow(fugu_gemma)*0.05)
fugu_gemma_sort_cut <- head(fugu_gemma %>% arrange(p_lrt), n=cutoff) 
#write_tsv(path="fugu_gemma_sort_cut.txt",fugu_gemma_sort_cut %>% select(scaf,base),col_names = F)
gemma_window_count_fugu <- read_tsv("gemma_window_count.txt") %>% arrange(-count) %>% rownames_to_column() %>% rename(GWAScount_rank=rowname) %>% mutate(GWAScount_rank = as.numeric(GWAScount_rank))

# combine SNP density and GWAS 
SNP_gemma_fugu <- full_join(fugu_SNP_filter_arrange,gemma_window_count_fugu) %>% select(scaf,base,GWAScount_rank,SNPdensity_rank)
top_SNP_gemma_fugu <- SNP_gemma_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100)

# now we want Fst
Fstugu <- read_tsv("biallelic_fst.weir.fst") %>% replace_na(list(WEIR_AND_COCKERHAM_FST=0)) %>% rename(scaf = CHROM, base=POS) %>% mutate(base=as.numeric(base)) %>% filter(WEIR_AND_COCKERHAM_FST >=0) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1")
cutoff <- round(nrow(Fstugu)*0.05)
Fstugu_sort_cut <- head(Fstugu %>% arrange(-WEIR_AND_COCKERHAM_FST), n=cutoff) 
#write_tsv(path="Fstugu_sort_cut.txt",Fstugu_sort_cut %>% select(scaf,base),col_names = F)

Fstugu_window_count <- read_tsv("Fstugu_window_count.txt") %>% arrange(-count) %>% rownames_to_column() %>% rename(Fstcount_rank=rowname) %>% mutate(Fstcount_rank = as.numeric(Fstcount_rank))

# combine SNP_gemma with Fstugu_window_count
all_SNP_based_fugu <- full_join(SNP_gemma_fugu,Fstugu_window_count %>% select(-count)) %>% mutate(total_rank = GWAScount_rank+SNPdensity_rank+Fstcount_rank)
top_SNP_based_fugu <- all_SNP_based_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100,Fstcount_rank <=100) %>% mutate(total_rank = GWAScount_rank+SNPdensity_rank+Fstcount_rank)
top_SNP_based_fugu %>% select(-total_rank) %>% write_tsv("Fugu_SexFindR_candidate_windows.txt")

# finally, the permuted depth
depth_fugu <- read_tsv("depth_p005_window_count.txt") %>% arrange(-count) %>% rownames_to_column() %>% rename(Depthcount_rank=rowname) %>% mutate(Depthcount_rank = as.numeric(Depthcount_rank)) %>% select(-count)

# all together now!
all_data_fugu <- full_join(all_SNP_based_fugu,depth_fugu)
top_all_data_fugu <- all_data_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100,Fstcount_rank <=100,Depthcount_rank<=10000)

top_all_data_fugu <- all_data_fugu %>% filter(GWAScount_rank<=100,SNPdensity_rank<=100,Fstcount_rank <=100,Depthcount_rank<=10000)

# some correlation tests:
# so top_all_data has nothing of note out to 10k....let's look at correlations:
cor.test(all_data_fugu$GWAScount_rank,all_data_fugu$SNPdensity_rank)
#all_data %>% ggplot(aes(x=GWAScount_rank,y=SNPdensity_rank)) + geom_point()
cor.test(all_data_fugu$GWAScount_rank,all_data_fugu$Fstcount_rank)
all_data_fugu %>% ggplot(aes(x=GWAScount_rank,y=Fstcount_rank)) + geom_point()
cor.test(all_data_fugu$SNPdensity_rank,all_data_fugu$Fstcount_rank)
#all_data %>% ggplot(aes(x=SNPdensity_rank,y=Fstcount_rank)) + geom_point()
# not correlated
cor.test(all_data_fugu$Depthcount_rank,all_data_fugu$Fstcount_rank)
cor.test(all_data_fugu$Depthcount_rank,all_data_fugu$SNPdensity_rank)
cor.test(all_data_fugu$Depthcount_rank,all_data_fugu$GWAScount_rank)


# want to make a plot for the top region

##### --- Fugu only NC_042303.1 --- #####
scaffold_lengths <- read_tsv("~/Desktop/fugu/GCF_901000725.2_fTakRub1.2_genomic.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("NC_",scaf)) %>% filter(scaf != "NC_004299.1") 
# 95% and 99% outliers
quantile(Fstugu$`WEIR_AND_COCKERHAM_FST`, c(0.95, 0.99), na.rm = T)

test <- left_join(scaffold_lengths,Fstugu) %>% filter(WEIR_AND_COCKERHAM_FST>=0) %>% filter(scaf == "NC_042303.1")

a <- test %>% ggplot(aes(x=base, y=WEIR_AND_COCKERHAM_FST)) +geom_point(alpha = 0.75,color="dodgerblue") + theme(axis.text.x = element_text(angle = 90)) + labs(x="",color="") + geom_hline(yintercept = 0.158070, linetype="dotted", color = "black", size=0.75) + geom_hline(yintercept = 0.250575, linetype="dotted", color = "black", size=1.5)  + geom_vline(xintercept = 12708100,linetype="dotted",color="red",size=1) + labs(title = "A. Fst")+ geom_vline(xintercept = 9270000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 2990000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 9290000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 10600000,linetype="dotted",color="coral4",size=1)
#+ labs(title="Fst for fugu with 5% and 1% cutoff and the FuguSDR as dotted lines")


b <- fugu_SNP%>% filter(scaf == "NC_042303.1") %>% ggplot(aes(x=base,y=mean_MvF_dif)) +geom_point(alpha = 0.75,color="mediumpurple")+ geom_vline(xintercept = 12708100,linetype="dotted",color="red",size=1)  + labs(x="",color="") + labs(y="male mean - female mean",color="")+ labs(title = "B. SNP Density")+ geom_vline(xintercept = 9270000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 2990000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 9290000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 10600000,linetype="dotted",color="coral4",size=1)

c <- gemma_window_count_fugu %>% filter(scaf == "NC_042303.1") %>% ggplot(aes(x=base,y=count)) +geom_point(alpha = 0.75,color="orchid2")+ geom_vline(xintercept = 12708100,linetype="dotted",color="red",size=1)  + labs(x="bp",color="") + labs(y="number of GWAS hits in top 5%",color="") + labs(title = "C. Gemma GWAS") + geom_vline(xintercept = 9270000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 2990000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 9290000,linetype="dotted",color="coral4",size=1) + geom_vline(xintercept = 10600000,linetype="dotted",color="coral4",size=1)

a/b/c
ggsave("NC_042303_fugu_sexFindR_results.pdf")
