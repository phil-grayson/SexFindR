library(tidyverse)
library(ggpubr)

# This script contains the information required to re-create the data that was provided in Table 1 of the Genomic_Basis_draft_dec2020.docx

setwd("~/Desktop/DifCover_Dec2020_Paper/")

##### --- Fugu --- #####
# M98 and F99 were run as male over female
# positive values are male enriched
# negative values are female enriched

# want to see if we can load in just the main file and filter or if we should go from the two filtered final files
fugu_up <- read_tsv(file = "fugu/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout.up0.7369656",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
fugu_down <- read_tsv(file = "fugu/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout.down-0.7369656",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
fugu_all <- read_tsv(file = "fugu/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
fugu_all_up <- fugu_all %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656)
fugu_all_down <- fugu_all %>%  filter(`log2(Male coverage/Female coverage)` <= -0.7369656)

all.equal(fugu_down,fugu_all_down)
all.equal(fugu_up,fugu_all_up)
# Given that these are both true, we can go ahead with just filtereing from the .DNAcopyout file as we have done with fugu_all for the rest of these species.

# number of male and female enriched regions is simply the number of rows in the DF.  for the number of bases, we need to do a quick sum

sum(fugu_all_down$bases_spanned)
sum(fugu_all_up$bases_spanned)
sum(fugu_all$bases_spanned)
# look for the enrichment at the region (non-existant)
fugu_all %>% filter(scaf == "NC_042303.1") %>% filter(base < 12708100,stop > 12708100)

##### --- Mosquito --- #####
# Male99 and Female96 were run with male over female
# up is male enriched. down is female enriched
mosquito_all <- read_tsv(file = "mosquito/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
mosquito_all_up <- mosquito_all %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656)
mosquito_all_down <- mosquito_all %>%  filter(`log2(Male coverage/Female coverage)` <= -0.7369656)

sum(mosquito_all$bases_spanned)
sum(mosquito_all_up$bases_spanned)
sum(mosquito_all_down$bases_spanned)

mosquito_all %>% filter(scaf == "NC_035107.1") %>% filter(base > 151680000, stop < 152960000)

##### --- Cannabis --- #####
# F23 and M15 were analyzed
# up is female enriched, down is male enriched

cannabis_all <- read_tsv(file = "cannabis/sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_0.558967673.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
cannabis_all_up <- cannabis_all %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656)
cannabis_all_down <- cannabis_all %>%  filter(`log2(Male coverage/Female coverage)` <= -0.7369656)

sum(cannabis_all$bases_spanned)
sum(cannabis_all_up$bases_spanned)
sum(cannabis_all_down$bases_spanned)

##### --- Rumex TX (XY) --- #####
# currently re-running
# up is female enrichment, down is male enrichment (column names wrong)
rumexTX_all <- read_tsv(file = "rumex_TX/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.2.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
rumexTX_all_up <- rumexTX_all %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656)
rumexTX_all_down <- rumexTX_all %>%  filter(`log2(Male coverage/Female coverage)` <= -0.7369656)

sum(rumexTX_all$bases_spanned)
sum(rumexTX_all_up$bases_spanned)
sum(rumexTX_all_down$bases_spanned)
##### --- Rumex NC (XYY) --- #####
# ran F20 over M21
# up is female enrichment, down is male enrichment (column names wrong)
rumexNC_all <- read_tsv(file = "rumex_NC/sample1_sample2.ratio_per_w_CC0_a10_A219_b10_B240_v1000_l500.log2adj_1.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Male coverage/Female coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
rumexNC_all_up <- rumexNC_all %>% filter(`log2(Male coverage/Female coverage)` >= 0.7369656)
rumexNC_all_down <- rumexNC_all %>%  filter(`log2(Male coverage/Female coverage)` <= -0.7369656)

sum(rumexNC_all$bases_spanned)
sum(rumexNC_all_up$bases_spanned)
sum(rumexNC_all_down$bases_spanned)

##### --- Chicken --- #####
# ran F634 over M635
# up is female enrichment, down is male enrichment
#
chicken_all <- read_tsv(file = "chicken/sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_1.222222222.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
chicken_all_up <- chicken_all %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656)
chicken_all_down <- chicken_all %>%  filter(`log2(Female coverage/Male coverage)` <= -0.7369656)

sum(chicken_all$bases_spanned)
sum(chicken_all_up$bases_spanned)
sum(chicken_all_down$bases_spanned)

chicken_all %>% filter(scaf == "NC_006126.5") #W
chicken_all %>% filter(scaf == "NC_006127.5") #Z
##### --- Lamprey --- #####
# F2 over M6
# up is female enrichment, down is male enrichment

lamprey_all <- read_tsv(file = "lamprey/sample1_sample2.ratio_per_w_CC0_a10_A500_b10_B500_v1000_l500.log2adj_1.080840086.DNAcopyout",col_names = F) %>% rename(scaf=X1,base=X2,stop=X3,windows=X4,"log2(Female coverage/Male coverage)"=X5) %>% mutate("bases_spanned" = stop-base)
lamprey_all_up <- lamprey_all %>% filter(`log2(Female coverage/Male coverage)` >= 0.7369656)
lamprey_all_down <- lamprey_all %>%  filter(`log2(Female coverage/Male coverage)` <= -0.7369656)

sum(lamprey_all$bases_spanned)
sum(lamprey_all_up$bases_spanned)
sum(lamprey_all_down$bases_spanned)


# Also, loading in the basic data from the 2/2 DifCover for chicken and lamprey
##### --- Chicken 2 Females versus 2 Males --- #####
# female over male means down is male enriched and up is female enriched
# annotation file:
chicken_annotate <- read_tsv(file="chicken/sorted_chick_10kb_windowed_difcover.txt")

chicken_annotate_controlled_down <- chicken_annotate %>% filter(F634M635_down >0, F634M637_down > 0, F634F707_down == 0, M635M637_down==0, F707M635_down > 0,F707M637_down >0, F634F707_up==0, M635M637_up==0)
table(chicken_annotate_controlled_down$scaf)
length(table(chicken_annotate_controlled_down$scaf))

chicken_annotate_controlled_up <- chicken_annotate %>% filter(F634M635_up >0, F634M637_up > 0, F634F707_up == 0, M635M637_up==0, F707M635_up > 0,F707M637_up >0, F634F707_down==0, M635M637_down==0)
table(chicken_annotate_controlled_up$scaf)
length(table(chicken_annotate_controlled_up$scaf))

##### --- Lamprey 2 Females versus 2 Males --- #####
# female over male, down is male enriched, up is female enriched
lamprey_annotate <- read_tsv(file="lamprey/sorted_vgp_10kb_windowed_difcover.txt")
lamprey_annotate_controlled_down <- lamprey_annotate %>% filter(F10M6_down >0,F10M8_down>0,F2F10_down==0,F2M6_down>0,F2M8_down>0,M6M8_down==0,F2F10_up==0,M6M8_up==0)
table(lamprey_annotate_controlled_down$scaf)
length(table(lamprey_annotate_controlled_down$scaf))

lamprey_annotate_controlled_up <- lamprey_annotate %>% filter(F10M6_up >0,F10M8_up>0,F2F10_up==0,F2M6_up>0,F2M8_up>0,M6M8_up==0,F2F10_down==0,M6M8_down==0)
table(lamprey_annotate_controlled_up$scaf)
length(table(lamprey_annotate_controlled_up$scaf))

# low coverage 
lamprey_annotate_low <- read_tsv(file="lamprey/vgp_10kb_windowed_difcover_lowcov_merged.txt") 

lamprey_annotate_together <- left_join(lamprey_annotate,lamprey_annotate_low)
lamprey_annotate_together_controlled_down <- lamprey_annotate_together %>% filter(F10M6_down >0,F10M8_down>0,F2F10_down==0,F2M6_down>0,F2M8_down>0,M6M8_down==0,F2F10_up==0,M6M8_up==0, down>0,up==0)
table(lamprey_annotate_together_controlled_down$scaf)
length(table(lamprey_annotate_together_controlled_down$scaf))

lamprey_annotate_together_controlled_up <- lamprey_annotate_together %>% filter(F10M6_up >0,F10M8_up>0,F2F10_up==0,F2M6_up>0,F2M8_up>0,M6M8_up==0,F2F10_down==0,M6M8_down==0,down==0,up>0)
table(lamprey_annotate_together_controlled_up$scaf)
length(table(lamprey_annotate_together_controlled_up$scaf))