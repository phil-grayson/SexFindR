library(tidyverse)
library(ggpubr)

# given the command run for SNP density, we have an easy way to isolate males and females

# load in the male data:
setwd("/Users/phil/Desktop/fugu/snpdensity/")
myFiles <- list.files(pattern="Male*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  #firstsplit <- strsplit(myFiles[i], "SNP_snpdensity_sorted_no_dups_sorted_trim_")
  firstsplit <- strsplit(myFiles[i], "_")
  column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}

SNPdensity.males <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

# load the female data:

myFiles <- list.files(pattern="Female*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  #firstsplit <- strsplit(myFiles[i], "SNP_snpdensity_sorted_no_dups_sorted_trim_")
  firstsplit <- strsplit(myFiles[i], "_")
  column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][2],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}


SNPdensity.females <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

#join them
SNPdensity <- full_join(SNPdensity.males,SNPdensity.females)

#make new tibble, change type for base, and replace NAs
SNPdensity.rows <- SNPdensity
SNPdensity.rows$base <- as.numeric(as.character(SNPdensity.rows$base))
SNPdensity.rows <- SNPdensity.rows %>% replace_na(list(Males = 0, Females = 0)) %>% replace(is.na(.), 0)

#write_tsv(SNPdensity.rows,"SNPdensity_rows_fugu.txt")
#SNPdensity.rows <- read_tsv("SNPdensity_rows_fugu.txt")

# use a subsetter to go back to get the male and female average densities (could have done it before, but will now)
male_n <- ncol(SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ])
males_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

female_n <- ncol(SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ])
females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity <- full_join(males_true,females_true) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity$mean_MvF_dif > 0) #17163
table(true_SNPdensity$mean_MvF_dif < 0) #15626
table(true_SNPdensity$mean_MvF_dif == 0) #4351

# below is a check that the same could be done with positional information so that we can do so in our permutations:

# given male_n and female_n we should be able to re-calculate the means with position information:
#males_true_position <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
#identical(males_true_position,males_true)
#females_true_position <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
#true_SNPdensity_position <- full_join(males_true_position,females_true_position) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
#table(true_SNPdensity_position$mean_MvF_dif > 0) 
#table(true_SNPdensity_position$mean_MvF_dif < 0) 
#table(true_SNPdensity_position$mean_MvF_dif == 0) 


# so, now we would like to permute!
# this means shuffling - we can shuffle the data columns and then use the code above.

# need an initial run of the permutation:

perm <- sample(3:length(SNPdensity.rows))
SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)

males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity_position_perm <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)


perm_backbone <- true_SNPdensity_position_perm %>% select(scaf,base,mean_MvF_dif) %>% rename(p1=mean_MvF_dif)

for(i in 2:1000){
  perm <- sample(3:length(SNPdensity.rows))
  SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)
  males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
  females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
  column_ID <- paste0("p",i)
  perm_upgrade <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females) %>% select(scaf,base,mean_MvF_dif) %>% dplyr::rename(!!column_ID := "mean_MvF_dif")
  perm_backbone <- full_join(perm_backbone,perm_upgrade)
}

# never write this again - needs to be saved because no see set write_tsv(perm_backbone,"perm_backbone_fugu.txt")
perm_backbone <- read_tsv("perm_backbone_fugu.txt")
perm_with_true <- full_join(true_SNPdensity %>% select(scaf,base,mean_MvF_dif),perm_backbone)

perm_with_true_p <- perm_with_true %>% 
  mutate(Pvalue = case_when(mean_MvF_dif > 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                             mean_MvF_dif < 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                             mean_MvF_dif == 0 ~ 1)) # here we are using a conditional mutate where if the MvF difference is > 0 we look for outliers above, if the MvF difference is < 0  we look for outliers below, and if the MvF difference is ==0, we say p value is 1.

#write_tsv(perm_with_true_p,"SNPdensity_perm_with_true_p_fugu.txt")


# want a look at proportion of scaffold:
fugu_index <- read_tsv("GCF_901000725.2_fTakRub1.2_genomic.fna.fai",col_names = F) %>% rename(scaf=X1,length=X2)

perm_with_true_p001_snp <- perm_with_true_p %>% filter(Pvalue <= 0.001) 
count_snp_p001_scaffolds <- perm_with_true_p001_snp %>% select(scaf) %>% count(scaf)
proportion_count_snp_p001_scaffolds <- left_join(count_snp_p001_scaffolds,fugu_index) %>% mutate(proportion = (n*10000)/length) %>% filter(scaf!="NC_004299.1")
proportion_count_snp_p001_scaffolds %>% filter(length > 10000000) %>% ggplot(aes(x=scaf,y=proportion)) + geom_col() + theme(axis.text.x = element_text(angle = 90))
#ggsave("fugu_scaf_proportion_main22chr_p001.pdf")

