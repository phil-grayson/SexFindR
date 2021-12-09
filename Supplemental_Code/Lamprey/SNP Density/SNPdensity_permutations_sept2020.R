library(tidyverse)
library(ggpubr)

# test.df <- data.frame(chrom,start,m1,m2,f1,f2) %>% tbl_df()
# 
# constant <- test.df[1:2]
# shuffle <- test.df[sample(c(3,4,5,6))]
# perm <- bind_cols(constant,shuffle)


# the data is already split in to males and females

# load in the male data:
setwd("/Users/phil/Desktop/SNPdensity_individual_306/")
myFiles <- list.files(pattern="Male*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][6],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  #firstsplit <- strsplit(myFiles[i], "SNP_snpdensity_sorted_no_dups_sorted_trim_")
  firstsplit <- strsplit(myFiles[i], "_")
  if(firstsplit[[1]][6]=="CTR"){ # deals with CTR
    temp <- paste0(firstsplit[[1]][6],firstsplit[[1]][7])
    column_ID <- paste(firstsplit[[1]][1],temp,sep="_")
  } else if (firstsplit[[1]][5]=="M"){ # deals with my samples
    temp <- paste0(firstsplit[[1]][5],firstsplit[[1]][6])
    column_ID <- paste(firstsplit[[1]][1],temp,sep="_")
  } else {
    column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][6],sep="_")
  }
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}

SNPdensity.males <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")
SNPdensity.males <- backbone.upgrade
# load the female data:

myFiles <- list.files(pattern="Female*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][6],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  #firstsplit <- strsplit(myFiles[i], "SNP_snpdensity_sorted_no_dups_sorted_trim_")
  firstsplit <- strsplit(myFiles[i], "_")
  if(firstsplit[[1]][6]=="CTR"){ # deals with CTR
    temp <- paste0(firstsplit[[1]][6],firstsplit[[1]][7])
    column_ID <- paste(firstsplit[[1]][1],temp,sep="_")
  } else if (firstsplit[[1]][5]=="F"){ # deals with my samples
    temp <- paste0(firstsplit[[1]][5],firstsplit[[1]][6])
    column_ID <- paste(firstsplit[[1]][1],temp,sep="_")
  } else {
    column_ID <- paste(firstsplit[[1]][1],firstsplit[[1]][6],sep="_")
  }
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}

SNPdensity.females <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")


SNPdensity <- full_join(SNPdensity.males,SNPdensity.females)

SNPdensity.rows <- SNPdensity
SNPdensity.rows$base <- as.numeric(as.character(SNPdensity.rows$base))

# want to replace all NA with 0
SNPdensity.rows <- SNPdensity.rows %>% replace_na(list(Males = 0, Females = 0)) %>% replace(is.na(.), 0)
#write_tsv(SNPdensity.rows,"SNPdensity_rows")
#write_tsv(SNPdensity.rows,"SNPdensity_rows_location.txt") # ran through without splitting the base and scaf

# use a subsetter to go back to get the male and female average densities (could have done it before, but will now)
male_n <- ncol(SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ])
males_true <- bind_cols(SNPdensity.rows %>% select(1),SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(LOCATION,mean_Males)
#males_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Male" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

female_n <- ncol(SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ])
females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(LOCATION,mean_Females)
#females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "Female" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity <- full_join(males_true,females_true) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity$mean_MvF_dif > 0) #54k
table(true_SNPdensity$mean_MvF_dif < 0) #36k
table(true_SNPdensity$mean_MvF_dif == 0) #10k

# given male_n and female_n we should be able to re-calculate the means with position information:
males_true_position <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
#identical(males_true_position,males_true)
females_true_position <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity_position <- full_join(males_true_position,females_true_position) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity_position$mean_MvF_dif > 0) #54k
table(true_SNPdensity_position$mean_MvF_dif < 0) #36k
table(true_SNPdensity_position$mean_MvF_dif == 0) #10k

# Great!

# so, now we would like to permute!
# this means shuffling - we can shuffle the data columns and then use the code above.

# here's the key!
perm <- sample(3:length(SNPdensity.rows))
SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)

males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity_position_perm <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity_position_perm$mean_MvF_dif > 0) 
table(true_SNPdensity_position_perm$mean_MvF_dif < 0) 
table(true_SNPdensity_position_perm$mean_MvF_dif == 0) 

# note that a seed was not set, and so, this would not be exactly replicated if run again (but we did save the output as perm_backbone.txt.  we can do this another time with a seed when we ramp up.
# need an initial run of the permutation:
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

# never write this again - needs to be saved because no see set write_tsv(perm_backbone,"perm_backbone.txt")
perm_backbone <- read_tsv("perm_backbone.txt")
perm_with_true <- full_join(true_SNPdensity %>% select(scaf,base,mean_MvF_dif),perm_backbone)

#%>% rowwise() %>% mutate(pvalue = mean(abs(c(p1,p2,p3)) > abs(mean_MvF_dif)))

#perm_with_true$pvalue <- rowSums((perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)/1000)

#perm_with_true$pvalue.abs <- rowSums(abs((perm_with_true[,grep("p", names(perm_with_true))]) > abs(perm_with_true$mean_MvF_dif))/1000)

# perm_with_true_p <- perm_with_true %>% 
#   mutate(p.value = case_when(mean_MvF_dif > 0 ~ rowSums((perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)/1000),
#                              mean_MvF_dif < 0 ~ rowSums((perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)/1000),
#                              mean_MvF_dif == 0 ~ 1)) # here we are using a conditional mutate where if the MvF difference is > 0 we look for outliers above, if the MvF difference is < 0  we look for outliers below, and if the MvF difference is ==0, we say p value is 1.


perm_with_true_p <- perm_with_true %>% 
  mutate(Pvalue = case_when(mean_MvF_dif > 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                             mean_MvF_dif < 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                             mean_MvF_dif == 0 ~ 1)) # here we are using a conditional mutate where if the MvF difference is > 0 we look for outliers above, if the MvF difference is < 0  we look for outliers below, and if the MvF difference is ==0, we say p value is 1.

#write_tsv(perm_with_true_p,"SNPdensity_perm_with_true_p.txt")
perm_with_true_p <- read_tsv("SNPdensity_perm_with_true_p.txt")

perm_with_true_p05_snp <- perm_with_true_p %>% filter(Pvalue < 0.05) %>% select(scaf,base,mean_MvF_dif,Pvalue)
perm_with_true_p001_snp <- perm_with_true_p05_snp %>% filter(Pvalue < 0.001) 

barplot(table(perm_with_true_p05_snp$scaf))
count_snp_p05_scaffolds <- perm_with_true_p05_snp %>% select(scaf) %>% count(scaf)
# original
table(true_SNPdensity_position$mean_MvF_dif > 0) #54k
table(true_SNPdensity_position$mean_MvF_dif < 0) #36k
table(true_SNPdensity_position$mean_MvF_dif == 0) #10k

# permuted and p05
table(perm_with_true_p05_snp$mean_MvF_dif>0) # 7.8k
table(perm_with_true_p05_snp$mean_MvF_dif<0) # 3.1 k
table(perm_with_true_p05_snp$mean_MvF_dif==0)

# remaining versus original windows
# male excess 
7794/54215 # 0.143761
# female excess
3180/36221 # 0.08779437

# maybe a slight excess in Male bias.

# look at top candidate
NC_046072_cand <- perm_with_true_p05_snp %>% filter(scaf=="NC_046072.1",base ==21620000) %>% gather(key="perm",value="MvF_perm",-scaf,-base,-mean_MvF_dif) %>% select(-perm)

hist(NC_046072_cand$MvF_perm,breaks = 50) 
abline(v=NC_046072_cand$mean_MvF_dif, lwd=2, col="purple")

true_SNPdensity_position %>% filter(scaf=="NC_046072.1",base ==21620000)
View(SNPdensity %>% filter(scaf=="NC_046072.1",base ==21620000))

# stacked bar plots generated in GWAS_comparison

# Nov 26 - ramping up for python permutations:
python <- read_tsv("~/Desktop/trans_test_permuttttter.txt", col_names = T)
test <- full_join(true_SNPdensity %>% select(LOCATION,mean_MvF_dif),python) %>% mutate(dif = mean_MvF_dif-true_m_v_f)

test <- read_tsv("~/Desktop/test.txt")
#test2 <- test %>% remove_rownames %>% column_to_rownames(var="LOCATION")
#test3 <- t(test2)

test4 <- t(test) # option 1
output <- as.data.frame(test4) %>% rownames_to_column()
write_tsv(output,"~/Desktop/test_pivot.txt",col_names = F)

#test %>% spread()

# Dec 1st

# coming back from the permutations:
perm_100k_p001 <- read_tsv("most_significant_permutations.txt",col_names =c("LOCATION","p_value")) %>% mutate(set = "100k")

perm_with_true_1k_p001 <- perm_with_true_p %>% filter(Pvalue <= 0.001) %>% select(scaf,base,mean_MvF_dif,Pvalue) %>% unite(LOCATION, c("scaf","base"),sep=":") %>% mutate(set="1k")

perms_consistency <- full_join(perm_with_true_1k_p001,perm_100k_p001,by="LOCATION") %>% mutate(same = !is.na(set.x) & !is.na(set.y))

table(perms_consistency$same)

one_thousand_checker <- perm_with_true_p %>% select(scaf,base,mean_MvF_dif,Pvalue) %>% unite(LOCATION, c("scaf","base"),sep=":")

View(perms_consistency %>% filter(set.y=="100k") %>% filter(mean_MvF_dif > 0.5 | mean_MvF_dif < -0.5))

gff <- read_table2("~/Desktop/kPetmar_gff_grep_gene_protein-coding.txt",col_names = F) %>% separate(X4,c("garbage","temp"),sep = "gene=") %>% separate(temp, c("gene",sep=";")) %>% select(X1,X2,X3,gene)

write_tsv(gff,"final_kPetmar_gff_grep_gene_protein-coding.txt",col_names = F)
gff %>% select(gene)
