#This script imports count data, tidies it and then runs GENIE3
#counts are available in the raw_counts folder on Google Drive

#set working directory
setwd('/Users/aideenmccabe/Documents/College Documents/BSPP2020/BSPP_2020')

#import necessary packages
library(readr)
library(dplyr)
library(tidyr)
library(doParallel)
library(doRNG)
library(plyr)
library(GENIE3)

#importing raw count data
ERP013829_count <- read_delim("Data/ERP013829_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ERP003465_count <- read_delim("Data/ERP003465_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP060670_count <- read_delim("Data/SRP060670_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ERP009837_count <- read_delim("Data/ERP009837_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP022869_count <- read_delim("Data/SRP022869_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP041017_count <- read_delim("Data/SRP041017_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP048912_count <- read_delim("Data/SRP048912_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP068165_count <- read_delim("Data/SRP068165_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP045409_count <- read_delim("Data/SRP045409_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

#function to clean the data
cleancounts <- function(counts){
  test <- counts %>% 
    separate(transcript, c("Gene", "bin"), remove = FALSE) %>% 
    select(-bin)
  
  test_agg <- aggregate(test, by = list(test$Gene), FUN = mean)
  names(test_agg[names(test_agg)=="Group.1"]) <- "GeneID"
  final_df <- test_agg %>% 
    select(-c(transcript, Gene))
  
  return(final_df)
}

#run cleaning function on all count data
ERP013829_count <- cleancounts(ERP013829_count)
ERP003465_count <- cleancounts(ERP003465_count)
SRP060670_count <- cleancounts(SRP060670_count)
ERP009837_count <- cleancounts(ERP009837_count)
SRP022869_count <- cleancounts(SRP022869_count)
SRP041017_count <- cleancounts(SRP041017_count)
SRP048912_count <- cleancounts(SRP048912_count)
SRP068165_count <- cleancounts(SRP068165_count)
SRP045409_count <- cleancounts(SRP045409_count)

#join by transcript ID
all_counts <- join_all(list(ERP003465_count, SRP060670_count, ERP009837_count, SRP022869_count,
                            SRP041017_count, SRP048912_count, SRP068165_count, SRP045409_count, ERP013829_count), 
                       by = 'Group.1', type = 'left')

#filter for hotspot genes and iTAK transcription factors
#import genes and tfs
all_hotspot_genes_and_tfs <- read_csv("all_hotspot_genes_and_tfs.csv")

#keep only genes that are in gene + tf list
all_counts_filtered <- all_counts %>% 
  filter(Group.1 %in% all_hotspot_genes_and_tfs$GeneID)

#change name of the column to GeneID
names(all_counts_filtered)[names(all_counts_filtered)=="Group.1"] <- "GeneID"

#check how many genes and TFs are present in final counts
#get genes only 
all_hotspot_genes <- read_csv("all_hotspot_genes.csv")

#see how many genes are present
gene_check <- all_counts_filtered %>% 
  filter(GeneID %in% all_hotspot_genes$GeneID)

#get tfs from this list
tfs <- all_hotspot_genes_and_tfs %>% 
  filter(!GeneID %in% all_hotspot_genes$GeneID)

#see how many tfs are present
tf_check <- all_counts_filtered %>% 
  filter(GeneID %in% tfs$GeneID)

#change to matrix and assign rownames
all_counts_filtered_nogene <- all_counts_filtered %>% 
  select(-GeneID)

count_matrix <- as.matrix(all_counts_filtered_nogene)

rownames(count_matrix) <- all_counts_filtered$GeneID 

#filter for tfs in matrix
tf_keep <- tfs %>% 
  filter(GeneID %in% all_counts_filtered$GeneID) %>% 
  distinct(GeneID, .keep_all = TRUE) %>% 
  pull(GeneID)

#run GENIE, specifying iTAK tfs in matrix as regulators
weightMat <- GENIE3(count_matrix, regulators = tf_keep)

#extract top 5000 links
linklist <- getLinkList(weightMat, reportMax = 5000)


#filter link list for hotspot genes 
filtered_link <- linklist %>% 
  filter(targetGene %in% all_hotspot_genes$GeneID)

#write to text file
write.table(filtered_link, "top5000interactions.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#write to csv file 
write.csv(filtered_link, "top5000interactions.csv", row.names = FALSE)

