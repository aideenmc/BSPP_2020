#big genie
library(readr)
library(dplyr)
library(tidyr)
library(doParallel)
library(doRNG)
library(plyr)
library(GENIE3)

#importing data
ERP013829_count <- read_delim("Data/ERP013829_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ERP003465_count <- read_delim("Data/ERP003465_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP060670_count <- read_delim("Data/SRP060670_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ERP009837_count <- read_delim("Data/ERP009837_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP022869_count <- read_delim("Data/SRP022869_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP041017_count <- read_delim("Data/SRP041017_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP048912_count <- read_delim("Data/SRP048912_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP068165_count <- read_delim("Data/SRP068165_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
SRP045409_count <- read_delim("Data/SRP045409_count.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


#join by transcript ID
all_counts <- join_all(list(ERP003465_count, SRP060670_count, ERP009837_count, SRP022869_count,
                            SRP041017_count, SRP048912_count, SRP068165_count, SRP045409_count, ERP013829_count), 
                       by = 'transcript', type = 'left')


#remove trailing decimal
all_counts_gene <- all_counts %>% 
  separate(transcript, c("GeneID", "bin")) %>% 
  select(-bin)
  
#Get unique values
all_counts_gene_unique <- all_counts_gene %>% 
  distinct(GeneID, .keep_all = TRUE)

#filter for hotspot genes
all_hotspot_genes_and_tfs <- read_csv("all_hotspot_genes_and_tfs.csv")

all_counts_filtered <- all_counts_gene_unique %>% 
  filter(GeneID %in% all_hotspot_genes_and_tfs$GeneID)

#check how many genes and TFs are present in final counts
all_hotspot_genes <- read_csv("all_hotspot_genes.csv")

gene_check <- all_counts_filtered %>% 
  filter(GeneID %in% all_hotspot_genes$GeneID)

tfs <- all_hotspot_genes_and_tfs %>% 
  filter(!GeneID %in% all_hotspot_genes$GeneID)

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

#run GENIE
weightMat <- GENIE3(count_matrix, regulators = tf_keep, nCores = 3)

linklist <- getLinkList(weightMat, reportMax = 1000)

#filter link list for hotspot genes 
filtered_link <- linklist %>% 
  filter(targetGene %in% all_hotspot_genes$GeneID)

write.table(filtered_link, "top1000interactions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(filtered_link, "top1000interactions.csv", row.names = FALSE)

#check tfs 
tf_check <- filtered_link %>% 
  filter(regulatoryGene %in% tfs$GeneID)

