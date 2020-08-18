#script to take input list of gene IDs and return a fasta file with sequences
#importing relevant packages
library(dplyr)
library(readr)
library(tidyr)

#import list of geneIDs to get fasta sequences for 
#column name must be geneID
#if using transcript ids, can trim the trailing number with the following code;
#list$geneid <- substr(list$transcriptid, 1, nchar(list$transcriptid)-2)
example_list <- read_csv("hotspot_id.txt")
colnames(example_list) <- c("geneID")

#import all fasta sequences from IWGSC
full_fasta <- read_csv("full_fasta.csv")

#function get_fasta takes two inputs
# 1. id_list: .txt file:  list of geneIDs to get sequences of
# 2. full_fasta: full fasta file (provided)
# returns: : a fasta file with sequences of query genes

get_fasta <- function(id_list, full_fasta){
  sequences <- full_fasta %>% 
    filter(geneID %in% id_list$geneID) %>% 
    distinct(geneID, .keep_all = TRUE)
  
  sequences$geneID <- paste(">", sequences$geneID, sep = "")
  
  return(sequences)
}

#run function
to_fasta <- get_fasta(example_list, full_fasta)

#write to file, specify the name of file
write.table(to_fasta, "example_fasta.txt", sep='\n', quote = FALSE, row.names = FALSE)
