#Hub analysis
library(dplyr)
library(readr)
library(Biostrings)
library(tidyr)

#importing file
final_hotspots <- read_csv("Data/final_hotspots.csv")

#extracting geneID to run through miRNA database
hotspot_id <- final_hotspots %>% 
  select(GeneID) %>% 
  unique()

#write to file
write.table(hotspot_id, "hotspot_id.txt", quote = FALSE, row.names = FALSE)

#this returns a csv file with miRNA interactions from PmiREN database
hotspot_miRNA <- read_csv("hotspot_miRNA.csv")

#now to get TFs
#first need sequence of genes
#fasta_hotspot <- read.fasta("/Users/aideenmccabe/Documents/College Documents/BSPP2020/BSPP_2020/Data/wheat_all.fa")
fasta_hotspot <- readDNAStringSet("/Users/aideenmccabe/Documents/College Documents/BSPP2020/BSPP_2020/Data/wheat_all.fa")
fastadf<-as.data.frame(fasta_hotspot)
fastadf$names <- rownames(fastadf)
rownames(fastadf) <- NULL
colnames(fastadf) <- c("sequence", "info") 

expanded <- fastadf %>% 
  separate(col = info, into = c("geneID", "cds", 'chromosome', 'gene', 'biotype', 'tbiotype'), sep = " ")

full_fasta <- expanded %>% 
  select(sequence, geneID)

#write fasta to csv format
write.csv(full_fasta, "full_fasta.csv", row.names = FALSE)

#cut off the trailing number
full_fasta$geneID = substr(full_fasta$geneID, 1, nchar(full_fasta$geneID)-2)

#crossover
crossover <-full_fasta %>% 
  filter(geneID %in% hotspot_id$GeneID) %>% 
  distinct(geneID, .keep_all = TRUE)

#reorder
to_fasta <- crossover[, c(2, 1)]

#now to convert to fasta file
#adding a >
to_fasta$geneID <- paste(">", to_fasta$geneID, sep="")
#write to file
write.table(to_fasta, "fasta_test.txt", sep='\n', quote = FALSE, row.names = FALSE)
