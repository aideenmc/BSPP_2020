#Hub analysis
library(dplyr)
library(readr)

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

#fasta_hotspot <- read.fasta(file = "/Users/aideenmccabe/Documents/College Documents/BSPP2020/BSPP_2020/Data/wheat_all.fa")

