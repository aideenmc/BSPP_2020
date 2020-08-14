library(readr)
library(tidyr)
library(dplyr)

new_ID_and_Type <- read_delim("new_ID_and_Type.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)


#1. Genes
final_hotspots <- read_csv("Data/final_hotspots.csv")

genes <- new_ID_and_Type[1:364,]
colnames(genes) <- c("GeneID", "Type")

GENE_ID_TYPE_DESC <- left_join(genes, final_hotspots, by = "GeneID")
GENE_ID_TYPE_DESC <- GENE_ID_TYPE_DESC %>% 
  select(GeneID, Type, `Gene Description`)
colnames(GENE_ID_TYPE_DESC) <- c("ID", "Type", "Desc")

#2. TFs
tf_info_full <- read_csv("tf_info_full.csv")

tfs <- new_ID_and_Type[365:408,]
colnames(tfs) <- c("TF", "Type")

TF_ID_TYPE_DESC <- left_join(tfs, tf_info_full, by = "TF")
TF_ID_TYPE_DESC <- TF_ID_TYPE_DESC %>% 
  select(TF, Type, Family) %>% 
  unique()
colnames(TF_ID_TYPE_DESC) <- c("ID", "Type", "Desc")
TF_ID_TYPE_DESC[is.na(TF_ID_TYPE_DESC)] <- "None"

#3. miRNAs
hotspot_miRNA <- read_csv("hotspot_miRNA.csv")
hotspot_miRNA$Desc <- hotspot_miRNA$miRNA
colnames(hotspot_miRNA) <- c("ID", "Type", "Desc")
hotspot_miRNA$Type <- "miRNA"

#combining all three
ID_TYPE_DESC_ALL <- bind_rows(GENE_ID_TYPE_DESC, TF_ID_TYPE_DESC, hotspot_miRNA)
#write to file for cytoscape 
write.table(ID_TYPE_DESC_ALL, "ID_Type_Desc.txt", sep = "\t", quote = FALSE, row.names = FALSE)
