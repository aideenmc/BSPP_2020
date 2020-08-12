#Transcription factor regulation analysis
library(readr)
library(dplyr)

#importing data
#1. Hotspot genes
final_hotspots <- read_csv("Data/final_hotspots.csv")

#2. Transcription factors
scan_final <- read_delim("scan_final.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
colnames(scan_final) <- c("TF", "Gene Target", "Start", "End", "Strand", "pvalue", "Sequence", "Species")

#there are 244? regulations between 44 TFs and 68 genes
length(unique(scan_final$`Gene Target`))
length(unique(scan_final$TF))

#3. miRNAs
hotspot_miRNA <- read_csv("hotspot_miRNA.csv")
#remove transcript id from the end of the gene id
hotspot_miRNA$target_gene = substr(hotspot_miRNA$target_gene, 1, nchar(hotspot_miRNA$target_gene)-2)

#now to manipulate into a format readable by cytoscape (including miRNAs)
#we need 3 files
#1. Interactor and Gene
#2. ID and type
#3. GeneID and Hotspot

#1. Interactor and Gene
tfs <- scan_final %>% 
  select(TF, `Gene Target`)
colnames(tfs) <- c('Interactor', 'GeneID')

colnames(hotspot_miRNA) <- c('Interactor','GeneID')

interactor_and_gene <- bind_rows(tfs, hotspot_miRNA)

interactor_and_gene <- interactor_and_gene %>% 
  unique()

write.table(interactor_and_gene, "new_interactions.txt", row.names= FALSE, quote = FALSE, sep = "\t")

#2. ID and type
#a.Genes
hotspot_id <- final_hotspots %>% 
  select(GeneID)
hotspot_id$Type <- "Gene"
colnames(hotspot_id) <- c("ID", "Type")

#b.TFs 
TF_id <- tfs%>% 
  select(Interactor) %>% 
  unique()
TF_id$Type <- "TF"
colnames(TF_id) <- c("ID", "Type")

#c.miRNA
hotspot_miRNA_id <- hotspot_miRNA %>% 
  select(Interactor)
hotspot_miRNA_id$Type <- "miRNA"
colnames(hotspot_miRNA_id) <- c("ID", "Type")


ID_and_Type <- bind_rows(hotspot_id, TF_id, hotspot_miRNA_id)

ID_and_Type <- ID_and_Type %>% 
  unique()

write.table(ID_and_Type, "new_ID_and_Type.txt", row.names = FALSE, quote = FALSE , sep = "\t")

#3. GeneID and hotspots
#see transcription_factor_binding.R for code 
#here it is if want to have a look
gene_hotspot_ids <- read_delim("gene_hotspot_ids.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)


