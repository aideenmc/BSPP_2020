#Transcription factor binding site analysis
library(readr)
library(dplyr)

#importing hotspot
final_hotspots <- read_csv("Data/final_hotspots.csv")

#import TF data
fimo <- read_delim("fimo.txt", "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

#import miRNA data
hotspot_miRNA <- read_csv("hotspot_miRNA.csv")

#Binding sites in 142 genes
length(unique(fimo$`sequence name`))

#9 different transcription factor families
length(unique(fimo$Family))

#need 3 files for plotting regulatory network
#1. Interactor and Gene
#2. ID and type
#3. GeneID and Hotspot

#1. Interactor and Gene
hotspot_miRNA$target_gene = substr(hotspot_miRNA$target_gene, 1, nchar(hotspot_miRNA$target_gene)-2)

fimo_id_fam <- fimo %>% 
  select(Family, `sequence name`)

colnames(fimo_id_fam) <- c('Interactor', 'GeneID')

colnames(hotspot_miRNA) <- c('Interactor', 'GeneID')

interactor_and_gene <- bind_rows(fimo_id_fam, hotspot_miRNA)
interactor_and_gene <- interactor_and_gene %>% 
  unique()

write.table(interactor_and_gene, "hotspot_interactions.txt", row.names= FALSE, quote = FALSE, sep = "\t")

#2. ID and type
#a.Genes
hotspot_id <- final_hotspots %>% 
  select(GeneID)
hotspot_id$Type <- "Gene"
colnames(hotspot_id) <- c("ID", "Type")

#b.TFs 
TF_id <- fimo %>% 
  select(Family) %>% 
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

write.table(ID_and_Type, "ID_and_Type.txt", row.names = FALSE, quote = FALSE , sep = "\t")

#3 Gene ID and Hotspot
gene_hotspot_ids <- final_hotspots %>% 
  select(GeneID, Hotspot)

#Add a column for group (non-par, non_met etc)
nonh_nonmet <- c(2,3,6,12,16,17,19,20,21,23,40,24,26,28,31,35,39,41,42,43)
nonh_met <- c(8,13,15,44,11,27,34,36)
par_met<- c(7,9,10,32)
par_nonmet<- c(1,4,5,14,18,22,25,29,30,33,37,38)

gene_hotspot_ids$Group <- ifelse(gene_hotspot_ids$Hotspot %in% nonh_nonmet == TRUE, "Non_Homologous_Non_Metabolic", 
                                 ifelse(gene_hotspot_ids$Hotspot %in% nonh_met == TRUE, "Non_Homologous_Metabolic",
                                        ifelse(gene_hotspot_ids$Hotspot %in% par_met==TRUE, "Paralagous_Metabolic",
                                               ifelse(gene_hotspot_ids$Hotspot %in% par_nonmet==TRUE, "Paralagous_Non_Metabolic", "x"))))

gene_hotspot_ids <- gene_hotspot_ids %>% 
  unique()

#write to table
write.table(gene_hotspot_ids, "gene_hotspot_ids.txt", row.names = FALSE, quote = FALSE, sep = "\t")


