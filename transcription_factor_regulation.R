#Transcription factor regulation analysis
library(readr)
library(dplyr)
library(ggplot2)

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

#investigating what hotspot tf interactions are located
gene_hotspot <- gene_hotspot_ids %>% 
  select(-Group)

gene_hotspot_tf <- left_join(gene_hotspot, tfs, by="GeneID")
nona_gene_hotspot_tf <- as.data.frame(na.omit(gene_hotspot_tf))
nona_gene_hotspot_tf <- unique(nona_gene_hotspot_tf)

#how many DIFFERENT TF interactions per hotspot
tf_hotspot_counts <- as.data.frame(table(nona_gene_hotspot_tf$Hotspot))

p <- ggplot(data=tf_hotspot_counts, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="deeppink") +
  theme_minimal()

p

#grouping by cluster groups
colnames(tf_hotspot_counts) <- c("Hotspot", "Frequency")
nonh_nonmet <- c(2,3,6,12,16,17,19,20,21,23,40,24,26,28,31,35,39,41,42,43)
nonh_met <- c(8,13,15,44,11,27,34,36)
par_met<- c(7,9,10,32)
par_nonmet<- c(1,4,5,14,18,22,25,29,30,33,37,38)

tf_hotspot_counts$Group <- ifelse(tf_hotspot_counts$Hotspot %in% nonh_nonmet == TRUE, "Nonhomologous Non-Metabolic", 
                                 ifelse(tf_hotspot_counts$Hotspot %in% nonh_met == TRUE, "Nonhomologous Metabolic",
                                        ifelse(tf_hotspot_counts$Hotspot %in% par_met==TRUE, "Paralagous Metabolic",
                                               ifelse(tf_hotspot_counts$Hotspot %in% par_nonmet==TRUE, "Paralagous Non-Metabolic", "x"))))

#factor hotspots to make visualising nicer
tf_hotspot_counts$Hotspot <- factor(tf_hotspot_counts$Hotspot, 
                                    levels = c(2,3,6,12,16,17,19,20,21,23,24,26,28,31,35,39,40,41,42,43,
                                               8,11,13,15,27,34,36,44,
                                               7,9,10,32,
                                               1,4,5,14,18,22,25,29,30,33,37,38))


q <- ggplot(data=tf_hotspot_counts, aes(x=Hotspot, y=Frequency, fill = Group)) +
  ggtitle("Number of Transcription Factor Interactions per Hotspot") +
  xlab("Hotspot") +
  ylab("Number of Transcription Factor Interactions") +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold.italic", family = "Helvetica", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold", family = "Helvetica"),
        axis.title.y = element_text(size = 14, face = "bold", family = "Helvetica"))

q

#plotting number of transcription factor interactions per group
tf_group_counts <- as.data.frame(table(gene_hotspot_ids$Group))
colnames(tf_group_counts) <- c("Group", "Freq")

r <- q <- ggplot(data=tf_group_counts, aes(x=Group, y=Freq, fill = Group)) +
  ggtitle("Number of Transcription Factor Interactions per Group") +
  xlab("Gene Cluster Group") +
  ylab("Number of Transcription Factor Interactions") +
  geom_bar(stat="identity") +
  theme_minimal()

r
