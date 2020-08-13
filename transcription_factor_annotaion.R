#transcription factor annotation
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

#from itak.com
Bread_wheat_transcription_factor <- read_delim("Bread_wheat-transcription factor.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(Bread_wheat_transcription_factor) <- c("TF_ID", "TF_Type", "Protein_Seq", "Nt_Seq")

#remove transcript id
wheat_tf_annotation <- data.frame(do.call('rbind', strsplit(as.character(Bread_wheat_transcription_factor$TF_ID), '.', fixed = TRUE)))
Bread_wheat_transcription_factor$TF_ID <- wheat_tf_annotation$X1
tf_info <- Bread_wheat_transcription_factor %>% 
  select(TF_ID, TF_Type) %>% 
  unique()
#plot transcription factor type by group
#need 3 columns
#a. Group
#b. tf_type
#c. number of interactions
#add gene interactor column
scan_final <- read_delim("scan_final.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
colnames(scan_final) <- c("TF", "Gene Target", "Start", "End", "Strand", "pvalue", "Sequence", "Species")

#get gene info
gene_tf <- scan_final %>% 
  filter(TF %in% Bread_wheat_transcription_factor$TF_ID) %>% 
  unique() %>% 
  select(TF, `Gene Target`)
colnames(gene_tf) <- c("TF", "GeneID")

#hotspot info
gene_hotspot_ids <- read_delim("gene_hotspot_ids.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

tf_hotspot_group <- left_join(gene_tf, gene_hotspot_ids, by = "GeneID")
colnames(tf_info)<- c("TF", "Family")
tf_hotspot_group <- left_join(tf_hotspot_group, tf_info, by="TF")

tf_hotspot_group <- unique(tf_hotspot_group)
#write to file
write.csv(tf_hotspot_group, "tf_info_full.csv", row.names = FALSE)


#extract group and family only
group_and_fam <- tf_hotspot_group %>% 
  select(Group, Family)

#counting tf fam by group
to_plot <- rename(count(group_and_fam, Group, Family), Freq = n)

s <- ggplot(data=to_plot, aes(x=Family, y=Freq, fill = Family)) +
  facet_wrap(~Group) +
  ggtitle("Number of Transcription Factor Interactions per Hotspot") +
  xlab("Hotspot") +
  ylab("Number of Transcription Factor Interactions") +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold.italic", family = "Helvetica", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold", family = "Helvetica"),
        axis.title.y = element_text(size = 14, face = "bold", family = "Helvetica")) +
        rotate_x_text(angle = 45)

s

##### 
#from plantregmap
scan_final <- read_delim("scan_final.txt", 
                         +     "\t", escape_double = FALSE, col_names = FALSE, 
                         +     trim_ws = TRUE)
colnames(scan_final) <- c("TF", "Gene Target", "Start", "End", "Strand", "pvalue", "Sequence", "Species")

tf_ids <- scan_final %>% 
  select(TF) %>% 
  unique()


#get info for our tfs
full_info <- tf_info %>% 
  filter(TF_ID %in% tf_ids$TF) %>% 
  unique()

#write to file for cytoscape
write.table(full_info, "40_tf_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
