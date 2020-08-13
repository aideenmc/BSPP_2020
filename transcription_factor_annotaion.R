#transcription factor annotation
library(dplyr)
library(readr)

#from itak.com
Bread_wheat_transcription_factor <- read_delim("Bread_wheat-transcription factor.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(Bread_wheat_transcription_factor) <- c("TF_ID", "TF_Type", "Protein_Seq", "Nt_Seq")

#remove transcript id
wheat_tf_annotation <- data.frame(do.call('rbind', strsplit(as.character(Bread_wheat_transcription_factor$TF_ID), '.', fixed = TRUE)))
Bread_wheat_transcription_factor$TF_ID <- wheat_tf_annotation$X1
tf_info <- Bread_wheat_transcription_factor %>% 
  select(TF_ID, TF_Type) %>% 
  unique()

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
