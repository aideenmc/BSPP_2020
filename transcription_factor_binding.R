#Transcription factor binding site analysis
library(readr)
library(dplyr)

#import data
fimo <- read_delim("fimo.txt", "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

#Binding sites in 142 genes
length(unique(fimo$`sequence name`))

#