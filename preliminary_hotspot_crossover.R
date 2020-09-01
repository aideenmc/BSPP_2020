#testing hotspots in common with seven datasets
#three fusarium, two abiotic, septoria, powdery mildew + stripe rust

library(dplyr)
library(readr)

#importing hotspot data
fusarium_og_hotspots <- read_csv("Gene_cluster_analysis/fusarium_processed.csv")
septoria_hotspots <- read_csv("Gene_cluster_analysis/densityERP009837_processed.csv")
powdery_rust_hotspots  <- read_csv("Gene_cluster_analysis/densitySRP041017_processed.csv")
peg_hotspots <- read_csv("Gene_cluster_analysis/densitySRP068165_processed.csv")
fusarium_erp_hotspots <-  read_csv("Gene_cluster_analysis/densityERP003465_processed.csv")
fusarium_srp_hotspots <- read_csv("Gene_cluster_analysis/densitySRP060670_processed.csv")
heat_drought_hotspots <- read_csv("Gene_cluster_analysis/densitySRP045409_processed.csv")


#adding dataset identifier to the hotspot column
fusarium_og_hotspots$Hotspot_Number <- sub("^", "Fusarium_og", fusarium_og_hotspots$Hotspot_Number)
septoria_hotspots$Hotspot_Number <- sub("^", "Septoria", septoria_hotspots$Hotspot_Number)
powdery_rust_hotspots$Hotspot_Number <- sub("^", "Powdery_and_Rust", powdery_rust_hotspots$Hotspot_Number)
peg_hotspots$Hotspot_Number <- sub("^", "PEG", peg_hotspots$Hotspot_Number)
fusarium_erp_hotspots$Hotspot_Number <- sub("^","Fusarium_erp", fusarium_erp_hotspots$Hotspot_Number)
fusarium_srp_hotspots$Hotspot_Number <- sub("^", "Fusarium_srp", fusarium_srp_hotspots$Hotspot_Number)
heat_drought_hotspots$Hotspot_Number <- sub("^", "Heat_Drought", heat_drought_hotspots$Hotspot_Number)

#concatanate, leaving out fusairum for now
all_data <- bind_rows(fusarium_og_hotspots, powdery_rust_hotspots, septoria_hotspots, peg_hotspots, fusarium_erp_hotspots, fusarium_srp_hotspots, heat_drought_hotspots)
fusarium_only <- bind_rows(fusarium_og_hotspots, fusarium_erp_hotspots, fusarium_srp_hotspots)
abiotic_only <- bind_rows(peg_hotspots, heat_drought_hotspots)
diseases_only <- bind_rows(fusarium_og_hotspots, fusarium_erp_hotspots, fusarium_srp_hotspots, powdery_rust_hotspots, septoria_hotspots)

#matchy function takes two lists of hotspots and returns a list with hotspots that are >80% similar
matchy <- function(out, out2){
  
  #create an empty list
  matches <- list()
  
  #iterate through list/dataset number 1
  for (i in 1:length(out)){
    
    #iterate through list/dataset number 2
    for (j in 1:length(out2)) {
      #if hotspots are the same, skip the comparison code      
      if (out[[i]][['Hotspot_Number']] == out2[[j]][['Hotspot_Number']]){
       next()
      }
      #need to determine which hotspot has most genes, assign this to max, and the hotspot with least genes is min
      max <- if (nrow(out[[i]]) > nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[i]] else out2[[j]]
      min <- if (nrow(out[[i]]) < nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out2[[j]] else out2[[j]]
      #determine which genes are in common
      crosscheck <- max %>% 
        filter(GeneID %in% min$GeneID)
      #classify as a match if the number of genes in common represents at least 80% of smaller hotspot
      #this allows for similar yet non identical hotspots to be classified as a match
      match <- if (nrow(crosscheck) > nrow(min)*0.8) TRUE else FALSE
      
      #if a match, then add a column specifying which datasets these genes are matched from then append to a list
      if (match == TRUE){
        #need to make sure name of hotspot in original inputs contains name of dataset
        crosscheck$comparison <- paste(out[[i]][[1]][[1]], out2[[j]][[1]][[1]], sep = "_")
        matches[[j]] <- as.data.frame(crosscheck)
        #print statement to see progress
        print(match)
      }
      
      else
        #if not a match, just print statement to see progress
        print(match)
    }
    
  }
  #once each hotspot has been analysed, return list of matches
  return(matches)
}


#split disease dataset into a list of smaller dfs based on hotspot
disease_dataset1 <- split(diseases_only , f = diseases_only$Hotspot_Number)

disease_dataset2 <- split(diseases_only , f = diseases_only$Hotspot_Number)

#running the function and turning the list into one big dataframe
disease_test <- matchy(disease_dataset1, disease_dataset2)
disease_testdf <- bind_rows(disease_test, .id = "column_label")

#fusarium
fusarium_dataset1 <- split(fusarium_only , f = fusarium_only$Hotspot_Number)

fusarium_dataset2 <- split(fusarium_only , f = fusarium_only$Hotspot_Number)

#running the function and turning the list into one big dataframe
fusarium_test <- matchy(fusarium_dataset1, fusarium_dataset2)
fusarium_testdf <- bind_rows(fusarium_test, .id = "column_label")


