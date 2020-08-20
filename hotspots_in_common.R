#this script compares hotspots from two datasets and isolates those with >80% similarity

#importing necessary packages
library(dplyr)
library(readr)

#importing example file
final_hotspots <- read_csv("Data/final_hotspots.csv")

#adding dataset identifier to the hotspot column
final_hotspots$Hotspot <- sub("^", "Fusarium", final_hotspots$Hotspot)

#split dataset into a list of smaller dfs based on hotspot
#for this example, using the same dataset
dataset1 <- split(final_hotspots , f = final_hotspots$Hotspot)

dataset2 <- split(final_hotspots , f = final_hotspots$Hotspot)


#matchy function takes two lists of hotspots and returns a list with hotspots that are >80% similar
matchy <- function(out, out2){
  
  #create an empty list
  matches <- list()
  
  #iterate through list/dataset number 1
  for (i in 1:length(out)){
  
    #iterate through list/dataset number 2
    for (j in 1:length(out2)) {
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

#running the function and turning the list into one big dataframe
test <- matchy(dataset1, dataset2)
testdf <- bind_rows(test, .id = "column_label")

#need to expand this to take more than two datasets as input 


