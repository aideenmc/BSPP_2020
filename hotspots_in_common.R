library(dplyr)
library(readr)

#importing example file
final_hotspots <- read_csv("Data/final_hotspots.csv")
final_hotspots$Hotspot <- sub("^", "Fusarium", final_hotspots$Hotspot)

#split dataset into a list of smaller dfs based on hotspot
dataset1 <- split(final_hotspots , f = final_hotspots$Hotspot)

dataset2 <- split(final_hotspots , f = final_hotspots$Hotspot)


#######this works!
matchy <- function(out, out2){
  
  matches <- list()
  
  for (i in 1:length(out)){
  
    for (j in 1:length(out2)) {
        max <- if (nrow(out[[i]]) > nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[i]] else out2[[j]]
        min <- if (nrow(out[[i]]) < nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out2[[j]] else out2[[j]]
        crosscheck <- max %>% 
          filter(GeneID %in% min$GeneID)
        match <- if (nrow(crosscheck) > nrow(min)*0.8) TRUE else FALSE
        
        if (match == TRUE){
          #need to make sure name of hotspot in original inputs contains name of dataset
          crosscheck$comparison <- paste(out[[i]][[1]][[1]], out2[[j]][[1]][[1]], sep = "_")
          matches[[j]] <- as.data.frame(crosscheck)
          print(match)
        }
        
        else
          print(match)
      }
  
  }

  return(matches)
}

#running the function and turning the list into one big dataframe
test <- matchy(dataset1, dataset2)
testdf <- bind_rows(test, .id = "column_label")

#need to expand this to take more than two datasets as input 


