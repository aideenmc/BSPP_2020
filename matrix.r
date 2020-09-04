#importing necessary packages
library(dplyr)
library(readr)
library(plyr)

#all data
fusarium_og_hotspots <- read_csv("Gene_cluster_analysis/fusarium_processed.csv")
septoria_erp_hotspots <- read_csv("Gene_cluster_analysis/densityERP009837_processed.csv")
powdery_rust_hotspots  <- read_csv("Gene_cluster_analysis/densitySRP041017_processed.csv")
peg_hotspots <- read_csv("Gene_cluster_analysis/densitySRP068165_processed.csv")
fusarium_erp_hotspots <-  read_csv("Gene_cluster_analysis/densityERP003465_processed.csv")
fusarium_srp_hotspots <- read_csv("Gene_cluster_analysis/densitySRP060670_processed.csv")
heat_drought_hotspots <- read_csv("Gene_cluster_analysis/densitySRP045409_processed.csv")
fusarium_pseudo_hotspots <- read_csv("Gene_cluster_analysis/densitySRP048912_processed.csv")
septoria_srp_hotspots <- read_csv("Gene_cluster_analysis/densitySRP022869_processed.csv")

#adding dataset identifier to the hotspot column
fusarium_og_hotspots$Hotspot_Number <- sub("^", "Fusarium_og", fusarium_og_hotspots$Hotspot_Number)
septoria_erp_hotspots$Hotspot_Number <- sub("^", "Septoria_erp", septoria_erp_hotspots$Hotspot_Number)
powdery_rust_hotspots$Hotspot_Number <- sub("^", "Powdery_and_Rust", powdery_rust_hotspots$Hotspot_Number)
peg_hotspots$Hotspot_Number <- sub("^", "PEG", peg_hotspots$Hotspot_Number)
fusarium_erp_hotspots$Hotspot_Number <- sub("^","Fusarium_erp", fusarium_erp_hotspots$Hotspot_Number)
fusarium_srp_hotspots$Hotspot_Number <- sub("^", "Fusarium_srp", fusarium_srp_hotspots$Hotspot_Number)
heat_drought_hotspots$Hotspot_Number <- sub("^", "Heat_Drought", heat_drought_hotspots$Hotspot_Number)
fusarium_pseudo_hotspots$Hotspot_Number <- sub("^", "Fusarium_pseudo", fusarium_pseudo_hotspots$Hotspot_Number)
septoria_srp_hotspots$Hotspot_Number <- sub("^", "Septoria_srp", septoria_srp_hotspots$Hotspot_Number)

#concatanate into one large dataframe
all_data <- bind_rows(fusarium_og_hotspots, septoria_erp_hotspots, powdery_rust_hotspots,
                      peg_hotspots, fusarium_erp_hotspots, fusarium_srp_hotspots,
                      heat_drought_hotspots, fusarium_pseudo_hotspots, septoria_srp_hotspots)

#####
#split by stress/disease
all_data_split <- split(all_data , f = all_data$`Disease/Stress`)

all_data_split[[1]][["Disease/Stress"]]
#create an empty matrix with rownames and column names = unique hotspot identifiers
m <- matrix(0, nrow = length(all_data_split), ncol = length(all_data_split))


#create an empty vector to be populated by hotspot IDs
unique_ids <- vector(mode = "character", length = length(all_data_split))

for (a in 1:length(all_data_split)){
  unique_ids[[a]] <- all_data_split[[a]][["Disease/Stress"]][[1]]
}


#assign these identifiers as row and column names
rownames(m) <- unique_ids
colnames(m) <- unique_ids

#####
#matchy function takes two lists of hotspots and returns a list with hotspots that are >80% similar
matchy <- function(out, out2){
  
  #create an empty list
  matches <- list()
  
  s = 1
  #iterate through list/dataset number 1
  for (i in 1:length(out)){
    #iterate through list/dataset number 2
    #for (j in i:length(out2)) {
    for (j in 1:length(out2)) {
      #if hotspots are the same, skip the comparison code      
      if (out[[i]][['Hotspot_Number']][[1]] == out2[[j]][['Hotspot_Number']][[1]]){
        next()
      }
      #need to determine which hotspot has most genes, assign this to max, and the hotspot with least genes is min
      max <- if (nrow(out[[i]]) > nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[i]] else out2[[j]]
      min <- if (nrow(out[[i]]) < nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[j]] else out2[[j]]
      
      #determine which genes are in common
      crosscheck <- max %>% 
        filter(GeneID %in% min$GeneID)
      
      #classify as a match if the number of genes in common represents at least 80% of smaller hotspot
      #this allows for similar yet non identical hotspots to be classified as a match
      match <- if (nrow(crosscheck) >= nrow(min)*0.8) TRUE else FALSE
      
      #if a match, then add a column specifying which datasets these genes are matched from then append to a list
      if (match == TRUE){
        #need to make sure name of hotspot in original inputs contains name of dataset
        crosscheck$stress1 <- paste(out[[i]][[10]][[1]])
        crosscheck$stress2 <- paste(out2[[j]][[10]][[1]])
        crosscheck$comparison1 <- paste(out[[i]][[1]][[1]])
        crosscheck$comparison2 <- paste(out2[[j]][[1]][[1]])
        #matches[[j]] <- as.data.frame(crosscheck)
        matches[[s]] <- as.data.frame(crosscheck)
        #print statement to see progress
        #print(match)
        s = s + 1
      }
      
      print(match)
      
    }
    
  }
  #once each hotspot has been analysed, return list of matches
  #return(m)
  return(matches)
}


#split dataset into a list of smaller dfs based on hotspot
#duplicate for pairwise comparison
dataset1 <- split(all_data , f = all_data$Hotspot_Number)

dataset2 <- split(all_data , f = all_data$Hotspot_Number)


#running the function and turning the list into one big dataframe
hotspots_in_common <- matchy(dataset1, dataset2)

testdf <- bind_rows(hotspots_in_common)

#need to remove duplicates
#unique_df <- testdf[!duplicated(testdf[,c('stress1','stress2')]),]

#now to add 1s to matrix 
#extract stress columns to use as i and j values
i_val <- testdf %>% 
  pull(stress1)

j_val <- testdf %>% 
  pull(stress2)

#assign 1 to these positions in the matrix
coordinates <-testdf %>% 
  select(stress1, stress2)

colnames(coordinates) <- c("x", "y")



for (i in 1:nrow(coordinates)){
  x <- coordinates[i,1]
  y <- coordinates[i,2]
  for (d in 1:nrow(m)){
    if (x == rownames(m)[d]){
      for(e in 1:ncol(m)){
        if(y==colnames(m)[e]){
          m[d,e]<-1
        }
      }
    }
  }
}

#now for all hotspot results



