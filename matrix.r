#importing necessary packages
library(dplyr)
library(readr)
library(plyr)

#all hotspots from permutation test
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

#concatanate all hotspots into one large dataframe
all_data <- bind_rows(fusarium_og_hotspots, septoria_erp_hotspots, powdery_rust_hotspots,
                      peg_hotspots, fusarium_erp_hotspots, fusarium_srp_hotspots,
                      heat_drought_hotspots, fusarium_pseudo_hotspots, septoria_srp_hotspots)


#####Finding hotspots in common
#matchy function takes two lists of hotspots and returns a list with hotspots that are >=80% similar
matchy <- function(out, out2){
  
  #create an empty list
  matches <- list()
  #variable s for iterating
  s = 1
  #iterate through list/dataset number 1
  for (i in 1:length(out)){
    
    #iterate through list/dataset number 2 (starting from the ith value)
    for (j in i:length(out2)) {
    
      #if hotspots are the same, skip the comparison code      
      if (out[[i]][['Hotspot_Number']][[1]] == out2[[j]][['Hotspot_Number']][[1]]){
        next()
      }
      
      #need to determine which hotspot has most genes, assign this to max, and the hotspot with least genes is min
      max <- if (nrow(out[[i]]) > nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[i]] else out2[[j]]
      min <- if (nrow(out[[i]]) < nrow(out2[[j]])) out[[i]] else if (nrow(out[[i]]) == nrow(out2[[j]])) out[[j]] else out2[[j]]
      
      #determine which genes are in common and store in variable crosscheck
      crosscheck <- max %>% 
        filter(GeneID %in% min$GeneID)
      
      #classify as a match if the number of genes in common represents at least 80% of smaller hotspot
      #this allows for similar yet non identical hotspots to be classified as a match
      match <- if (nrow(crosscheck) >= nrow(min)*0.8) TRUE else FALSE
      
      #if a match, then add 4 columns; specifying which hotspots and which stresses these genes are matched from then append to the list matches
      if (match == TRUE){
       
        crosscheck$stress1 <- paste(out[[i]][[10]][[1]])
        crosscheck$stress2 <- paste(out2[[j]][[10]][[1]])
        crosscheck$comparison1 <- paste(out[[i]][[1]][[1]])
        crosscheck$comparison2 <- paste(out2[[j]][[1]][[1]])
        
        #add to list
        matches[[s]] <- as.data.frame(crosscheck)
        #next s
        s = s + 1
      }
      #print match to see progress
      print(match)
      
    }
    
  }
  #once each hotspot has been analysed, return list of matches
  return(matches)
}


#split all data into a list of smaller dfs based on hotspot
dataset1 <- split(all_data , f = all_data$Hotspot_Number)

#duplicate for pairwise comparison
dataset2 <- split(all_data , f = all_data$Hotspot_Number)

#split dataframe into a list by stress/disease (for matrix)
all_data_split <- split(all_data , f = all_data$`Disease/Stress`)

#running the function  #turning the list into one big dataframe
hotspots_in_common <- matchy(dataset1, dataset2)

#turning the list into one big dataframe
testdf <- bind_rows(hotspots_in_common)

#check which hotspots are expressed across multiple diseases
#create an empty matrix with rownames = hotspots and column names = stress
n <- matrix(0, nrow = length(dataset1), ncol = length(all_data_split))


#create two empty vectors to be populated by hotspot and stress IDs
stress_ids <- vector(mode = "character", length = length(all_data_split))
hotspot_ids <- vector(mode = "character", length = length(dataset1))

#populate stress vector with names of stresses
for (a in 1:length(all_data_split)){
  stress_ids[[a]] <- all_data_split[[a]][["Disease/Stress"]][[1]]
}

#populate hotspot vector with names of hotspots
for (b in 1:length(dataset1)){
  hotspot_ids[[b]] <- dataset1[[b]][["Hotspot_Number"]][[1]]
}

#assign these identifiers as row and column names
colnames(n) <- stress_ids
rownames(n) <- hotspot_ids

#now to add 1s to matrix 
#extract hotspot results to use as i (row) values and stress results as corresponding j (column) values
coordinates <-testdf %>% 
  select(comparison1, stress2)

#change col names to x and y
colnames(coordinates) <- c("x", "y")


#loop to add 1 to matrix for specified coordinates
for (i in 1:nrow(coordinates)){
  x <- coordinates[i,1]
  y <- coordinates[i,2]
  for (d in 1:nrow(n)){
    if (x == rownames(n)[d]){
      for(e in 1:ncol(n)){
        if(y==colnames(n)[e]){
          n[d,e]<-1
        }
      }
    }
  }
}

#convert to daraframe
n_df <- as.data.frame(n)

#add a column specifying the number of diseases the hotspot is present in
n_df$sum <- rowSums(n[,1:9])+1

#table showing the count distribution
count(n_df$sum)
