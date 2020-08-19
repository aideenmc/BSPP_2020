library(dplyr)
library(readr)


matchy <- function(hotspot1, hotspot2){
  
  max <- if (length(hotspot1) > length(hotspot2)) hotspot1 else hotspot2
  min <- if (length(hotspot1) < length(hotspot2)) hotspot1 else hotspot2
  crosscheck <- max$GeneID %in% min$GeneID
  match <- if (length(crosscheck[crosscheck==TRUE]) > length(min)*0.8) TRUE else FALSE
  
  if (match == TRUE) {
    return(as.data.frame(bind_rows(max, min)))
  }
  else {
    return()
  }
  
}

test <- matchy(final_hotspots, hotspots2)

