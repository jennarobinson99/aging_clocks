# This file contains all functions needed for overlap analysis. 

overlap <- function(A, B){
  # A and B should be Ranges objects which can be reduced and overlapped
}


chr_string <- function(number) {
  # function to convert numeric chromosome number to .bed file chromosome string (e.g.:"chr1")
  
  if (is.integer(number)){
    char <- as.character(number)
    chromosome = paste("chr", char, sep = "")
  }
  else if (str_detect(number, "chr")[1]) 
    chromosome = number
  else{
    chromosome=NaN
  }
  return(chromosome)
}