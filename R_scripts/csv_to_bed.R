# R script
# Usage: Rscript --vanilla csv_to_bed.R input_file.csv output_file.csv window_size

# IMPORTS
library(dplyr)

# FUNCTIONS
# function to convert numeric chromosome number to .bed file chromosome string (e.g.:"chr1")
chr_string <- function(number) {
  if (is.integer(number)){
    char <- as.character(number)
    chromosome = paste("chr", char, sep = "")
  }
  else 
    chromosome = NaN
  return(chromosome)
}

#MAIN

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
# test for correct arguments
if (length(args)!=2) {
  stop("Two arguments need to be supplied: input_file output_file", call.=TRUE)
} else {
  # define all file names and parameters
  input_file <- args[1]
  output_file <- args[2]
}

# load data
horvath_data <- read.csv(input_file, header = TRUE, sep = ";")
# get chromosomes vector in .bed file format
chr_list <- list(horvath_data[,"Chr"])
chromosomes <- lapply(chr_list, chr_string)
# get start end end positions
starts <-  horvath_data[, "MapInfo"]
ends <-  horvath_data[, "MapInfo"] + 1
coefficients <- horvath_data[,"CoefficientTraining"]
# create data frame with .bed file appropriate format
bed_data <-  cbind(data.frame(chromosomes), starts, ends, coefficients)
# write .bed file
write.table(bed_data, file=output_file, sep="\t", col.names = F, row.names = F, quote = F)

