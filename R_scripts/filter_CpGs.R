# This script filters CpGs into positive and negatively correlated CpGs and
# extends the base window

# R script
# Usage: Fixed window: Rscript --vanilla filter_CpGs.R input_file output_file window_size 
#        Variable window: Rscript --vanilla filter_CpGs.R input_file output_file start_size end_size interval



#MAIN

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
# test for correct arguments
# error message
message <- "Use script in one of the two versions: \n 
  Fixed window: Rscript --vanilla filter_CpGs.R input_file output_file window_size \n
  Variable window: Rscript --vanilla filter_CpGs.R input_file output_file start_size end_size interval \n
  Arguments must be integers."

if (length(args) == 3 ){
 window_size<- as.integer(args[3])
 input_file <- args[1]
 output_file <- args[2]
 range <- FALSE
  
} else if (length(args) == 5){
  start_size <- as.integer(args[3])
  end_size <- as.integer(args[4])
  interval <- as.integer(args[5])
  input_file <- args[1]
  output_file <- args[2]
  range <- TRUE
} else {
  stop(message, call.=TRUE)
}

# IMPORTS
library(dplyr)

# define positive and negative output paths
pos_output_file <- paste(substr(output_file, 0, regexpr(pattern = ".", 
                                                        text = output_file, fixed = TRUE)[1]-1),
                         "_pos.bed", sep="")
neg_output_file <- paste(substr(output_file, 0, regexpr(pattern = ".", 
                                                        text = output_file, fixed = TRUE)[1]-1),
                         "_neg.bed", sep="")

# load data
CpG_locs <- read.csv(input_file, header=F, col.names=c("chr", "start", "end", "weight"), sep='\t')

# extend by window size
if (range){
  
} else {
  extension_bases <- round(window_size/2)
  CpG_locs[,"start"] <- CpG_locs[,"start"] - extension_bases
  CpG_locs[,"end"] <- CpG_locs[,"end"] + extension_bases
  # write data to file
  write.table(CpG_locs, file=output_file, sep="\t", col.names = F, row.names = F, quote = F)
}

# get positively CpGs in separate dataframe
CpGs_pos <- CpG_locs %>% filter(weight > 0)
# write .bed file
write.table(CpGs_pos, file=pos_output_file, sep="\t", col.names = F, row.names = F, quote = F)

# get positively CpGs in separate dataframe
CpGs_neg <- CpG_locs %>% filter(weight < 0)
# write .bed file
write.table(CpGs_neg, file=neg_output_file, sep="\t", col.names = F, row.names = F, quote = F)


