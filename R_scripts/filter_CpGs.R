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
  Fixed window: Rscript --vanilla filter_CpGs.R input_file output_file CGI_map_file window_size \n
  Arguments must be integers."

if (length(args) == 4 ){
 window_size<- as.integer(args[4])
 input_file <- args[1]
 output_file <- args[2]
 CGI_map_file <- args[3]
 range <- FALSE
  
} else {
  stop(message, call.=TRUE)
}
# for interactive mode:
# input_file <- "temp/CpGs_lifted.bed"
# output_file <- "temp/CpGs_ext.bed"
# range <- F
# window_size <- 100
# CGI_map_file <- "-"

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

# get CpGs split according to CGI context 
if (CGI_map_file=="-"){
  print("No CGI map file supplied")
} else {
  intersected_file_name <- "temp/CpG_CGI_intersect.bed"
  command <- paste("bedtools intersect -wa -a", output_file, "-b", CGI_map_file, ">", intersected_file_name, sep=" ")
  cat(command, "\n")
  try(system(command))
  CpGs_in_CGI <- read.csv(intersected_file_name, header=F, sep="\t", col.names=c("chr", "start", "end", "weight", "correlation"))
  CpGs_unmethylated <- CpG_locs %>% filter(weight %in% CpGs_in_CGI[,"weight"]) %>% mutate(CGI = "within CGI")
  CpGs_methylated <- CpG_locs %>% filter(!(weight %in% CpGs_in_CGI[,"weight"])) %>% mutate(CGI = "outside CGI")
  CpGs_CGI <- rbind(CpGs_unmethylated, CpGs_methylated) 
  write.table(CpGs_CGI, file="temp/all_CpGs_CGI_context.bed", sep="\t",col.names=F,row.names=F,quote=F)
  write.table(CpGs_unmethylated, file="temp/CpGs_within.bed", sep="\t",col.names=F,row.names=F, quote=F)
  write.table(CpGs_methylated, file="temp/CpGs_outside.bed", sep="\t",col.names=F,row.names=F,quote=F)
}
  
