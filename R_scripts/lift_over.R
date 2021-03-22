#Usage: Rscript --vanilla lift_over.R input_file output_file chain_file

library(rtracklayer)

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

# MAIN 

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
# test for correct arguments
if (length(args)!=3) {
  stop("Three arguments need to be supplied: input_file output_file chain_file", call.=TRUE)
} else {
  # define all file names and parameters
  input_file <- args[1]
  output_file <- args[2]
  chain_file <- args[3]
}

#input_file <- "temp/CpGs_locs.bed"
#output_file <- "CpGs_lifted.bed"
#chain_file <- "data/chain_files/hg18ToHg38.over.chain"

chain <- import.chain(chain_file)

#load the bed file containing genome coordinates 
coordinates <- read.csv(input_file, header=F, col.names=c("chr", "start", "end", "weight"), sep='\t')

#coerce dataframe into GRanges
ranges <- makeGRangesFromDataFrame(coordinates)

# use R implementation of Lift Over tool to convert ranges to hg38 coordinate system
hg38_ranges <- liftOver(ranges, chain)
# get the object in the right format and write to file
hg38_coordinates <- data.frame(iranges=hg38_ranges)
hg38_coordinates <- data.frame(cbind.data.frame(coordinates$chr, hg38_coordinates$iranges.start, 
                          hg38_coordinates$iranges.end, coordinates$weight))
names(hg38_coordinates) <- c("chr", "start", "end", "weight")
write.table(hg38_coordinates, file=output_file, sep="\t", col.names = F, row.names = F, quote = F)
