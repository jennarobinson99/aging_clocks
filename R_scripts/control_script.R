# This is a refactored version of all the bash scripts used before. 

# IMPORTS 
library(dplyr)
library(rjson)
library(stringr)

# FUNCTIONS
source("R_scripts/functions.R")

# load config.json file 
params <- fromJSON(file="config.json")
CpG_file <- params[["CpG_file"]]
G4_file_plus <- params[["G4_file_plus"]]
G4_file_minus <- params[["G4_file_minus"]]               
chain_file <- params[["chain_file"]]
genome_file <- params[["genome_file"]]
window_size_CpG_dist_plot <- params[["window_size_CpG_dist_plot"]]
all_CpGs_file <- params[["all_CpGs_file"]]
window_sizes <- params[["window_sizes"]]
window_size <- window_sizes[1]

# 1) Convert CpG data from .csv to .bed format
# load data
CpG_locs <- read.csv(CpG_file, header = TRUE, sep = ";")
# get chromosomes vector in .bed file format
chr_list <- list(CpG_locs[,"Chr"])
chromosomes <- lapply(chr_list, chr_string)
# get start end end positions
starts <-  CpG_locs[, "MapInfo"]
ends <-  CpG_locs[, "MapInfo"] + 1
coefficients <- CpG_locs[,"CoefficientTraining"]
# create data frame with .bed file appropriate format
CpG_locs <-  cbind(data.frame(chromosomes), starts, ends, coefficients)
names(CpG_locs) <- c("chromosome", "start", "end", "coefficient")


# 2) Lift over coordinates 
if (chain_file=="-" | chain_file=="" | chain_file=="/" | chain_file=="\\") {
  #load the bed file containing genome coordinates 
  print("No lift over requested, proceeding without it.")
} else {
  library(rtracklayer)
  #import chain file
  chain <- import.chain(chain_file)
  #coerce dataframe into GRanges
  ranges <- makeGRangesFromDataFrame(CpG_locs)
  # use R implementation of Lift Over tool to convert ranges to hg38 coordinate system
  lifted_ranges <- liftOver(ranges, chain)
  # get the object in the right format and write to file
  temp_df <- data.frame(iranges=lifted_ranges)
  CpG_locs <- data.frame(cbind.data.frame(CpG_locs$chromosome, temp_df$iranges.start, 
                                                  temp_df$iranges.end, CpG_locs$coefficient))
  names(CpG_locs) <- c("chromosome", "start", "end", "coefficient")
}


# 3) Extend bases
CpG_locs <- CpG_locs %>% mutate(start=start-(window_size/2), end=end+(window_size/2)-1)


# 4) Load G4 data and catenate both strands
G4_plus <- read.csv(G4_file_plus, header=F, sep="\t")
G4_minus <- read.csv(G4_file_minus, header=F, sep="\t")
G4_locs <- rbind.data.frame(G4_plus, G4_minus)
names(G4_locs) <- c("chromosome", "start", "end", "score")
G4_locs <- G4_locs %>% arrange(chromosome, start, end)
