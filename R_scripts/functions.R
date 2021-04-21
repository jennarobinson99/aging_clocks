# This file contains all functions needed for overlap analysis. 


overlap <- function(query=NULL, search_set=NULL){
  ### A and B should be Ranges objects which can be reduced and overlapped (subset overlap, bedtools -wa mode)
  library(GenomicRanges)
  ranges_A <- makeGRangesFromDataFrame(query, keep.extra.columns = T)
  ranges_B <- makeGRangesFromDataFrame(search_set, keep.extra.columns = T)

  # subset overlap mode: giving subset of query that is overlapping
  overlap_ranges <- subsetByOverlaps(ranges_A, ranges_B) #returns subset of ranges_A which has overlap to ranges_B
  red_overlap_ranges <- GenomicRanges::reduce(overlap_ranges) # merges potential duplicates
  # TODO: need to keep extra columns (e.g. coefficient / score) -> for duplicates, take average
  
  # convert to data frame, return resulting dataframe
  overlap_coor <- GRanges_to_df(ranges=red_overlap_ranges)
  
  return(overlap_coor)
}


analyse_overlap <- function(query = NULL, search_set = NULL, genome_file=NULL, random_trials=5) {
  ### Function overlaps two files and analyses the overlap
  
  library(regioneR)
  library(IRanges)
  library(Biostrings)
  library(stats)
  
  # Overlap the file and get query entries that are in the search set
  overlap_coor <- overlap(query=query, search_set = search_set)
  overlap_ranges <- makeGRangesFromDataFrame(overlap_coor)
  num_overlaps <- dim(overlap_coor)[1]
  
  # Load and process genome file
  genome <- read.csv(genome_file, header = F,sep="\t", col.names = c("chromosome", "length", "unknown"))
  genome <- cbind.data.frame(genome$chromosome, 1, genome$length)
  names(genome) <- c("chromosome", "start", "end")
  
  # Random shuffling overlap 
  rand_overlaps <- vector()
  for(i in 1:random_trials){
    set.seed(i)
    shuffled_coor <- GRanges_to_df(randomizeRegions(query, genome=genome))
    shuffled_overlap <- overlap(query=shuffled_coor, search_set = search_set)
    rand_overlaps[i] <- dim(shuffled_overlap)[1]
  }
  average_rand_overlaps <- mean(rand_overlaps)
  
  # Calculate fold enrichment
  if(average_rand_overlaps == 0) {
    print("Unable to calculate enrichment, zero random overlaps found. ")
  } else {
    enrichment <- num_overlaps / average_rand_overlaps
  }
  
  # Calculate p-value against the whole genome (using hypergeometric dist)
  num_bases <- sum(genome$end)
  num_search_set_bases <- sum(search_set$end - search_set$start)
  num_query_bases <- sum(query$end - query$start)
  num_overlap_bases <- sum(overlap_coor$end - overlap_coor$start)
  p_value <- phyper(q = num_overlap_bases, m = num_search_set_bases, n = (num_bases-num_search_set_bases), k = num_query_bases, lower.tail = F)
  
  return(list("overlap_coordinates" = overlap_coor, "num_overlaps" = num_overlaps, "fold_enrichment" = enrichment, "p_value" = p_value))
}


analyse_window_size <- function(query=NULL, search_set=NULL, genome_file=NULL, window_sizes=NULL){
  ### Function loops through all given window sizes and calculates statistical metrics about distribution and enrichment for all. Outputs dataframe containing all results and produces graphs. 
  
  library(dplyr)
  i <- 1
  num_overlaps <- numeric()
  fold_enrichment <- numeric()
  p_value <- numeric()
  for(window_size in window_sizes){
    query_ext <- query %>% mutate(start=start-round(window_size/2), end=end+round(window_size/2)-1)
    temp_results <- analyse_overlap(query=query_ext, search_set = search_set, genome_file=genome_file)
    num_overlaps[i] <- temp_results["num_overlaps"]
    fold_enrichment[i] <- temp_results["fold_enrichment"]
    p_value[i] <-  temp_results["p_value"]
    i <- i + 1
  }
  
  results <- cbind.data.frame(window_sizes, as.numeric(num_overlaps), as.numeric(fold_enrichment), as.numeric(p_value))
  names(results) <- c("window_size", "num_overlaps", "fold_enrichment", "p_value")
  return(results)
}


plot_results <- function(overlap_results=NULL, query_name="CpGs", search_set_name="G4s", figure_name_window_size="out/window_size_vs_enrichment.pdf") {
  ### Function visualises the results of overlap
  
  library(tidyverse)
  library(scales)
  ylabel <- paste("Fold enrichment of", query_name, "in", search_set_name, sep=" ")
  xlabel <- paste("Window size in bp")
  # plot G4s
  ggplot(overlap_results) +
    geom_point(aes(x = window_size, y = fold_enrichment)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    ylab(ylabel) + 
    xlab(xlabel)
  ggsave(figure_name_window_size)
}


load_CpG_data <- function(CpG_file=NULL){
  ### Function loads csv data into dataframe and returns it
  
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
  return(CpG_locs)
}


lift_over <- function(coordinates=NULL, chain_file=NULL) {
  ### Function lifts over genome coordinates from one assembly to another given a chain file
  
  if (chain_file=="-" | chain_file=="" | chain_file=="/" | chain_file=="\\" | is.null(chain_file) ) {
    #load the bed file containing genome coordinates 
    print("No lift over requested, proceeding without it.")
  } else {
    library(rtracklayer)
    #import chain file
    chain <- import.chain(chain_file)
    #coerce dataframe into GRanges
    ranges <- makeGRangesFromDataFrame(coordinates)
    # use R implementation of Lift Over tool to convert ranges to hg38 coordinate system
    lifted_ranges <- liftOver(ranges, chain)
    # get the object in the right format and write to file
    temp_df <- data.frame(iranges=lifted_ranges)
    lifted_coor <- data.frame(cbind.data.frame(coordinates$chromosome, temp_df$iranges.start, 
                                            temp_df$iranges.end, coordinates$coefficient))
    names(lifted_coor) <- c("chromosome", "start", "end", "coefficient")
    return(lifted_coor)
  }
}


load_G4_data <- function(G4_file_plus=NULL, G4_file_minus=NULL){
  ### Function loads in given G4 data
  library(dplyr)
  G4_plus <- read.csv(G4_file_plus, header=F, sep="\t")
  G4_minus <- read.csv(G4_file_minus, header=F, sep="\t")
  G4_locs <- rbind.data.frame(G4_plus, G4_minus)
  names(G4_locs) <- c("chromosome", "start", "end", "score")
  G4_locs <- G4_locs %>% arrange(chromosome, start, end) 
  G4_locs <- G4_locs %>% filter(chromosome!="chrX" & chromosome!="chrY" & chromosome!="chrM")
  return(G4_locs)
}


GRanges_to_df <- function(ranges=NULL){
  ### Function converts GRanges object to data frame
  df <- data.frame(iranges=ranges)
  df <- data.frame(cbind.data.frame(df$iranges.seqnames, df$iranges.start, df$iranges.end))
  names(df) <- c("chr", "start", "end")
  return(df)
}


chr_string <- function(number) {
  ## Function to convert numeric chromosome number to .bed file chromosome string (e.g.:"chr1")
  library(stringr)
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