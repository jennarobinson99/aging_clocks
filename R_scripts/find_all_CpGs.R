setwd("D:/proj_epigen/aging_clocks")
library(tidyverse)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

# load genome
genome <- BSgenome.Hsapiens.UCSC.hg38
chromosomes <-seqnames(genome)[0:24]
pattern <- "CG"

# get number of CpGs in all chromosomes and find genomic coordinates
CpG_occurence <- vector(mode='list', length=25)
names(CpG_occurence) <- c(chromosomes, 'total')
total <- 0
CpG_coordinates <- tibble()

# loop through chromosomes to get CpGs for each one
for (chr_name in chromosomes) {
  current_chr <- genome[[chr_name]]
  #count CpG occurence and save in named list
  n <- countPattern(pattern, current_chr)
  CpG_occurence[chr_name] <- n
  total <- total + n
  #match pattern to get genomic coordinates of each CpG
  seqs <- matchPattern(pattern, current_chr) #gets XStringView object of coordinates
  seqs <- as(seqs, "IRanges") #convert object to IRanges object
  # get start and end coordinates, save in temporary tibble, then combine with
  # data from other chromosomes
  starts <- start(seqs) 
  ends <- end(seqs)
  cur_chr_tibble <- tibble(chr=chr_name, start=starts, end=ends)
  CpG_coordinates <- rbind(CpG_coordinates, cur_chr_tibble)
}
CpG_occurence[[25]] <- total

#apply window to CpG locations and save in extended CpG coordinates tibble
extension_bases <- 50
CpG_coordinates_ext <- CpG_coordinates %>% mutate(end=start+extension_bases, start=start-extension_bases)

# write data to file
# write coordinates to .bed file
write.table(CpG_coordinates, file="CpG_lists/bed_files/global_CpGs.bed", sep="\t", col.names = F, row.names = F, quote = F)
# write CpG occurances into .csv file
write.table(CpG_occurence, file="CpG_lists/csv_files/global_CpGs_occurence.csv", sep=",", col.names = T, row.names = F, quote = F)
# write extended coordinates table (where window function was applied)
write.table(CpG_coordinates_ext, file="CpG_lists/bed_files/global_CpGs_ext.bed", sep="\t", col.names = F, row.names = F, quote = F)
