library(tidyverse)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)
setwd("D:/proj_epigen/aging_clocks")

# load genome
genome <- BSgenome.Hsapiens.UCSC.hg19
chromosomes <-seqnames(genome)[0:22]
pattern <- "CG"

# get number of CpGs in all chromosomes and find genomic coordinates
CpG_occurence <- vector(mode='list', length=23)
names(CpG_occurence) <- c(chromosomes, 'total')
total <- 0
CpG_coordinates <- tibble()
sequences <- tibble()
seqs_chr <- vector(mode='list', length=22)
#apply window to CpG locations and save in extended CpG coordinates tibble
extension_bases <- 25
first <- 1

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
  starts <- start(seqs) - extension_bases 
  ends <- end(seqs) + extension_bases
  cur_chr_tibble <- tibble(chr=chr_name, start=starts, end=ends)
  CpG_coordinates <- rbind(CpG_coordinates, cur_chr_tibble)
  
  #get sequence
  sequence <- as.list(unname(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr_name, starts, ends, as.character=T)))
  
  # create unique identifiers for each sequence
  final <- first + length(sequence) - 1
  ids <- c(first:final)
  id_string <- paste("Sequence_", as.character(ids), sep="")
  identifiers <- paste(id_string,chr_name,sep="|")
  
  # construct new fasta file with sequences
  write.fasta(sequence, identifiers, "sequences/all_CpG_seq.fasta", open="a", as.string=T)
  
}
CpG_occurence[[23]] <- total

write.table(CpG_coordinates, file="data/CpG_lists/all_CpGs_hg19.bed", sep="\t", col.names = F, row.names = F, quote = F)
