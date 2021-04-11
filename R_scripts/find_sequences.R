setwd("D:/proj_epigen/aging_clocks")

#load the bed file containing genome coordinates
input_file <- "data/CpG_lists/all_CpGs_hg19_ws_50_merged.bed"
coordinates <- read.csv(input_file, header=F, col.names=c("chr", "start", "end"), sep='\t')

#import biostrings and genome
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)
genome <- BSgenome.Hsapiens.UCSC.hg19

# random subsampling to reduce amount of sequences
length <-  dim(coordinates)[[1]]
rand_idx <- sample(1:length, 500000)
rand_sample <- coordinates[rand_idx,]

# reducing doesnt seem to effect the sample -> already non-overlapping
# sample_ranges <- makeGRangesFromDataFrame(rand_sample, keep.extra.columns = T)
# sample_ranges <- reduce(sample_ranges)

# get sequence from genome
chromosomes <-rand_sample[,"chr"]
starts <- rand_sample[,"start"]
ends <- rand_sample[,"end"]
sequences <- as.list(unname(Biostrings::getSeq(genome, chromosomes, starts, ends, as.character=T)))
# create unique identifiers for each sequence
ids <- c(1:length(chromosomes))
id_string <- paste("Sequence_", as.character(ids), sep="")
identifiers <- paste(id_string,chromosomes,sep="|")
  
# construct new fasta file with sequences
write.fasta(sequences, identifiers, "sequences/all_CpGs_human_subsampled.fasta", open="w", as.string=T)
