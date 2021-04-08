setwd("D:/proj_epigen/aging_clocks")

#load the bed file containing genome coordinates
input_file <- "data/CpG_lists/all_CpGs_mm10_ws_1000_merged.bed"
coordinates <- read.csv(input_file, header=F, col.names=c("chr", "start", "end"), sep='\t')

#import biostrings and genome
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(rGADEM)

# get sequence from genome
chromosomes <-coordinates[,"chr"]
starts <- coordinates[,"start"]
ends <- coordinates[,"end"]
sequences <- as.list(unname(Biostrings::getSeq(genome, chromosomes, starts, ends, as.character=T)))
# create unique identifiers for each sequence
ids <- c(1:length(chromosomes))
id_string <- paste("Sequence_", as.character(ids), sep="")
identifiers <- paste(id_string,chromosomes,sep="|")

# construct new fasta file with sequences
write.fasta(sequences, identifiers, "sequences/all_CpGs_human.fasta", open="w", as.string=T)

# run motif discovery
ranges <- makeGRangesFromDataFrame(coordinates)
gadem <- GADEM(ranges, verbose=1, genome=Mmusculus)
