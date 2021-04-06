setwd("D:/proj_epigen/aging_clocks")

#load the bed file containing genome coordinates 
coordinates <- read.csv("out/overlapped_clock_CpGs_pooled.bed", header=F, col.names=c("chr", "start", "ends", "weights"), sep='\t')

#import biostrings and genome
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)
human_genome <- BSgenome.Hsapiens.UCSC.hg19

# get sequence from genome 
chromosomes <-coordinates[,"chr"]
starts <- coordinates[,"start"]
ends <- coordinates[,"ends"]
sequences <- as.list(unname(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chromosomes, starts, ends, as.character=T)))
# create unique identifiers for each sequence
ids <- c(1:length(chromosomes))
id_string <- paste("Sequence_", as.character(ids), sep="")
identifiers <- paste(id_string,chromosomes,sep="|")

# construct new fasta file with sequences
write.fasta(sequences, identifiers, "sequences/overlapped_clock_CpGs_pooled_seqs.fasta", open="w", as.string=T)
  
