library(dplyr)

setwd("D:/proj_epigen/aging_clocks")

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
# define number of bases to include left and right of CpG into range
extension_bases <- 50
# load data
horvath_data <- read.csv("./CpG_lists/csv_files/Horvath_clock_CpGs.csv", header = TRUE, sep = ";")
# drop first row (intercept)
horvath_data <- horvath_data[-1,]
# get chromosomes vector in .bed file format
chr_list <- list(horvath_data[,"Chr"])
chromosomes <- lapply(chr_list, chr_string)
# get start end end positions
starts <-  horvath_data[, "MapInfo"] - extension_bases
ends <-  horvath_data[, "MapInfo"] + extension_bases
# create data frame with .bed file appropriate format
bed_data <-  cbind(data.frame(chromosomes), starts, ends)
# write .bed file
write.table(bed_data, file="CpG_lists/bed_files/horvath_CpGs.bed", sep="\t", col.names = F, row.names = F, quote = F)


# get positively CpGs in separate dataframe
horvath_pos <- horvath_data %>% filter(CoefficientTraining > 0)
# get start end end positions
starts <-  horvath_pos[, "MapInfo"] - extension_bases
ends <-  horvath_pos[, "MapInfo"] + extension_bases
# get chromosomes vector in .bed file format
chr_list <- list(horvath_pos[,"Chr"])
chromosomes <- lapply(chr_list, chr_string)
# create data frame with .bed file appropriate format
bed_data <-  cbind(data.frame(chromosomes), starts, ends)
# write .bed file
write.table(bed_data, file="CpG_lists/bed_files/horvath_CpGs_positive.bed", sep="\t", col.names = F, row.names = F, quote = F)


horvath_neg <- horvath_data %>% filter(CoefficientTraining < 0)
# get start end end positions
starts <-  horvath_neg[, "MapInfo"] - extension_bases
ends <-  horvath_neg[, "MapInfo"] + extension_bases
# get chromosomes vector in .bed file format
chr_list <- list(horvath_neg[,"Chr"])
chromosomes <- lapply(chr_list, chr_string)
# create data frame with .bed file appropriate format
bed_data <-  cbind(data.frame(chromosomes), starts, ends)
# write .bed file
write.table(bed_data, file="CpG_lists/bed_files/horvath_CpGs_negative.bed", sep="\t", col.names = F, row.names = F, quote = F)