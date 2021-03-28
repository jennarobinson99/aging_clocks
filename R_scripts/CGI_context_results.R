library(tidyverse)
library(scales)

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
# test for correct arguments
if (length(args)!=8) {
  stop("Eight arguments need to be supplied: FE_file CpG_within_file G4_outside_file CGI_map_file figure_name_within figure_name_outside figure_name_CpG_distribution figure_name_CGI_context", call.=TRUE)
} else {
  # define all file names and parameters
  FE_file <- args[1]
  CpG_within_file <- args[2]
  CpG_outside_file <- args[3]
  CGI_map_file <- args[4]
  figure_name_within <- args[5]
  figure_name_outside <- args[6]
  figure_name_chr_dist <- args[7]
  figure_name_CGI_context <- args[8]
}

# # Read in the overlapped CpG and G4 data 
# FE_file <- "out/Levine_FE.csv"
# CpG_bed_file <- "out/CpGs_in_G4s_ws_50.bed"
# G4_bed_file <- "out/G4s_in_CpGs_ws_50.bed"
# figure_name_G4 <- "out/FE_G4s.pdf"
# figure_name_CpGs <- "out/FE_CpGs.pdf"
# figure_name_chr_dist <- "out/CpG_distribution.pdf"

# 1. plot fold enrichment vs. window size 
# load data
overlap_data <- read.csv(FE_file, header = FALSE, sep = ";", col.names = c("input_file","window_size","FE_within", "FE_outside"))
# plot G4s
ggplot(overlap_data) +
  geom_point(aes(x = window_size, y= FE_within)) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of CGI-associated CpGs at G4s") + 
  xlab("Basepair window size around CpGs")
ggsave(figure_name_G4)
# plot CpGs
ggplot(overlap_data) + 
  geom_point(aes(x=window_size,y=FE_outside)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of non-CGI CpGs at G4s") + 
  xlab("Basepair window size around CpGs")
ggsave(figure_name_CpGs)

# Read in CpG data
CpGs_within <- read.csv(CpG_within_file, header=FALSE, sep="\t", col.names = c("chr", "start", "end", "weight"))
CpGs_within <- CpGs_within %>% mutate(CGI_context="inside CGI")
CpGs_outside <- read.csv(CpG_outside_file, header=FALSE, sep="\t", col.names = c("chr", "start", "end", "weight"))
CpGs_outside <- CpGs_outside %>% mutate(CGI_context="outside CGI")
CpGs_oI <- rbind(CpGs_within, CpGs_outside)
# 2. Analyse which CpGs overlap --> find commonalities?
CpGs_pos <- CpGs_oI %>% filter(weight > 0) %>% mutate(correlation = "Positive")
CpGs_neg <- CpGs_oI %>% filter(weight < 0) %>% mutate(correlation = "Negative")
CpGs_oI <- rbind(CpGs_pos, CpGs_neg)
sprintf("%.2f percent of CpGs correlate positively.", 100*nrow(CpGs_pos)/nrow(CpGs_oI))
sprintf("%.2f percent of CpGs correlate negatively.", 100*nrow(CpGs_neg)/nrow(CpGs_oI))
#plot the distribution
ggplot(CpGs_oI) +
  geom_bar(aes(x=chr)) + 
  facet_grid(correlation~.)
ggsave(figure_name_chr_dist) 

# 3. Check which CpGs lie within CGIs
ggplot(CpGs_oI) + 
  geom_bar(aes(x=chr)) +
  facet_grid(CGI_context~.)
ggsave(figure_name_CGI_context)
print("Results:")
sprintf("%.2f percent of CpGs lie within CGIs.", 100* nrow(CpGs_unmethylated)/nrow(CpGs_oI))
sprintf("%.2f percent of CpGs lie outside of CGIs.", 100 * nrow(CpGs_methylated)/nrow(CpGs_oI)) 
