library(tidyverse)
library(scales)
# Read in the overlapped CpG and G4 data 
FE_file <- "out/Levine_FE.csv"
CpG_bed_file <- "out/CpGs_in_G4s_ws_50.bed"
G4_bed_file <- "out/G4s_in_CpGs_ws_50.bed"

# 1. plot fold enrichment vs. window size 
# load data
overlap_data <- read.csv(FE_file, header = FALSE, sep = ";", col.names = c("input_file","window_size","FE_G4s", "FE_CpGs"))
# plot G4s
ggplot(overlap_data) +
  geom_point(aes(x = window_size, y= FE_G4s)) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of G4s in AC CpG regions") + 
  xlab("Basepair window size around CpGs")
ggsave("out/FE_G4s.pdf")
# plot CpGs
ggplot(overlap_data) + 
  geom_point(aes(x=window_size,y=FE_CpGs)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of AC CpGs in G4 regions") + 
  xlab("Basepair window size around CpGs")
ggsave("out/FE_CpGs.pdf")

# 2. Analyse which G4s and CpGs overlap --> find commonalities? 
CpGs_oI <- read.csv(CpG_bed_file, header=FALSE, sep="\t", col.names = c("chr", "start", "end", "weight"))
G4s_oI <- read.csv(G4_bed_file, header=FALSE, sep="\t", col.names=c("chr", "start", "end", "score"))
# 2.1 CpGs
CpGs_pos <- CpGs_oI %>% filter(weight > 0) %>% mutate(correlation = "Positive")
CpGs_neg <- CpGs_oI %>% filter(weight < 0) %>% mutate(correlation = "Negative")
CpGs_oI <- rbind(CpGs_pos, CpGs_neg)
# plot CpG distribution 
ggplot(CpGs_oI) + 
  geom_bar(aes(x=chr)) + 
  facet_grid(correlation~.)
ggsave("out/CpG_distribution.pdf")
