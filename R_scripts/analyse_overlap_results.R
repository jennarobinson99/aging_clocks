library(tidyverse)
library(scales)

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
# test for correct arguments
if (length(args)!=8) {
  stop("Eight arguments need to be supplied: FE_file CpG_bed_file G4_bed_file CGI_map_file figure_name_g4 figure_name_CpGs figure_name_CpG_distribution figure_name_CGI_context", call.=TRUE)
} else {
  # define all file names and parameters
  FE_file <- args[1]
  CpG_bed_file <- args[2]
  G4_bed_file <- args[3]
  CGI_map_file <- args[4]
  figure_name_G4 <- args[5]
  figure_name_CpGs <- args[6]
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
overlap_data <- read.csv(FE_file, header = TRUE, sep = ";")
# plot G4s
ggplot(overlap_data) +
  geom_point(aes(x = window_size, y= FE_G4s)) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of G4s in AC CpG regions") + 
  xlab("Basepair window size around CpGs")
ggsave(figure_name_G4)
# plot CpGs
ggplot(overlap_data) + 
  geom_point(aes(x=window_size,y=FE_CpGs)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylab("Fold enrichment of AC CpGs in G4 regions") + 
  xlab("Basepair window size around CpGs")
ggsave(figure_name_CpGs)

# 2. Analyse which G4s and CpGs overlap --> find commonalities? 
CpGs_oI <- read.csv(CpG_bed_file, header=FALSE, sep="\t", col.names = c("chr", "start", "end", "weight"))
G4s_oI <- read.csv(G4_bed_file, header=FALSE, sep="\t", col.names=c("chr", "start", "end", "score"))
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
if (CGI_map_file=="-"){
	print("No CGI map file supplied")
} else {
	intersected_file_name <- "temp/CpG_CGI_intersect.bed"
	command <- paste("bedtools intersect -wa -a", CpG_bed_file, "-b", CGI_map_file, ">", intersected_file_name, sep=" ")
	cat(command, "\n")
	try(system(command))
	CpGs_in_CGI <- read.csv(intersected_file_name, header=F, sep="\t", col.names=c("chr", "start", "end", "weight", "correlation"))
	CpGs_unmethylated <- CpGs_oI %>% filter(weight %in% CpGs_in_CGI[,"weight"]) %>% mutate(CGI = "within CGI")
	CpGs_methylated <- CpGs_oI %>% filter(!(weight %in% CpGs_in_CGI[,"weight"])) %>% mutate(CGI = "outside CGI")
	CpGs_oI <- rbind(CpGs_unmethylated, CpGs_methylated) 
	write.table(CpGs_oI, file="out/CpG_info.bed", sep="\t",quote=F)
	# plot CGI context
	ggplot(CpGs_oI) + 
		geom_bar(aes(x=chr)) +
		facet_grid(CGI~.)
	ggsave(figure_name_CGI_context)
	print("Results:")
	sprintf("%.2f percent of CpGs lie within CGIs.", 100* nrow(CpGs_unmethylated)/nrow(CpGs_oI))
	sprintf("%.2f percent of CpGs lie outside of CGIs.", 100 * nrow(CpGs_methylated)/nrow(CpGs_oI)) 
}
