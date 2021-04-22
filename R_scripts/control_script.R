# This is a refactored version of all the bash scripts used before. 

# IMPORTS 
library(dplyr)
library(rjson)

# FUNCTIONS
source("R_scripts/functions.R")

# load config.json file 
params <- fromJSON(file="config.json")
CpG_file <- params[["CpG_file"]]
G4_file_plus <- params[["G4_file_plus"]]
G4_file_minus <- params[["G4_file_minus"]]               
chain_file <- params[["chain_file"]]
genome_file <- params[["genome_file"]]
window_size_CpG_dist_plot <- params[["window_size_CpG_dist_plot"]]
all_CpGs_file <- params[["all_CpGs_file"]]
window_sizes <- params[["window_sizes"]]


### Steps 1 to 5 show an overview over the complete pipeline: 

# 1) Load G4 data and catenate both strands
G4_locs <- load_G4_data(G4_file_plus)

# 2) Load CpG data and convert from .csv to .bed format
CpG_locs <- load_CpG_data(CpG_file = CpG_file)

# 3) Lift over coordinates 
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)

# 3.5a) Get number of all CpGs and overlaps with G4s genome-wide (not aging clock, control case) (LONG execution time)
results_all_CpGs <- global_CpG_overlap(all_CpGs_file=all_CpGs_file, G4_locs=G4_locs,genome_file=genome_file, reduce=T)
# 3.5b) Or read data from file (SHORT execution time)
results_all_CpGs <- read.csv(file="out/all_CpGs_props_all_G4s.csv", header=T, sep=";")

# 4) Analyse window size
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)

# 5) Plot results 
plot_results(overlap_results=results, query_name="CpGs", search_set_name="G4s", figure_name_window_size = "out/test_window_size_vs_enrichtment.pdf")


### Miscellaneous specific analyses
# Analyse separately positive and negative correlation coefficients and compare their enrichments
CpGs_pos <- CpG_locs %>% filter(coefficient > 0)
CpGs_neg <- CpG_locs %>% filter(coefficient < 0)
results_pos <- analyse_window_size(query=CpGs_pos, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes)
results_neg <- analyse_window_size(query=CpGs_neg, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes)
write.table(results_pos, file="out/pos_vs_neg/Levine_pos.csv", sep=";", quote=F, row.names = F)
write.table(results_neg, file="out/pos_vs_neg/Levine_neg.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_pos, query_name="positive CpGs", search_set_name="G4s", figure_name_window_size="out/pos_vs_neg/Levine_pos.pdf")
plot_results(overlap_results = results_neg, query_name="negative CpGs", search_set_name="G4s", figure_name_window_size="out/pos_vs_neg/Levine_neg.pdf")

# Analyse BG4 chip seq data from papers (2016 Nat Gen, 2018 Nat S&M Bio)
# 2016 Nat Gen
G4_locs <- load_G4_data(G4_file_plus = G4_file_plus)
CpG_locs <- load_CpG_data(CpG_file)
CpG_locs <- lift_over(coordinates = CpG_locs, chain_file = chain_file)
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes)



plot_results(overlap_results = results, query_name="CpGs", search_set_name="BG4 ChIP-seq data (2016)", figure_name_window_size = "out/BG4_2016_output_R/Horvath.pdf")
write.table(results, file="out/BG4_2018_output_R/Horvath.csv", sep=";", quote=F, row.names = F)
