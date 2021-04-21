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

# 1) Load G4 data and catenate both strands
G4_locs <- load_G4_data(G4_file_plus, G4_file_minus)

# 2) Load CpG data and convert from .csv to .bed format
CpG_locs <- load_CpG_data(CpG_file = CpG_file)

# 3) Lift over coordinates 
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)

# 4) Analyse window size
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes)

# 5) Plot results 
plot_results(overlap_results=results, query_name="CpGs", search_set_name="G4s", figure_name_window_size = "out/test_window_size_vs_enrichtment.pdf")

