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
ATAC_file <- params[["ATAC_file"]]
chain_file <- params[["chain_file"]]
genome_file <- params[["genome_file"]]
window_size_CpG_dist_plot <- params[["window_size_CpG_dist_plot"]]
figure_name_window_size <- params[["figure_name_window_size"]]
figure_name_coefficient_dist <- params[["figure_name_coefficient_dist"]]
figure_name_CGI_dist <- params[["figure_name_CGI_dist"]]
all_CpGs_file <- params[["all_CpGs_file"]]
window_sizes <- params[["window_sizes"]]
CGI_map_file <- params[["CGI_map_file"]]
all_CpGs_props_file <- params[["all_CpGs_props_file"]]
output_file <- params[["output_file"]]


### Steps 1 to 5 show an overview over the complete pipeline: 

# 1) Load G4 data and catenate both strands
G4_locs <- load_G4_data(G4_file_plus, narrow_peak = T)

# 2) Load CpG data and convert from .csv to .bed format
CpG_locs <- load_CpG_data(CpG_file = CpG_file)

# 3) Lift over coordinates 
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)

# 4) Annotate CpGs according to CGI context 
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file)

# 4.5a) Get number of all CpGs and overlaps with G4s genome-wide (not aging clock, control case) (LONG execution time)
# results_all_CpGs <- global_CpG_overlap(all_CpGs_file=all_CpGs_file, G4_locs=G4_locs,genome_file=genome_file, reduce=T)
# 4.5b) Or read data from file (SHORT execution time)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")

# 5) Analyse window size
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)


# 6) Plot results 
plot_results(overlap_results=results, query_name="CpGs", search_set_name="G4s", 
             window_size=window_size_CpG_dist_plot, 
             figure_name_window_size = figure_name_window_size,
             figure_name_coeff_dist = figure_name_coefficient_dist,
             figure_name_CGI_dist = figure_name_CGI_dist)

# 7) Write result statistics to file
write.table(results$stats, file=output_file, sep=";", row.names = F)

################################################################################

### Miscellaneous specific analyses


# Analyse separately positive and negative correlation coefficients and compare their enrichments
# load config.json file 
params <- fromJSON(file="config.json")
CpG_file <- params[["CpG_file"]]
G4_file_plus <- params[["G4_file_plus"]]
G4_file_minus <- params[["G4_file_minus"]]               
chain_file <- params[["chain_file"]]
genome_file <- params[["genome_file"]]
CGI_map_file <- params[["CGI_map_file"]]
window_size_CpG_dist_plot <- params[["window_size_CpG_dist_plot"]]
window_sizes <- params[["window_sizes"]]
all_CpGs_file <- params[["all_CpGs_file"]]
all_CpGs_props_file <- params[["all_CpGs_props_file"]]

G4_locs <- load_G4_data(G4_file_plus, narrow_peak = T)
CpG_locs <- load_CpG_data(CpG_file = CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file)
CpGs_pos <- CpG_locs %>% filter(coefficient > 0)
CpGs_neg <- CpG_locs %>% filter(coefficient < 0)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
results_pos <- analyse_window_size(query=CpGs_pos, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
results_neg <- analyse_window_size(query=CpGs_neg, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
write.table(results_pos$stats, file="out/pos_vs_neg/BG4s_K562/Levine_pos_stats.csv", sep=";", quote=F, row.names = F)
write.table(results_neg$stats, file="out/pos_vs_neg/BG4s_K562/Levine_neg_stats.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_pos, query_name="positive CpGs", search_set_name="BG4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/pos_vs_neg/BG4s_K562/Levine_pos_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/pos_vs_neg/BG4s_K562/redundant.pdf",
             figure_name_CGI_dist="out/pos_vs_neg/BG4s_K562/Levine_pos_CGI_context_dist.pdf")
plot_results(overlap_results = results_neg, query_name="negative CpGs", search_set_name="BG4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/pos_vs_neg/BG4s_K562/Levine_neg_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/pos_vs_neg/BG4s_K562/redundant.pdf",
             figure_name_CGI_dist="out/pos_vs_neg/BG4s_K562/Levine_neg_CGI_context_dist.pdf")

# Analyse BG4 chip seq data from papers (2016 Nat Gen, 2018 Nat S&M Bio)
# 2016 Nat Gen
G4_locs <- load_G4_data(G4_file_plus = G4_file_plus, narrow_peak = T)
CpG_locs <- load_CpG_data(CpG_file)
CpG_locs <- lift_over(coordinates = CpG_locs, chain_file = chain_file)
results_all_CpGs <- read.csv(file="out/all_CpGs_props/all_CpGs_props_BG4_HEK_rep1.csv", header=T, sep=";")
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props = results_all_CpGs)
plot_results(overlap_results = results, query_name="CpGs", search_set_name="BG4 ChIP-seq data (HEK)", figure_name_window_size = "out/BG4_HEK_output/rep1_Horvath.pdf")
write.table(results, file="out/BG4_HEK_output/rep1_Horvath.csv", sep=";", quote=F, row.names = F)


# Analyse ATAC-seq data for HaCaT and HEK cell lines
ATAC_data <- load_ATAC_data(ATAC_file=ATAC_file)
CpG_locs <- load_CpG_data(CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file = CGI_map_file)
G4_locs <- load_G4_data(G4_file_plus=G4_file_plus, narrow_peak=T)
all_CpGs_props <- read.csv(all_CpGs_props_file, sep=";")
# Overlap G4 and ATAC-seq data
overlap_G4_ATAC <- analyse_overlap(query=G4_locs, search_set = ATAC_data, genome_file = genome_file)
overlap_CpG_ATAC <- analyse_window_size(query = CpG_locs, search_set=ATAC_data, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props = all_CpGs_props)
overlap_G4_ATAC_hacat <- overlap_G4_ATAC
overlap_CpG_ATAC_hacat <- overlap_CpG_ATAC
