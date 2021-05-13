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
figure_name_chromatin_dist <- params[["figure_name_chromatin_dist"]]
all_CpGs_file <- params[["all_CpGs_file"]]
window_sizes <- params[["window_sizes"]]
CGI_map_file <- params[["CGI_map_file"]]
all_CpGs_props_file <- params[["all_CpGs_props_file"]]
output_file <- params[["output_file"]]


### Steps 1 to 5 show an overview over the complete pipeline: 

# 1) Load G4 data and catenate both strands
G4_locs <- load_G4_data(G4_file_plus=G4_file_plus, narrow_peak = T)

# 2) Load ATAC-seq data if desired
ATAC_data <- load_ATAC_data(ATAC_file = ATAC_file)

# 2) Load CpG data and convert from .csv to .bed format
CpG_locs <- load_CpG_data(CpG_file = CpG_file)

# 3) Lift over coordinates 
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)

# 4) Annotate CpGs according to CGI context 
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file, ATAC_map_file = ATAC_file)

# 4.5a) Get number of all CpGs and overlaps with G4s genome-wide (not aging clock, control case) (LONG execution time)
#results_all_CpGs <- global_CpG_overlap(all_CpGs_file=all_CpGs_file, G4_locs=G4_locs,genome_file=genome_file, reduce=F, window_sizes=window_sizes)
#write.table(results_all_CpGs, file = "out/all_CpGs_props/all_CpGs_props_all_G4s_mouse.csv", sep=";", row.names=F)

# 4.5b) Or read data from file (SHORT execution time)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")

# 5) Analyse window size
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)


# 6) Plot results 
plot_results(overlap_results=results, query_name="CpGs", search_set_name="G4s", 
             window_size=window_size_CpG_dist_plot, 
             #figure_name_window_size = figure_name_window_size,
             #figure_name_coeff_dist = figure_name_coefficient_dist,
             #figure_name_CGI_dist = figure_name_CGI_dist, 
             figure_name_chromatin_dist = figure_name_chromatin_dist)

# 7) Write result statistics to file
write.table(results$stats, file=output_file, sep=";", row.names = F)

################################################################################
################################################################################
### Miscellaneous specific analyses

# Check CGI enrichment at G4s and reverse
CGIs <- read.csv("data/CGI_maps/hg19_CGI_map.bed", sep="\t", header = F)
names(CGIs) <- c("chromosome", "start", "end", "id", "CG_contenct", "score")
CGIs <- CGIs %>% select(chromosome, start, end)
CGIs <- CGIs %>% filter(!(chromosome=="chr11_gl000202_random"))
G4s <- load_G4_data(G4_file_plus="data/G4_maps/G4s_human_plus.bed", G4_file_minus="data/G4_maps/G4s_human_minus.bed")
#CGIs in G4s
res_CGIs <- analyse_overlap(query=CGIs, search_set=G4s, genome_file="data/genome_files/hg19_chromInfo.txt")
stats_CGI <- list("direction"="CGIs_in_G4s", "num_overlap"=res_CGIs$num_overlaps, "fold_enrichment"=res_CGIs$fold_enrichment, "p_value_bp"=res_CGIs$p_value_bp)
stats_CGI <- data.frame(stats_CGI)
# G4s in CGIs
res_G4s <- analyse_overlap(query=G4s, search_set=CGIs, genome_file="data/genome_files/hg19_chromInfo.txt")
stats_G4s <- list("direction"="G4s_in_CGIs", "num_overlap"=res_G4s$num_overlaps, "fold_enrichment"=res_G4s$fold_enrichment, "p_value_bp"=res_G4s$p_value_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_CGI, stats_G4s)
# write to file
write.table(stats_all, file="out/CGIs_G4s_stats.csv", sep=";", row.names=F)

################################################################################
# Enrichment of G4s in all CpGs genome-wide vs. G4s in AC CpGs 
G4_locs <- load_G4_data(G4_file_plus = G4_file_plus, G4_file_minus=G4_file_minus)
CpG_locs <- read.csv(all_CpGs_file, header=F, sep="\t", col.names = c("chromosome", "start", "end"))
#CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
# search_set needs to be CpG_locs, cant shuffle such large files
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, shuffle="s")
write.table(results$stats, file=output_file, sep=";", quote=F, row.names = F)

################################################################################
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

G4_locs <- load_G4_data(G4_file_plus, narrow_peak = F)
CpG_locs <- load_CpG_data(CpG_file = CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file)
CpGs_pos <- CpG_locs %>% filter(coefficient > 0)
CpGs_neg <- CpG_locs %>% filter(coefficient < 0)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
results_pos <- analyse_window_size(query=CpGs_pos, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
results_neg <- analyse_window_size(query=CpGs_neg, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
write.table(results_pos$stats, file="out/pos_vs_neg/gw_G4s/Levine_pos_stats.csv", sep=";", quote=F, row.names = F)
write.table(results_neg$stats, file="out/pos_vs_neg/gw_G4s/Levine_neg_stats.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_pos, query_name="positive CpGs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/pos_vs_neg/gw_G4s/Levine_pos_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/pos_vs_neg/gw_G4s/redundant.pdf",
             figure_name_CGI_dist="out/pos_vs_neg/gw_G4s/Levine_pos_CGI_context_dist.pdf")
plot_results(overlap_results = results_neg, query_name="negative CpGs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/pos_vs_neg/gw_G4s/Levine_neg_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/pos_vs_neg/gw_G4s/redundant.pdf",
             figure_name_CGI_dist="out/pos_vs_neg/gw_G4s/Levine_neg_CGI_context_dist.pdf")

################################################################################
# Analyse CGI contexts separately (human DNAm clocks)
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

G4_locs <- load_G4_data(G4_file_plus, narrow_peak = F)
CpG_locs <- load_CpG_data(CpG_file = CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file)
CpGs_in <- CpG_locs %>% filter(CGI_context=="inside")
CpGs_out <- CpG_locs %>% filter(CGI_context=="outside")
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
results_in <- analyse_window_size(query=CpGs_in, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
results_out <- analyse_window_size(query=CpGs_out, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs)
write.table(results_in$stats, file="out/CGI_in_vs_out/BG4s_K562/Levine_in_stats.csv", sep=";", quote=F, row.names = F)
write.table(results_out$stats, file="out/CGI_in_vs_out/BG4s_K562/Levine_out_stats.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_in, query_name="CpGs in CGIs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/CGI_in_vs_out/BG4s_K562/Levine_in_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/CGI_in_vs_out/BG4s_K562/Levine_in_coeff_dist.pdf",
             figure_name_CGI_dist="out/CGI_in_vs_out/BG4s_K562/redundant.pdf")
plot_results(overlap_results = results_out, query_name="CpGs outside CGIs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_window_size="out/CGI_in_vs_out/BG4s_K562/Levine_out_ws_vs_enrichment.pdf",
             figure_name_coeff_dist="out/CGI_in_vs_out/BG4s_K562/Levine_out_coeff_dist.pdf",
             figure_name_CGI_dist="out/CGI_in_vs_out/BG4s_K562/redundant.pdf")

################################################################################
# Analyse ATAC-seq data for cell lines
ATAC_data <- load_ATAC_data(ATAC_file=ATAC_file)
CpG_locs <- load_CpG_data(CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file = CGI_map_file, ATAC_map_file = ATAC_file)
# split CpGs in open and closed chromatin groups 
CpG_open <- CpG_locs %>% filter(chromatin=="open")
CpG_closed <- CpG_locs %>% filter(chromatin=="closed")
G4_locs <- load_G4_data(G4_file_plus=G4_file_plus, narrow_peak=F)
all_CpGs_props <- read.csv(all_CpGs_props_file, sep=";")
# Overlap clock CpGs with G4s
results_open_CpGs <- analyse_window_size(query=CpG_open, search_set=G4_locs, genome_file=genome_file,all_CpGs_props = all_CpGs_props, window_sizes=window_sizes)
plot_results(overlap_results=results_open_CpGs, query_name="CpGs", search_set_name="G4s",
             window_size = window_size_CpG_dist_plot, 
             figure_name_window_size = "out/chromatin_open_vs_closed/BG4_K562/Levine_open_enrichment_vs_ws.pdf", 
             figure_name_CGI_dist = "out/test_CGI.pdf", 
             figure_name_coeff_dist = "out/test_coeff.pdf", 
             figure_name_chromatin_dist = "out/test_chromatin.pdf")
stats_open <- results_open_CpGs$stats
write.table(stats_open, file="out/chromatin_open_vs_closed/BG4_K562/Levine_open_stats.csv", sep=";", quote=F, row.names=F)

results_closed_CpGs <- analyse_window_size(query=CpG_closed, search_set=G4_locs, genome_file=genome_file,all_CpGs_props = all_CpGs_props, window_sizes=window_sizes)
plot_results(overlap_results=results_closed_CpGs, query_name="CpGs", search_set_name="G4s",
             window_size = window_size_CpG_dist_plot, 
             figure_name_window_size = "out/chromatin_open_vs_closed/BG4_K562/Levine_closed_enrichment_vs_ws.pdf", 
             figure_name_CGI_dist = "out/test_CGI.pdf", 
             figure_name_coeff_dist = "out/test_coeff.pdf", 
             figure_name_chromatin_dist = "out/test_chromatin.pdf")
stats_closed <- results_closed_CpGs$stats
write.table(stats_closed, file="out/chromatin_open_vs_closed/BG4_K562/Levine_closed_stats.csv", sep=";", quote=F, row.names=F)


# Overlap G4 and ATAC-seq data
overlap_G4_ATAC <- analyse_overlap(query=G4_locs, search_set = ATAC_data, genome_file = genome_file)
overlap_CpG_ATAC <- analyse_window_size(query = CpG_locs, search_set=ATAC_data, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props = all_CpGs_props)
overlap_G4_ATAC_hacat <- overlap_G4_ATAC
overlap_CpG_ATAC_hacat <- overlap_CpG_ATAC

################################################################################
# Analyse TET enzyme overlap with G4s in human and mouse (genome-wide G4s)
# load TET data 
TET_data <- read.csv("data/TET_data/TET2_mouse_peaks.bed", sep="\t", header=F)
names(TET_data) <- c("chromosome", "start", "end", ".", "..", "strand")
TET_data <- TET_data %>% select(chromosome, start, end) %>% filter(chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")) 
TET_data <- lift_over(coordinates=TET_data, chain_file="data/chain_files/mm9ToMm10.over.chain")
# load G4s 
G4_locs <- load_G4_data(G4_file_plus="data/G4_maps/G4s_mouse_plus.bed", G4_file_minus="data/G4_maps/G4s_mouse_minus.bed", autosomes=T)
# filter chrM out
G4_locs <- G4_locs %>% filter(!(chromosome=="chrM"))
# Overlap TET peaks and G4s 
TETs_in_G4s <- analyse_overlap(query=TET_data, search_set=G4_locs, genome_file="data/genome_files/mm10_chromInfo.txt", shuffle = "s")
G4s_in_TETs <- analyse_overlap(query=G4_locs, search_set=TET_data, genome_file="data/genome_files/mm10_chromInfo.txt", shuffle="q")
# construct output dataframe to write to file
stats_TETs <- list("direction"="TETs_in_G4s", "num_overlap"=TETs_in_G4s$num_overlaps, "fold_enrichment"=TETs_in_G4s$fold_enrichment, "p_value_bp"=TETs_in_G4s$p_value_bp)
stats_TETs <- data.frame(stats_TETs)
stats_G4s <- list("direction"="G4s_in_TETs", "num_overlap"=G4s_in_TETs$num_overlaps, "fold_enrichment"=G4s_in_TETs$fold_enrichment, "p_value_bp"=G4s_in_TETs$p_value_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_TETs, stats_G4s)
write.table(stats_all, file="out/TET_analysis/mouse_TET2_all_G4_stats.csv", sep=";", row.names = F)

################################################################################
# DNMT1 / DNMT3b analysis
DNMT_data <- read.csv("data/DNMT_data/DNMT1_human_HepG2.bed", sep="\t", header=F)
names(DNMT_data) <- c("chromosome", "start", "end", ".", "..", "strand")
DNMT_data <- DNMT_data %>% select(chromosome, start, end) %>% filter(chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")) 
# load G4s 
G4_locs <- load_G4_data(G4_file_plus="data/G4_maps/G4s_human_plus.bed", G4_file_minus="data/G4_maps/G4s_human_minus.bed", autosomes=F)
G4_locs <- lift_over(coordinates=G4_locs, chain_file="data/chain_files/hg19ToHg38.over.chain")
# filter chrM out
G4_locs <- G4_locs %>% filter(!(chromosome %in% c("chrX", "chrY", "chrM")))
# Overlap TET peaks and G4s 
DNMTs_in_G4s <- analyse_overlap(query=DNMT_data, search_set=G4_locs, genome_file="data/genome_files/hg19_chromInfo.txt", shuffle = "s")
G4s_in_DNMTs <- analyse_overlap(query=G4_locs, search_set=DNMT_data, genome_file="data/genome_files/hg19_chromInfo.txt", shuffle="q")
# construct output dataframe to write to file
stats_DNMTs <- list("direction"="DNMTs_in_G4s", "num_overlap"=DNMTs_in_G4s$num_overlaps, "fold_enrichment"=DNMTs_in_G4s$fold_enrichment, "p_value_bp"=DNMTs_in_G4s$p_value_bp)
stats_DNMTs <- data.frame(stats_DNMTs)
stats_G4s <- list("direction"="G4s_in_DNMTs", "num_overlap"=G4s_in_DNMTs$num_overlaps, "fold_enrichment"=G4s_in_DNMTs$fold_enrichment, "p_value_bp"=G4s_in_DNMTs$p_value_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_DNMTs, stats_G4s)
write.table(stats_all, file="out/DNMT_analysis/DNMT1_HepG2_all_G4_stats_shuffled_G4s.csv", sep=";", row.names = F)
