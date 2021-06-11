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
figure_name_enrichment <- params[["figure_name_enrichment"]]
figure_name_p_value <- params[["figure_name_p_value"]]
figure_name_coefficient_dist <- params[["figure_name_coefficient_dist"]]
figure_name_CGI_dist <- params[["figure_name_CGI_dist"]]
figure_name_chromatin_dist <- params[["figure_name_chromatin_dist"]]
all_CpGs_file <- params[["all_CpGs_file"]]
window_sizes <- params[["window_sizes"]]
CGI_map_file <- params[["CGI_map_file"]]
mask_file <- params[["mask_file"]]
all_CpGs_props_file <- params[["all_CpGs_props_file"]]
output_file <- params[["output_file"]]


### Steps 1 to 5 show an overview over the complete pipeline: 

# 1) Load G4 data and catenate both strands
G4_locs <- load_G4_data(G4_file_plus=G4_file_plus, narrow_peak = F)

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

# 5) Load mask file
mask <- read.csv(mask_file, header=F, sep="\t", col.names = c("chromosome", "start", "end"))

# 5) Analyse window size
results <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30, shuffle = "q",mask = mask)

# 6) Plot results 
plot_results(overlap_results=results, query_name="CpGs", search_set_name="G4s", 
             window_size=window_size_CpG_dist_plot, 
             figure_name_enrichment = figure_name_enrichment,
             figure_name_p_value = figure_name_p_value,
             figure_name_coeff_dist = figure_name_coefficient_dist,
             figure_name_CGI_dist = figure_name_CGI_dist, 
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
res_CGIs <- analyse_overlap(query=CGIs, search_set=G4s, genome_file="data/genome_files/hg19_chromInfo.txt", random_trials = 50, shuffle = "q")
stats_CGI <- list("direction"="CGIs_in_G4s", "num_overlap_subset"=res_CGIs$num_overlap_subset, "num_overlap_bases"=res_CGIs$num_overlap_bases, "FE_MC_subset"=res_CGIs$FE_MC_subset, "p_value_MC_subset"=res_CGIs$p_value_MC_subset, "FE_MC_bp"=res_CGIs$FE_MC_bp, "p_value_MC_bp"=res_CGIs$p_value_MC_bp)
stats_CGI <- data.frame(stats_CGI)
# G4s in CGIs
res_G4s <- analyse_overlap(query=G4s, search_set=CGIs, genome_file="data/genome_files/hg19_chromInfo.txt", random_trials=50, shuffle="s")
stats_G4s<- list("direction"="G4s_in_CGIs", "num_overlap_subset"=res_G4s$num_overlap_subset, "num_overlap_bases"=res_G4s$num_overlap_bases, "FE_MC_subset"=res_G4s$FE_MC_subset, "p_value_MC_subset"=res_G4s$p_value_MC_subset, "FE_MC_bp"=res_G4s$FE_MC_bp, "p_value_MC_bp"=res_G4s$p_value_MC_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_CGI, stats_G4s)
# write to file
write.table(stats_all, file="out/CGIs_G4s_stats.csv", sep=";", row.names=F)

################################################################################
# Enrichment of G4s in all CpGs genome-wide vs. G4s in AC CpGs 
G4_locs <- load_G4_data(G4_file_plus = "data/G4_maps/G4s_human_plus.bed", G4_file_minus="data/G4_maps/G4s_human_minus.bed", narrow_peak = F)
CpG_locs <- read.csv("data/CpG_lists/all_CpGs_hg19.bed", header=F, sep="\t", col.names = c("chromosome", "start", "end"))
results_all_CpGs <- read.csv(file="out/all_CpGs_props/all_CpGs_props_all_G4s.csv", header=T, sep=";")
# search_set needs to be CpG_locs, cant shuffle such large files
results_G4 <- analyse_window_size(query=G4_locs, search_set=CpG_locs, genome_file="data/genome_files/hg19_chromInfo.txt", window_sizes=c(10),  shuffle="q", random_trials = 3, all_CpGs_props = results_all_CpGs)
results_CpG <- analyse_window_size(query=CpG_locs, search_set=G4_locs, genome_file="data/genome_files/hg19_chromInfo.txt", window_sizes=c(10),  shuffle="s", random_trials = 3, all_CpGs_props = results_all_CpGs)
write.table(results$stats, file="out/global_CpGs_in_G4s/all_G4s_in_gw_CpGs.csv", sep=";", quote=F, row.names = F)

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
ATAC_file <- params[["ATAC_file"]]
window_size_CpG_dist_plot <- params[["window_size_CpG_dist_plot"]]
window_sizes <- params[["window_sizes"]]
all_CpGs_file <- params[["all_CpGs_file"]]
all_CpGs_props_file <- params[["all_CpGs_props_file"]]

G4_locs <- load_G4_data(G4_file_plus,G4_file_minus, narrow_peak = F)
CpG_locs <- load_CpG_data(CpG_file = CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file, ATAC_map_file = ATAC_file)
CpGs_pos <- CpG_locs %>% filter(coefficient > 0)
CpGs_neg <- CpG_locs %>% filter(coefficient < 0)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
results_pos <- analyse_window_size(query=CpGs_pos, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30, shuffle = "q")
results_neg <- analyse_window_size(query=CpGs_neg, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30, shuffle="q")
write.table(results_pos$stats, file="out/pos_vs_neg/gw_G4s/Horvath_pos_stats.csv", sep=";", quote=F, row.names = F)
write.table(results_neg$stats, file="out/pos_vs_neg/gw_G4s/Horvath_neg_stats.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_pos, query_name="positive CpGs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_enrichment="out/pos_vs_neg/gw_G4s/Horvath_pos_enrichment.pdf",
             figure_name_p_value = "out/pos_vs_neg/gw_G4s/Horvath_pos_p_value.pdf")
plot_results(overlap_results = results_neg, query_name="negative CpGs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_enrichment="out/pos_vs_neg/gw_G4s/Horvath_neg_enrichment.pdf",
             figure_name_p_value = "out/pos_vs_neg/gw_G4s/Horvath_neg_p_value.pdf")

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
ATAC_map_file <- params[["ATAC_file"]]

G4_locs <- load_G4_data(G4_file_plus, narrow_peak = F)
CpG_locs <- load_CpG_data(CpG_file = CpG_file)
CpG_locs <- lift_over(coordinates=CpG_locs, chain_file=chain_file)
CpG_locs <- annotate_CpGs(CpG_coordinates=CpG_locs, CGI_map_file=CGI_map_file, ATAC_map_file = ATAC_map_file)
CpGs_in <- CpG_locs %>% filter(CGI_context=="inside")
CpGs_out <- CpG_locs %>% filter(CGI_context=="outside")
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
results_in <- analyse_window_size(query=CpGs_in, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30)
results_out <- analyse_window_size(query=CpGs_out, search_set=G4_locs, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props=results_all_CpGs, random_trials = 30)
write.table(results_in$stats, file="out/CGI_in_vs_out/BG4s_K562/Horvath_in_stats.csv", sep=";", quote=F, row.names = F)
write.table(results_out$stats, file="out/CGI_in_vs_out/BG4s_K562/Horvath_out_stats.csv", sep=";", quote=F, row.names = F)
plot_results(overlap_results = results_in, query_name="CpGs in CGIs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_enrichment ="out/CGI_in_vs_out/BG4s_K562/Horvath_in_enrichment.pdf",
             figure_name_p_value = "out/CGI_in_vs_out/BG4s_K562/Horvath_in_p_value.pdf")
plot_results(overlap_results = results_out, query_name="CpGs outside CGIs", search_set_name="G4s",
             window_size=window_size_CpG_dist_plot,
             figure_name_enrichment="out/CGI_in_vs_out/BG4s_K562/Horvath_out_enrichment.pdf",
             figure_name_p_value ="out/CGI_in_vs_out/BG4s_K562/Horvath_out_p_value.pdf")

################################################################################
# Analyse ATAC-seq data for cell lines
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
ATAC_file <- params[["ATAC_file"]]

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
results_open_CpGs <- analyse_window_size(query=CpG_open, search_set=G4_locs, genome_file=genome_file,all_CpGs_props = all_CpGs_props, window_sizes=window_sizes, random_trials = 30)
plot_results(overlap_results=results_open_CpGs, query_name="CpGs", search_set_name="G4s",
             window_size = window_size_CpG_dist_plot, 
             figure_name_enrichment = "out/chromatin_open_vs_closed/BG4_K562/Horvath_open_enrichment.pdf", 
             figure_name_p_value = "out/chromatin_open_vs_closed/BG4_K562/Horvath_open_p_value.pdf")
stats_open <- results_open_CpGs$stats
write.table(stats_open, file="out/chromatin_open_vs_closed/BG4_K562/Horvath_open_stats.csv", sep=";", quote=F, row.names=F)

results_closed_CpGs <- analyse_window_size(query=CpG_closed, search_set=G4_locs, genome_file=genome_file,all_CpGs_props = all_CpGs_props, window_sizes=window_sizes, random_trials = 30)
plot_results(overlap_results=results_closed_CpGs, query_name="CpGs", search_set_name="G4s",
             window_size = window_size_CpG_dist_plot, 
             figure_name_enrichment ="out/chromatin_open_vs_closed/BG4_K562/Horvath_closed_enrichment.pdf",
             figure_name_p_value = "out/chromatin_open_vs_closed/BG4_K562/Horvath_closed_p_value.pdf")
stats_closed <- results_closed_CpGs$stats
write.table(stats_closed, file="out/chromatin_open_vs_closed/BG4_K562/Horvath_closed_stats.csv", sep=";", quote=F, row.names=F)


# Overlap ATAC-seq data with G4s -> enrichment of accessible sites with G4s
overlap_accessible_sites_with_G4s <- analyse_overlap(query=G4_locs, search_set = ATAC_data, genome_file = genome_file, shuffle = "q", random_trials = 30)
stats_accessible_sites <- list("direction"="G4s_in_K562_accessible_sites", 
                               "num_overlap_subset"=overlap_accessible_sites_with_G4s$num_overlap_subset,
                               "num_overlap_bases"=overlap_accessible_sites_with_G4s$num_overlap_bases,
                               "p_value_MC_subset"=overlap_accessible_sites_with_G4s$p_value_MC_subset,
                               "FE_MC_subset"=overlap_accessible_sites_with_G4s$FE_MC_subset,
                               "p_value_MC_bp"=overlap_accessible_sites_with_G4s$p_value_MC_bp,
                               "FE_MC_bp"=overlap_accessible_sites_with_G4s$FE_MC_bp)
stats_accessible_sites <- data.frame(stats_accessible_sites)
stats_K562 <- stats_accessible_sites# write to file
stats_all <- rbind(stats_HaCaT, stats_HEK, stats_K562)
write.table(stats_all, file="out/G4s_in_accessible_sites_all_cell_lines_stats.csv", sep=";", row.names=F)

overlap_CpG_ATAC <- analyse_window_size(query = CpG_locs, search_set=ATAC_data, genome_file=genome_file, window_sizes=window_sizes, all_CpGs_props = all_CpGs_props)
overlap_G4_ATAC_hacat <- overlap_G4_ATAC
overlap_CpG_ATAC_hacat <- overlap_CpG_ATAC

################################################################################
# Analyse TET enzyme overlap with G4s in human and mouse (genome-wide G4s)
# load TET data 
TET_data <- read.csv("data/TET_data/TET1_human_peaks.bed", sep="\t", header=F)
names(TET_data) <- c("chromosome", "start", "end", ".", "strand")
TET_data <- TET_data %>% select(chromosome, start, end) %>% filter(chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) 
#TET_data <- lift_over(coordinates=TET_data, chain_file="data/chain_files/mm9ToMm10.over.chain")
# load G4s 
G4_locs <- load_G4_data(G4_file_plus="data/G4_maps/G4s_human_plus.bed", G4_file_minus="data/G4_maps/G4s_human_minus.bed", autosomes=T)
# filter chrM out
G4_locs <- G4_locs %>% filter(!(chromosome=="chrM"))
# Overlap TET peaks and G4s 
TETs_in_G4s <- analyse_overlap(query=TET_data, search_set=G4_locs, genome_file="data/genome_files/hg19_chromInfo.txt", shuffle = "q", random_trials = 30, mask=mask)
 G4s_in_TETs <- analyse_overlap(query=G4_locs, search_set=TET_data, genome_file="data/genome_files/hg19_chromInfo.txt", shuffle="s", random_trials = 30, mask=mask)
# construct output dataframe to write to file
stats_TETs <- list("direction"="TETs_in_G4s", "num_overlap_subset"=TETs_in_G4s$num_overlap_subset, "num_overlap_bp"=TETs_in_G4s$num_overlap_bases, "p_value_MC_subset"=TETs_in_G4s$p_value_MC_subset, "FE_MC_subset"=TETs_in_G4s$FE_MC_subset, "p_value_MC_bp"=TETs_in_G4s$p_value_MC_bp, "FE_MC_bp"=TETs_in_G4s$FE_MC_bp)
stats_TETs <- data.frame(stats_TETs)
stats_G4s <- list("direction"="G4s_in_TETs", "num_overlap_subset"=G4s_in_TETs$num_overlap_subset, "num_overlap_bp"=G4s_in_TETs$num_overlap_bases, "p_value_MC_subset"=G4s_in_TETs$p_value_MC_subset, "FE_MC_subset"=G4s_in_TETs$FE_MC_subset, "p_value_MC_bp"=G4s_in_TETs$p_value_MC_bp, "FE_MC_bp"=G4s_in_TETs$FE_MC_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_TETs, stats_G4s)
write.table(stats_all, file="out/TET_analysis/mouse_TET1_C_all_G4_stats.csv", sep=";", row.names = F)

################################################################################
# DNMT1 / DNMT3b analysis
DNMT_data <- read.csv("data/DNMT_data/DNMT1_human_K562.bed", sep="\t", header=F)
names(DNMT_data) <- c("chromosome", "start", "end", ".", "..", "strand", "...", "....", ".....", "......")
DNMT_data <- DNMT_data %>% select(chromosome, start, end) %>% filter(chromosome %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) 
# load G4s 
G4_locs <- load_G4_data(G4_file_plus="data/G4_maps/BG4_K562.bed", narrow_peak = F)
#G4_locs <- lift_over(coordinates=G4_locs, chain_file="data/chain_files/hg19ToHg38.over.chain")
# filter chrM out
G4_locs <- G4_locs %>% filter(!(chromosome %in% c("chrX", "chrY", "chrM")))
# Overlap TET peaks and G4s 
DNMTs_in_G4s <- analyse_overlap(query=DNMT_data, search_set=G4_locs, genome_file="data/genome_files/hg19_chromInfo.txt", shuffle = "q", random_trials = 30, mask=mask)
G4s_in_DNMTs <- analyse_overlap(query=G4_locs, search_set=DNMT_data, genome_file="data/genome_files/hg38_chromInfo.txt", shuffle="s", random_trials = 30)
# construct output dataframe to write to file
stats_DNMTs <- list("direction"="DNMTs_in_G4s", "num_overlap_subset"=DNMTs_in_G4s$num_overlap_subset, "num_overlap_bp"=DNMTs_in_G4s$num_overlap_bases, "p_value_MC_subset"=DNMTs_in_G4s$p_value_MC_subset, "FE_MC_subset"=DNMTs_in_G4s$FE_MC_subset, "p_value_MC_bp"=DNMTs_in_G4s$p_value_MC_bp, "FE_MC_bp"=DNMTs_in_G4s$FE_MC_bp)
stats_DNMTs <- data.frame(stats_DNMTs)
stats_G4s <- list("direction"="G4s_in_DNMTs", "num_overlap_subset"=G4s_in_DNMTs$num_overlap_subset, "num_overlap_bp"=G4s_in_DNMTs$num_overlap_bases, "p_value_MC_subset"=G4s_in_DNMTs$p_value_MC_subset, "FE_MC_subset"=G4s_in_DNMTs$FE_MC_subset, "p_value_MC_bp"=G4s_in_DNMTs$p_value_MC_bp, "FE_MC_bp"=G4s_in_DNMTs$FE_MC_bp)
stats_G4s <- data.frame(stats_G4s)
stats_all <- rbind(stats_DNMTs, stats_G4s)
write.table(stats_all, file="out/DNMT_analysis/DNMT3b_HepG2_all_G4.csv", sep=";", row.names = F)

################################################################################
# Overlap CGIs and ATAC-seq data to get chromatin accessibility context of CGIs
CGIs <- read.csv("data/CGI_maps/hg19_CGI_map.bed", sep="\t", header = F)
CGIs <- CGIs %>% filter(!(chromosome=="chr11_gl000202_random"))
names(CGIs) <- c("chromosome", "start", "end", "id", "GC_content", "score")
ATAC_data <- load_ATAC_data(ATAC_file=ATAC_file)
results_all_CpGs <- read.csv(file=all_CpGs_props_file, header=T, sep=";")
# search_set needs to be CpG_locs, cant shuffle such large files
results_CGIs <- analyse_overlap(query=CGIs, search_set=ATAC_data, genome_file=genome_file, all_CpGs_props=results_all_CpGs, shuffle="q")
stats_ATAC <- analyse_overlap(query=ATAC_data, search_set=CGIs, genome_file=genome_file, all_CpGs_props=results_all_CpGs, shuffle="q")
stats_CGIs <- list("direction"="CGIs_in_ATAC_peaks", "num_overlap"=results_CGIs$num_overlaps, "fold_enrichment"=results_CGIs$fold_enrichment, "p_value_bp"=results_CGIs$p_value_bp)
stats_CGIs <- data.frame(stats_CGIs)
stats_ATAC <- list("direction"="ATAC_peaks_in_CGIs", "num_overlap"=stats_ATAC$num_overlaps, "fold_enrichment"=stats_ATAC$fold_enrichment, "p_value_bp"=stats_ATAC$p_value_bp)
stats_ATAC <- data.frame(stats_ATAC)
stats_all <- rbind(stats_CGIs, stats_ATAC)
write.table(stats_all, file=output_file, sep=";", row.names = F)

################################################################################
# Overlap BG4s of different cell lines to see commonalities
hacat <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_hacat.narrowPeak", narrow_peak = T)
HEK <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_HEK_merged.narrowPeak", narrow_peak = T)
K562 <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_K562.bed")
gw_G4s <- load_G4_data(G4_file_plus = "data/G4_maps/G4s_human_plus.bed", G4_file_minus="data/G4_maps/G4s_human_minus.bed")
data <- list(gw_G4s, hacat, HEK, K562)

# loop through all data set and calculate percentage of overlaps between datasets 
# in terms of G4s and basepairs 
matrix_G4s <- matrix(nrow=length(data), ncol=length(data))
matrix_bps <- matrix(nrow=length(data), ncol=length(data))
for(i in 1:length(data)){
  query <- data[[i]]
  for(j in 1:length(data)){
    search_set <- data[[j]]
    results_G4s <- analyse_overlap(query = query, search_set=search_set, 
                               genome_file = "data/genome_files/hg19_chromInfo.txt", 
                               shuffle="q", random_trials=30)
    # G4s
    matrix_G4s[i,j] <- results_G4s$num_overlap_subset / dim(query)[1]
    # Basepairs
    matrix_bps[i,j] <- results_G4s$num_overlap_bases / (query %>% mutate(length=end-start) %>% 
                                       summarise(num_bases=sum(length)))[["num_bases"]]
  }
}

shared_G4s <- data.frame(matrix_G4s, row.names=c("gw_G4s", "HaCaT", "HEK", "K562"))
names(shared_G4s) <- c("gw_G4s", "HaCaT", "HEK", "K562")
write.table(x=shared_G4s, file="out/shared_BG4_peaks_.csv", sep=";", row.names = T, col.names = T)
shared_bps <- data.frame(matrix_bps, row.names=c("gw_G4s", "HaCaT", "HEK",  "K562"))
names(shared_bps) <- c("gw_G4s", "HaCaT", "HEK", "K562")
write.table(x=shared_bps, file="out/shared_basepairs_in_BG4_peaks_.csv", sep=";", row.names = T, col.names = T)


###############################################################################
# Overlap ATAC-seq data from different cell lines to see commonalitites
hacat <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_hacat_merged.narrowPeak")
BG4_hacat <- load_G4_data(G4_file_plus="data/G4_maps/BG4_hacat.narrowPeak", narrow_peak = T)
HEK <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_HEK_merged.narrowPeak")
BG4_HEK_1 <- load_G4_data(G4_file_plus="data/G4_maps/BG4_HEK_rep1.narrowPeak", narrow_peak = T)
BG4_HEK_2 <- load_G4_data(G4_file_plus="data/G4_maps/BG4_HEK_rep2.narrowPeak", narrow_peak = T)
K562 <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_K562.bed")
BG4_K562 <- load_G4_data(G4_file_plus="data/G4_maps/BG4_K562.bed", narrow_peak = F)
data <- list(hacat, HEK, K562)

# loop through all data set and calculate percentage of overlaps between datasets 
# in terms of G4s and basepairs 
matrix_ATAC_peaks <- matrix(nrow=length(data), ncol=length(data))
matrix_bps <- matrix(nrow=length(data), ncol=length(data))
for(i in 1:length(data)){
  query <- data[[i]]
  for(j in 1:length(data)){
    search_set <- data[[j]]
    results_ATAC_peaks <- analyse_overlap(query = query, search_set=search_set, 
                                   genome_file = "data/genome_files/hg19_chromInfo.txt", shuffle="q",
                                   random_trials=30)
    # G4s
    matrix_ATAC_peaks[i,j] <- results_ATAC_peaks$num_overlap_subset / dim(query)[1]
    # Basepairs
    matrix_bps[i,j] <- results_ATAC_peaks$num_overlap_bases / (query %>% mutate(length=end-start+1) %>% 
                                      summarise(num_bases=sum(length)))[["num_bases"]]
  }
}

shared_ATAC_peaks <- data.frame(matrix_ATAC_peaks, row.names=c("HaCaT", "HEK", "K562"))
names(shared_ATAC_peaks) <- c("HaCaT", "HEK", "K562")
write.table(x=shared_ATAC_peaks, file="out/shared_ATAC_peaks.csv", sep=";", row.names = T, col.names = T)
shared_bps <- data.frame(matrix_bps, row.names=c("HaCaT", "HEK", "K562"))
names(shared_bps) <- c("HaCaT", "HEK", "K562")
write.table(x=shared_bps, file="out/shared_basepairs_in_ATAC_peaks.csv", sep=";", row.names = T, col.names = T)


################################################################################
# data for shared annotations Venn diagram

# G4s
hacat <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_hacat.narrowPeak", narrow_peak = T)
HEK <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_HEK_merged.narrowPeak", narrow_peak = T)
K562 <- load_G4_data(G4_file_plus = "data/G4_maps/BG4_K562.bed")
all <- rbind(hacat, HEK, K562)
HEK_hacat <-sum_bases(overlap(query = HEK, search_set=hacat)$overlapping_bases)
HEK_K562 <- sum_bases(overlap(query = HEK, search_set=K562)$overlapping_bases)
K562_hacat <- sum_bases(overlap(query = K562, search_set=hacat)$overlapping_bases)
HEK_hacat_K562 <- sum_bases(overlap(query = HEK, search_set=(overlap(query = K562, search_set = hacat)$overlapping_bases))$overlapping_bases) 

#absolut base pairs
all <- HEK_hacat_K562
HEK_hacat_only <- HEK_hacat - all
HEK_K562_only <- HEK_K562 - all
K562_hacat_only <- K562_hacat - all
HEK_only <- sum_bases(HEK) - HEK_hacat_only - HEK_K562_only - all
hacat_only <- sum_bases(hacat) - HEK_hacat_only - K562_hacat_only - all
K562_only <- sum_bases(K562) - HEK_K562_only - K562_hacat_only - all
total <- sum(all, HEK_hacat_only, HEK_K562_only, K562_hacat_only, HEK_only, hacat_only, K562_only)

# percentages 
all_p <- all / total
HEK_hacat_only_p <- HEK_hacat_only / total
HEK_K562_only_p <- HEK_K562_only / total 
K562_hacat_only_p <- K562_hacat_only / total
HEK_only_p <- HEK_only / total 
hacat_only_p <- hacat_only / total
K562_only_p <- K562_only/total


# ATAC 
hacat <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_hacat_merged.narrowPeak")
HEK <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_HEK_merged.narrowPeak")
K562 <- load_ATAC_data(ATAC_file = "data/ATAC_maps/ATAC_K562.bed")
all <- rbind(hacat, HEK, K562)

HEK_hacat <-sum_bases(overlap(query = HEK, search_set=hacat)$overlapping_bases)
HEK_K562 <- sum_bases(overlap(query = HEK, search_set=K562)$overlapping_bases)
K562_hacat <- sum_bases(overlap(query = K562, search_set=hacat)$overlapping_bases)
HEK_hacat_K562 <- sum_bases(overlap(query = HEK, search_set=(overlap(query = K562, search_set = hacat)$overlapping_bases))$overlapping_bases) 

#absolut base pairs
all <- HEK_hacat_K562
HEK_hacat_only <- HEK_hacat - all
HEK_K562_only <- HEK_K562 - all
K562_hacat_only <- K562_hacat - all
HEK_only <- sum_bases(HEK) - HEK_hacat_only - HEK_K562_only - all
hacat_only <- sum_bases(hacat) - HEK_hacat_only - K562_hacat_only - all
K562_only <- sum_bases(K562) - HEK_K562_only - K562_hacat_only - all
total <- sum(all, HEK_hacat_only, HEK_K562_only, K562_hacat_only, HEK_only, hacat_only, K562_only)

# percentages 
all_p <- all / total
HEK_hacat_only_p <- HEK_hacat_only / total
HEK_K562_only_p <- HEK_K562_only / total 
K562_hacat_only_p <- K562_hacat_only / total
HEK_only_p <- HEK_only / total 
hacat_only_p <- hacat_only / total
K562_only_p <- K562_only/total

################################################################################
# Merge data sets
merged <- merge_replicates(BG4_HEK_1, BG4_HEK_2)
merged <- mutate(merged, dummy1=".", dummy2="..")
write.table(merged, file="data/G4_maps/BG4_HEK_merged.narrowPeak", quote=F, col.names = F, row.names=F,sep="\t")


################################################################################
# Construct mask where to constrict shuffling
genome <- read.csv("data/genome_files/hg19_chromInfo.txt", sep="\t", 
                   col.names = c("chromosome","length","source"), header=F)
lengths <- genome$length
CGIs <- read.csv("data/CGI_maps/hg19_CGI_map.bed", sep="\t", header = F, col.names = c("chromosome", "start", "end", "id", "GC_content", "score"))
CGIs <- CGIs %>% filter(!(chromosome=="chr11_gl000202_random"))
CGIs_r <- makeGRangesFromDataFrame(CGIs)
seqlengths(CGIs_r) <- c(chr=lengths)
CGI_gaps <- gaps(CGIs_r)
gaps <- GRanges_to_df(CGI_gaps)
gaps <- gaps %>% filter(!(start==1 & end %in% lengths))
write.table(gaps, file="data/CGI_maps/hg19_CGI_opposite_mask.bed", quote=F, sep="\t", row.names=F, col.names = F)


################################################################################
# Plot the results 
stats <- read.csv("out/BG4_K562_output/Horvath_stats.csv", header=T, sep=";")
plot_results(overlap_results=NULL, window_size=NULL, stats=stats,
             query_name="CpGs", search_set_name="G4s",
             legend_position_fe = c(0.79, 0.9), legend_position_p = c(0.79, 0.8),
             figure_name_enrichment = "out/BG4_K562_output/Horvath_enrichment.pdf", 
             figure_name_p_value = "out/BG4_K562_output/Horvath_p_value.pdf" )
stats_1 <- read.csv("out/chromatin_open_vs_closed/BG4_HEK/rep1_Levine_closed_stats.csv", header=T, sep=";")
stats_2 <- read.csv("out/chromatin_open_vs_closed/BG4_HEK/rep1_Levine_open_stats.csv", header=T, sep=";")
plot_grouped_results(overlap_results=NULL, window_size=NULL, 
                     stats_1=stats_1, stats_2 = stats_2, 
                     group_1 = "closed", group_2 = "open",
                     query_name="CpGs", search_set_name="G4s",
                     legend_position_fe = NULL, legend_position_p = NULL,
                     figure_name_enrichment = "out/chromatin_open_vs_closed/BG4_HEK/Levine_enrichment.pdf", 
                     figure_name_p_value = "out/chromatin_open_vs_closed/BG4_HEK/Levine_p_value.pdf" )





