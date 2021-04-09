CpG_file="data/CpG_lists/Horvath_clock_CpGs.csv"
G4_file_minus="data/G4_maps/BG4_peaks_2016.bed"
G4_file_plus="data/G4_maps/empty_file.txt"
chain_file="data/chain_files/hg18ToHg19.over.chain"
genome_file="data/genome_files/hg19_chromInfo.txt"
output_file="out/Horvath_FE.csv"
name="Horvath"
CGI_map_file="-"
window_size_CpG_dist_plot="50"
organism="Human"
CpG_vs_CGI="0" #enrichment test: only set to 1 if wanna test G4 enrichment at CpGs vs CGIs, otherwise leave blank
all_CpGs_file="data/CpG_lists/all_CpGs_mm10_ws_1000.bed" #only important for when running script AC_CpG_vs_global_CpG.sh
