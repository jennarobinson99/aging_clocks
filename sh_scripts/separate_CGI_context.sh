#!/bin/sh

# script takes in range of values for window size and input file
# to execute contextual_overlap.sh once for each case
# USAGE: sh sh_scripts/separate_CGI_context (window sizes...]
#         e.g. sh sh_script/separate_CGI_context.sh 10 100 1000 10000

# read in configuration data from sh_scripts/config.sh
config_file="sh_scripts/config.sh"
source $config_file
#echo "CpG file: $CpG_file"
#echo "G4 plus file: $G4_file_plus"
#echo "G4 minus file: $G4_file_minus"
#echo "Chain file: $chain_file"
#echo "Genome file: $genome_file"
#echo "Output file: $output_file"
#echo $name
#echo $window_size_CpG_dist_plot


for window_size in "$@"; do
	echo "Analysing overlap with window size = $window_size"
	sh sh_scripts/contextual_overlap.sh $CpG_file $G4_file_plus $G4_file_minus $chain_file $genome_file $output_file $window_size $name
done
echo "Producing figures..."
CpG_within_file="out/${name}_CpGs_within_CGI_ws_${window_size_CpG_dist_plot}.bed"
CpG_outside_file="out/${name}_CpGs_outside_CGI_ws_${window_size_CpG_dist_plot}.bed"
figure_within="out/${name}_FE_G4s.pdf"
figure_outside="out/${name}_FE_CpGs.pdf"
figure_dist="out/${name}_CpG_distribution.pdf"
figure_CGI="out/${name}_CGI_context.pdf"
Rscript --vanilla R_scripts/CGI_context_results.R $output_file $CpG_within_file $CpG_outside_file $CGI_map_file $figure_within $figure_outside $figure_dist $figure_CGI
echo "... saved figures."
