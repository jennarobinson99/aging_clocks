#!/bin/sh

# USAGE: sh sh_scripts/CGI_FE.sh
# Wrapper for CGI_G4_overlap using input from config.sh
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

echo "Analysing G4 enrichment in CpG islands (CGIs)..."
sh sh_scripts/CGI_G4_overlap.sh $CGI_map_file $G4_file_plus $G4_file_minus $genome_file $output_file $organism
