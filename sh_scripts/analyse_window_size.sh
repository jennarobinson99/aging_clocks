#!/bin/sh

# script takes in range of values for window size and input file
# to execute CpG_G4_overlap.sh once for each case 
# USAGE: sh sh_scripts/analyse_window_size.sh input_file output_file [window sizes...]
#         e.g. sh sh_script/analyse_window_size.sh data/CpG_lists/Horvath_CpGs out/fold_enrichment.csv 10 100 1000 10000
echo "$*"
input_file=$1 
output_file=$2
shift
shift 

for window_size in "$@"; do
	echo "Analysing overlap with window size = $window_size"
	sh sh_scripts/CpG_G4_overlap.sh $input_file $output_file $window_size
done


