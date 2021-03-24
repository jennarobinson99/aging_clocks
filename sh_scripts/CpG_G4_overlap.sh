#!/bin/sh
#Need to execute this script from project root directory aging_clocks/
# Usage: sh sh_scripts/CpG_G4_overlap.sh CpG_file G4_file_plus G4_file_minus chain_file genome_file output_file window_size name

# first transform standardised CpG .csv file into .bed file
CpG_file=$1 # file containing the CpG coordinates
G4_file_plus=$2 # file containing plus strand of G4 coordinates
G4_file_minus=$3 # file containing minus strand of G4 coordinates
chain_file=$4 # chain file for lifting over between genome builds
genome_file=$5 # bed file specifying shuffleling ranges for overlap
output_file=$6 # desired output file name for fold enrichment results
window_size=$7 # window size around CpGs 
name=$8 #name extension of output files
echo -e "Processing files..."
echo -e "CpG file: \t $CpG_file"
echo -e "G4 files: \t $G4_file_plus (plus), $G4_file_minus (minus)"
echo -e "Chain file: \t $chain_file"
echo -e "Genome file: \t $genome_file"
echo -e "Output file: \t $output_file"
echo -e "Window size: \t $window_size"
echo -e "Name: \t $name"
echo "Transforming CpG format from .csv to .bed ..."
Rscript --vanilla R_scripts/csv_to_bed.R $CpG_file temp/CpGs_locs.bed
echo "... finished transforming."

# next lift over CpG coordinates into horvath_CpGs_hg38
echo "Lifting over CpGs into hg38..."
Rscript --vanilla R_scripts/lift_over.R temp/CpGs_locs.bed temp/CpGs_lifted.bed $chain_file
echo "... lifted CpGs."

# Next filter coordinates (positive, negative correlation) and extend bases by window length n
echo "Filtering and extending bases..."
Rscript --vanilla R_scripts/filter_CpGs.R temp/CpGs_lifted.bed temp/CpGs_ext.bed $window_size
echo "... filtered bases."

# OVERLAPPING G4 AND CpG DATA

# concatenate both strands of G4 data
echo "Concatenating both strands of G4 file ..."
cat $G4_file_plus $G4_file_minus > temp/G4s.bed
echo "... concatenated G4 files."

# Intersecting files
echo "Finding G4s with overlap to CpG windows..."
name_G4_out_file="out/${name}_G4s_in_CpGs_ws_${window_size}.bed"
bedtools intersect -wa -a temp/G4s.bed -b temp/CpGs_ext.bed > $name_G4_out_file
n_G4s=$(wc -l $name_G4_out_file | awk '{print $1}')
echo "Number of G4s found to lie within CpGs: $n_G4s"

echo "Finding CpGs with overlap to G4s..."
name_CpG_out_file="out/${name}_CpGs_in_G4s_ws_${window_size}.bed"
bedtools intersect -wa -a temp/CpGs_ext.bed -b temp/G4s.bed > $name_CpG_out_file
n_CpGs=$(wc -l $name_CpG_out_file | awk '{print $1}')
echo "Number of CpGs found to lie within G4s: $n_CpGs"

# Assess Fold enrichment 
echo "Control case: shuffle G4 and CpGs coordinates and overlap again..."
# loop through multiple random seeds to get good value for control
# get mean for CpGs and G4s
mean_CpGs=0
mean_G4s=0
for i in 1 2 3
do
echo "Loop at i = $i"
#Shuffle the data
bedtools shuffle -i temp/G4s.bed -g $genome_file -seed $i > temp/G4s_shuffled.bed
bedtools shuffle -i temp/CpGs_ext.bed -g $genome_file -seed $i > temp/CpGs_ext_shuffled.bed

#calculate number of G4s that lie within shuffled CpGs
bedtools intersect -wa -a temp/G4s.bed -b temp/CpGs_ext_shuffled.bed > temp/G4s_in_shuffled_CpGs_temp.bed
n_G4s_temp=$(wc -l temp/G4s_in_shuffled_CpGs_temp.bed | awk '{print $1}')
# calculate running average with bc program
mean_G4s=$(echo "scale=6;$mean_G4s+1/3*$n_G4s_temp" | bc)

#calculate number of CpGs that lie within G4s
bedtools intersect -wa -a temp/CpGs_ext.bed -b temp/G4s_shuffled.bed > temp/CpGs_in_shuffled_G4s_temp.bed
n_CpGs_temp=$(wc -l temp/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
# calcuate running average with bc program
mean_CpGs=$(echo "scale=6;$mean_CpGs+1/3*$n_CpGs_temp" | bc)
done

#Calculate fold enrichment 
FE_G4s=$(echo "scale=6;$n_G4s/$mean_G4s" | bc)
FE_CpGs=$(echo "scale=6;$n_CpGs/$mean_CpGs" | bc)
echo "G4s are enriched at CpG sites to a fold enrichment value of: $FE_G4s"
echo "CpGs are enriched at G4 sites to a fold enrichment value of: $FE_CpGs"

#write input_file, FE and window size to file 
echo -e "$CpG_file; $window_size; $FE_G4s; $FE_CpGs" >> $output_file


