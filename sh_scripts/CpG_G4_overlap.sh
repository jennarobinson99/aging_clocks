#!/bin/sh
#Need to execute this script from project root directory aging_clocks/
# Usage: sh sh_scripts/CpG_G4_overlap.sh input_file output_file window_size
# first transform standardised CpG .csv file into .bed file
input_file=$1
output_file=$2
window_size=$3
echo $1
echo $2
echo "Transforming CpG format from .csv to .bed ..."
Rscript --vanilla R_scripts/csv_to_bed.R $input_file  temp/CpGs_locs.bed
echo "... finished transforming."

# next lift over CpG coordinates into horvath_CpGs_hg38
echo "Lifting over CpGs into hg38..."
Rscript --vanilla R_scripts/lift_over.R temp/CpGs_locs.bed temp/CpGs_lifted.bed data/chain_files/hg18ToHg38.over.chain
echo "... lifted CpGs."

# Next filter coordinates (positive, negative correlation) and extend bases by window length n
echo "Filtering and extending bases..."
Rscript --vanilla R_scripts/filter_CpGs.R temp/CpGs_lifted.bed temp/CpGs_ext.bed $window_size
echo "... filtered bases."

# OVERLAPPING G4 AND CpG DATA

# concatenate both strands of G4 data
echo "Concatenating both strands of G4 file ..."
cat data/G4_maps/G4s_human_minus.bed data/G4_maps/G4s_human_plus.bed > temp/G4s_human.bed
echo "... concatenated G4 files."

# Intersecting files
echo "Finding G4s with overlap to CpG windows..."
bedtools intersect -wa -a temp/G4s_human.bed -b temp/CpGs_ext.bed > out/G4s_in_CpGs.bed
n_G4s=$(wc -l out/G4s_in_CpGs.bed | awk '{print $1}')
echo "Number of G4s found to lie within CpGs: $n_G4s"

echo "Finding CpGs with overlap to G4s..."
bedtools intersect -wa -a temp/CpGs_ext.bed -b temp/G4s_human.bed > out/CpGs_in_G4s.bed
n_CpGs=$(wc -l out/CpGs_in_G4s.bed | awk '{print $1}')
echo "Number of CpGs found to lie within G4s: $n_CpGs"

# Write overlap results to file

echo "Control case: shuffle G4 and CpGs coordinates and overlap again..."
# loop through multiple random seeds to get good value for control
# get mean for CpGs and G4s
mean_CpGs=0
mean_G4s=0
for i in 1 2 3
do
echo "Loop at i = $i"
#Shuffle the data
bedtools shuffle -i temp/G4s_human.bed -g data/genome_files/chromInfo.txt -seed $i > temp/G4s_human_shuffled.bed
bedtools shuffle -i temp/CpGs_ext.bed -g data/genome_files/chromInfo.txt -seed $i > temp/CpGs_ext_shuffled.bed

#calculate number of G4s that lie within shuffled CpGs
bedtools intersect -wa -a temp/G4s_human.bed -b temp/CpGs_ext_shuffled.bed > out/G4s_in_shuffled_CpGs_temp.bed
n_G4s_temp=$(wc -l out/G4s_in_shuffled_CpGs_temp.bed | awk '{print $1}')
# calculate running average with bc program
mean_G4s=$(echo "scale=6;$mean_G4s+1/3*$n_G4s_temp" | bc)

#calculate number of CpGs that lie within G4s
bedtools intersect -wa -a temp/CpGs_ext.bed -b temp/G4s_human_shuffled.bed > out/CpGs_in_shuffled_G4s_temp.bed
n_CpGs_temp=$(wc -l out/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
# calcuate running average with bc program
mean_CpGs=$(echo "scale=6;$mean_CpGs+1/3*$n_CpGs_temp" | bc)
done

#Calculate fold enrichment 
FE_G4s=$(echo "scale=6;$n_G4s/$mean_G4s" | bc)
FE_CpGs=$(echo "scale=6;$n_CpGs/$mean_CpGs" | bc)
echo "G4s are enriched at CpG sites to a fold enrichment value of: $FE_G4s"
echo "CpGs are enriched at G4 sites to a fold enrichment value of: $FE_CpGs"

#write input_file, FE and window size to file 
echo -e "$input_file; $window_size; $FE_G4s; $FE_CpGs" >> $output_file


