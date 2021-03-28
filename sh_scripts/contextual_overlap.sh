#!/bin/sh
#Need to execute this script from project root directory aging_clocks/
# Usage: sh sh_scripts/contextual_overlap.sh  $CpG_file $G4_file_plus $G4_file_minus $chain_file $genome_file $output_file $window_size $name

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
Rscript --vanilla R_scripts/filter_CpGs.R temp/CpGs_lifted.bed temp/CpGs_ext.bed data/CGI_maps/hg19_CGI_map.bed $window_size
echo "... filtered bases."

# OVERLAPPING G4 AND CpG DATA

# concatenate both strands of G4 data
echo "Concatenating both strands of G4 file ..."
cat $G4_file_plus $G4_file_minus > temp/G4s.bed
echo "... concatenated G4 files."

# Intersecting files
# within CGI
echo "Finding CGI-associated CpGs that lie in G4s..."
name_CpGs_within_file="out/${name}_CpGs_within_CGI_ws_${window_size}.bed"
bedtools intersect -wa -a temp/CpGs_within.bed -b temp/G4s.bed > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > $name_CpGs_within_file
n_within_CpGs=$(wc -l $name_CpGs_within_file | awk '{print $1}')
echo "Number of CGI-associated CpGs found to lie in G4s: $n_within_CpGs"

# outside CGI
echo "Finding non-CGI CpGs that lie in G4s..."
name_CpGs_outside_file="out/${name}_CpGs_outside_CGI_ws_${window_size}.bed"
bedtools intersect -wa -a temp/CpGs_outside.bed -b temp/G4s.bed > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed >temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > $name_CpGs_outside_file
n_outside_CpGs=$(wc -l $name_CpGs_outside_file | awk '{print $1}')
echo "Number of non-CGI CpGs found to lie within G4s: $n_outside_CpGs"

# Assess Fold enrichment
echo "Control case: shuffle G4 coordinates and overlap again..."
# loop through multiple random seeds to get good value for control
# get mean for CpGs and G4s

# when measuring enrichment globally
echo "Measuring global enrichment of G4s and CpGs..."
mean_CpGs_within=0
mean_CpGs_outside=0
for i in 1 2 3
do
    echo "Loop at i = $i"

    #Shuffle the data
    bedtools shuffle -i temp/G4s.bed -g $genome_file -seed $i > temp/G4s_shuffled.bed

    #calculate number of CGI-CpGs that lie within shuffled G4s
    bedtools intersect -wa -a temp/CpGs_within.bed -b temp/G4s_shuffled.bed > temp/overlap.bed
    sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
    bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > temp/CpGs_in_shuffled_G4s_temp.bed
    n_CpGs_within_temp=$(wc -l temp/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
    # calculate running average with bc program
    mean_CpGs_within=$(echo "scale=6;$mean_CpGs_within+1/3*$n_CpGs_within_temp" | bc)

    #calculate number of non-CGI-CpGs that lie within G4s   
bedtools intersect -wa -a temp/CpGs_outside.bed -b temp/G4s_shuffled.bed > temp/overlap.bed
    sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
    bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > temp/CpGs_in_shuffled_G4s_temp.bed
    n_CpGs_outside_temp=$(wc -l temp/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
    # calcuate running average with bc program
    mean_CpGs_outside=$(echo "scale=6;$mean_CpGs_outside+1/3*$n_CpGs_outside_temp" | bc)
done

#Calculate fold enrichment
FE_inside_CpGs=$(echo "scale=6;$n_within_CpGs/$mean_CpGs_within" | bc)
FE_outside_CpGs=$(echo "scale=6;$n_outside_CpGs/$mean_CpGs_outside" | bc)
echo "CGI-associated CpGs are enriched at G4s to a value of: $FE_inside_CpGs"
echo "Non-CGI CpGs are enriched at G4s to a  value of: $FE_outside_CpGs"

#write input_file, FE and window size to file
echo -e "$CpG_file; $window_size; $FE_inside_CpGs; $FE_outside_CpGs" >> $output_file
