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
CpG_vs_CGI=$9 #test enrichment vs global or vs CGI sites?
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
echo "Lifting over CpGs into hg19..."
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
echo "Finding G4s with overlap to AC CpG windows..."
name_G4_out_file="out/${name}_G4s_in_CpGs_ws_${window_size}.bed"
bedtools intersect -wa -a temp/G4s.bed -b temp/CpGs_ext.bed > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > $name_G4_out_file
n_G4s_AC=$(wc -l $name_G4_out_file | awk '{print $1}')
echo "Number of G4s found to lie within aging clock CpGs: $n_G4s_AC"

echo "Finding G4s with overlap to all CpGs..."
name_CpG_out_file="out/global_CpGs_in_G4s_ws_${window_size}.bed"
bedtools intersect -wa -a temp/G4s.bed -b data/CpG_lists/all_CpGs_hg19_ext.bed > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed >temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > $name_CpG_out_file
n_G4s_global=$(wc -l $name_CpG_out_file | awk '{print $1}')
echo "Number of G4s found to lie within all CpGs in the genome: $n_G4s_global"

# Assess Fold enrichment
echo "Control case: shuffle G4 coordinates and overlap again..."
# loop through multiple random seeds to get good value for control
# get mean for CpGs and G4s
# when measuring enrichment globally
echo "Measuring enrichment of G4s in AC CpGs vs all CpGs"
mean_G4s_AC=0
mean_G4s_global=0
for i in 1 2 3
do
    echo "Loop at i = $i"

    #Shuffle the data
    bedtools shuffle -i temp/G4s.bed -g $genome_file -seed $i > temp/G4s_shuffled.bed

    #calculate number of G4s that lie within shuffled aging clock CpGs
    bedtools intersect -wa -a temp/G4s_shuffled.bed -b temp/CpGs_ext.bed > temp/overlap.bed
    sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
    bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > temp/G4s_in_shuffled_CpGs_temp.bed
    n_G4s_temp=$(wc -l temp/G4s_in_shuffled_CpGs_temp.bed | awk '{print $1}')
    # calculate running average with bc program
    mean_G4s_AC=$(echo "scale=6;$mean_G4s_AC+1/3*$n_G4s_temp" | bc)

    #calculate number of G4s that lie within all CpGs in the genome
    bedtools intersect -wa -a temp/G4s_shuffled.bed -b data/CpG_lists/all_CpGs_hg19_ext.bed > temp/overlap.bed
    sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
    bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > temp/CpGs_in_shuffled_G4s_temp.bed
    n_G4s_temp=$(wc -l temp/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
    # calcuate running average with bc program
    mean_CpGs_global=$(echo "scale=6;$mean_CpGs_global+1/3*$n_G4s_temp" | bc)
done

#Calculate fold enrichment
FE_G4s_AC=$(echo "scale=6;$n_G4s_AC/$mean_G4s_AC" | bc)
FE_G4s_global=$(echo "scale=6;$n_G4s_global/$mean_G4s_global" | bc)
echo "G4s are enriched at aging clock CpG sites to a fold enrichment value of: $FE_G4s_AC"
echo "G4s are enriched at genome-wide CpGs to a fold enrichment value of: $FE_G4s_global"

#write input_file, FE and window size to file
echo -e "$CpG_file; $window_size; $FE_G4s_AC; $FE_G4s_global" >> $output_file
