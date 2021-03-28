#!/bin/sh
#Need to execute this script from project root directory aging_clocks/
# Usage: sh sh_scripts/CpG_G4_overlap.sh CGI_file G4_file_plus G4_file_minus chain_file genome_file output_file name

# first transform standardised CpG .csv file into .bed file
CGI_file=$1 # file containing the CpG coordinates
G4_file_plus=$2 # file containing plus strand of G4 coordinates
G4_file_minus=$3 # file containing minus strand of G4 coordinates
genome_file=$4 # bed file specifying shuffleling ranges for overlap
output_file=$5 # desired output file name for fold enrichment results
organism=$6 #name extension of output files
echo -e "Processing files..."
echo -e "CGI file: \t $CGI_file"
echo -e "G4 files: \t $G4_file_plus (plus), $G4_file_minus (minus)"
echo -e "Genome file: \t $genome_file"
echo -e "Output file: \t $output_file"
echo -e "Organism: \t $organism"

# OVERLAPPING G4 AND CGI DATA

# concatenate both strands of G4 data
echo "Concatenating both strands of G4 file ..."
cat $G4_file_plus $G4_file_minus > temp/G4s.bed
echo "... concatenated G4 files."

# Intersecting files
echo "Finding G4s with overlap to CGIs"
name_G4_out_file="out/G4s_in_CGIs_${organism}.bed"
bedtools intersect -wa -a temp/G4s.bed -b $CGI_file > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > $name_G4_out_file
n_G4s=$(wc -l $name_G4_out_file | awk '{print $1}')
echo "Number of G4s found to lie within CGIs: $n_G4s"

name_CGI_out_file="out/CGIs_in_G4s_${organism}.bed"
echo "Finding CGIs with overlap to G4s"
bedtools intersect -wa -a $CGI_file -b temp/G4s.bed > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed 
bedtools merge -i temp/overlap_sorted.bed -c 4,5,6 -o distinct,mean,mean > $name_CGI_out_file
n_CGIs=$(wc -l $name_CGI_out_file | awk '{print $1}')
echo "Number of CpGs found to lie within G4s: $n_CGIs"

# Assess Fold enrichment
echo "Control case: shuffle G4 coordinates and overlap again..."
# loop through multiple random seeds to get good value for control
# get mean for CpGs and G4s
mean_G4s=0
for i in 1 2 3
do
    echo "Loop at i = $i"
    #Shuffle the data
    bedtools shuffle -i temp/G4s.bed -g $genome_file -seed $i > temp/G4s_shuffled.bed

    #calculate number of G4s that lie within shuffled CpGs
    bedtools intersect -wa -a temp/G4s_shuffled.bed -b $CGI_file > temp/overlap.bed
sort -k1,1 -k2,2n temp/overlap.bed > temp/overlap_sorted.bed
bedtools merge -i temp/overlap_sorted.bed -c 4 -o mean > temp/shuffled_G4s_in_CGIs.bed
    n_G4s_temp=$(wc -l temp/shuffled_G4s_in_CGIs.bed | awk '{print $1}')
    # calculate running average with bc program
    mean_G4s=$(echo "scale=6;$mean_G4s+1/3*$n_G4s_temp" | bc)
done

#Calculate fold enrichment
FE_G4s=$(echo "scale=6;$n_G4s/$mean_G4s" | bc)
echo "G4s are enriched at CGI sites to a fold enrichment value of: $FE_G4s"

#write input_file, FE and window size to file
echo -e "$CGI_file; $organism; $FE_G4s" > $output_file
