#!/bin/sh

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

echo "Control case: shuffle G4 and CpGs coordinates and overlap again..."
# loop through multiple random seeds to get good value for control 
for i in 1 2 3 4 5
do
echo "Loop at i = $i"
bedtools shuffle -i temp/G4s_human.bed -g data/genome_files/chromInfo.txt -seed $i > temp/G4s_human_shuffled.bed
bedtools shuffle -i temp/CpGs_ext.bed -g data/genome_files/chromInfo.txt -seed $i > temp/CpGs_ext_shuffled.bed

bedtools intersect -wa -a temp/G4s_human.bed -b temp/CpGs_ext_shuffled.bed > out/G4s_in_shuffled_CpGs_temp.bed
n_G4s=$(wc -l out/G4s_in_shuffled_CpGs_temp.bed | awk '{print $1}')
echo "G4s in CpG regions: $n_G4s"

bedtools intersect -wa -a temp/CpGs_ext.bed -b temp/G4s_human_shuffled.bed > out/CpGs_in_shuffled_G4s_temp.bed
n_CpGs=$(wc -l out/CpGs_in_shuffled_G4s_temp.bed | awk '{print $1}')
echo "CpG regions in G4s: $n_CpGs"

done
