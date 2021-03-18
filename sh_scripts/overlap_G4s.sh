#!/bin/sh

# concatenate both strands of G4 data
echo "Concatenating both strands of G4 file ..."
cat data/G4_maps/G4s_human_minus.bed data/G4_maps/G4s_human_plus.bed > temp/G4s_human.bed
echo "... concatenated G4 files."

#echo "Intersecting G4s and CpGs..."

#cd D:/proj_epigen/aging_clocks/
#bedtools intersect -a G4_maps/Human_G4s.bed -b CpG_lists/bed_files/horvath_CpGs_hg38.bed > out/intersect_G4_horvath.bed

#echo "Successfully overlapped."
