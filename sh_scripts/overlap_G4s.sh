#!/bin/sh

echo "Intersecting G4s and CpGs..."

cd D:/proj_epigen/aging_clocks/
bedtools intersect -a G4_maps/Human_G4s.bed -b CpG_lists/bed_files/horvath_CpGs_hg38.bed > out/intersect_G4_horvath.bed

echo "Successfully overlapped."

