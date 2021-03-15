#!/bin/sh
#Need to execute this script from project root directory aging_clocks/

# first transform standardised CpG .csv file into .bed file
echo "Transforming CpG format from .csv to .bed ..."
Rscript --vanilla R_scripts/csv_to_bed.R data/CpG_lists/Horvath_CpGs.csv temp/Horvath_CpGs.bed
echo "... finished transforming."

# next lift over CpG coordinates into horvath_CpGs_hg38
echo "Lifting over CpGs into hg38..."
Rscript --vanilla R_scripts/lift_over.R temp/Horvath_CpGs.bed temp/Horvath_CpGs_lifted.bed data/chain_files/hg18ToHg38.over.chain
echo "... lifted CpGs."

# Next filter coordinates (positive, negative correlation) and extend bases by window length n
echo "Filtering and extending bases..."
