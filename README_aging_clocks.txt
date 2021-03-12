Workflow for DNAm clock CpG analysis

MEME analysis:
1) Obtain list of CpG sites (genomic coordinates, chromosome number) in form of
    .csv or similar
2) Extend single CpGs to genomic ranges of n bases around each CpG and convert
    to .bed file (script csv_to_bed.R).
3) If using older genome assembly data, convert genome coordinates to
    appropriate (e.g. most recent) assembly using LiftOver
    (https://genome.ucsc.edu/cgi-bin/hgLiftOver) on .bed files.
4) Obtain respective DNA sequences from genome coordinates data (.bed file)
    using script find_sequences.R and save as .fasta file (using BSgenome, Biostrings
    and seqinr packages in R).
5) Queue up MEME job with your .fasta file on website (https://meme-suite.org/meme/info/status?service=MEME&id=appMEME_5.3.316142543884241624959992)

G4 overlap analysis:
1) Get data from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110582)
2) Intersect CpGs with G4s using bedtools intersect -wa
3) Explore overlap data: Correlation coeff of overlapping CpGs, number of overlaps, Differential enrichment against control case (all CpGs). 

Control case: 
1) Using BSgenome, get hg38 and extract all CpGs in sequence
Note: Only perfect matches, excludes SNPs or other mutations.
