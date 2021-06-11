Directory structure:
    data/ -> Contains all the data needed for analysis: aging clock CpG lists and global CpG lists (CpG_lists/), whole-genome G4 and BG4-ChIP data (G4_maps/), chain files for lift over of genome assemblies (chain_files/), maps of CpG islands in the genome (CGI_maps/) and genome files specifying the size of each chromosome which is required for bedtools shuffle (genome_files/).
    out/ -> Output directory containing all final outputs of scripts.
    temp/ -> Directory containing all temporarily produced data, NOT final results.
    sequences/ -> Contains all sequences (fasta format) that need to be analysed by other tools (e.g. MEME analysis).
    R_scripts/ -> Contains all scripts written in R.
    sh_scripts/ -> Contains all scripts written for Bash.

How to use this repo:
    - The control_script.R gives an overview over the workflow for colocalization. Further, it contains miscellaneous experimental setups that vary slightly from the standard workflow. Refer to comments in the script to see the details.
    - The functions.R script defines all functions needed to conduct colocalization analysis. This includes functions for overlapping features, loading data, data wrangling and plotting. Refer to comments in the script to see the details.
    - The config.json file defines important parameters for the colocalization analysis. This is loaded into the control_script.R to define the experimental setup.
    - Other scripts not classified as main scripts are DEPRECATED and can be ignored.

Main scripts:
    - control_script.R: controls the colocalization workflow
    - functions.R: defines all important functions for the colocalization analysis. This script is called from control_script.R

Environment:
    The following packages are required. The utilised versions are given in brackets.
    -

Description of other scripts for Bash:
*** NOTE: The Bash scripts were originally used to automate the workflow between scripts (Bash and R scripts). They are usually written in pairs, one iterating through various window sizes (control script c) and the other one executing the overlap an analysis for a given window size (overlap script o). ***
    0) config.sh -> This file is read by the control script to specify all inputs (e.g. input file name, output file name, chain file, ...) except the window sizes.
    1) analyse_window_size.sh (c) & CpG_G4_overlap.sh (o) -> Executes standard overlap analysis of a set of aging clock CpGs with a set of G4s. Fold enrichment can be calculated as overlaps between G4s and AC CpGs divided by shuffled version of either within the whole genome (e.g. AC CpGs within G4s vs. random distribution of AC CpGs throughout the genome). Alternatively if the parameter CpG_vs_CGI in config.sh is set to 1, the enrichment is calculated by dividing through G4 overlaps with CpGs shuffled within known CGI islands (only available for human clocks). This may elsewhere be referred to as the "pseudo-CGIs" method. Output file contains fold enrichment as well as number of G4s and CpGs overlapping.
    2) AC_CpG_vs_global_CpG.sh (c) & AC_CpG_vs_global_CpG_overlap.sh (o) -> Overlaps G4s with AC CpGs and G4s with all CpGs in the genome. It outputs both the fold enrichment of G4s at AC CpG sites and the fold enrichment of G4s at all CpGs in the genome.
    3) global_CpG_vs_G4.sh (c) & global_CpG_vs_G4_overlap.sh (o) -> Overlaps AC CpGs with G4s and all CpGs in the genome with G4s. Outputs fold enrichment of AC CpGs at G4s and all CpGs at G4s (reverse overlap than 2).00
    4) separate_CGI_context (c) & contextual_overlap.sh (o) -> Separates CpGs according to CGI context and overlaps both separately with G4s. Outputs the fold enrichment of these CpG sub-lists calculated with respect to random shuffled CpGs in G4s.
    5) CGI_FE.sh (c) & CGI_G4_overlap.sh (o) -> Overlaps G4s with CpG islands and gives the fold enrichment of G4s in CpG islands with respect to random distribution of G4s throughout the genome.

Description of other scripts in R:
*** NOTE: The following scripts are designed to be used from the command line and usually are interfaced by Bash scripts which control the higher level workflow (see "Description of scripts for Bash") ***
    Preprocessing for overlap analysis:
        csv_to_bed.R -> converts an input file from .csv format to .bed format.
        lift_over.R -> lifts over genome coordinates of input file, using a supplied chain file (e.g. needed to go from hg18 -> hg19). hg19 is the standard assembly used in this project.
        filter_CpGs.R -> splits up the supplied CpG lists according to positive / negative correlation and CGI context (inside / outside CpG island) and saves all generated sub-lists into temp/ directory.
    Analysis of overlap results:
        analyse_overlap_results.R -> takes in output file of overlap analysis (e.g. fold enrichment file, CpG and G4 lists) and plots the enrichment against all window sizes tested. Also generates bar charts showing the distribution of overlapping CpGs between positive/negative correlation and inside/outside CGIs.
        CGI_context_results.R -> Plots the fold enrichment at given window size and CpG distributions for the specific case of analysing CpGs within CGIs and CpGs outside CGIs separately.

*** NOTE: The following scripts are designed to be used in interactive mode (e.g. RStudio) and are independent of other scripts. ***
    Miscellaneous:
        G4_map_exploration.R -> generates some key statistics on the given G4 data. Used for initial data exploration.
        CpG_distance.R -> generates some key statistics on the given CpG lists. Used for initial data exploration.
        find_all_CpGs.R -> finds all CpGs within a genome and generates a file containing their genome coordinates extended by a chosen window size.
        find_sequences.R -> takes in a list of genome coordinates (e.g. generated by find_all_CpGs.R) and finds their corresponding sequences. The sequences are saved in .fasta format in the sequences/ directory.


Virtual environment:
    - The following packages are necessary to conduct the analysis: tidyverse, GenomicRanges, regioneR, rtracklayer, scales, IRanges, Biostrings, stringr and stats.
    - For detailed version numbers, refer to the below output of sessionInfo():

    >sessionInfo()
        R version 4.0.3 (2020-10-10)
        Platform: x86_64-w64-mingw32/x64 (64-bit)
        Running under: Windows 10 x64 (build 19041)

        Matrix products: default

        locale:
        [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252
        [4] LC_NUMERIC=C                    LC_TIME=German_Germany.1252

        attached base packages:
        [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base

        other attached packages:
         [1] scales_1.1.1         forcats_0.5.1        stringr_1.4.0        purrr_0.3.4          readr_1.4.0
         [6] tidyr_1.1.3          tibble_3.1.0         ggplot2_3.3.3        tidyverse_1.3.0      Biostrings_2.58.0
        [11] XVector_0.30.0       regioneR_1.22.0      GenomicRanges_1.42.0 GenomeInfoDb_1.26.5  IRanges_2.24.1
        [16] S4Vectors_0.28.1     BiocGenerics_0.36.0  rjson_0.2.20         dplyr_1.0.5

        loaded via a namespace (and not attached):
         [1] MatrixGenerics_1.2.1        Biobase_2.50.0              httr_1.4.2
         [4] jsonlite_1.7.2              modelr_0.1.8                assertthat_0.2.1
         [7] BSgenome_1.58.0             GenomeInfoDbData_1.2.4      cellranger_1.1.0
        [10] Rsamtools_2.6.0             pillar_1.5.1                backports_1.2.1
        [13] lattice_0.20-41             glue_1.4.2                  rvest_1.0.0
        [16] colorspace_2.0-0            Matrix_1.2-18               XML_3.99-0.6
        [19] pkgconfig_2.0.3             broom_0.7.6                 haven_2.3.1
        [22] zlibbioc_1.36.0             BiocParallel_1.24.1         generics_0.1.0
        [25] ellipsis_0.3.1              cachem_1.0.4                withr_2.4.1
        [28] SummarizedExperiment_1.20.0 cli_2.4.0                   magrittr_2.0.1
        [31] crayon_1.4.1                readxl_1.3.1                memoise_2.0.0
        [34] fs_1.5.0                    fansi_0.4.2                 xml2_1.3.2
        [37] tools_4.0.3                 hms_1.0.0                   lifecycle_1.0.0
        [40] matrixStats_0.58.0          munsell_0.5.0               reprex_2.0.0
        [43] DelayedArray_0.16.3         compiler_4.0.3              tinytex_0.31
        [46] rlang_0.4.10                grid_4.0.3                  RCurl_1.98-1.3
        [49] rstudioapi_0.13             bitops_1.0-6                gtable_0.3.0
        [52] DBI_1.1.1                   R6_2.5.0                    GenomicAlignments_1.26.0
        [55] lubridate_1.7.10            rtracklayer_1.49.5          fastmap_1.1.0
        [58] utf8_1.2.1                  stringi_1.5.3               Rcpp_1.0.6
        [61] vctrs_0.3.7                 dbplyr_2.1.1                tidyselect_1.1.0
        [64] xfun_0.20

Miscellaneous:
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

    global CpGs source:
    1) Using BSgenome, get hg19 and extract all CpGs in sequence
    Note: Only perfect matches, excludes SNPs or other mutations.
