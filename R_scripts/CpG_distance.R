setwd("D:/proj_epigen/aging_clocks")
library(dplyr)
#load the bed file containing genome coordinates 
coordinates <- read.csv("Cpg_lists/bed_files/global_CpGs_ext.bed", header=F, col.names=c("chr", "start", "end"), sep='\t')
#add row with CpG location
coordinates <-mutate(coordinates, cpg_loc = (end-start)/2+start)
#get mean, standard deviation, min and maximum gap of CpG distances per chromosome 
# plus other interesting statistics
window_length <- 100
gaps <- coordinates %>% 
  arrange(desc(cpg_loc)) %>% 
  group_by(chr) %>% 
  summarise(mean_gap = mean(lag(cpg_loc)-cpg_loc,na.rm=T), 
            std_gap = sd(lag(cpg_loc)-cpg_loc, na.rm=T), 
            min_gap = min(lag(cpg_loc)-cpg_loc, na.rm=T), 
            max_gap = max(lag(cpg_loc)-cpg_loc, na.rm=T), 
            number_of_CpGs=n(), 
            #neighbouring CpGs are ones that lie within a distance from each other 
            #that is smaller than half the base pair window
            neighbouring_CpGs=sum((lag(cpg_loc)-cpg_loc) < window_length/2, na.rm=T) 
            ) %>%
  mutate(chromosome = as.integer(substr(chr, 4, 6))) %>%
  arrange(chromosome)

#record gap data in csv file
write.table(gaps, file="CpG_lists/csv_files/global_CpGs_gaps.csv", sep=",", col.names = T, row.names = F, quote = F)
