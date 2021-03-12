setwd("D:/proj_epigen/aging_clocks")
library(dplyr)
coordinates <- read.csv("G4_maps/Human_G4s.bed", header=F, col.names=c("chr", "start", "ends", "score"), sep='\t')
coordinates <- coordinates %>% mutate(len = ends - start)
G4_summary <- coordinates %>% summarise(average_length = mean(len),
                                        std_length = sd(len),
                                        min_length = min(len),
                                        max_length = max(len),
                                        number_G4s = n())
#record gap data in csv file
write.table(G4_summary, file="G4_maps/Human_G4s_summary.csv", sep=",", col.names = T, row.names = F, quote = F)
