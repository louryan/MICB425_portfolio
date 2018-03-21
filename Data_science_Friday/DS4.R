#Load tidyverse
library(tidyverse)

#Load data
md<- read.table(file="Saanich.metadata.txt", header=TRUE, row.names = 1, sep="\t", na.strings=c("NAN","NA","."))
otu<- read.table(file="Saanich.OTU.txt", header=TRUE, row.names = 1, sep="\t", na.strings=c("NAN","NA","."))

#Filter out methane above 100 nM and a temperature below 10 degrees C
md %>%
  filter(CH4_nM >100) %>%
  filter(Temperature_C <10) %>%
  select(Depth_m, CH4_nM, Temperature_C)

md %>%
  filter(CH4_nM >100, Temperature_C <10) %>%
  mutate(CH4_uM = CH4_nM/1000) %>%
  select(Depth_m, CH4_nM, CH4_uM, Temperature_C)

md %>%
  mutate(CH4_uM = CH4_nM/1000) %>%
  mutate(N2O_uM = N2O_nM/1000) %>%
  select(Depth_m, CH4_nM,CH4_uM,N2O_nM,N2O_uM)

ggplot(md, aes(x=O2_uM, y=Depth_m)) +
  geom_p


#Convert all variables that are in nM to uM. Output a table showing only the original nM and converted uM variables.

