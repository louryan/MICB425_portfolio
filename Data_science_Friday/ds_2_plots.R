library(tidyverse)
library(phyloseq)
metadata <- read.delim("Saanich.metadata.txt")
OTU <- read.delim("Saanich.OTU.txt")
load("phyloseq_object.RData")


#Exercise 1
ggplot(metadata, aes(x=NO3_uM, y=Depth_m))+
  geom_point(color = "purple", shape=17, size = 1)

#Exercise 2
#Change from °C to °F
metadata %>%
  mutate(Temperature_F = (Temperature_C*9/5)+32) %>%
  select(Temperature_F, Depth_m) %>%
#Plot data
  ggplot(., aes(x=Temperature_F, y=Depth_m))+
  geom_point()

#Exercise 3
#convert physeq data to % and plot
physeq_percent = transform_sample_counts(physeq, function(x) 100*x/sum(x))
plot_bar(physeq_percent, fill="Phylum")+
  geom_bar(aes(fill=Phylum), stat = "identity")+
#Add title and axis labels
  ggtitle("Phyla from 10 to 200 m in Saanich Inlet")+
  xlab("Sample depth") + ylab("Percent relative abundance")

#Exercise 4
#Gather data in 1 column
metadata_nutrients <- metadata %>%
  gather(Nutrient, Concentration, O2_uM, PO4_uM, SiO2_uM, NO3_uM, NH4_uM, NO2_uM)%>%
  select(Depth_m, Nutrient, Concentration)
#Plot facet plot
ggplot(metadata_nutrients, aes(x=Depth_m, y=Concentration))+
  geom_line(color="grey")+
  geom_point()+
  facet_wrap(~Nutrient, scales="free_y")+
  theme(legend.position = "none")+
  ggtitle("Nutrient distribution across a depth of 10-200m in Saanich Inlet")+
  xlab("Depth (m)")+
  ylab("Nutrient Concentration (uM)")
