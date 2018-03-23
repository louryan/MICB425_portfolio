library(phyloseq)
library(knitr)
library(kableExtra)
library(tidyverse)

subsample_bag3 = data.frame(
  Species_number = c(seq(1,32, by=1)),
  Species_type = c("MMs","MMs","MMs","MMs","MMs","MMs",
                   "Skittles","Skittles","Skittles","Skittles",
                  "Skittles","Gummy_Bears","Gummy_Bears","Gummy_Bears",
                   "Gummy_Bears","Gummy_Bears","Gummy_Bears","Rods","Rods",
                   "Rods","Rods","Rods","Circular","Filamentous",
                   "Filamentous","Filamentous", "Brick","Brick",
                   "Round","Round","Round","Round"),
  Species_name = c("Orange", "Yellow","Green","Red", "Brown","Blue",
                   "Yellow","Orange","Red","Brown","Green",
                   "Orange", "Pink","Yellow","White","Red","Green",
                   "Yellow","Orange","Pink","Green","Red",
                   "Blue_Swirl", "String_like","Spider",
                   "Coke_Bottle","Big_blue","Small_blue","Green", "Yellow",
                   "Red","Purple"),
  Abundance = c(13,8,4,3,4,12,6,7,4,10,8,1,3,2,4,3,4,5,8,5,2,4,
                1,1,1,1,1,1,1,1,1,1)
)

subsample_bag3 %>%
  kable("html") %>%
  kable_styling(bootstrap_options = "striped", font_size = 10,
                full_width = F)

collectors_curve = data.frame(
  x = c(seq(1, 130,by = 1)),
  y = c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
        3,3,3,3,4,4,4,5,5,5,5,6,6,6,6,6,6,6,6,6,6,
        6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,10,10,10,
        10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,12,
        13,13,13,14,14,15,15,15,15,16,16,16,17,17,17,17,
        18,18,18,18,18,19,19,19,19,19,19,19,19,20,20,20,20,
        20,21,21,22,22,22,22,23,24,25,26,27,28,29,30,31,32)
)

ggplot(collectors_curve, aes(x = x, y = y)) +
         geom_point() +
         geom_smooth()+
  labs(x="Cumulative number of individuals classified", 
       y="Cumulative number of species observed", title = "Collector's Curve of a subsample of Bag #3")

invsimpson_MMs = 1/((13/44)^2+(8/44)^2+(4/44)^2+(3/44)^2+
                      (4/44)^2+(12/44)^2)
invsimpson_Skittles = 1/((6/35)^2+(7/35)^2+(4/35)^2+(10/35)^2+
                           (8/35)^2)
invsimpson_Gummy_Bears = 1/((1/17)^2+(3/17)^2+(2/17)^2+(4/17)^2+
                              (3/17)^2+(4/17)^2)
invsimpson_Rods = 1/((5/24)^2+(8/24)^2+(5/24)^2+(2/24)^2+
                       (4/24)^2)
invsimpson_Circular = 1
invsimpson_Filamentous = 1/((1/3)^2+(1/3)^2+(1/3)^2)
invsimpson_Brick = 1/((1/2)^2+(1/2)^2)
invsimpson_Round = 1/((1/4)^2+(1/4)^2+(1/4)^2+(1/4)^2)

invsimpson_subsample = 1/((13/130)^2+(8/130)^2+(4/130)^2+(3/130)^2+
                            (4/130)^2+(12/130)^2+(6/130)^2+(7/130)^2+
                            (4/130)^2+(10/130)^2+(8/130)^2+(1/130)^2+
                            (3/130)^2+(2/130)^2+(4/130)^2+(3/130)^2+(4/130)^2+
                            (5/130)^2+(8/130)^2+(5/130)^2+(2/130)^2+(4/130)^2+
                            (1/130)^2+(1/130)^2+(1/130)^2+(1/130)^2+(1/130)^2+
                            (1/130)^2+(1/130)^2+(1/130)^2+(1/130)^2)
## Singletons = 11; Species observed >2 times: 21 
chao1 = 32 + ((11^2)/(2*21))

#Using Vegan
library(vegan)

subsample_bag3_diversity = 
  subsample_bag3 %>%
  select(Species_number, Abundance) %>%
  spread(Species_number, Abundance)

diversity(subsample_bag3_diversity, index ="invsimpson")

specpool(subsample_bag3_diversity)