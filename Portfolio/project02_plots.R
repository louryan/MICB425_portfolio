library(tidyverse)
library(cowplot)


#Load DNA data
nirS.DNA.10m = read_tsv("nirS_DNA_10m_marker_contig_map.tsv") %>%
  select(Tax.DNA.10 = Confident_Taxonomy, Abund.DNA.10 = Abundance, Query)
nirS.DNA.100m = read_tsv("nirS_DNA_100m_marker_contig_map.tsv") %>%
  select(Tax.DNA.100 = Confident_Taxonomy, Abund.DNA.100 = Abundance, Query)
nirS.DNA.120m = read_tsv("nirS_DNA_120m_marker_contig_map.tsv") %>%
  select(Tax.DNA.120 = Confident_Taxonomy, Abund.DNA.120 = Abundance, Query)
nirS.DNA.135m = read_tsv("nirS_DNA_135m_marker_contig_map.tsv") %>%
  select(Tax.DNA.135 = Confident_Taxonomy, Abund.DNA.135 = Abundance, Query)
nirS.DNA.150m = read_tsv("nirS_DNA_150m_marker_contig_map.tsv") %>%
  select(Tax.DNA.150 = Confident_Taxonomy, Abund.DNA.150 = Abundance, Query)
nirS.DNA.165m = read_tsv("nirS_DNA_165m_marker_contig_map.tsv") %>%
  select(Tax.DNA.165 = Confident_Taxonomy, Abund.DNA.165 = Abundance, Query)
nirS.DNA.200m = read_tsv("nirS_DNA_200m_marker_contig_map.tsv") %>%
  select(Tax.DNA.200 = Confident_Taxonomy, Abund.DNA.200 = Abundance, Query)

#Load RNA data
nirS.RNA.10m = read_tsv("nirS_RNA_10m_marker_contig_map.tsv") %>%
  select(Tax.RNA.10 = Confident_Taxonomy, Abund.RNA.10 = Abundance, Query)
nirS.RNA.100m = read_tsv("nirS_RNA_100m_marker_contig_map.tsv") %>%
  select(Tax.RNA.100 = Confident_Taxonomy, Abund.RNA.100 = Abundance, Query)
nirS.RNA.120m = read_tsv("nirS_RNA_120m_marker_contig_map.tsv") %>%
  select(Tax.RNA.120 = Confident_Taxonomy, Abund.RNA.120 = Abundance, Query)
nirS.RNA.135m = read_tsv("nirS_RNA_135m_marker_contig_map.tsv") %>%
  select(Tax.RNA.135 = Confident_Taxonomy, Abund.RNA.135 = Abundance, Query)
nirS.RNA.150m = read_tsv("nirS_RNA_150m_marker_contig_map.tsv") %>%
  select(Tax.RNA.150 = Confident_Taxonomy, Abund.RNA.150 = Abundance, Query)
nirS.RNA.165m = read_tsv("nirS_RNA_165m_marker_contig_map.tsv") %>%
  select(Tax.RNA.165 = Confident_Taxonomy, Abund.RNA.165 = Abundance, Query)
nirS.RNA.200m = read_tsv("nirS_RNA_200m_marker_contig_map.tsv") %>%
  select(Tax.RNA.200 = Confident_Taxonomy, Abund.RNA.200 = Abundance, Query)

#Combine all data frames - 10m is excluded as Relative Abundance = 0
nirS.all = nirS.DNA.100m %>%
  full_join(nirS.DNA.120m, by = "Query") %>%
  full_join(nirS.DNA.135m, by = "Query") %>%  
  full_join(nirS.DNA.150m, by = "Query") %>%  
  full_join(nirS.DNA.165m, by = "Query") %>%
  full_join(nirS.DNA.200m, by = "Query") %>%
  full_join(nirS.RNA.100m, by = "Query") %>%
  full_join(nirS.RNA.120m, by = "Query") %>%
  full_join(nirS.RNA.135m, by = "Query") %>%
  full_join(nirS.RNA.150m, by = "Query") %>%
  full_join(nirS.RNA.165m, by = "Query") %>%
  full_join(nirS.RNA.200m, by = "Query") %>%
  mutate(Taxonomy = ifelse(!is.na(Tax.DNA.100), Tax.DNA.100,
                    ifelse(!is.na(Tax.DNA.120), Tax.DNA.120,
                    ifelse(!is.na(Tax.DNA.135), Tax.DNA.135,
                    ifelse(!is.na(Tax.DNA.150), Tax.DNA.150,
                    ifelse(!is.na(Tax.DNA.165), Tax.DNA.165,
                    ifelse(!is.na(Tax.DNA.200), Tax.DNA.200,
                    ifelse(!is.na(Tax.RNA.100), Tax.RNA.100,
                    ifelse(!is.na(Tax.RNA.120), Tax.RNA.120,                           
                    ifelse(!is.na(Tax.RNA.135), Tax.RNA.135,
                    ifelse(!is.na(Tax.RNA.150), Tax.RNA.150,
                    ifelse(!is.na(Tax.RNA.165), Tax.RNA.165,
                    ifelse(!is.na(Tax.RNA.200), Tax.RNA.200, 
             "unclassified"))))))))))))) %>%
  select(-starts_with("Tax.")) %>%
  gather("Key", "Abundance", starts_with("Abund")) %>%
  separate(Key, c("Key","Type","Depth_m"), by = ".") %>%
  select(Depth_m, Type, Abundance, Taxonomy, Query) %>%
  mutate(Depth_m = as.numeric(Depth_m)) %>%
  separate(Taxonomy, into = c("Domain", "Phylum", "Class", "Order", 
                              "Family","Genus", "Species"), sep=";")

#Combine and sum all abundances at a certain depth
nirS.abund = aggregate(Abundance ~ Type + Depth_m, data = nirS.all, sum)

#Plot DNA vs. RNA Abundance
ggplot(nirS.abund, aes(x = Abundance, y = Depth_m)) +
  geom_path(aes(colour = Type), size = 1.5) +
  scale_y_reverse(lim=c(200,0)) +
  ggtitle("Abundance of nirS gene (DNA vs. RNA) at different depths") + 
  theme(plot.title = element_text(face = "bold", hjust=0.5)) +
  xlab("Relative Abundance") +
  ylab("Depth (m)") 

#Load phyloseq data
load("mothur_phyloseq.RData")
metadata = data.frame(mothur@sam_data)
#Select NO3, NO2, O2 + set line colours
mn <- metadata %>% 
  gather(Nutrient, Concentration, NO3_uM, NO2_uM, O2_uM) %>% 
  select(Depth_m, Nutrient, Concentration)
colours = c("NO3_uM" = "green", "NO2_uM" = "red", "O2_uM" = "blue") 

#Plot by Order
plot1 = nirS.all %>% 
  mutate(Order = ifelse(is.na(Order), "unclassified", Order)) %>% 
  ggplot(aes(x = Order, y = Depth_m)) +
  geom_point(aes(size = ifelse(Abundance == 0, NA, Abundance), colour = Type), 
             position = position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 12, hjust=0.5))+
  theme(axis.text.x = element_text(angle=-90,hjust=1,vjust=0.5))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_y_reverse(lim=c(200,10)) +
  scale_size_continuous(name = "Abundance")+
  labs(y = "Depth (m)",
       x = "Order")

#Plot phyloseq metadata
plot2 = ggplot()+
  geom_point(data = mn, aes(x = log10(Concentration), y = Depth_m, colour = Nutrient))+
  geom_path(data = mn,aes(x = log10(Concentration), y = Depth_m, colour = Nutrient)) +
  scale_color_manual(values=colours, 
                     labels=c(expression(paste("NO"[2])),
                              expression(paste("NO"[3])), 
                              expression(paste("O"[2]))))+
  scale_y_reverse(lim=c(200,10))+
  labs(y = "Depth (m)",
       x = expression(paste("Log"[10], " Concentration", " (μM)")))

plot_grid(plot2, plot1, labels=c("A","B"), rel_widths = c(1/4, 3/4), align = 'h')


#Plot by Class
plot3 = nirS.all %>% 
  mutate(Class = ifelse(is.na(Class), "unclassified", Class)) %>% 
  ggplot(aes(x = Class, y = Depth_m)) +
  geom_point(aes(size = ifelse(Abundance == 0, NA, Abundance), colour = Type), 
             position = position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 12, hjust=0.5))+
  theme(axis.text.x = element_text(angle=-90,hjust=1,vjust=0.5))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_y_reverse(lim=c(200,10)) +
  scale_size_continuous(name = "Abundance")+
  labs(y = "Depth (m)",
       x = "Class")

plot_grid(plot2, plot3, labels=c("A","B"), rel_widths = c(1/4, 3/4), align = 'h')

#Plot by Family
plot4 = nirS.all %>% 
  mutate(Class = ifelse(is.na(Family), "unclassified", Family)) %>% 
  ggplot(aes(x = Family, y = Depth_m)) +
  geom_point(aes(size = ifelse(Abundance == 0, NA, Abundance), colour = Type), 
             position = position_dodge(0.5)) +
  theme(axis.text.x = element_text(size = 12, hjust=0.5))+
  theme(axis.text.x = element_text(angle=-90,hjust=1,vjust=0.5))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  scale_y_reverse(lim=c(200,10)) +
  scale_size_continuous(name = "Abundance")+
  labs(y = "Depth (m)",
       x = "Family")

plot_grid(plot2, plot4, labels=c("A","B"), rel_widths = c(1/4, 3/4), align = 'h')
