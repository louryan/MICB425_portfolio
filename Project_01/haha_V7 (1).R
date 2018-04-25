# Libraries
library(tidyverse)
library(phyloseq)
library(kableExtra)
library(vegan)
#library(magrittr) pretty sure it's unnecessary, it's as a part of the tidyverse no?

# Loading in Data
load("mothur_phyloseq.RData")
load("qiime2_phyloseq.RData")

# Things
mothur_sam = 
  mothur@sam_data %>% 
  rownames_to_column("Sample") %>% 
  data.frame() %>%
  na.omit()

mothur_percent = 	mothur %>% transform_sample_counts(function(x) 100 * x/sum(x)) 

mothur_percent %>% 
	plot_bar(fill="Phylum") + 
	geom_bar(aes(fill=Phylum), stat="identity") +
	labs(x = "Sample depth", y = "Relative abundance (%)", 
			 title = "Phyla from 10 to 200 m in Saanich Inlet")

sample_data(mothur) %>% 
	select(matches("O2_uM|depth"),-matches("NO2|SiO2")) %>%
	gather(key = "Nutrient", value = "Concentration", -Depth_m) %>%
	ggplot(., aes(x = Depth_m, y = Concentration)) +
	geom_point() +
	geom_line() +
	facet_wrap( ~ Nutrient, scales = "free") +
	theme(legend.position = "none") +
	labs(x = "Depth (m)", y = expression(paste("Concentration (", mu, "M)")))

qiime2_percent = 	qiime2 %>% transform_sample_counts(function(x) 100 * x/sum(x)) 

qiime2_percent %>% 
	plot_bar(fill="Phylum") + 
	geom_bar(aes(fill=Phylum), stat="identity") +
	labs(x = "Sample depth", y = "Relative abundance (%)", 
			 title = "Phyla from 10 to 200 m in Saanich Inlet")

mothur %>% 
	plot_heatmap(max.label = 30, sample.order = mothur_sam[,1]) +
	theme(axis.text.y = element_blank(), 
				axis.ticks.y = element_blank())

qiime2 %>% 
	plot_heatmap(max.label = 30, sample.order = mothur_sam[,1]) +
	theme(axis.text.y = element_blank(), 
				axis.ticks.y = element_blank())

#mothur_fr =
#	mothur %>% 
#	prune_taxa(as.vector(!grepl("uncultured|unclassified|Candidatus",
#															tax_table(.)[,"Genus"])), .) %>% 
#	filter_taxa(., function(x) sum (x >= 1) >= 3, TRUE)

#qiime2_fr = 
#	qiime2 %>% 
#	prune_taxa(as.vector(substr(tax_table(.)[,"Genus"], 0, 8)!="D_5__unc"), .) %>% 
#	prune_taxa(as.vector(tax_table(.)[,"Family"]!="D_4__Unknown Family"), .) %>% 
#	prune_taxa(as.vector(tax_table(.)[,"Genus"]!="D_5__"), .) %>% 
#	filter_taxa(., function(x) sum (x >= 1) >= 3, TRUE)

qiime2 %>% 
	subset_taxa(Genus=="D_5__Sulfurimonas") %>% 
	prune_samples(sample_sums(.)>=1, .) %>% 
	plot_heatmap(taxa.label="Genus")

mothur %>% 
	subset_taxa(Genus=="Sulfurimonas") %>% 
	prune_samples(sample_sums(.)>=1, .) %>% 
	plot_heatmap(taxa.label="Genus")

# !!!
# Must rarify data before estimating and plotting richness 
# Unless this is automatically done by phyloseq (unlikely)

mothur_alpha = 
	mothur %>% 
	prune_taxa(taxa_sums(.) > 0, .) %>% 
	estimate_richness(., measures = c("Chao1", "Shannon", "InvSimpson")) %>% 
	rownames_to_column(., var = "Sample") 

qiime2_alpha = 
	qiime2 %>% 
	prune_taxa(taxa_sums(.) > 0, .) %>% 
	estimate_richness(., measures = c("Chao1", "Shannon", "InvSimpson")) %>% 
	rownames_to_column(., var = "Sample") 

# Alternate colour scheme that isn't R/G color-blind ambiguous
#	scale_colour_continuous(low="yellow", high="blue") +

mothur %>% 
	prune_taxa(taxa_sums(.) > 0, .) %>% 
	plot_richness(x = "Depth_m", measures=c("Chao1", "InvSimpson")) +
	geom_point(aes(colour = select(data.frame(mothur@sam_data), 
																 starts_with("O2_uM"))), na.rm = TRUE) +
	scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
	labs(x = "Depth (m)", y = "Alpha Diversity Indicator Value", 
			 colour = expression(paste("[O"[2], "] (", mu, "M)")))

qiime2 %>% 
	prune_taxa(taxa_sums(.) > 0, .) %>% 
	plot_richness(x = "Depth_m", measures=c("Chao1", "InvSimpson")) +
	geom_point(aes(colour = select(data.frame(qiime2@sam_data), 
																 starts_with("O2_uM"))), na.rm = TRUE) +
	scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
	labs(x = "Depth (m)", y = "Alpha Diversity Indicator Value", 
			 colour = expression(paste("[O"[2], "] (", mu, "M)")))

mothur_otus = 
	mothur %>% 
	subset_taxa(Genus=="Sulfurimonas") %>% 
	.@otu_table %>% 
	data.frame() 

qiime2_asvs = 
	qiime2 %>% 
	subset_taxa(Genus=="D_5__Sulfurimonas") %>% 
	.@otu_table %>% 
	data.frame() %>%
	t() %>% 
	data.frame()

mothur_totals =
	mothur %>% 
	.@otu_table %>% 
	data.frame() %>% 
	rowSums() %>% 
	data.frame() %>% 
	rename(., TotAbundance = ".") %>% 
	rownames_to_column(., "Sample") 

qiime2_totals =
	qiime2 %>% 
	.@otu_table %>% 
	data.frame() %>% 
	t() %>% 
	data.frame() %>% 
	rowSums() %>% 
	data.frame() %>% 
	rename(., TotAbundance = ".") %>% 
	rownames_to_column(., "Sample") 


mothur_taxon =
	mothur_otus %>% 
	rowSums(.) %>% 
	data.frame() %>% 
	rename(., TaxAbundance = ".") %>% 
	rownames_to_column(., var = "Sample") %>% 
	merge(., mothur_sam, "Sample") %>% 
	merge(., mothur_totals, "Sample") 

qiime2_taxon =
	qiime2_asvs %>% 
	rowSums(.) %>% 
	data.frame() %>% 
	rename(., TaxAbundance = ".") %>% 
	rownames_to_column(., "Sample") %>% 
	merge(., mothur_sam, "Sample") %>% 
	merge(., qiime2_totals, "Sample") 

#linear model  
# !!! needs to be FDR-corrected
summary(lm(formula = TaxAbundance ~ O2_uM, data = mothur_taxon))
summary(lm(formula = TaxAbundance ~ O2_uM, data = qiime2_taxon))
summary(lm(formula = TaxAbundance ~ Depth_m, data = mothur_taxon)) # verge of significance
summary(lm(formula = TaxAbundance ~ Depth_m, data = qiime2_taxon)) # significant
summary(lm(formula = TaxAbundance ~ H2S_uM, data = mothur_taxon)) # significant
summary(lm(formula = TaxAbundance ~ H2S_uM, data = qiime2_taxon)) # significant

#LM graphs of how taxon differs in abundance with depth and oxygen concentration

ggplot(mothur_taxon, aes(x=O2_uM, y=TaxAbundance)) + 
    geom_point(aes(colour = select(data.frame(mothur_taxon), 
          starts_with("Depth_m"))), na.rm = TRUE) +
    geom_smooth(method=lm)+
    scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
    labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "Taxon Abundance", 
         colour = "Depth (m)") 
  
ggplot(qiime2_taxon, aes(x=O2_uM, y=TaxAbundance)) + 
  geom_point(aes(colour = select(data.frame(qiime2_taxon), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "Taxon Abundance", 
       colour = "Depth (m)") 

ggplot(mothur_taxon, aes(x=Depth_m, y=TaxAbundance)) + 
  geom_point(aes(colour = select(data.frame(mothur_taxon), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "Taxon Abundance", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_taxon, aes(x=Depth_m, y=TaxAbundance)) + 
  geom_point(aes(colour = select(data.frame(qiime2_taxon), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "Taxon Abundance", 
       colour = expression(paste("[O"[2], "] (", mu, "M)")))  

ggplot(mothur_taxon, aes(x=Depth_m, y=TaxAbundance)) + 
	geom_point(aes(colour = select(data.frame(mothur_taxon), 
																 starts_with("H2S_uM"))), na.rm = TRUE) +
	geom_smooth(method=lm)+
	scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
	labs(x = "Depth (m)", y = "Taxon Abundance", 
			 colour = expression(paste("[H"[2], "S] (", mu, "M)"))) 

ggplot(qiime2_taxon, aes(x=H2S_uM, y=TaxAbundance)) + 
	geom_point(aes(colour = select(data.frame(qiime2_taxon), 
																 starts_with("H2S_uM"))), na.rm = TRUE) +
	geom_smooth(method=lm)+
	scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
	labs(x = "Depth (m)", y = "Taxon Abundance", 
			 colour = expression(paste("[H"[2], "S] (", mu, "M)")))  

#set up for LM

mothur_otus_total =
  mothur_otus %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Sample") %>% 
  merge(., mothur_taxon, "Sample")

qiime2_asvs_total =
  qiime2_asvs %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Sample") %>% 
  merge(., qiime2_taxon, "Sample")

#linear model
#everything is no longer significant once corrected for multiple-comparisons testing
  #mothur
summary(lm(formula = Otu0308 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu0308 ~ Depth_m, data = mothur_otus_total)) #verge of significance
summary(lm(formula = Otu0666 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu0666 ~ Depth_m, data = mothur_otus_total)) #significant
summary(lm(formula = Otu0704 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu0704 ~ Depth_m, data = mothur_otus_total))
summary(lm(formula = Otu0751 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu0751 ~ Depth_m, data = mothur_otus_total))
summary(lm(formula = Otu1315 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu1315 ~ Depth_m, data = mothur_otus_total))
summary(lm(formula = Otu2793 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu2793 ~ Depth_m, data = mothur_otus_total))
summary(lm(formula = Otu3512 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu3512 ~ Depth_m, data = mothur_otus_total))
summary(lm(formula = Otu3610 ~ O2_uM, data = mothur_otus_total))
summary(lm(formula = Otu3610 ~ Depth_m, data = mothur_otus_total))

  #qiime2
summary(lm(formula = Asv277 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv277 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv561 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv561 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv578 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv578 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv1153 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv1153 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv1250 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv1250 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv1620 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv1620 ~ Depth_m, data = qiime2_asvs_total))
summary(lm(formula = Asv2216 ~ O2_uM, data = qiime2_asvs_total))
summary(lm(formula = Asv2216 ~ Depth_m, data = qiime2_asvs_total))


# Results 4: graphs of how each OTU or ASV differs in abundance with depth and oxygen concentration, using colourful dot graph

mothur_percent %>% 
  subset_taxa(Genus=="Sulfurimonas") %>%
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=Depth_m, y=OTU, size=Abundance, color=OTU)) + 
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of OTUs within Sulfurimonas genus across depth")

mothur_percent %>% 
  subset_taxa(Genus=="Sulfurimonas") %>%
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=O2_uM, y=OTU, size=Abundance, color=OTU)) + 
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of OTUs within Sulfurimonas genus across oxygen concentration")

qiime2_percent %>% 
  subset_taxa(Genus=="D_5__Sulfurimonas") %>%
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=Depth_m, y=OTU, size=Abundance, color=OTU)) + 
  scale_size_continuous(range = c(0,5)) +
  labs(y="ASV",title="Abundance of ASVs within Sulfurimonas genus across depth", color="ASV")

qiime2_percent %>% 
  subset_taxa(Genus=="D_5__Sulfurimonas") %>%
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=O2_uM, y=OTU, size=Abundance, color=OTU)) + 
  scale_size_continuous(range = c(0,5)) +
  labs(y="ASV",title="Abundance of ASVs within Sulfurimonas genus across oxygen concentration", color="ASV")

# Results 4: graphs of how each OTU or ASV differs in abundance with depth and oxygen concentration, displaying LM

mothur_percent %>% 
  subset_taxa(Genus=="Sulfurimonas") %>% 
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance)) +
  geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of OTUs within Sulfurimonas genus across depth")

mothur_percent %>% 
  subset_taxa(Genus=="Sulfurimonas") %>% 
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=O2_uM, y=Abundance)) +
  geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of OTUs within Sulfurimonas genus across oxygen concentration")

qiime2_percent %>% 
  subset_taxa(Genus=="D_5__Sulfurimonas") %>% 
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance)) +
  geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of ASVs within Sulfurimonas genus across depth")

qiime2_percent %>% 
  subset_taxa(Genus=="D_5__Sulfurimonas") %>% 
  psmelt() %>% 
  ggplot() +
  geom_point(aes(x=O2_uM, y=Abundance)) +
  geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of ASVs within Sulfurimonas genus across oxygen concentration")

#Result 4 supplemental: use this to look at detailed LM graphs of how each OTU or ASV differs in abundance with depth and oxygen concentration
  #mothur OTUs

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu0308)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU0308",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu0308)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU0308", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu0666)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU0666",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu0666)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU0666", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu0704)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU0704",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu0704)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU0704", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu0751)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU0751",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu0751)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU0751", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu1315)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU1315",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu1315)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU1315", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu2793)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU2793",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu2793)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU2793", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu3512)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU3512",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu3512)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU3512", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(mothur_otus_total, aes(x=O2_uM, y=Otu3610)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "OTU Abundance", title ="OTU3610",
       colour = "Depth (m)") 

ggplot(mothur_otus_total, aes(x=Depth_m, y=Otu3610)) + 
  geom_point(aes(colour = select(data.frame(mothur_otus_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "OTU Abundance",title ="OTU3610", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

#qiime2 ASVs

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv277)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV277",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv277)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV277", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv561)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV561",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv561)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV561", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv578)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV578",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv578)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV578", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv1153)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV1153",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv1153)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV1153", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv1250)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV1250",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv1250)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV1250", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv1620)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV1620",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv1620)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV1620", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 

ggplot(qiime2_asvs_total, aes(x=O2_uM, y=Asv2216)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("Depth_m"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = expression(paste("[O"[2], "] (", mu, "M)")), y = "ASV Abundance", title ="ASV2216",
       colour = "Depth (m)") 

ggplot(qiime2_asvs_total, aes(x=Depth_m, y=Asv2216)) + 
  geom_point(aes(colour = select(data.frame(qiime2_asvs_total), 
                                 starts_with("O2_uM"))), na.rm = TRUE) +
  geom_smooth(method=lm)+
  scale_colour_distiller(palette = "RdYlGn") + # replace if using alt colours
  labs(x = "Depth (m)", y = "ASV Abundance",title ="ASV2216", 
       colour = expression(paste("[O"[2], "] (", mu, "M)"))) 
