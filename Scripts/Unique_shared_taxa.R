
# Calculate ASVs, Species, Genera and Phyla overlap in HCB and MPYG plates #

#### Package setup ####

library(phyloseq)
library(microViz)
library(ggvenn)
library(gridExtra)
library(openxlsx)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

load("Results/RData/1.physeq.original.RData")  #unfiltered
load("Results/RData/2.phyloseq.filtered.RData") #filtered

# Create directories for results
dir.create("Results/6.Taxonomy", recursive = T)

#### End ####

#### Divide phyloseq object according to condition ####

phy_HCB <- prune_samples(sample_data(physeq)$Plate=="HCB", physeq)
phy_HCB <- subset_taxa(phy_HCB, rowSums(otu_table(phy_HCB)) > 0)

phy_MPYG <- prune_samples(sample_data(physeq)$Plate=="MPYG", physeq)
phy_MPYG <- subset_taxa(phy_MPYG, rowSums(otu_table(phy_MPYG)) > 0)

# Subset not NA taxa per condition

Sp_HCB <- as.data.frame(tax_table(phy_HCB))$Species[!is.na(as.data.frame(tax_table(phy_HCB))$Species)]
Sp_MPYG <- as.data.frame(tax_table(phy_MPYG))$Species[!is.na(as.data.frame(tax_table(phy_MPYG))$Species)]

Gen_HCB <- as.data.frame(tax_table(phy_HCB))$Genus[!is.na(as.data.frame(tax_table(phy_HCB))$Genus)]
Gen_MPYG <- as.data.frame(tax_table(phy_MPYG))$Genus[!is.na(as.data.frame(tax_table(phy_MPYG))$Genus)]

Phy_HCB <- as.data.frame(tax_table(phy_HCB))$Phylum[!is.na(as.data.frame(tax_table(phy_HCB))$Phylum)]
Phy_MPYG <- as.data.frame(tax_table(phy_MPYG))$Phylum[!is.na(as.data.frame(tax_table(phy_MPYG))$Phylum)]

#### End ####

#### Venn Diagrams ####

# Shared and unique ASVs
ASV <- ggvenn(
  list(HCB=rownames(tax_table(phy_HCB)),
       MPYG=rownames(tax_table(phy_MPYG))), 
  fill_color = c("#F8766D", "#00BFC4"),
  fill_alpha = 0.4,
  stroke_color = "gray",
  stroke_size = 0.5, set_name_size = 6) +
  ggtitle("ASV") +
  theme(plot.title = element_text(size= 20, hjust = 0.5, vjust = -9)) +
  annotate("rect",xmin=-1.8,xmax=1.8,ymin=1.5,ymax=2,alpha=.2,color="black",  fill = "grey1") 

# Shared and unique Species
Sp <- ggvenn(
  list(HCB=unique(Sp_HCB), MPYG=unique(Sp_MPYG)), 
  fill_color = c("#F8766D", "#00BFC4"),
  fill_alpha = 0.4,
  stroke_color = "gray",
  stroke_size = 0.5, set_name_size = 6) +
  ggtitle("Species") +
  theme(plot.title = element_text(size= 20, hjust = 0.5, vjust = -9)) +
  annotate("rect",xmin=-1.8,xmax=1.8,ymin=1.5,ymax=2,alpha=.2,color="black",  fill = "grey1")

# Shared and unique Genera
Gen <- ggvenn(
  list(HCB=unique(Gen_HCB), MPYG=unique(Gen_MPYG)), 
  fill_color = c("#F8766D", "#00BFC4"),
  fill_alpha = 0.4,
  stroke_color = "gray",
  stroke_size = 0.5, set_name_size = 6) +
  ggtitle("Genus") +
  theme(plot.title = element_text(size= 20, hjust = 0.5, vjust = -9)) +
  annotate("rect",xmin=-1.8,xmax=1.8,ymin=1.5,ymax=2,alpha=.2,color="black",  fill = "grey1")

# Shared and unique Phyla
Phy <- ggvenn(
  list(HCB=unique(Phy_HCB), MPYG=unique(Phy_MPYG)), 
  fill_color = c("#F8766D", "#00BFC4"),
  fill_alpha = 0.4,
  stroke_color = "gray",
  stroke_size = 0.5, set_name_size = 6) +
  ggtitle("Phylum") +
  theme(plot.title = element_text(size= 20, hjust = 0.5, vjust = -9)) +
  annotate("rect",xmin=-1.8,xmax=1.8,ymin=1.5,ymax=2,alpha=.2,color="black",  fill = "grey1")

# Common plot
grid.arrange(ASV, Sp, Gen, Phy, nrow=2)

#### End ####

#### Unique genera and species ####

# Unique Genera Relative
Gens <- gplots::venn(list(HCB=unique(Gen_HCB), MPYG=unique(Gen_MPYG)), show.plot = F)
Gens <- attr(Gens,"intersections")

# HCB
HCB_per <- transform_sample_counts(phy_HCB, function(OTU) OTU/sum(OTU))
HCB_per <- tax_fix(HCB_per)

Gen_HCB <- subset_taxa(HCB_per, tax_table(HCB_per)[,"Genus"] %in% Gens$HCB)

plot_bar(Gen_HCB, y =  "Abundance", title = "Relative abundance based on Genus", fill="Genus")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

sum = tapply(taxa_sums(Gen_HCB), tax_table(Gen_HCB)[, "Genus"], sum, na.rm=TRUE)
genusH = sort(sum, TRUE)

# MPYG
MPYG_per <- transform_sample_counts(phy_MPYG, function(OTU) OTU/sum(OTU))
MPYG_per <- tax_fix(MPYG_per)

Gen_MPYG <- subset_taxa(MPYG_per, tax_table(MPYG_per)[,"Genus"] %in% Gens$MPYG)

sum = tapply(taxa_sums(Gen_MPYG), tax_table(Gen_MPYG)[, "Genus"], sum, na.rm=TRUE)
genusM = sort(sum, TRUE)

plot_bar(Gen_MPYG, y =  "Abundance", title = "Abundance based on Genus", fill="Genus")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

# Unique Species Relative
Sps <- gplots::venn(list(HCB=unique(Sp_HCB), MPYG=unique(Sp_MPYG)), show.plot = F)
Sps <- attr(Sps,"intersections")

# HCB
HCB_per <- transform_sample_counts(phy_HCB, function(OTU) OTU/sum(OTU))
HCB_per <- tax_fix(HCB_per)

Sp_HCB <- subset_taxa(HCB_per, tax_table(HCB_per)[,"Species"] %in% Sps$HCB)

plot_bar(Sp_HCB, y =  "Abundance", title = "Relative abundance based on Species", fill="Species")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack")

sum = tapply(taxa_sums(Sp_HCB), tax_table(Sp_HCB)[, "Species"], sum, na.rm=TRUE)
topspH = sort(sum, TRUE)

# MPYG
MPYG_per <- transform_sample_counts(phy_MPYG, function(OTU) OTU/sum(OTU))
MPYG_per <- tax_fix(MPYG_per)

Sp_MPYG <- subset_taxa(MPYG_per, tax_table(MPYG_per)[,"Species"] %in% Sps$MPYG)

sum = tapply(taxa_sums(Sp_MPYG), tax_table(Sp_MPYG)[, "Species"], sum, na.rm=TRUE)
topspM = sort(sum, TRUE)

plot_bar(Sp_MPYG, y =  "Abundance", title = "Abundance based on Species", fill="Species")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack")

# Combine in excel
list(Genus_HCB=cbind(genusH), Genus_MPYG=genusM,
     Species_HCB=topspH, Species_MPYG=topspM)
write.xlsx(list(Genus_HCB=cbind(genusH), Genus_MPYG=cbind(genusM),
                Species_HCB=cbind(topspH), Species_MPYG=cbind(topspM)),
           "Results/6.Taxonomy/Unique_Gen_Sp.xlsx", rowNames = T)

#### End ####
