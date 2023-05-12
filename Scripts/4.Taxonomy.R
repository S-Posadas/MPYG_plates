
# 6.Taxonomy analysis #

# Statistics, bar plots and bubble plots

#### Package setup ####

library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(gridExtra)
library(cowplot)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

load("Results/RData/2.phyloseq.filtered.RData")

# Create directories for results
dir.create("Results/6.Taxonomy/Bar_plots", recursive = T)
dir.create("Results/6.Taxonomy/Bubble_plots")
dir.create("Results/6.Taxonomy/Stats")

#### End ####

#### Bar plots ####

# Assign upper taxonomic ranks to unassigned taxa

physeq <- tax_fix(physeq)

# Change Species names for visualization

tax <- as.data.frame(tax_table(physeq))
tax[,"Species"] <- tax[,"Species"] %>% gsub("\\(RS_.*", "", .) %>% gsub("\\(GB_.*", "", .) %>% gsub("_RS_GCF_.*", "", .)
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq) <- tax; rm(tax)

# Normalize to get relative abundance
physeq_per <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

# Absolute abundance

plot_bar(physeq, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")
plot_bar(physeq, y =  "Abundance", title = "Abundance based on Class", fill="Class")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Class, fill=Class), stat = "Identity", position = "stack")

# Relative abundance

plot_bar(physeq_per, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(~Plate, scales="free_x") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")
plot_bar(physeq_per, y =  "Abundance", title = "Abundance based on Class", fill="Class")+ facet_wrap(~Plate, scales="free_x") +
  geom_bar(aes(color=Class, fill=Class), stat = "Identity", position = "stack")
plot_bar(physeq_per, y =  "Abundance", title = "Abundance based on Order", fill="Order")+ facet_wrap(~Plate, scales="free_x") +
  geom_bar(aes(color=Order, fill=Order), stat = "Identity", position = "stack")
plot_bar(physeq_per, y =  "Abundance", title = "Abundance based on Genus", fill="Genus")+ facet_wrap(~Plate, scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

#### End ####

#### Plot abundance from top ASVs ####
# Top 20/50/100 ASVs

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]  # adjust number to wished top ones
top50 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:50]
top100 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:100]

ps.top20_per <- prune_taxa(top20, physeq_per)
ps.top20 <- prune_taxa(top20, physeq)

ps.top50_per <- prune_taxa(top50, physeq_per)
ps.top50 <- prune_taxa(top50, physeq)

ps.top100_per <- prune_taxa(top100, physeq_per)
ps.top100 <- prune_taxa(top100, physeq)

# Relative abundance

plot_bar(ps.top100_per, y =  "Abundance", title = "Abundance based on top 100 ASV", fill="Species")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack")

plot_bar(ps.top100_per, y =  "Abundance", title = "Abundance based on top 100 ASV", fill="Genus")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

plot_bar(ps.top100_per, y =  "Abundance", title = "Abundance based on top 100 ASV", fill="Phylum")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")


# Absolute abundance

plot_bar(ps.top100, y =  "Abundance", title = "Abundance based on top 100 ASV", fill="Genus")+ facet_wrap(~Plate, scales="free_x") +
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

plot_bar(ps.top100, y =  "Abundance", title = "Abundance based on top 100 ASV", fill="Phylum")+ facet_wrap(~Plate, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")

#### End ####

#### Plot abundance per condition ####

# Transform all variables to factors just in case...
df <- as.data.frame(lapply(sample_data(physeq_per),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq_per)
sample_data(physeq_per) <- sample_data(df)

ps.Plate <- merge_samples(physeq_per, "Plate")  # Merging 
ps.Plate_per <- transform_sample_counts(ps.Plate, function(OTU) OTU/sum(OTU)) #2nd transformation to make it again in percentage
sample_data(ps.Plate_per)$Plate <- rownames(sample_data(ps.Plate_per))

ps.Plate_plot = plot_bar(ps.Plate_per, y =  "Abundance", fill="Phylum") +
  facet_wrap(~Plate, scales="free_x")

# Bar plot paper 1
BP1 <- ps.Plate_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 


# Plot top 20 ASVs per condition

ps.Plate_per.20 <- names(sort(taxa_sums(ps.Plate_per), decreasing=TRUE))[1:20]  # adjust number to wished top ones
ps.Plate_per.20 <- prune_taxa(ps.Plate_per.20, ps.Plate_per)

ps.Plate_plot = plot_bar(ps.Plate_per.20, y =  "Abundance", title = "Abundance based on top 20 ASV", fill="Species")+ facet_wrap(~Plate, scales="free_x")

ps.Plate_plot + geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + theme_bw()

#### End ####

#### Plot abundance of different conditions per taxa ####

# Relative abundance

# Bar plot paper 2
BP2 <- plot_bar(ps.Plate_per, "Plate", fill="Plate", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(text = element_text(size=18),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 26)) 

plot_bar(ps.Plate_per, "Phylum", fill="Phylum", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack") +
  theme(text = element_text(size=18)) 

plot_bar(ps.top20_per, "Plate", fill="Plate", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack")  + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30)) +
  theme(text = element_text(size=18)) 
plot_bar(ps.top20_per, "Genus", fill="Genus", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30)) +
  theme(text = element_text(size=18)) 

# Plot top 15 genera
sum = tapply(taxa_sums(ps.Plate_per), tax_table(physeq_per)[, "Genus"], sum, na.rm=TRUE)
top10 = names(sort(sum, TRUE))[1:15]
topgenus = prune_taxa((tax_table(ps.Plate_per)[, "Genus"] %in% top10), ps.Plate_per)

BP3 <- plot_bar(topgenus, "Genus", fill="Genus", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30),
        text = element_text(size=18)) 

plot_bar(topgenus, "Plate", fill="Plate", facet_grid=~Genus) + geom_bar(stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(text = element_text(size=18),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 26)) 

# Absolute abundance

plot_bar(physeq, "Phylum", fill="Plate", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack")
plot_bar(ps.top20, "Plate", fill="Plate", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack")  + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))
plot_bar(ps.top20, "Phylum", fill="Phylum", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack")
plot_bar(ps.top20, "Genus", fill="Genus", facet_grid=~Plate) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

#### End ####

#### Taxonomy statistis ####
# Select phyloseq object to plot/perform stats                                  # Edit

phy_stats <- physeq_per             #Relative abundance phyloseq object

# Function for creating data frame for stats

fseltax <- function(taxon){

  # Get taxa to plot that are present in the phyloseq object
  sel_tax <- unique(tax_table(phy_stats)[,taxon])
  sel_tax <- as.vector(sel_tax)
  
  # Create relative abundance matrix for each group and specific taxon level
  tax_sums <- matrix(nrow=length(sample_names(phy_stats)), ncol=length(sel_tax)) # define matrix dimensions like this
  rownames(tax_sums) <- sample_names(phy_stats) # specify row names
  colnames(tax_sums) <- sel_tax # specify column names
  
  for (i in sample_names(phy_stats)) { 
    for (j in sel_tax) {
      tax_sums[[i,j]] <- sum(otu_table(phy_stats)[which(tax_table(phy_stats)[,taxon]==j), sample_names(phy_stats)==i])
      # this is calculating the average percent abundance of each taxon
    }
  }
  
  # Merge with metadata
  seltax <- merge(as.data.frame(tax_sums), sample_data(phy_stats), by=0)
  return(seltax)
  
}

seltax_phy <- fseltax("Phylum")
seltax_ord <- fseltax("Order")
seltax_fam <- fseltax("Family")
seltax_gen <- fseltax("Genus")
seltax_sp <- fseltax("Species")

# Save tables used for statistics

write.xlsx(seltax_phy, "Results/6.Taxonomy/Stats/Phylum_total_abundance.xlsx")
write.xlsx(seltax_ord, "Results/6.Taxonomy/Stats/Order_total_abundance.xlsx")
write.xlsx(seltax_fam, "Results/6.Taxonomy/Stats/Family_total_abundance.xlsx")
write.xlsx(seltax_gen, "Results/6.Taxonomy/Stats/Genus_total_abundance.xlsx")
write.xlsx(seltax_sp, "Results/6.Taxonomy/Stats/Species_total_abundance.xlsx")

# Function for calculating Wilcoxon and Mann-Whitney-U-test

stats_tax <- function(taxon, seltax_stats){
  
taxa_stats <- colnames(seltax_stats)[!(colnames(seltax_stats) %in%
                                         c("Row.names", colnames(sample_data(phy_stats))))]
  
  
  # Mann-Whitney-U-Test -----------------------------------------------------
  
  # U-Test 
  mwut <-list()
  for (g in taxa_stats){
    U_test <- wilcox.test(seltax_stats[,g] ~ seltax_stats[,"Plate"], data = seltax_stats, paired = F, exact = T,
                          p.adjust.methods = "BH")
    z <- abs(qnorm(U_test$p.value/2))
    r <- z/sqrt(nrow(seltax_stats))
    tab <- c(U_test$method, paste("data:", g, "and Plate"), U_test$statistic, U_test$p.value, r)
    mwut[[g]] <- tab
  }
  mwut <- as.data.frame(do.call(rbind, mwut))
  colnames(mwut) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")
  mwut$p.adj <- p.adjust(mwut$`p-value`, method = "BH")
  
  mwut <- mwut[order(mwut$`p-value`),]
  
  write.xlsx(mwut, paste0("Results/6.Taxonomy/Stats/", taxon, "_U-test_rel_ab.xlsx"))
  
  # Wilcoxon Test -----------------------------------------------------------
  
  wt <-list()
  for (g in taxa_stats){
    W_test <- wilcox.test(seltax_stats[,g] ~ seltax_stats[,"Plate"], data = seltax_stats, paired = T, exact = T,
                          p.adjust.methods = "BH")
    z <- abs(qnorm(W_test$p.value/2))
    r <- z/sqrt(nrow(seltax_stats))
    tab <- c(W_test$method, paste("data:", g, "and Plate"), W_test$statistic, W_test$p.value, r)
    wt[[g]] <- tab
  }
  wt <- as.data.frame(do.call(rbind, wt))
  colnames(wt) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")
  wt$p.adj <- p.adjust(wt$`p-value`, method = "BH")
  wt <- wt[order(wt$`p-value`),]
  
  write.xlsx(wt, paste0("Results/6.Taxonomy/Stats/", taxon, "_Wilcoxon_rel_ab.xlsx"))
  return(list(mwut, wt))
  
}

stats_tax("Phylum", seltax_phy)
stats_tax("Order", seltax_ord)
stats_tax("Family", seltax_fam)
stats_tax("Genus", seltax_gen)
stats_tax("Species", seltax_sp)

#### End ####

#### Taxonomy per condition ####
# Same as in prepare for stats but only to get a table of rel abundance per condition

# Select phyloseq object to plot/perform stats                                  # Edit

phy_stats <- ps.Plate_per             #Relative abundance phyloseq object

# Transpose rows and columns for the function to work

otu_table(phy_stats) <- t(otu_table(phy_stats))

# Function for creating data frame 
fseltax_cond <- function(taxon){
  
  # Apply fseltax function
  seltax <- fseltax(taxon)
  
  # Modify the data frame to look how we want
  rownames(seltax) <- seltax$Row.names 
  seltax <- seltax[, !(colnames(seltax) %in% c("Row.names", "Plate", "Sample_ID"))]
  seltax <- as.data.frame(t(seltax))
  seltax$Total <- rowSums(seltax)
  seltax <- seltax[order(seltax$Total, decreasing = T),]
  
  return(seltax)
  
}

# Create data frame
seltax_phy_cond <- fseltax_cond("Phylum")
seltax_ord_cond <- fseltax_cond("Order")
seltax_fam_cond <- fseltax_cond("Family")
seltax_gen_cond <- fseltax_cond("Genus")
seltax_sp_cond <- fseltax_cond("Species")

write.xlsx(seltax_phy_cond, "Results/6.Taxonomy/Stats/Phylum_rel_abundance_Plate.xlsx", rowNames=T)
write.xlsx(seltax_ord_cond, "Results/6.Taxonomy/Stats/Order_rel_abundance_Plate.xlsx", rowNames=T)
write.xlsx(seltax_ord_cond, "Results/6.Taxonomy/Stats/Family_rel_abundance_Plate.xlsx", rowNames=T)
write.xlsx(seltax_gen_cond, "Results/6.Taxonomy/Stats/Genus_rel_abundance_Plate.xlsx", rowNames=T)
write.xlsx(seltax_sp_cond, "Results/6.Taxonomy/Stats/Species_rel_abundance_Plate.xlsx", rowNames=T)

#### End ####

#### Bubble plot ####

# Select specific taxa to plot

ord <- c("Bacteroidales", "Peptostreptococcales", "Lactobacillales",
         "Oscillospirales", "Erysipelotrichales", "Coriobacteriales")

genera <- c("Lactobacillus", "Ligilactobacillus", "Limosilactobacillus",
            "Merdibacter", "Bacteroides", "Butyricimonas", "Clostridioides",
            "Odoribacter", "Lawsonibacter", "Erysipelatoclostridium")

sp <- c("Bacteroides thetaiotaomicron", "Ligilactobacillus hayakitensis",
        "Limosilactobacillus fermentum","Merdibacter massiliensis",
        "Lawsonibacter sp900066645")

# Prepare data for bubble plot

# Function for creating long data frames for plotting
long_df <- function(seltax){
  
  seltax$tax <- rownames(seltax)
  seltax$Plate <- NA
  
  seltax1 <- seltax[,c("HCB", "Plate", "tax")]
  seltax1$Plate <- "HCB"
  seltax2 <- seltax[,c("MPYG", "Plate", "tax")]
  seltax2$Plate <- "MPYG"
  
  colnames(seltax1) <- colnames(seltax2) <- c("Abundance", "Plate", "tax")
  longdf <- rbind(seltax1, seltax2)
  
  return(longdf)
  
}

# Create long data frame for Phylum
longdf_phy <- long_df(seltax_phy_cond)

# Create long data frame for Genus
longdf_ord <- long_df(seltax_ord_cond)

# Create long data frame for Family
longdf_fam <- long_df(seltax_fam_cond)

# Create long data frame for Genus
longdf_gen <- long_df(seltax_gen_cond)

# Create long data frame for Species
longdf_sp <- long_df(seltax_sp_cond)


# Function for creating bubble plots 

bubble_plot <- function(longdf){
  ggplot(longdf,aes(Plate,tax)) +
    geom_point(aes(size=Abundance, fill=tax, color=tax),shape=21) +
    facet_wrap(~Plate, nrow = 1, scales = "free_x", strip.position = "bottom") +
    theme(text = element_text(size=15),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_blank())+
    ylab("") +
    xlab(NULL) +
    guides(fill = "none", color = "none") +
    scale_size(name = "Relative\nabundance")
}

# Bubble plots

bubble_plot(longdf_phy)

bubble_plot(longdf_ord)

# Select specific Orders to plot

longdf_ord <- longdf_ord[longdf_ord$tax %in% ord,]
bubble_plot(longdf_ord) 

# Select specific Genera to plot

longdf_gen <- longdf_gen[longdf_gen$tax %in% genera,]
bubble_plot(longdf_gen)

# Select specific Species to plot
longdf_sp <- longdf_sp[longdf_sp$tax %in% sp,]
bubble_plot(longdf_sp)

# Plot only results from LEfSe

lefse <- c("Lactobacillales", "Ligilactobacillus hayakitensis", "Ligilactobacillus",
           "Limosilactobacillus fermentum", "Merdibacter massiliensis", "Merdibacter",
           "Butyricimonas", "Marifinilaceae", "Lawsonibacter sp900066645",
           "Peptostreptococcales", "Bacteroides thetaiotaomicron")

# Create long data frames
longdf_ord <- long_df(seltax_ord_cond)
longdf_fam <- long_df(seltax_fam_cond)
longdf_gen <- long_df(seltax_gen_cond)
longdf_sp <- long_df(seltax_sp_cond)

longdf_lefse <- rbind(longdf_ord, longdf_fam, longdf_gen, longdf_sp)
longdf_lefse <- longdf_lefse[longdf_lefse$tax %in% lefse,]
longdf_lefse$tax <- factor(longdf_lefse$tax, levels=rev(lefse))

bubble_plot(longdf_lefse)

#### End ####

#### Violin plots of specific bacteria ####

#Abundance value transformation function

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Class",
                          Color = "Phylum",
                          n = NULL){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "Plate", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
    scale_y_log10()
}


violin = subset_taxa(physeq_re, Order == "Lactobacillales")
violin <- subset_taxa(violin, !is.na(Genus))
plot_abundance(violin, Facet = "Genus", Color = "Plate") 

# Plot top 15 most abundant Species
sum = tapply(taxa_sums(ps.Plate_per), tax_table(physeq_per)[, "Species"], sum, na.rm=TRUE)
top10 = names(sort(sum, TRUE))[1:15]
violin = subset_taxa(physeq_per, Species %in% top10)
BP4 <- plot_abundance(violin, Facet = "Species", Color = "Plate", n = 5) 

#### End ####

#### Bar plots paper ####

library(cowplot)

# Figure 3

as_ggplot(grid.arrange(arrangeGrob(BP1 + theme(legend.position = "none",
                                     axis.title.y = element_text(size = 25)),
                         get_legend(BP1 + theme(text = element_text(size = 30))), nrow=1), 
             arrangeGrob(BP2), 
             nrow=2)) + draw_plot_label(label = c("A", "B"), size = 20,
                                     x = c(0, 0), y = c(1, 0.5)) 

#Supplementary Figure 5


as_ggplot(grid.arrange(arrangeGrob(BP3 + theme(text =element_text(size = 25),
                                               axis.title.y = element_text(size = 20),
                                               axis.text.x.bottom = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               legend.text = element_text(size = 16))),
                       arrangeGrob(BP4 + theme(text = element_text(size = 25),
                                               axis.text.y = element_text(size = 15),
                                               strip.text.x = element_text(size = 15))),
                       heights = c(1,2.5))) + draw_plot_label(label = c("A", "B"), size = 20,
                                                    x = c(0, 0), y = c(1, 0.72)) 

#### End ####