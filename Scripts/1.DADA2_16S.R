####Introduction####

# This pipeline was used to evaluate the impact of growth media on the culture
# of anaerobic bacteria from human stool samples.
# First, we recovered the ASVs following the DADA2 pipeline with some modifications

#### End ####

#### Package setup ####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
#remotes::install_github("kstagaman/phyloseqCompanion")

library(dada2); packageVersion("dada2")
library(vegan)
library(phyloseq)
library(msa)
library(phangorn)
library(ggplot2)
library(gridExtra)

#### End ####

#### Set Directory ####
setwd("path/to/working/directory") 
path <- "Data/1-Fastq" # Change to the directory containing the fastq files
list.files(path)
theme_set(theme_bw())

# Create directory to save results of the analysis and subdirectory for quality control

dir.create("Results/1.QC", recursive = T)

#### End ####

#### Quality Control and Filtering####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

plotQualityProfile(fnFs)

plotQualityProfile(fnRs)

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim 

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE
# and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read,
# which is a better filter than simply averaging quality scores.
# Watch out with trunclen, reads have to overlap at the end, you have to try out.
# maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# Check quality after filtering

plotQualityProfile(filtFs)

plotQualityProfile(filtRs)

#### End ####

####Learn the Error Rates####

errF <- learnErrors(filtFs, multithread=FALSE)

errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#### End ####

#### Dada and Merge####
# Apply the core sample inference algorithm to the the filtered and trimmed sequence data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

# Inspecting the returned dada-class object

dadaFs[[1]]

#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) # min overlap is 12 as default, but can be adjusted

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Most of your reads should successfully merge. If that is not the case upstream parameters may need
#to be revisited: Did you trim away the overlap between your reads?

#Construct sequence table

seqtab <- makeSequenceTable(mergers)  

dim(seqtab)

#If you wish to filter very short or long merged reads:-
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]
#dim(seqtab2)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)

# Percentage of not chimeric reads

sum(seqtab.nochim)/sum(seqtab)

#### End ####

#### Track Pipeline ####
#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

# Save the track

write.csv(track, file = "Results/1.QC/Pipeline_Track_final.csv")
write.csv(sample.names, file = "Results/1.QC/sample.names.csv")

#### End ####

#### Assign taxonomy for Bacteria with GTDB####

ref_fasta = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz"

#tryRC = TRUE -> reverse-complement orientation
set.seed(123456789)
taxa_GTDB <- assignTaxonomy(seqtab.nochim, ref_fasta,tryRC = TRUE )
unname(taxa_GTDB)

taxa.print_GTDB  <- taxa_GTDB # Removing sequence rownames for display only
rownames(taxa.print_GTDB) <- NULL
head(taxa.print_GTDB)

#### End ####

#### Make Count and Taxa Tables ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab <- c(rbind(asv_headers, asv_seqs))
write(fasta_tab, "Results/fasta_tab.fa")

# count table:

count_tab <- t(seqtab.nochim)
row.names(count_tab) <- sub(">", "", asv_headers)
write.table(count_tab, "Results/count_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:

tax_tab <- taxa_GTDB
row.names(tax_tab) <- sub(">", "", asv_headers)
write.table(tax_tab, "Results/tax_tab_GTDB.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

#### Add samples metadata and match to count table####

sample_info_tab <- read.table("Data/Sample_Info.txt",header=T, 
                              row.names = 1, check.names=F)
# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

# Examine consistancy in order between count_tab colnames and coldata rownames 
# (They have to be in the same order or Deseq2 won't work out)

all(rownames(sample_info_tab) %in% colnames(count_tab))
all(rownames(sample_info_tab) == colnames(count_tab)) # Mock community is not in sample_info_tab


all(rownames(tax_tab) %in% rownames(count_tab))
all(rownames(tax_tab) == rownames(count_tab))

gplots::venn(list(taxonomy=rownames(tax_tab), featuretable=rownames(count_tab)))

#### End ####

#### Rarefaction curves ####

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

raremax <- min(rowSums(t(count_tab)))
raremax
rarefy(t(count_tab), raremax)

# Option 1
rarecurve(t(count_tab), step = 100,  col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)
abline(v=(raremax))

# Option 2
rarecurve(t(count_tab), step = 100, sample = raremax, col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)

#### End ####

#### Phylogenetic tree ####

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

# The phangorn package is then used to construct a phylogenetic tree. 
# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum
# likelihood tree using the neighbor-joining tree as a starting point.

phang.align <- as.phyDat(mult, type="dna", names=getSequence(seqtab.nochim))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# Negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

#### End ####

#### Making our phyloseq object ####
# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

physeq <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                   sample_data(sample_info_tab), 
                   tax_table(tax_tab))

# Adding ASV Fasta sequences and Phylogenetic tree to phyloseq object

dna <- Biostrings::DNAStringSet(asv_seqs)  # Making ASV Fasta sequences 
names(dna) <- taxa_names(physeq)

ph_tree = phy_tree(fitGTR$tree)            # Making Phylogenetic tree
taxa_names(ph_tree) = taxa_names(dna)

physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

taxa_names(physeq_dna_tree) <- paste0("ASV", seq(ntaxa(physeq_dna_tree)))

physeq
physeq_dna_tree
physeq = physeq_dna_tree

#### End ####

## Phyloseq Object Filtering ## First filtering step for low count samples and NA phyla ##

#### Remove samples with less than 100 total reads ####
sample_sums(physeq) #Nr of reads per Sample

#77-H  77-M  78-H  78-M  80-H  80-M  81-H  81-M  82-H  82-M  83-H  83-M  84-H  84-M  85-H  85-M 
#11860 20420  2807 40444 51321 32495 25008 13653 39722 15412 12287  9456 21606  1723 10909 31215 

#No samples below 100 reads, thus no filtering necessary

#physeq_above100reads = prune_samples(sample_sums(physeq)>=100, physeq)

#physeq = physeq_above100reads

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

dir.create("Results/RData")

save.image("Results/RData/1.physeq.original.RData")

#### End ####

## Filter low Abundance Taxa and count table normalization ##

#### Remove NA Phyla####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

physeq = physeq_oNA

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

#### End ####

####Define prevalence of each taxa (in how many samples did each taxa appear at least once)####

prev0 = apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))

#save ASV Prevalence and Abundance table before filtering

write.table(prevdf, "Results/asv_prevdf.tsv", sep="\t", quote=F, col.names=NA)

#Plot Taxa prevalence v. total counts. Each point is a different taxa. 

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#### End ####

#### Remove taxa not seen more than 3 times in at least 5% of the samples #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperphyla = 3
Samplepercentage = 0.05

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)
physeq = physeq_filtered

#### End ####

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq, standf)

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

#### End ####

#### Exploratory plots after filtering and normalization ####
# Check individual phylum Abundance
# Abundance value transformation function

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


#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image("Results/RData/2.phyloseq.filtered.RData")

#### End ####