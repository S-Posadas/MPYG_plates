#Save phyloseq object as input file for lefse

#### Package setup ####

library(phyloseq)
library(phyloseqCompanion)
library(dplyr)

#### End ####

#### Load RData with phyloseq object ####

#load("Results/RData/1.physeq.original.RData")          #unfiltered
load("Results/RData/2.phyloseq.filtered.RData")         #filtered 

#### End ####

#### Save as input for LEfSe #### 

#Select phyloseq object to save
#phy2lefse <- physeq               #total
phy2lefse <- physeq_re            #relative

#Changing Species names
tax <- as.data.frame(tax_table(phy2lefse))
tax[,"Species"] <- tax[,"Species"] %>% gsub("\\(RS_.*", "", .) %>% gsub("\\(GB_.*", "", .) %>% gsub("_RS_GCF_.*", "", .) 

# Obtaining right percentage in taxa with _letter. Does not make a difference
#for (n in rank_names(phy2lefse)) {
#  rows <- grep("_", tax[,n])
#  subs <- tax[rows,n]
#  root <- unique(gsub("_.*", "", subs))
#  for(i in root){
#    tax[,n] <- sub(i, paste0(i, "."), tax[,n])
#  }
#  tax[rows,n] <- subs 
#}

tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(phy2lefse) <- tax; rm(tax)


phyloseq2lefse(phy2lefse, c("Plate", "Sample_ID"), file.name = "Results/lefse/physeq_re_lfse.txt",
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)

#### End ####