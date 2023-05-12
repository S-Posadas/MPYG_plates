
# Alpha diversity #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(grid)
library(openxlsx)

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

load("Results/RData/1.physeq.original.RData")

# Create directories for results
dir.create("Results/2.Alpha_stats")
dir.create("Results/3.Alpha_plots")

#### End ####

#### Alpha-diversity (Richness and diversity estimates) ####
# Should be done in non filtered physeq
sample_variables(physeq)

# Visualize alpha-diversity on unfiltered phyloseq object

P1 = plot_richness(physeq, x="Plate", color = "Plate", title = "Alpha Diversity", measures=c("Observed", "Chao1", "Shannon", "InvSimpson"))

P1.1 = P1+
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  theme(strip.background = element_blank(),
        text = element_text(size=20),
        axis.text.y = element_text(size = 13),
        axis.text.x.bottom = element_text(angle = -30)) 

P1.1

# Violin plot

P1 + geom_violin(alpha = 0.5) +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90)) + geom_point(position = position_dodge(width = 0.75))

#### End ####

#### Richness statistics####
# Has to be counts not relative abundances

rich = estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon","InvSimpson"))
rownames(rich) <- gsub("X", "", rownames(rich)) %>% gsub("\\.", "-", .)

# Merge with sample_data for easy plotting

rich <- merge(rich, sample_data(physeq), by = 0)

# Save table with alpha diversity

write.xlsx(rich, "Results/2.Alpha_stats/richness.xlsx")

#Check if the indices are normally distributed with "Normality_check.R" script and come back to this one

#Shannon and Invsimpson are not normally distributed and paired -> paired Wilcoxon 
#Observed and Chao1 are normally distributed and paired -> paired t test
#However, we performed all tests in all indexes

# Mann-Whitney-U-Test 

mwut <- list()

for (n in c("Observed", "Chao1", "Shannon", "InvSimpson")) {
  U_test <- wilcox.test(rich[, n] ~ rich[, "Plate"], data = rich, paired = F, exact = T)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(U_test$method, paste("data:", n, "and Plate"), U_test$statistic, U_test$p.value, r)
  mwut[[paste0(n)]] <- tab
  
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")

mwut$p.adj <- p.adjust(mwut$`p-value`, method = "BH")

write.xlsx(mwut, "Results/2.Alpha_stats/Mann_Whitney_u_test.xlsx", rowNames = T)

# Wilcoxon-Test 

wt <- list()

for (n in c("Observed", "Chao1", "Shannon", "InvSimpson")) {
  W_test <- wilcox.test(rich[, n] ~ rich[, "Plate"], data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(W_test$method, paste("data:", n, "and Plate"), W_test$statistic, W_test$p.value, r)
  wt[[paste0(n)]] <- tab
  
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")

wt$p.adj <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj_notnorm <- c(NA, NA, p.adjust(wt$`p-value`[3:4], method = "BH"))

write.xlsx(wt, "Results/2.Alpha_stats/Wilcoxon_test.xlsx", rowNames = T)

# One-Way Anova

#https://scienceparkstudygroup.github.io/microbiome-lesson/aio/index.html

aov <- list()

for(n in c("Observed", "Chao1", "Shannon", "InvSimpson")){
  aov_test <- aov(rich[,n] ~ rich[, "Plate"], data = rich)
  aov_sum <- summary(aov_test)
  aov_test$df.residual
  tab <- c(n, aov_sum[[1]]$`Sum Sq`, aov_sum[[1]]$`Mean Sq`, aov_sum[[1]]$`F value`, aov_sum[[1]]$`Pr(>F)`)
  aov[[n]] <- tab
}
aov <- as.data.frame(do.call(rbind, aov)) ; aov <- aov[ , colSums(is.na(aov)) == 0]
colnames(aov) <- c("Variable", "Sum Sq [i]", "Sum Sq Residuals", "Mean Sq [i]", "Mean Sq Residuals", "F value", "Pr(>F)")

aov$p.adj <- p.adjust(aov$`Pr(>F)`, method = "BH")

write.xlsx(aov, "Results/2.Alpha_stats/One_way_Anova.xlsx")

# Paired t-test

t <- list()
for(n in c("Observed", "Chao1", "Shannon", "InvSimpson")){
  t_test <- t.test(rich[rich$Plate == "HCB", n], rich[rich$Plate == "MPYG", n], paired = T)
  tab <- c(t_test$method, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
  t[[n]] <- tab
}

t <- as.data.frame(do.call(rbind, t)) 
colnames(t) <- c("Method", colnames(t)[2:4], "p-value")

t$p.adj <- p.adjust(t$`p-value`, method = "BH")
t$p.adj_norm <- c(p.adjust(t$`p-value`[1:2], method = "BH"), NA, NA)

write.xlsx(t, "Results/2.Alpha_stats/Paired_t_test.xlsx", rowNames = T)

#### End ####

#### Alpha-diversity (Richness and diversity estimates) ####

# Create plot with the right statistics

# Visualize alpha-diversity on unfiltered phyloseq object

P2 = plot_richness(physeq, x="Plate", color = "Plate", title = "Alpha Diversity", measures=c("Observed", "Chao1"))

P2 <- P2+
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  xlab(NULL) +
  theme(strip.background = element_blank(),
        text = element_text(size=20),
        axis.text.y = element_text(size = 13),
        axis.text.x.bottom = element_text(angle = -30)) +
  stat_compare_means(method = "t.test", paired = T, label.y = 1.8) 


P3 = plot_richness(physeq, x="Plate", color = "Plate", title = "Alpha Diversity", measures=c("Shannon", "InvSimpson"))

P3 <- P3+
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  xlab(NULL) + 
  theme(strip.background = element_blank(),
        text = element_text(size=20),
        axis.text.y = element_text(size = 13),
        axis.text.x.bottom = element_text(angle = -30)) +
  stat_compare_means(method = "wilcox.test", paired = T, label.y = 1.8) 


# Combine all plots
grid.arrange(nrow = 1,
             P2 + ggtitle(NULL) + theme(legend.position = "none"),
             P3 + ggtitle(NULL) + theme(legend.position = "none") + ylab(NULL),
             get_legend(P1), widths = c(4,4,1),
             top = textGrob("Alpha Diversity",gp=gpar(fontsize=25)),
             bottom = textGrob("Plate",gp=gpar(fontsize=20)))

#### End ####