
# Alpha diversity #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
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

# Wilcoxon-Test 

wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  W_test <- wilcox.test(rich[, n] ~ rich[, "Plate"], data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(W_test$method, paste("data:", n, "and Plate"), W_test$statistic, W_test$p.value, r)
  wt[[paste0(n)]] <- tab
  
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variables", "Statistic", "p-value", "Effect size")

wt$p.adj <- p.adjust(wt$`p-value`, method = "BH")

write.xlsx(wt, "Results/2.Alpha_stats/Wilcoxon_test.xlsx", rowNames = T)

# Paired t-test

t <- list()
for(n in c("Observed", "Chao1")){
  t_test <- t.test(rich[rich$Plate == "HCB", n], rich[rich$Plate == "MPYG", n], paired = T)
  tab <- c(t_test$method, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
  t[[n]] <- tab
}

t <- as.data.frame(do.call(rbind, t)) 
colnames(t) <- c("Test", "Statistic", "df", "mean difference", "p-value")

t$p.adj <- p.adjust(t$`p-value`, method = "BH")

write.xlsx(t, "Results/2.Alpha_stats/Paired_t_test.xlsx", rowNames = T)

# Adjust p-values

all <- rbind(wt[,c("Test", "Statistic", "p-value")], t[,c("Test", "Statistic", "p-value")])
all$p.adj <- p.adjust(all$`p-value`, method = "BH")

write.xlsx(t, "Results/2.Alpha_stats/Alpha_tests.xlsx", rowNames = T)

#### End ####

#### Alpha-diversity (Richness and diversity estimates) ####

# Customize plot for figure 1

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
P4 <- grid.arrange(nrow = 1,
             P2 + ggtitle(NULL) + theme(legend.position = "none"),
             P3 + ggtitle(NULL) + theme(legend.position = "none") + ylab(NULL),
             get_legend(P1.1), widths = c(4.1,4,1),
           #  top = textGrob("Alpha Diversity",gp=gpar(fontsize=25)),
             bottom = textGrob("Plate",gp=gpar(fontsize=20)))

ggsave(plot = P4,"Results/3.Alpha_plots/Figure_1.png", width = 11, height = 6.4, dpi = 300)
ggsave(plot = P4,"Results/3.Alpha_plots/Figure_1.svg", width = 11, height = 6.4, dpi = 300)

#### End ####
