library("phyloseq")
library("dplyr")
library("ggplot2")
library("qqplotr")
library("gridExtra")
library("DescTools")


#https://rstudio-pubs-static.s3.amazonaws.com/713954_d40760746cd3402fb0d1012ff67e5ab9.html

#rich.test <- read.xlsx("Results/richness.xlsx", rowNames = T)
rich.test <- rich
rownames(rich.test) <- gsub("X", "", rownames(rich.test)) %>% gsub("\\.", "-", .)
seq.depth <- as.data.frame(sample_sums(physeq)) ; colnames(seq.depth) = c("seq.depth")
rich.test <- as.data.frame(merge(rich.test, seq.depth, by = 0))


# Normality check for sequencing depth ------------------------------------

hist(rich.test$seq.depth)

ggplot(data = rich.test, mapping = aes(sample = seq.depth)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

shapiro.test(rich.test$seq.depth)


# Normality check for alpha diversity -------------------------------------
rich.test$Ob <- "Observed"
rich.test$Ch <- "Chao1"
rich.test$Sh <- "Shannon"
rich.test$Si <- "InvSimpson"

#Create histograms 
P_Ob <- ggplot(data = rich.test, mapping = aes(x = Observed)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Ob) 

P_Ch <- ggplot(data = rich.test, mapping = aes(x = Chao1)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Ch)

P_Sh <- ggplot(data = rich.test, mapping = aes(x = Shannon)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Sh) 

P_Si <- ggplot(data = rich.test, mapping = aes(x = InvSimpson)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Si) 

grid.arrange(P_Ob, P_Ch, P_Sh, P_Si, nrow = 1)

#Produce descriptive statistics

alpha <- list()
for(i in c("Observed", "Chao1", "Shannon", "InvSimpson")){
  alpha[[i]] <- cbind("alpha.index" = i,
    "n" = length(rich.test[,i]), 
    "mean" = mean(rich.test[,i], na.rm = TRUE), 
    "sd "= sd(rich.test[,i], na.rm = TRUE),
    "stderr" = sd(rich.test[,i], na.rm = TRUE)/sqrt(length(rich.test[,i])),
    "LCL" = mean(rich.test[,i], na.rm = TRUE) - qt(1 - (0.05 / 2), length(rich.test[,i]) - 1) * sd(rich.test[,i], na.rm = TRUE)/sqrt(length(rich.test[,i])),
    "UCL" = mean(rich.test[,i], na.rm = TRUE) + qt(1 - (0.05 / 2), length(rich.test[,i]) - 1) * sd(rich.test[,i], na.rm = TRUE)/sqrt(length(rich.test[,i])),
    "median" = median(rich.test[,i], na.rm = TRUE),
    "min" = min(rich.test[,i], na.rm = TRUE), 
    "max" = max(rich.test[,i], na.rm = TRUE),
    "IQR" = IQR(rich.test[,i], na.rm = TRUE),
    "LCLmed" = MedianCI(rich.test[,i], na.rm=TRUE)[2],
    "UCLmed" = MedianCI(rich.test[,i], na.rm=TRUE)[3],
    "shapiro.W.Stat" = shapiro.test(rich.test[,i])$statistic,
    "shapiro.p.value" = shapiro.test(rich.test[,i])$p.value)
  
}
alpha <- as.data.frame(do.call(rbind, alpha))
alpha[,c("alpha.index", "shapiro.W.Stat", "shapiro.p.value")]

write.xlsx(alpha, "Results/Alpha_stats/shapiro.xlsx")

#Perform QQ plots
P_Ob <- ggplot(data = rich.test, mapping = aes(sample = Observed)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Ob) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

P_Ch <- ggplot(data = rich.test, mapping = aes(sample = Chao1)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Ch) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

P_Sh <- ggplot(data = rich.test, mapping = aes(sample = Shannon)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Sh) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

P_Si <- ggplot(data = rich.test, mapping = aes(sample = InvSimpson)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Si) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

grid.arrange(P_Ob, P_Ch, P_Sh, P_Si, nrow = 1)


# Normality check for alpha diversity by groups ---------------------------

#Change each variable
#Designate Plate as a categorical factor
rich.test$Plate<-as.factor(gsub("*..-", "", rich.test$Row.names))

# Produce descriptive statistics by group
rich.test %>% select("Plate", "Shannon") %>% group_by(Plate) %>% 
  dplyr::summarise(n = n(), 
            mean = mean(Shannon, na.rm = TRUE), 
            sd = sd(Shannon, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(Shannon, na.rm = TRUE),
            min = min(Shannon, na.rm = TRUE), 
            max = max(Shannon, na.rm = TRUE),
            IQR = IQR(Shannon, na.rm = TRUE),
            LCLmed = MedianCI(Shannon, na.rm=TRUE)[2],
            UCLmed = MedianCI(Shannon, na.rm=TRUE)[3],
            "shapiro.W.Stat" = shapiro.test(rich.test[,i])$statistic,
            "shapiro.p.value" = shapiro.test(rich.test[,i])$p.value)


# Produce Boxplots and visually check for outliers
ggplot(rich.test, aes(x = Plate, y = Shannon, fill = Plate)) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(fill = "light blue") + 
  stat_summary(fun.y=mean, geom="point", shape=10, size=3.5, color="black") + 
  ggtitle("Boxplot of Treatments C and D") + 
  theme_bw() + theme(legend.position="none")

# Perform QQ plots by group
ggplot(data = rich.test, mapping = aes(sample = Shannon, color = Plate, fill = Plate)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_wrap(~ Plate, scales = "free") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()


