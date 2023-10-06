#### Average SNP depth by sample ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/05_Analyses/sambamba-depth/")


#### LOAD LIBRARIES ####

library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)
library(ggthemes)
devtools::install_github("thomasp85/ggforce")
library(ggforce)

# 1. No missing data (maxmiss100), ld001 ----
## 1.1. Load data -----
maxmiss100_ld001 <- read.csv("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001_sambamba_depth.txt", header=TRUE, sep = "\t")
maxmiss100_ld001$NumMites <- as.factor(maxmiss100_ld001$NumMites)

# double check everything
str(maxmiss100_ld001)

mean(maxmiss100_ld001$meanCoverage)
median(maxmiss100_ld001$meanCoverage)
sd(maxmiss100_ld001$meanCoverage)
sd(maxmiss100_ld001$meanCoverage)/(sqrt(length(maxmiss100_ld001$meanCoverage)))

facet_label_num_mites <- as_labeller(
    c(`1` = "1 mite", `5` = "5 mites",`20` = "20 mites"))

levels(factor(maxmiss100_ld001$sampleName_short))

ggplot(maxmiss100_ld001, aes(x = factor(sampleName_short, level=c("PROW_953_R1","PROW_954_R1","PROW_981_R1","PROW_984_R1","PROW_953_R5","PROW_954_R5","PROW_981_R5","PROW_984_R5","PROW_953_R20","PROW_954_R20","PROW_981_R20","PROW_984_R20")), y = meanCoverage, fill = sample_ID)) +
    geom_violin() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "SNP Coverage (Number of Reads)") +
    xlab(label = "Sample") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    labs(fill = "Bird ID")

# save as landscape 6x8


maxmiss100_ld001_box<-ggplot(maxmiss100_ld001, aes(x = factor(sampleName_short, level=c("PROW_953_R1","PROW_954_R1","PROW_981_R1","PROW_984_R1","PROW_953_R5","PROW_954_R5","PROW_981_R5","PROW_984_R5","PROW_953_R20","PROW_954_R20","PROW_981_R20","PROW_984_R20")), y = meanCoverage, fill = sample_ID)) +
    geom_boxplot() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "SNP Coverage (Number of Reads)") +
    xlab(label = "Sample") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    labs(fill = "Bird ID")

# save as landscape 6x8


# add dotted lines at sample breaks
maxmiss100_ld001_box_vline_cb <- maxmiss100_ld001_box +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dotted") +
    scale_fill_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c"))



# viridis
maxmiss100_ld001_box_vline_viridis <- maxmiss100_ld001_box +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dotted") +
    scale_fill_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    theme(axis.title.x = element_blank())




ggplot(maxmiss100_ld001, aes(x = sampleName_short, y = meanCoverage, fill = sample_ID)) +
    geom_boxplot() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "Mean Read Coverage") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20))




ggplot(maxmiss100_ld001, aes(x = sampleName_short, y = meanCoverage, fill = sample_ID)) +
    geom_violin() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "Mean Read Coverage") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20))


























# 2. No missing data (maxmiss100), ld99 ----
## 2.1. Load data -----
maxmiss100_ld99 <- read.csv("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99_sambamba_depth.txt", header=TRUE, sep = "\t")
maxmiss100_ld99$NumMites <- as.factor(maxmiss100_ld99$NumMites)

# double check everything
str(maxmiss100_ld99)

mean(maxmiss100_ld99$meanCoverage)
median(maxmiss100_ld99$meanCoverage)
sd(maxmiss100_ld99$meanCoverage)
sd(maxmiss100_ld99$meanCoverage)/(sqrt(length(maxmiss100_ld99$meanCoverage)))

facet_label_num_mites <- as_labeller(
    c(`1` = "1 mite", `5` = "5 mites",`20` = "20 mites"))



ggplot(maxmiss100_ld99, aes(x = factor(sampleName_short, level=c("PROW_953_R1","PROW_954_R1","PROW_981_R1","PROW_984_R1","PROW_953_R5","PROW_954_R5","PROW_981_R5","PROW_984_R5","PROW_953_R20","PROW_954_R20","PROW_981_R20","PROW_984_R20")), y = meanCoverage, fill = sample_ID)) +
    geom_violin() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "SNP Coverage (Number of Reads)") +
    xlab(label = "Sample") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    labs(fill = "Bird ID")

# save as landscape 6x8




maxmiss100_ld99_box<-ggplot(maxmiss100_ld99, aes(x = factor(sampleName_short, level=c("PROW_953_R1","PROW_954_R1","PROW_981_R1","PROW_984_R1","PROW_953_R5","PROW_954_R5","PROW_981_R5","PROW_984_R5","PROW_953_R20","PROW_954_R20","PROW_981_R20","PROW_984_R20")), y = meanCoverage, fill = sample_ID)) +
    geom_boxplot() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "SNP Coverage (Number of Reads)") +
    xlab(label = "Sample") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold")) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    labs(fill = "Bird ID")

# save as landscape 6x8

maxmiss100_ld99_box_vline_cb <- maxmiss100_ld99_box +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dotted") +
    scale_fill_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c"))


# viridis
maxmiss100_ld99_box_vline_viridis <- maxmiss100_ld99_box +
    geom_vline(xintercept = c(4.5, 8.5), linetype = "dotted") +
    scale_fill_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    theme(axis.title.x = element_blank())



ggplot(maxmiss100_ld99, aes(x = sampleName_short, y = meanCoverage, fill = sample_ID)) +
    geom_boxplot() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    facet_zoom(ylim = c(0, 100)) +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "Mean Read Coverage") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20))




ggplot(maxmiss100_ld99, aes(x = sampleName_short, y = meanCoverage, fill = sample_ID)) +
    geom_violin() +
    facet_grid(.~ NumMites, labeller = facet_label_num_mites,scales="free") +
    scale_fill_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    ylab(label = "Mean Read Coverage") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20))




















########## 3. all together now! ######

maxmiss100_ld001_ld99_box<-
    ggarrange(
        maxmiss100_ld001_box, maxmiss100_ld99_box,
        labels = c("A", "B"),
        ncol = 1, nrow = 2, common.legend = TRUE)

maxmiss100_ld001_ld99_box

# 12x8 portrait




maxmiss100_ld001_ld99_box_vline<-
    ggarrange(
        maxmiss100_ld001_box_vline, maxmiss100_ld99_box_vline,
        labels = c("A", "B"),
        ncol = 1, nrow = 2, common.legend = TRUE)

maxmiss100_ld001_ld99_box_vline

# 12x8 portrait



maxmiss100_ld001_ld99_box_vline_cb<-
    ggarrange(
        maxmiss100_ld001_box_vline_cb, maxmiss100_ld99_box_vline_cb,
        labels = c("A", "B"),
        ncol = 1, nrow = 2, common.legend = TRUE)

maxmiss100_ld001_ld99_box_vline_cb

# 12x8 portrait





maxmiss100_ld001_ld99_box_vline_viridis<-
    ggarrange(
        maxmiss100_ld001_box_vline_viridis, maxmiss100_ld99_box_vline_viridis,
        labels = c("A", "B"),
        font.label = list(size = 24),
        hjust = -0.1,
        ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

maxmiss100_ld001_ld99_box_vline_viridis

# 12x8 portrait
