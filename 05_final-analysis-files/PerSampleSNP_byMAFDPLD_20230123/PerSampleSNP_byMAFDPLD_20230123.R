#### Number of SNPs versus LD threshold ####
## Using the CORRECT values now, previous one was using [6]nHets+[7]nTransitions columns instead of [5]nNonRefHom+[6]nHets to calculate these total # of SNP numbers - shoot!
# updated 20230921 to jitter the points

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/04_SNP_20221013/PerSampleSNP/")

#### LOAD LIBRARIES ####

library(ggplot2)
library(scales)
library(ggpubr)

data <- read.csv("PerSampleSNP_byMAFDPLD_20230123.csv", header=TRUE, sep = ",")
str(data)

data$NumMites <- as.factor(data$NumMites)
data$MAF <- as.factor(data$MAF)
data$LD <- as.factor(data$LD)

str(data)



#### subset different datasets
# missing included (for comparison) vs. missing excluded (what we use for most analyses)

# missing excluded (I'm calling it 'nmd' for 'no missing data')
data_nmd<-subset(data, MissingExcluded == "Yes")
levels(data_nmd$MissingExcluded) # missingexcluded == "no" still included
data_nmd$MissingExcluded <- factor(data_nmd$MissingExcluded) # factor it again to remove the extra level
levels(data_nmd$MissingExcluded) # gone!

# missing excluded, two different ld levels
data_nmd_maf05dp10ld01<-subset(data_nmd, LD == "0.01")
data_nmd_maf05dp10ld01$LD <- factor(data_nmd_maf05dp10ld01$LD)
str(data_nmd_maf05dp10ld01)

data_nmd_maf05dp10ld99<-subset(data_nmd, LD == "0.99")
data_nmd_maf05dp10ld99$LD <- factor(data_nmd_maf05dp10ld99$LD)
str(data_nmd_maf05dp10ld99)

# missing included
data_md<-subset(data, MissingExcluded == "No")
levels(data_md$MissingExcluded) # missingexcluded == "yes" still included
data_md$MissingExcluded <- factor(data_md$MissingExcluded) # factor it again to remove the extra level
levels(data_md$MissingExcluded) # gone!

# missing included, two different ld levels
data_md_maf05dp10ld01<-subset(data_md, LD == "0.01")
data_md_maf05dp10ld01$LD <- factor(data_md_maf05dp10ld01$LD)
str(data_md_maf05dp10ld01)

data_md_maf05dp10ld99<-subset(data_md, LD == "0.99")
data_md_maf05dp10ld99$LD <- factor(data_md_maf05dp10ld99$LD)
str(data_md_maf05dp10ld99)




#### PLOTS ####

### missing excluded ####
plot_nmd_maf05dp10ld01<- ggplot(data_nmd_maf05dp10ld01, aes(x = Sample, y = NumSNPs)) +
    geom_jitter(aes(col=NumMites), size = 5, alpha = 0.9, width = 0.15) +
    scale_color_manual(values = c("lightcyan3", "steelblue3", "mediumpurple4")) +
    labs(color = "Number of Mites") +
    ylab("Number of SNPs") +
    xlab("Sample") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, vjust = 0.5),
          axis.title.x = element_blank())

plot_nmd_maf05dp10ld01


plot_nmd_maf05dp10ld99<- ggplot(data_nmd_maf05dp10ld99, aes(x = Sample, y = NumSNPs)) +
    geom_jitter(aes(col=NumMites), size = 5, alpha = 0.9, width = 0.15) +
    scale_color_manual(values = c("lightcyan3", "steelblue3", "mediumpurple4")) +
    labs(color = "Number of Mites") +
    ylab("Number of SNPs") +
    xlab("Sample") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, vjust = 0.5),
          axis.title.x = element_blank())

plot_nmd_maf05dp10ld99


plot_nmd_maf05_dp10<-
    ggarrange(
        plot_nmd_maf05dp10ld01, plot_nmd_maf05dp10ld99,
        labels = c("A", "B"),
        font.label = list(size = 24),
        ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

plot_nmd_maf05_dp10

# 4x7 landscape










### missing included ####
plot_md_maf05dp10ld01<- ggplot(data_md_maf05dp10ld01, aes(x = Sample, y = NumSNPs)) +
    geom_jitter(aes(col=NumMites), size = 5, alpha = 0.9, width = 0.15) +
    scale_color_manual(values = c("lightcyan3", "steelblue3", "mediumpurple4")) +
    labs(color = "Number of Mites") +
    ylab("Number of SNPs") +
    xlab("Sample") +
    scale_y_continuous(labels = scales::comma) +
    ggtitle(bquote('LD'~r^2:0.01)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, vjust = 0.5),
          axis.title.x = element_blank())

plot_md_maf05dp10ld01


plot_md_maf05dp10ld99<- ggplot(data_md_maf05dp10ld99, aes(x = Sample, y = NumSNPs)) +
    geom_jitter(aes(col=NumMites), size = 5, alpha = 0.9, width = 0.15) +
    scale_color_manual(values = c("lightcyan3", "steelblue3", "mediumpurple4")) +
    labs(color = "Number of Mites") +
    ylab("Number of SNPs") +
    xlab("Sample") +
    scale_y_continuous(labels = scales::comma) +
    ggtitle(bquote('LD'~r^2:0.99)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 60, vjust = 0.5),
          axis.title.x = element_blank())

plot_md_maf05dp10ld99


plot_md_maf05_dp10<-
    ggarrange(
        plot_md_maf05dp10ld01, plot_md_maf05dp10ld99,
        labels = c("A", "B"),
        font.label = list(size = 24),
        ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

plot_md_maf05_dp10

# 4x7 landscape

