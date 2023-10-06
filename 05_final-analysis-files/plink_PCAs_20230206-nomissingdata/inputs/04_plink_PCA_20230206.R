#### plink PCAs of different maf, dp, and ld thresholds using no missing data ####
#### 6 Feb 2023
# updated 20230921 with viridis colors


#### LOAD LIBRARIES ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/04_SNP_20221013/plink_PCAs_20230206/inputs")


# following this tutorial: https://speciationgenomics.github.io/pca/

library(tidyverse)
library(ggplot2)
library(ggpubr)




#################################################################################
########################### ld = 0.01 ##########################################
#################################################################################



#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld001 ####
prow_maf05_dp10_ld001_pca_nomiss <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld001.eigenvec", col_names = FALSE)
prow_maf05_dp10_ld001_eigenval_nomiss <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp10_ld001_pca_nomiss <- prow_maf05_dp10_ld001_pca_nomiss[,-1]
# set names
names(prow_maf05_dp10_ld001_pca_nomiss)[1] <- "ind"
names(prow_maf05_dp10_ld001_pca_nomiss)[2:ncol(prow_maf05_dp10_ld001_pca_nomiss)] <- paste0("PC", 1:(ncol(prow_maf05_dp10_ld001_pca_nomiss)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp10_ld001_bird_nomiss <- rep(NA, length(prow_maf05_dp10_ld001_pca_nomiss$ind))
prow_maf05_dp10_ld001_bird_nomiss[grep("PROW953", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "PROW953"
prow_maf05_dp10_ld001_bird_nomiss[grep("PROW954", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "PROW954"
prow_maf05_dp10_ld001_bird_nomiss[grep("PROW981", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "PROW981"
prow_maf05_dp10_ld001_bird_nomiss[grep("PROW984", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "PROW984"

# number of mites
prow_maf05_dp10_ld001_nummites_nomiss <- rep(NA, length(prow_maf05_dp10_ld001_pca_nomiss$ind))
prow_maf05_dp10_ld001_nummites_nomiss[grep("R1", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "R1"
prow_maf05_dp10_ld001_nummites_nomiss[grep("R5", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "R5"
prow_maf05_dp10_ld001_nummites_nomiss[grep("R20", prow_maf05_dp10_ld001_pca_nomiss$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp10_ld001_pca_nomiss <- as_tibble(data.frame(prow_maf05_dp10_ld001_pca_nomiss, prow_maf05_dp10_ld001_bird_nomiss, prow_maf05_dp10_ld001_nummites_nomiss))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp10_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp10_ld001_eigenval_nomiss/sum(prow_maf05_dp10_ld001_eigenval_nomiss)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp10_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp10_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp10_ld001_pca_nomiss_plot<-
    ggplot(prow_maf05_dp10_ld001_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld001_bird_nomiss)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_ld001_nummites_nomiss)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_ld001_pca_nomiss_plot



prow_maf05_dp10_ld001_pca_nomiss_plot_cb<-
    ggplot(prow_maf05_dp10_ld001_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld001_bird_nomiss)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_ld001_nummites_nomiss)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_ld001_pca_nomiss_plot_cb

prow_maf05_dp10_ld001_pca_nomiss_plot_viridis<-
    ggplot(prow_maf05_dp10_ld001_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld001_bird_nomiss)) +
    geom_point(size = 5, alpha = 0.8, aes(shape=prow_maf05_dp10_ld001_nummites_nomiss)) +
    scale_color_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld001_pve$pve[2], 3), "%)")) +
    theme(aspect.ratio = 1)

prow_maf05_dp10_ld001_pca_nomiss_plot_viridis
























#################################################################################
########################### ld = 0.99 ##########################################
#################################################################################



#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld99 ####
prow_maf05_dp10_ld99_pca_nomiss <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld99.eigenvec", col_names = FALSE)
prow_maf05_dp10_ld99_eigenval_nomiss <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp10_ld99_pca_nomiss <- prow_maf05_dp10_ld99_pca_nomiss[,-1]
# set names
names(prow_maf05_dp10_ld99_pca_nomiss)[1] <- "ind"
names(prow_maf05_dp10_ld99_pca_nomiss)[2:ncol(prow_maf05_dp10_ld99_pca_nomiss)] <- paste0("PC", 1:(ncol(prow_maf05_dp10_ld99_pca_nomiss)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp10_ld99_bird_nomiss <- rep(NA, length(prow_maf05_dp10_ld99_pca_nomiss$ind))
prow_maf05_dp10_ld99_bird_nomiss[grep("PROW953", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "PROW953"
prow_maf05_dp10_ld99_bird_nomiss[grep("PROW954", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "PROW954"
prow_maf05_dp10_ld99_bird_nomiss[grep("PROW981", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "PROW981"
prow_maf05_dp10_ld99_bird_nomiss[grep("PROW984", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "PROW984"

# number of mites
prow_maf05_dp10_ld99_nummites_nomiss <- rep(NA, length(prow_maf05_dp10_ld99_pca_nomiss$ind))
prow_maf05_dp10_ld99_nummites_nomiss[grep("R1", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "R1"
prow_maf05_dp10_ld99_nummites_nomiss[grep("R5", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "R5"
prow_maf05_dp10_ld99_nummites_nomiss[grep("R20", prow_maf05_dp10_ld99_pca_nomiss$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp10_ld99_pca_nomiss <- as_tibble(data.frame(prow_maf05_dp10_ld99_pca_nomiss, prow_maf05_dp10_ld99_bird_nomiss, prow_maf05_dp10_ld99_nummites_nomiss))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp10_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp10_ld99_eigenval_nomiss/sum(prow_maf05_dp10_ld99_eigenval_nomiss)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp10_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp10_ld99_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp10_ld99_pca_nomiss_plot<-
    ggplot(prow_maf05_dp10_ld99_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld99_bird_nomiss)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_ld99_nummites_nomiss)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_ld99_pca_nomiss_plot





prow_maf05_dp10_ld99_pca_nomiss_plot_cb<-
    ggplot(prow_maf05_dp10_ld99_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld99_bird_nomiss)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_ld99_nummites_nomiss)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_ld99_pca_nomiss_plot_cb






prow_maf05_dp10_ld99_pca_nomiss_plot_viridis<-
    ggplot(prow_maf05_dp10_ld99_pca_nomiss, aes(PC1, PC2, col = prow_maf05_dp10_ld99_bird_nomiss)) +
    geom_point(size = 5, alpha = 0.8, aes(shape=prow_maf05_dp10_ld99_nummites_nomiss)) +
    scale_color_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_ld99_pve$pve[2], 3), "%)")) +
    theme(aspect.ratio = 1)

prow_maf05_dp10_ld99_pca_nomiss_plot_viridis






prow_maf05_dp10_nomiss<-
    ggarrange(
        prow_maf05_dp10_ld001_pca_nomiss_plot, prow_maf05_dp10_ld99_pca_nomiss_plot,
        labels = c("A", "B"),
        ncol = 2, nrow = 1, common.legend = TRUE)

prow_maf05_dp10_nomiss

# 4x9 landscape




prow_maf05_dp10_nomiss_cb<-
    ggarrange(
        prow_maf05_dp10_ld001_pca_nomiss_plot_cb, prow_maf05_dp10_ld99_pca_nomiss_plot_cb,
        labels = c("A", "B"),
        ncol = 2, nrow = 1, common.legend = TRUE)

prow_maf05_dp10_nomiss_cb

# 4x9 landscape



prow_maf05_dp10_nomiss_viridis<-
    ggarrange(
        prow_maf05_dp10_ld001_pca_nomiss_plot_viridis, prow_maf05_dp10_ld99_pca_nomiss_plot_viridis,
        labels = c("A", "B"),
        font.label = list(size = 24),
        ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

prow_maf05_dp10_nomiss_viridis

# 4x9 landscape