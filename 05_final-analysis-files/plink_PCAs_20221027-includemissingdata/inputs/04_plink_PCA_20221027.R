#### plink PCAs of different maf, dp, and ld thresholds ####
#### 27 Oct 2022


#### LOAD LIBRARIES ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/04_SNP_20221013/plink_PCAs_20221027/inputs")


# following this tutorial: https://speciationgenomics.github.io/pca/

library(tidyverse)
library(ggplot2)
library(ggpubr)





#################################################################################
########################### maf = 0.01 ##########################################
#################################################################################


#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld001 ####
prow_maf01_dp05_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf01_dp05_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp05_maxmiss100_ld001_pca <- prow_maf01_dp05_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf01_dp05_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf01_dp05_maxmiss100_ld001_pca)[2:ncol(prow_maf01_dp05_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp05_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp05_maxmiss100_ld001_bird <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld001_pca$ind))
prow_maf01_dp05_maxmiss100_ld001_bird[grep("PROW953", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf01_dp05_maxmiss100_ld001_bird[grep("PROW954", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf01_dp05_maxmiss100_ld001_bird[grep("PROW981", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf01_dp05_maxmiss100_ld001_bird[grep("PROW984", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp05_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld001_pca$ind))
prow_maf01_dp05_maxmiss100_ld001_nummites[grep("R1", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf01_dp05_maxmiss100_ld001_nummites[grep("R5", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf01_dp05_maxmiss100_ld001_nummites[grep("R20", prow_maf01_dp05_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp05_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf01_dp05_maxmiss100_ld001_pca, prow_maf01_dp05_maxmiss100_ld001_bird, prow_maf01_dp05_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp05_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp05_maxmiss100_ld001_eigenval/sum(prow_maf01_dp05_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp05_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp05_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp05_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf01_dp05_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf01_dp05_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp05_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp05_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp05_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf01_dp05_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld05 ####
prow_maf01_dp05_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf01_dp05_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp05_maxmiss100_ld05_pca <- prow_maf01_dp05_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf01_dp05_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf01_dp05_maxmiss100_ld05_pca)[2:ncol(prow_maf01_dp05_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp05_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp05_maxmiss100_ld05_bird <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld05_pca$ind))
prow_maf01_dp05_maxmiss100_ld05_bird[grep("PROW953", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf01_dp05_maxmiss100_ld05_bird[grep("PROW954", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf01_dp05_maxmiss100_ld05_bird[grep("PROW981", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf01_dp05_maxmiss100_ld05_bird[grep("PROW984", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp05_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld05_pca$ind))
prow_maf01_dp05_maxmiss100_ld05_nummites[grep("R1", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf01_dp05_maxmiss100_ld05_nummites[grep("R5", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf01_dp05_maxmiss100_ld05_nummites[grep("R20", prow_maf01_dp05_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp05_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf01_dp05_maxmiss100_ld05_pca, prow_maf01_dp05_maxmiss100_ld05_bird, prow_maf01_dp05_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp05_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp05_maxmiss100_ld05_eigenval/sum(prow_maf01_dp05_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp05_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp05_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp05_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf01_dp05_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf01_dp05_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp05_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp05_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp05_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf01_dp05_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99 ####
prow_maf01_dp05_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf01_dp05_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp05_maxmiss100_ld99_pca <- prow_maf01_dp05_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf01_dp05_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf01_dp05_maxmiss100_ld99_pca)[2:ncol(prow_maf01_dp05_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp05_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp05_maxmiss100_ld99_bird <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld99_pca$ind))
prow_maf01_dp05_maxmiss100_ld99_bird[grep("PROW953", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf01_dp05_maxmiss100_ld99_bird[grep("PROW954", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf01_dp05_maxmiss100_ld99_bird[grep("PROW981", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf01_dp05_maxmiss100_ld99_bird[grep("PROW984", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp05_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf01_dp05_maxmiss100_ld99_pca$ind))
prow_maf01_dp05_maxmiss100_ld99_nummites[grep("R1", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf01_dp05_maxmiss100_ld99_nummites[grep("R5", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf01_dp05_maxmiss100_ld99_nummites[grep("R20", prow_maf01_dp05_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp05_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf01_dp05_maxmiss100_ld99_pca, prow_maf01_dp05_maxmiss100_ld99_bird, prow_maf01_dp05_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp05_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp05_maxmiss100_ld99_eigenval/sum(prow_maf01_dp05_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp05_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp05_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf01_dp05_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf01_dp05_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf01_dp05_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp05_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp05_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp05_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf01_dp05_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld001 ####
prow_maf01_dp10_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf01_dp10_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp10_maxmiss100_ld001_pca <- prow_maf01_dp10_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf01_dp10_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf01_dp10_maxmiss100_ld001_pca)[2:ncol(prow_maf01_dp10_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp10_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp10_maxmiss100_ld001_bird <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld001_pca$ind))
prow_maf01_dp10_maxmiss100_ld001_bird[grep("PROW953", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf01_dp10_maxmiss100_ld001_bird[grep("PROW954", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf01_dp10_maxmiss100_ld001_bird[grep("PROW981", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf01_dp10_maxmiss100_ld001_bird[grep("PROW984", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp10_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld001_pca$ind))
prow_maf01_dp10_maxmiss100_ld001_nummites[grep("R1", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf01_dp10_maxmiss100_ld001_nummites[grep("R5", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf01_dp10_maxmiss100_ld001_nummites[grep("R20", prow_maf01_dp10_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp10_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf01_dp10_maxmiss100_ld001_pca, prow_maf01_dp10_maxmiss100_ld001_bird, prow_maf01_dp10_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp10_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp10_maxmiss100_ld001_eigenval/sum(prow_maf01_dp10_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp10_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp10_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp10_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf01_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf01_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf01_dp10_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld05 ####
prow_maf01_dp10_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf01_dp10_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp10_maxmiss100_ld05_pca <- prow_maf01_dp10_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf01_dp10_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf01_dp10_maxmiss100_ld05_pca)[2:ncol(prow_maf01_dp10_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp10_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp10_maxmiss100_ld05_bird <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld05_pca$ind))
prow_maf01_dp10_maxmiss100_ld05_bird[grep("PROW953", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf01_dp10_maxmiss100_ld05_bird[grep("PROW954", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf01_dp10_maxmiss100_ld05_bird[grep("PROW981", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf01_dp10_maxmiss100_ld05_bird[grep("PROW984", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp10_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld05_pca$ind))
prow_maf01_dp10_maxmiss100_ld05_nummites[grep("R1", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf01_dp10_maxmiss100_ld05_nummites[grep("R5", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf01_dp10_maxmiss100_ld05_nummites[grep("R20", prow_maf01_dp10_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp10_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf01_dp10_maxmiss100_ld05_pca, prow_maf01_dp10_maxmiss100_ld05_bird, prow_maf01_dp10_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp10_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp10_maxmiss100_ld05_eigenval/sum(prow_maf01_dp10_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp10_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp10_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp10_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf01_dp10_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf01_dp10_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp10_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp10_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp10_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf01_dp10_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99 ####
prow_maf01_dp10_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf01_dp10_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp10_maxmiss100_ld99_pca <- prow_maf01_dp10_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf01_dp10_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf01_dp10_maxmiss100_ld99_pca)[2:ncol(prow_maf01_dp10_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp10_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp10_maxmiss100_ld99_bird <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld99_pca$ind))
prow_maf01_dp10_maxmiss100_ld99_bird[grep("PROW953", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf01_dp10_maxmiss100_ld99_bird[grep("PROW954", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf01_dp10_maxmiss100_ld99_bird[grep("PROW981", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf01_dp10_maxmiss100_ld99_bird[grep("PROW984", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp10_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf01_dp10_maxmiss100_ld99_pca$ind))
prow_maf01_dp10_maxmiss100_ld99_nummites[grep("R1", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf01_dp10_maxmiss100_ld99_nummites[grep("R5", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf01_dp10_maxmiss100_ld99_nummites[grep("R20", prow_maf01_dp10_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp10_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf01_dp10_maxmiss100_ld99_pca, prow_maf01_dp10_maxmiss100_ld99_bird, prow_maf01_dp10_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp10_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp10_maxmiss100_ld99_eigenval/sum(prow_maf01_dp10_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp10_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp10_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf01_dp10_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf01_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf01_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf01_dp10_maxmiss100_ld99_pca_plot



























#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld001 ####
prow_maf01_dp15_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf01_dp15_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp15_maxmiss100_ld001_pca <- prow_maf01_dp15_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf01_dp15_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf01_dp15_maxmiss100_ld001_pca)[2:ncol(prow_maf01_dp15_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp15_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp15_maxmiss100_ld001_bird <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld001_pca$ind))
prow_maf01_dp15_maxmiss100_ld001_bird[grep("PROW953", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf01_dp15_maxmiss100_ld001_bird[grep("PROW954", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf01_dp15_maxmiss100_ld001_bird[grep("PROW981", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf01_dp15_maxmiss100_ld001_bird[grep("PROW984", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp15_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld001_pca$ind))
prow_maf01_dp15_maxmiss100_ld001_nummites[grep("R1", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf01_dp15_maxmiss100_ld001_nummites[grep("R5", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf01_dp15_maxmiss100_ld001_nummites[grep("R20", prow_maf01_dp15_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp15_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf01_dp15_maxmiss100_ld001_pca, prow_maf01_dp15_maxmiss100_ld001_bird, prow_maf01_dp15_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp15_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp15_maxmiss100_ld001_eigenval/sum(prow_maf01_dp15_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp15_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp15_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp15_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf01_dp15_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf01_dp15_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp15_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp15_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp15_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf01_dp15_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld05 ####
prow_maf01_dp15_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf01_dp15_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp15_maxmiss100_ld05_pca <- prow_maf01_dp15_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf01_dp15_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf01_dp15_maxmiss100_ld05_pca)[2:ncol(prow_maf01_dp15_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp15_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp15_maxmiss100_ld05_bird <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld05_pca$ind))
prow_maf01_dp15_maxmiss100_ld05_bird[grep("PROW953", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf01_dp15_maxmiss100_ld05_bird[grep("PROW954", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf01_dp15_maxmiss100_ld05_bird[grep("PROW981", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf01_dp15_maxmiss100_ld05_bird[grep("PROW984", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp15_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld05_pca$ind))
prow_maf01_dp15_maxmiss100_ld05_nummites[grep("R1", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf01_dp15_maxmiss100_ld05_nummites[grep("R5", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf01_dp15_maxmiss100_ld05_nummites[grep("R20", prow_maf01_dp15_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp15_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf01_dp15_maxmiss100_ld05_pca, prow_maf01_dp15_maxmiss100_ld05_bird, prow_maf01_dp15_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp15_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp15_maxmiss100_ld05_eigenval/sum(prow_maf01_dp15_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp15_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp15_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp15_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf01_dp15_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf01_dp15_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp15_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp15_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp15_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf01_dp15_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld99 ####
prow_maf01_dp15_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf01_dp15_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp15_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp15_maxmiss100_ld99_pca <- prow_maf01_dp15_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf01_dp15_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf01_dp15_maxmiss100_ld99_pca)[2:ncol(prow_maf01_dp15_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp15_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp15_maxmiss100_ld99_bird <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld99_pca$ind))
prow_maf01_dp15_maxmiss100_ld99_bird[grep("PROW953", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf01_dp15_maxmiss100_ld99_bird[grep("PROW954", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf01_dp15_maxmiss100_ld99_bird[grep("PROW981", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf01_dp15_maxmiss100_ld99_bird[grep("PROW984", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp15_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf01_dp15_maxmiss100_ld99_pca$ind))
prow_maf01_dp15_maxmiss100_ld99_nummites[grep("R1", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf01_dp15_maxmiss100_ld99_nummites[grep("R5", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf01_dp15_maxmiss100_ld99_nummites[grep("R20", prow_maf01_dp15_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp15_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf01_dp15_maxmiss100_ld99_pca, prow_maf01_dp15_maxmiss100_ld99_bird, prow_maf01_dp15_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp15_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp15_maxmiss100_ld99_eigenval/sum(prow_maf01_dp15_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp15_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp15_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf01_dp15_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf01_dp15_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf01_dp15_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp15_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp15_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp15_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf01_dp15_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld001 ####
prow_maf01_dp20_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf01_dp20_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp20_maxmiss100_ld001_pca <- prow_maf01_dp20_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf01_dp20_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf01_dp20_maxmiss100_ld001_pca)[2:ncol(prow_maf01_dp20_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp20_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp20_maxmiss100_ld001_bird <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld001_pca$ind))
prow_maf01_dp20_maxmiss100_ld001_bird[grep("PROW953", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf01_dp20_maxmiss100_ld001_bird[grep("PROW954", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf01_dp20_maxmiss100_ld001_bird[grep("PROW981", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf01_dp20_maxmiss100_ld001_bird[grep("PROW984", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp20_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld001_pca$ind))
prow_maf01_dp20_maxmiss100_ld001_nummites[grep("R1", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf01_dp20_maxmiss100_ld001_nummites[grep("R5", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf01_dp20_maxmiss100_ld001_nummites[grep("R20", prow_maf01_dp20_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp20_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf01_dp20_maxmiss100_ld001_pca, prow_maf01_dp20_maxmiss100_ld001_bird, prow_maf01_dp20_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp20_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp20_maxmiss100_ld001_eigenval/sum(prow_maf01_dp20_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp20_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp20_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp20_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf01_dp20_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf01_dp20_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp20_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp20_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp20_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf01_dp20_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld05 ####
prow_maf01_dp20_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf01_dp20_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp20_maxmiss100_ld05_pca <- prow_maf01_dp20_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf01_dp20_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf01_dp20_maxmiss100_ld05_pca)[2:ncol(prow_maf01_dp20_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp20_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp20_maxmiss100_ld05_bird <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld05_pca$ind))
prow_maf01_dp20_maxmiss100_ld05_bird[grep("PROW953", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf01_dp20_maxmiss100_ld05_bird[grep("PROW954", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf01_dp20_maxmiss100_ld05_bird[grep("PROW981", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf01_dp20_maxmiss100_ld05_bird[grep("PROW984", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp20_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld05_pca$ind))
prow_maf01_dp20_maxmiss100_ld05_nummites[grep("R1", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf01_dp20_maxmiss100_ld05_nummites[grep("R5", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf01_dp20_maxmiss100_ld05_nummites[grep("R20", prow_maf01_dp20_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp20_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf01_dp20_maxmiss100_ld05_pca, prow_maf01_dp20_maxmiss100_ld05_bird, prow_maf01_dp20_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp20_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp20_maxmiss100_ld05_eigenval/sum(prow_maf01_dp20_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp20_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp20_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf01_dp20_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf01_dp20_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf01_dp20_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp20_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp20_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp20_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf01_dp20_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld99 ####
prow_maf01_dp20_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf01_dp20_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf01_dp20_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf01_dp20_maxmiss100_ld99_pca <- prow_maf01_dp20_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf01_dp20_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf01_dp20_maxmiss100_ld99_pca)[2:ncol(prow_maf01_dp20_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf01_dp20_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf01_dp20_maxmiss100_ld99_bird <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld99_pca$ind))
prow_maf01_dp20_maxmiss100_ld99_bird[grep("PROW953", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf01_dp20_maxmiss100_ld99_bird[grep("PROW954", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf01_dp20_maxmiss100_ld99_bird[grep("PROW981", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf01_dp20_maxmiss100_ld99_bird[grep("PROW984", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf01_dp20_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf01_dp20_maxmiss100_ld99_pca$ind))
prow_maf01_dp20_maxmiss100_ld99_nummites[grep("R1", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf01_dp20_maxmiss100_ld99_nummites[grep("R5", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf01_dp20_maxmiss100_ld99_nummites[grep("R20", prow_maf01_dp20_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf01_dp20_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf01_dp20_maxmiss100_ld99_pca, prow_maf01_dp20_maxmiss100_ld99_bird, prow_maf01_dp20_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf01_dp20_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf01_dp20_maxmiss100_ld99_eigenval/sum(prow_maf01_dp20_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf01_dp20_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf01_dp20_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf01_dp20_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf01_dp20_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf01_dp20_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf01_dp20_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf01_dp20_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf01_dp20_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf01_dp20_maxmiss100_ld99_pca_plot



























#### ++++ ALL PLOTS maf01, variable DPs, variable LD thresholds ####


prow_maf01<-
    ggarrange(
        prow_maf01_dp05_maxmiss100_ld001_pca_plot, prow_maf01_dp05_maxmiss100_ld05_pca_plot, prow_maf01_dp05_maxmiss100_ld99_pca_plot, prow_maf01_dp10_maxmiss100_ld001_pca_plot, prow_maf01_dp10_maxmiss100_ld05_pca_plot, prow_maf01_dp10_maxmiss100_ld99_pca_plot, prow_maf01_dp15_maxmiss100_ld001_pca_plot, prow_maf01_dp15_maxmiss100_ld05_pca_plot, prow_maf01_dp15_maxmiss100_ld99_pca_plot, prow_maf01_dp20_maxmiss100_ld001_pca_plot, prow_maf01_dp20_maxmiss100_ld05_pca_plot, prow_maf01_dp20_maxmiss100_ld99_pca_plot,
        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
        ncol = 3, nrow = 4, common.legend = TRUE)

annotate_figure(prow_maf01, top = text_grob("MAF = 0.01",
                                                 color = "black", face = "bold", size = 24))


# 10x14 portrait





















#################################################################################
########################### maf = 0.05 ##########################################
#################################################################################


#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001 ####
prow_maf05_dp05_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf05_dp05_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp05_maxmiss100_ld001_pca <- prow_maf05_dp05_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf05_dp05_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf05_dp05_maxmiss100_ld001_pca)[2:ncol(prow_maf05_dp05_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp05_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp05_maxmiss100_ld001_bird <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld001_pca$ind))
prow_maf05_dp05_maxmiss100_ld001_bird[grep("PROW953", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf05_dp05_maxmiss100_ld001_bird[grep("PROW954", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf05_dp05_maxmiss100_ld001_bird[grep("PROW981", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf05_dp05_maxmiss100_ld001_bird[grep("PROW984", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp05_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld001_pca$ind))
prow_maf05_dp05_maxmiss100_ld001_nummites[grep("R1", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf05_dp05_maxmiss100_ld001_nummites[grep("R5", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf05_dp05_maxmiss100_ld001_nummites[grep("R20", prow_maf05_dp05_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp05_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf05_dp05_maxmiss100_ld001_pca, prow_maf05_dp05_maxmiss100_ld001_bird, prow_maf05_dp05_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp05_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp05_maxmiss100_ld001_eigenval/sum(prow_maf05_dp05_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp05_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp05_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp05_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf05_dp05_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp05_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp05_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp05_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp05_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp05_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05 ####
prow_maf05_dp05_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf05_dp05_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp05_maxmiss100_ld05_pca <- prow_maf05_dp05_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf05_dp05_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf05_dp05_maxmiss100_ld05_pca)[2:ncol(prow_maf05_dp05_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp05_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp05_maxmiss100_ld05_bird <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld05_pca$ind))
prow_maf05_dp05_maxmiss100_ld05_bird[grep("PROW953", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf05_dp05_maxmiss100_ld05_bird[grep("PROW954", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf05_dp05_maxmiss100_ld05_bird[grep("PROW981", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf05_dp05_maxmiss100_ld05_bird[grep("PROW984", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp05_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld05_pca$ind))
prow_maf05_dp05_maxmiss100_ld05_nummites[grep("R1", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf05_dp05_maxmiss100_ld05_nummites[grep("R5", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf05_dp05_maxmiss100_ld05_nummites[grep("R20", prow_maf05_dp05_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp05_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf05_dp05_maxmiss100_ld05_pca, prow_maf05_dp05_maxmiss100_ld05_bird, prow_maf05_dp05_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp05_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp05_maxmiss100_ld05_eigenval/sum(prow_maf05_dp05_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp05_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp05_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp05_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf05_dp05_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf05_dp05_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp05_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp05_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp05_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf05_dp05_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99 ####
prow_maf05_dp05_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf05_dp05_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp05_maxmiss100_ld99_pca <- prow_maf05_dp05_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf05_dp05_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf05_dp05_maxmiss100_ld99_pca)[2:ncol(prow_maf05_dp05_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp05_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp05_maxmiss100_ld99_bird <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld99_pca$ind))
prow_maf05_dp05_maxmiss100_ld99_bird[grep("PROW953", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf05_dp05_maxmiss100_ld99_bird[grep("PROW954", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf05_dp05_maxmiss100_ld99_bird[grep("PROW981", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf05_dp05_maxmiss100_ld99_bird[grep("PROW984", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp05_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf05_dp05_maxmiss100_ld99_pca$ind))
prow_maf05_dp05_maxmiss100_ld99_nummites[grep("R1", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf05_dp05_maxmiss100_ld99_nummites[grep("R5", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf05_dp05_maxmiss100_ld99_nummites[grep("R20", prow_maf05_dp05_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp05_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf05_dp05_maxmiss100_ld99_pca, prow_maf05_dp05_maxmiss100_ld99_bird, prow_maf05_dp05_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp05_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp05_maxmiss100_ld99_eigenval/sum(prow_maf05_dp05_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp05_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp05_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf05_dp05_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf05_dp05_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp05_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp05_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp05_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp05_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp05_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001 ####
prow_maf05_dp10_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf05_dp10_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp10_maxmiss100_ld001_pca <- prow_maf05_dp10_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf05_dp10_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf05_dp10_maxmiss100_ld001_pca)[2:ncol(prow_maf05_dp10_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp10_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp10_maxmiss100_ld001_bird <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld001_pca$ind))
prow_maf05_dp10_maxmiss100_ld001_bird[grep("PROW953", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf05_dp10_maxmiss100_ld001_bird[grep("PROW954", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf05_dp10_maxmiss100_ld001_bird[grep("PROW981", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf05_dp10_maxmiss100_ld001_bird[grep("PROW984", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp10_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld001_pca$ind))
prow_maf05_dp10_maxmiss100_ld001_nummites[grep("R1", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf05_dp10_maxmiss100_ld001_nummites[grep("R5", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf05_dp10_maxmiss100_ld001_nummites[grep("R20", prow_maf05_dp10_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp10_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf05_dp10_maxmiss100_ld001_pca, prow_maf05_dp10_maxmiss100_ld001_bird, prow_maf05_dp10_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp10_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp10_maxmiss100_ld001_eigenval/sum(prow_maf05_dp10_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp10_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp10_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp10_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld001_pca_plot





prow_maf05_dp10_maxmiss100_ld001_pca_plot_cb<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld001_pca_plot_cb





prow_maf05_dp10_maxmiss100_ld001_pca_plot_viridis<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 5, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)")) +
    theme(aspect.ratio = 1)

prow_maf05_dp10_maxmiss100_ld001_pca_plot_viridis







#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05 ####
prow_maf05_dp10_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf05_dp10_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp10_maxmiss100_ld05_pca <- prow_maf05_dp10_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf05_dp10_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf05_dp10_maxmiss100_ld05_pca)[2:ncol(prow_maf05_dp10_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp10_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp10_maxmiss100_ld05_bird <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld05_pca$ind))
prow_maf05_dp10_maxmiss100_ld05_bird[grep("PROW953", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf05_dp10_maxmiss100_ld05_bird[grep("PROW954", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf05_dp10_maxmiss100_ld05_bird[grep("PROW981", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf05_dp10_maxmiss100_ld05_bird[grep("PROW984", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp10_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld05_pca$ind))
prow_maf05_dp10_maxmiss100_ld05_nummites[grep("R1", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf05_dp10_maxmiss100_ld05_nummites[grep("R5", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf05_dp10_maxmiss100_ld05_nummites[grep("R20", prow_maf05_dp10_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp10_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf05_dp10_maxmiss100_ld05_pca, prow_maf05_dp10_maxmiss100_ld05_bird, prow_maf05_dp10_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp10_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp10_maxmiss100_ld05_eigenval/sum(prow_maf05_dp10_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp10_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp10_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp10_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf05_dp10_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99 ####
prow_maf05_dp10_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf05_dp10_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp10_maxmiss100_ld99_pca <- prow_maf05_dp10_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf05_dp10_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf05_dp10_maxmiss100_ld99_pca)[2:ncol(prow_maf05_dp10_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp10_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp10_maxmiss100_ld99_bird <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld99_pca$ind))
prow_maf05_dp10_maxmiss100_ld99_bird[grep("PROW953", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf05_dp10_maxmiss100_ld99_bird[grep("PROW954", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf05_dp10_maxmiss100_ld99_bird[grep("PROW981", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf05_dp10_maxmiss100_ld99_bird[grep("PROW984", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp10_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf05_dp10_maxmiss100_ld99_pca$ind))
prow_maf05_dp10_maxmiss100_ld99_nummites[grep("R1", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf05_dp10_maxmiss100_ld99_nummites[grep("R5", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf05_dp10_maxmiss100_ld99_nummites[grep("R20", prow_maf05_dp10_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp10_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf05_dp10_maxmiss100_ld99_pca, prow_maf05_dp10_maxmiss100_ld99_bird, prow_maf05_dp10_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp10_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp10_maxmiss100_ld99_eigenval/sum(prow_maf05_dp10_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp10_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp10_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf05_dp10_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf05_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld99_pca_plot



prow_maf05_dp10_maxmiss100_ld99_pca_plot_cb<-
    ggplot(prow_maf05_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld99_pca_plot_cb























#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001 ####
prow_maf05_dp15_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf05_dp15_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp15_maxmiss100_ld001_pca <- prow_maf05_dp15_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf05_dp15_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf05_dp15_maxmiss100_ld001_pca)[2:ncol(prow_maf05_dp15_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp15_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp15_maxmiss100_ld001_bird <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld001_pca$ind))
prow_maf05_dp15_maxmiss100_ld001_bird[grep("PROW953", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf05_dp15_maxmiss100_ld001_bird[grep("PROW954", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf05_dp15_maxmiss100_ld001_bird[grep("PROW981", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf05_dp15_maxmiss100_ld001_bird[grep("PROW984", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp15_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld001_pca$ind))
prow_maf05_dp15_maxmiss100_ld001_nummites[grep("R1", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf05_dp15_maxmiss100_ld001_nummites[grep("R5", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf05_dp15_maxmiss100_ld001_nummites[grep("R20", prow_maf05_dp15_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp15_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf05_dp15_maxmiss100_ld001_pca, prow_maf05_dp15_maxmiss100_ld001_bird, prow_maf05_dp15_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp15_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp15_maxmiss100_ld001_eigenval/sum(prow_maf05_dp15_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp15_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp15_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp15_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf05_dp15_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp15_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp15_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp15_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp15_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp15_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05 ####
prow_maf05_dp15_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf05_dp15_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp15_maxmiss100_ld05_pca <- prow_maf05_dp15_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf05_dp15_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf05_dp15_maxmiss100_ld05_pca)[2:ncol(prow_maf05_dp15_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp15_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp15_maxmiss100_ld05_bird <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld05_pca$ind))
prow_maf05_dp15_maxmiss100_ld05_bird[grep("PROW953", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf05_dp15_maxmiss100_ld05_bird[grep("PROW954", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf05_dp15_maxmiss100_ld05_bird[grep("PROW981", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf05_dp15_maxmiss100_ld05_bird[grep("PROW984", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp15_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld05_pca$ind))
prow_maf05_dp15_maxmiss100_ld05_nummites[grep("R1", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf05_dp15_maxmiss100_ld05_nummites[grep("R5", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf05_dp15_maxmiss100_ld05_nummites[grep("R20", prow_maf05_dp15_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp15_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf05_dp15_maxmiss100_ld05_pca, prow_maf05_dp15_maxmiss100_ld05_bird, prow_maf05_dp15_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp15_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp15_maxmiss100_ld05_eigenval/sum(prow_maf05_dp15_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp15_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp15_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp15_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf05_dp15_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf05_dp15_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp15_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp15_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp15_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf05_dp15_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld99 ####
prow_maf05_dp15_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf05_dp15_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp15_maxmiss100_ld99_pca <- prow_maf05_dp15_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf05_dp15_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf05_dp15_maxmiss100_ld99_pca)[2:ncol(prow_maf05_dp15_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp15_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp15_maxmiss100_ld99_bird <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld99_pca$ind))
prow_maf05_dp15_maxmiss100_ld99_bird[grep("PROW953", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf05_dp15_maxmiss100_ld99_bird[grep("PROW954", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf05_dp15_maxmiss100_ld99_bird[grep("PROW981", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf05_dp15_maxmiss100_ld99_bird[grep("PROW984", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp15_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf05_dp15_maxmiss100_ld99_pca$ind))
prow_maf05_dp15_maxmiss100_ld99_nummites[grep("R1", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf05_dp15_maxmiss100_ld99_nummites[grep("R5", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf05_dp15_maxmiss100_ld99_nummites[grep("R20", prow_maf05_dp15_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp15_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf05_dp15_maxmiss100_ld99_pca, prow_maf05_dp15_maxmiss100_ld99_bird, prow_maf05_dp15_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp15_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp15_maxmiss100_ld99_eigenval/sum(prow_maf05_dp15_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp15_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp15_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf05_dp15_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf05_dp15_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp15_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp15_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp15_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp15_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp15_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001 ####
prow_maf05_dp20_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf05_dp20_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp20_maxmiss100_ld001_pca <- prow_maf05_dp20_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf05_dp20_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf05_dp20_maxmiss100_ld001_pca)[2:ncol(prow_maf05_dp20_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp20_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp20_maxmiss100_ld001_bird <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld001_pca$ind))
prow_maf05_dp20_maxmiss100_ld001_bird[grep("PROW953", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf05_dp20_maxmiss100_ld001_bird[grep("PROW954", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf05_dp20_maxmiss100_ld001_bird[grep("PROW981", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf05_dp20_maxmiss100_ld001_bird[grep("PROW984", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp20_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld001_pca$ind))
prow_maf05_dp20_maxmiss100_ld001_nummites[grep("R1", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf05_dp20_maxmiss100_ld001_nummites[grep("R5", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf05_dp20_maxmiss100_ld001_nummites[grep("R20", prow_maf05_dp20_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp20_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf05_dp20_maxmiss100_ld001_pca, prow_maf05_dp20_maxmiss100_ld001_bird, prow_maf05_dp20_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp20_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp20_maxmiss100_ld001_eigenval/sum(prow_maf05_dp20_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp20_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp20_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp20_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf05_dp20_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp20_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp20_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp20_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp20_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp20_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05 ####
prow_maf05_dp20_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf05_dp20_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp20_maxmiss100_ld05_pca <- prow_maf05_dp20_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf05_dp20_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf05_dp20_maxmiss100_ld05_pca)[2:ncol(prow_maf05_dp20_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp20_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp20_maxmiss100_ld05_bird <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld05_pca$ind))
prow_maf05_dp20_maxmiss100_ld05_bird[grep("PROW953", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf05_dp20_maxmiss100_ld05_bird[grep("PROW954", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf05_dp20_maxmiss100_ld05_bird[grep("PROW981", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf05_dp20_maxmiss100_ld05_bird[grep("PROW984", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp20_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld05_pca$ind))
prow_maf05_dp20_maxmiss100_ld05_nummites[grep("R1", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf05_dp20_maxmiss100_ld05_nummites[grep("R5", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf05_dp20_maxmiss100_ld05_nummites[grep("R20", prow_maf05_dp20_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp20_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf05_dp20_maxmiss100_ld05_pca, prow_maf05_dp20_maxmiss100_ld05_bird, prow_maf05_dp20_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp20_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp20_maxmiss100_ld05_eigenval/sum(prow_maf05_dp20_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp20_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp20_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf05_dp20_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf05_dp20_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf05_dp20_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp20_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp20_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp20_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf05_dp20_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld99 ####
prow_maf05_dp20_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf05_dp20_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf05_dp20_maxmiss100_ld99_pca <- prow_maf05_dp20_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf05_dp20_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf05_dp20_maxmiss100_ld99_pca)[2:ncol(prow_maf05_dp20_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf05_dp20_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf05_dp20_maxmiss100_ld99_bird <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld99_pca$ind))
prow_maf05_dp20_maxmiss100_ld99_bird[grep("PROW953", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf05_dp20_maxmiss100_ld99_bird[grep("PROW954", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf05_dp20_maxmiss100_ld99_bird[grep("PROW981", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf05_dp20_maxmiss100_ld99_bird[grep("PROW984", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf05_dp20_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf05_dp20_maxmiss100_ld99_pca$ind))
prow_maf05_dp20_maxmiss100_ld99_nummites[grep("R1", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf05_dp20_maxmiss100_ld99_nummites[grep("R5", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf05_dp20_maxmiss100_ld99_nummites[grep("R20", prow_maf05_dp20_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf05_dp20_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf05_dp20_maxmiss100_ld99_pca, prow_maf05_dp20_maxmiss100_ld99_bird, prow_maf05_dp20_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf05_dp20_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf05_dp20_maxmiss100_ld99_eigenval/sum(prow_maf05_dp20_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf05_dp20_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf05_dp20_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf05_dp20_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf05_dp20_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp20_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp20_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp20_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp20_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp20_maxmiss100_ld99_pca_plot



























#### ++++ ALL PLOTS maf05, variable DPs, variable LD thresholds ####


prow_maf05<-
    ggarrange(
        prow_maf05_dp05_maxmiss100_ld001_pca_plot, prow_maf05_dp05_maxmiss100_ld05_pca_plot, prow_maf05_dp05_maxmiss100_ld99_pca_plot, prow_maf05_dp10_maxmiss100_ld001_pca_plot, prow_maf05_dp10_maxmiss100_ld05_pca_plot, prow_maf05_dp10_maxmiss100_ld99_pca_plot, prow_maf05_dp15_maxmiss100_ld001_pca_plot, prow_maf05_dp15_maxmiss100_ld05_pca_plot, prow_maf05_dp15_maxmiss100_ld99_pca_plot, prow_maf05_dp20_maxmiss100_ld001_pca_plot, prow_maf05_dp20_maxmiss100_ld05_pca_plot, prow_maf05_dp20_maxmiss100_ld99_pca_plot,
        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
        ncol = 3, nrow = 4, common.legend = TRUE)

annotate_figure(prow_maf05, top = text_grob("MAF = 0.05",
                                            color = "black", face = "bold", size = 24))


# 10x14 portrait





prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2



prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_cb<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_cb


prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_viridis<-
    ggplot(prow_maf05_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 5, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)")) +
    theme(aspect.ratio = 1)

prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_viridis





prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2<-
    ggplot(prow_maf05_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2

prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_cb<-
    ggplot(prow_maf05_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("#90005d", "#b085b9", "#174417", "#7cbf7c")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_cb



prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_viridis<-
    ggplot(prow_maf05_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf05_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 5, alpha = 0.8, aes(shape=prow_maf05_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("#440154", "#414487", "#22a884", "#bddf26")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf05_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)")) +
    theme(aspect.ratio = 1)

prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_viridis





prow_maf05_dp10_v2<-
    ggarrange(
        prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2, prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2,
        labels = c("A", "B"),
        ncol = 2, nrow = 1, common.legend = TRUE)

# 4x9 landscape


prow_maf05_dp10_v2_cb<-
    ggarrange(
        prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_cb, prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_cb,
        labels = c("A", "B"),
        ncol = 2, nrow = 1, common.legend = TRUE)

# 4x9 landscape


prow_maf05_dp10_v2_viridis<-
    ggarrange(
        prow_maf05_dp10_maxmiss100_ld001_pca_plot_v2_viridis, prow_maf05_dp10_maxmiss100_ld99_pca_plot_v2_viridis,
        labels = c("A", "B"),
        font.label = list(size = 24),
        ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

# 4x9 landscape








































#################################################################################
########################### maf = 0.10 ##########################################
#################################################################################



#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld001 ####
prow_maf10_dp05_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf10_dp05_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp05_maxmiss100_ld001_pca <- prow_maf10_dp05_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf10_dp05_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf10_dp05_maxmiss100_ld001_pca)[2:ncol(prow_maf10_dp05_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp05_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp05_maxmiss100_ld001_bird <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld001_pca$ind))
prow_maf10_dp05_maxmiss100_ld001_bird[grep("PROW953", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf10_dp05_maxmiss100_ld001_bird[grep("PROW954", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf10_dp05_maxmiss100_ld001_bird[grep("PROW981", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf10_dp05_maxmiss100_ld001_bird[grep("PROW984", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp05_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld001_pca$ind))
prow_maf10_dp05_maxmiss100_ld001_nummites[grep("R1", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf10_dp05_maxmiss100_ld001_nummites[grep("R5", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf10_dp05_maxmiss100_ld001_nummites[grep("R20", prow_maf10_dp05_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp05_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf10_dp05_maxmiss100_ld001_pca, prow_maf10_dp05_maxmiss100_ld001_bird, prow_maf10_dp05_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp05_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp05_maxmiss100_ld001_eigenval/sum(prow_maf10_dp05_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp05_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp05_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp05_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf10_dp05_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf10_dp05_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp05_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp05_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp05_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf10_dp05_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld05 ####
prow_maf10_dp05_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf10_dp05_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp05_maxmiss100_ld05_pca <- prow_maf10_dp05_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf10_dp05_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf10_dp05_maxmiss100_ld05_pca)[2:ncol(prow_maf10_dp05_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp05_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp05_maxmiss100_ld05_bird <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld05_pca$ind))
prow_maf10_dp05_maxmiss100_ld05_bird[grep("PROW953", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf10_dp05_maxmiss100_ld05_bird[grep("PROW954", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf10_dp05_maxmiss100_ld05_bird[grep("PROW981", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf10_dp05_maxmiss100_ld05_bird[grep("PROW984", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp05_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld05_pca$ind))
prow_maf10_dp05_maxmiss100_ld05_nummites[grep("R1", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf10_dp05_maxmiss100_ld05_nummites[grep("R5", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf10_dp05_maxmiss100_ld05_nummites[grep("R20", prow_maf10_dp05_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp05_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf10_dp05_maxmiss100_ld05_pca, prow_maf10_dp05_maxmiss100_ld05_bird, prow_maf10_dp05_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp05_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp05_maxmiss100_ld05_eigenval/sum(prow_maf10_dp05_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp05_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp05_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp05_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf10_dp05_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf10_dp05_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp05_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp05_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp05_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf10_dp05_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99 ####
prow_maf10_dp05_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf10_dp05_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp05_maxmiss100_ld99_pca <- prow_maf10_dp05_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf10_dp05_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf10_dp05_maxmiss100_ld99_pca)[2:ncol(prow_maf10_dp05_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp05_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp05_maxmiss100_ld99_bird <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld99_pca$ind))
prow_maf10_dp05_maxmiss100_ld99_bird[grep("PROW953", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf10_dp05_maxmiss100_ld99_bird[grep("PROW954", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf10_dp05_maxmiss100_ld99_bird[grep("PROW981", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf10_dp05_maxmiss100_ld99_bird[grep("PROW984", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp05_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf10_dp05_maxmiss100_ld99_pca$ind))
prow_maf10_dp05_maxmiss100_ld99_nummites[grep("R1", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf10_dp05_maxmiss100_ld99_nummites[grep("R5", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf10_dp05_maxmiss100_ld99_nummites[grep("R20", prow_maf10_dp05_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp05_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf10_dp05_maxmiss100_ld99_pca, prow_maf10_dp05_maxmiss100_ld99_bird, prow_maf10_dp05_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp05_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp05_maxmiss100_ld99_eigenval/sum(prow_maf10_dp05_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp05_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp05_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf10_dp05_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf10_dp05_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf10_dp05_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp05_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 5X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp05_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp05_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf10_dp05_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld001 ####
prow_maf10_dp10_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf10_dp10_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp10_maxmiss100_ld001_pca <- prow_maf10_dp10_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf10_dp10_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf10_dp10_maxmiss100_ld001_pca)[2:ncol(prow_maf10_dp10_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp10_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp10_maxmiss100_ld001_bird <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld001_pca$ind))
prow_maf10_dp10_maxmiss100_ld001_bird[grep("PROW953", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf10_dp10_maxmiss100_ld001_bird[grep("PROW954", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf10_dp10_maxmiss100_ld001_bird[grep("PROW981", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf10_dp10_maxmiss100_ld001_bird[grep("PROW984", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp10_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld001_pca$ind))
prow_maf10_dp10_maxmiss100_ld001_nummites[grep("R1", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf10_dp10_maxmiss100_ld001_nummites[grep("R5", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf10_dp10_maxmiss100_ld001_nummites[grep("R20", prow_maf10_dp10_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp10_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf10_dp10_maxmiss100_ld001_pca, prow_maf10_dp10_maxmiss100_ld001_bird, prow_maf10_dp10_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp10_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp10_maxmiss100_ld001_eigenval/sum(prow_maf10_dp10_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp10_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp10_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp10_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf10_dp10_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf10_dp10_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp10_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp10_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp10_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf10_dp10_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld05 ####
prow_maf10_dp10_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf10_dp10_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp10_maxmiss100_ld05_pca <- prow_maf10_dp10_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf10_dp10_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf10_dp10_maxmiss100_ld05_pca)[2:ncol(prow_maf10_dp10_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp10_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp10_maxmiss100_ld05_bird <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld05_pca$ind))
prow_maf10_dp10_maxmiss100_ld05_bird[grep("PROW953", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf10_dp10_maxmiss100_ld05_bird[grep("PROW954", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf10_dp10_maxmiss100_ld05_bird[grep("PROW981", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf10_dp10_maxmiss100_ld05_bird[grep("PROW984", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp10_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld05_pca$ind))
prow_maf10_dp10_maxmiss100_ld05_nummites[grep("R1", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf10_dp10_maxmiss100_ld05_nummites[grep("R5", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf10_dp10_maxmiss100_ld05_nummites[grep("R20", prow_maf10_dp10_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp10_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf10_dp10_maxmiss100_ld05_pca, prow_maf10_dp10_maxmiss100_ld05_bird, prow_maf10_dp10_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp10_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp10_maxmiss100_ld05_eigenval/sum(prow_maf10_dp10_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp10_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp10_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp10_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf10_dp10_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf10_dp10_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp10_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp10_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp10_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf10_dp10_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99 ####
prow_maf10_dp10_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf10_dp10_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp10_maxmiss100_ld99_pca <- prow_maf10_dp10_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf10_dp10_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf10_dp10_maxmiss100_ld99_pca)[2:ncol(prow_maf10_dp10_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp10_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp10_maxmiss100_ld99_bird <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld99_pca$ind))
prow_maf10_dp10_maxmiss100_ld99_bird[grep("PROW953", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf10_dp10_maxmiss100_ld99_bird[grep("PROW954", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf10_dp10_maxmiss100_ld99_bird[grep("PROW981", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf10_dp10_maxmiss100_ld99_bird[grep("PROW984", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp10_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf10_dp10_maxmiss100_ld99_pca$ind))
prow_maf10_dp10_maxmiss100_ld99_nummites[grep("R1", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf10_dp10_maxmiss100_ld99_nummites[grep("R5", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf10_dp10_maxmiss100_ld99_nummites[grep("R20", prow_maf10_dp10_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp10_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf10_dp10_maxmiss100_ld99_pca, prow_maf10_dp10_maxmiss100_ld99_bird, prow_maf10_dp10_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp10_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp10_maxmiss100_ld99_eigenval/sum(prow_maf10_dp10_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp10_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp10_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf10_dp10_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf10_dp10_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf10_dp10_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp10_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 10X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp10_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp10_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf10_dp10_maxmiss100_ld99_pca_plot



























#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld001 ####
prow_maf10_dp15_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf10_dp15_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp15_maxmiss100_ld001_pca <- prow_maf10_dp15_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf10_dp15_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf10_dp15_maxmiss100_ld001_pca)[2:ncol(prow_maf10_dp15_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp15_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp15_maxmiss100_ld001_bird <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld001_pca$ind))
prow_maf10_dp15_maxmiss100_ld001_bird[grep("PROW953", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf10_dp15_maxmiss100_ld001_bird[grep("PROW954", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf10_dp15_maxmiss100_ld001_bird[grep("PROW981", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf10_dp15_maxmiss100_ld001_bird[grep("PROW984", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp15_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld001_pca$ind))
prow_maf10_dp15_maxmiss100_ld001_nummites[grep("R1", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf10_dp15_maxmiss100_ld001_nummites[grep("R5", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf10_dp15_maxmiss100_ld001_nummites[grep("R20", prow_maf10_dp15_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp15_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf10_dp15_maxmiss100_ld001_pca, prow_maf10_dp15_maxmiss100_ld001_bird, prow_maf10_dp15_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp15_maxmiss100_ld001_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp15_maxmiss100_ld001_eigenval/sum(prow_maf10_dp15_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp15_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp15_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp15_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf10_dp15_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf10_dp15_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp15_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp15_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp15_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf10_dp15_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld05 ####
prow_maf10_dp15_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf10_dp15_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp15_maxmiss100_ld05_pca <- prow_maf10_dp15_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf10_dp15_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf10_dp15_maxmiss100_ld05_pca)[2:ncol(prow_maf10_dp15_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp15_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp15_maxmiss100_ld05_bird <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld05_pca$ind))
prow_maf10_dp15_maxmiss100_ld05_bird[grep("PROW953", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf10_dp15_maxmiss100_ld05_bird[grep("PROW954", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf10_dp15_maxmiss100_ld05_bird[grep("PROW981", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf10_dp15_maxmiss100_ld05_bird[grep("PROW984", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp15_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld05_pca$ind))
prow_maf10_dp15_maxmiss100_ld05_nummites[grep("R1", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf10_dp15_maxmiss100_ld05_nummites[grep("R5", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf10_dp15_maxmiss100_ld05_nummites[grep("R20", prow_maf10_dp15_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp15_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf10_dp15_maxmiss100_ld05_pca, prow_maf10_dp15_maxmiss100_ld05_bird, prow_maf10_dp15_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp15_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp15_maxmiss100_ld05_eigenval/sum(prow_maf10_dp15_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp15_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp15_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp15_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf10_dp15_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf10_dp15_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp15_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp15_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp15_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf10_dp15_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld99 ####
prow_maf10_dp15_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf10_dp15_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp15_maxmiss100_ld99_pca <- prow_maf10_dp15_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf10_dp15_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf10_dp15_maxmiss100_ld99_pca)[2:ncol(prow_maf10_dp15_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp15_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp15_maxmiss100_ld99_bird <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld99_pca$ind))
prow_maf10_dp15_maxmiss100_ld99_bird[grep("PROW953", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf10_dp15_maxmiss100_ld99_bird[grep("PROW954", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf10_dp15_maxmiss100_ld99_bird[grep("PROW981", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf10_dp15_maxmiss100_ld99_bird[grep("PROW984", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp15_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf10_dp15_maxmiss100_ld99_pca$ind))
prow_maf10_dp15_maxmiss100_ld99_nummites[grep("R1", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf10_dp15_maxmiss100_ld99_nummites[grep("R5", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf10_dp15_maxmiss100_ld99_nummites[grep("R20", prow_maf10_dp15_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp15_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf10_dp15_maxmiss100_ld99_pca, prow_maf10_dp15_maxmiss100_ld99_bird, prow_maf10_dp15_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp15_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp15_maxmiss100_ld99_eigenval/sum(prow_maf10_dp15_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp15_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp15_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf10_dp15_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf10_dp15_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf10_dp15_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp15_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 15X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp15_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp15_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf10_dp15_maxmiss100_ld99_pca_plot





















#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld001 ####
prow_maf10_dp20_maxmiss100_ld001_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld001.eigenvec", col_names = FALSE)
prow_maf10_dp20_maxmiss100_ld001_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld001.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp20_maxmiss100_ld001_pca <- prow_maf10_dp20_maxmiss100_ld001_pca[,-1]
# set names
names(prow_maf10_dp20_maxmiss100_ld001_pca)[1] <- "ind"
names(prow_maf10_dp20_maxmiss100_ld001_pca)[2:ncol(prow_maf10_dp20_maxmiss100_ld001_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp20_maxmiss100_ld001_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp20_maxmiss100_ld001_bird <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld001_pca$ind))
prow_maf10_dp20_maxmiss100_ld001_bird[grep("PROW953", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "PROW953"
prow_maf10_dp20_maxmiss100_ld001_bird[grep("PROW954", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "PROW954"
prow_maf10_dp20_maxmiss100_ld001_bird[grep("PROW981", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "PROW981"
prow_maf10_dp20_maxmiss100_ld001_bird[grep("PROW984", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp20_maxmiss100_ld001_nummites <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld001_pca$ind))
prow_maf10_dp20_maxmiss100_ld001_nummites[grep("R1", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "R1"
prow_maf10_dp20_maxmiss100_ld001_nummites[grep("R5", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "R5"
prow_maf10_dp20_maxmiss100_ld001_nummites[grep("R20", prow_maf10_dp20_maxmiss100_ld001_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp20_maxmiss100_ld001_pca <- as_tibble(data.frame(prow_maf10_dp20_maxmiss100_ld001_pca, prow_maf10_dp20_maxmiss100_ld001_bird, prow_maf10_dp20_maxmiss100_ld001_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp20_maxmiss100_ld001_pve <- data.frame(PC = 1:11, pve = prow_maf10_dp20_maxmiss100_ld001_eigenval/sum(prow_maf10_dp20_maxmiss100_ld001_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp20_maxmiss100_ld001_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp20_maxmiss100_ld001_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp20_maxmiss100_ld001_pca_plot<-
    ggplot(prow_maf10_dp20_maxmiss100_ld001_pca, aes(PC1, PC2, col = prow_maf10_dp20_maxmiss100_ld001_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp20_maxmiss100_ld001_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.01)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp20_maxmiss100_ld001_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp20_maxmiss100_ld001_pve$pve[2], 3), "%)"))

prow_maf10_dp20_maxmiss100_ld001_pca_plot








#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld05 ####
prow_maf10_dp20_maxmiss100_ld05_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld05.eigenvec", col_names = FALSE)
prow_maf10_dp20_maxmiss100_ld05_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld05.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp20_maxmiss100_ld05_pca <- prow_maf10_dp20_maxmiss100_ld05_pca[,-1]
# set names
names(prow_maf10_dp20_maxmiss100_ld05_pca)[1] <- "ind"
names(prow_maf10_dp20_maxmiss100_ld05_pca)[2:ncol(prow_maf10_dp20_maxmiss100_ld05_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp20_maxmiss100_ld05_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp20_maxmiss100_ld05_bird <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld05_pca$ind))
prow_maf10_dp20_maxmiss100_ld05_bird[grep("PROW953", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "PROW953"
prow_maf10_dp20_maxmiss100_ld05_bird[grep("PROW954", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "PROW954"
prow_maf10_dp20_maxmiss100_ld05_bird[grep("PROW981", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "PROW981"
prow_maf10_dp20_maxmiss100_ld05_bird[grep("PROW984", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp20_maxmiss100_ld05_nummites <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld05_pca$ind))
prow_maf10_dp20_maxmiss100_ld05_nummites[grep("R1", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "R1"
prow_maf10_dp20_maxmiss100_ld05_nummites[grep("R5", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "R5"
prow_maf10_dp20_maxmiss100_ld05_nummites[grep("R20", prow_maf10_dp20_maxmiss100_ld05_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp20_maxmiss100_ld05_pca <- as_tibble(data.frame(prow_maf10_dp20_maxmiss100_ld05_pca, prow_maf10_dp20_maxmiss100_ld05_bird, prow_maf10_dp20_maxmiss100_ld05_nummites))



# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp20_maxmiss100_ld05_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp20_maxmiss100_ld05_eigenval/sum(prow_maf10_dp20_maxmiss100_ld05_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp20_maxmiss100_ld05_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp20_maxmiss100_ld05_pve$pve)

# plot pca by bird/num mites
prow_maf10_dp20_maxmiss100_ld05_pca_plot<-
    ggplot(prow_maf10_dp20_maxmiss100_ld05_pca, aes(PC1, PC2, col = prow_maf10_dp20_maxmiss100_ld05_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp20_maxmiss100_ld05_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.5)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp20_maxmiss100_ld05_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp20_maxmiss100_ld05_pve$pve[2], 3), "%)"))

prow_maf10_dp20_maxmiss100_ld05_pca_plot





#### ++ PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld99 ####
prow_maf10_dp20_maxmiss100_ld99_pca <- read_table("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld99.eigenvec", col_names = FALSE)
prow_maf10_dp20_maxmiss100_ld99_eigenval <- scan("./PROW_EXP_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100_ld99.eigenval")


# Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

# sort out the pca data
# remove nuisance column
prow_maf10_dp20_maxmiss100_ld99_pca <- prow_maf10_dp20_maxmiss100_ld99_pca[,-1]
# set names
names(prow_maf10_dp20_maxmiss100_ld99_pca)[1] <- "ind"
names(prow_maf10_dp20_maxmiss100_ld99_pca)[2:ncol(prow_maf10_dp20_maxmiss100_ld99_pca)] <- paste0("PC", 1:(ncol(prow_maf10_dp20_maxmiss100_ld99_pca)-1))

# Next we can add a species, location and if required, a species x location vector. We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual pops, two ways
# bird
prow_maf10_dp20_maxmiss100_ld99_bird <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld99_pca$ind))
prow_maf10_dp20_maxmiss100_ld99_bird[grep("PROW953", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "PROW953"
prow_maf10_dp20_maxmiss100_ld99_bird[grep("PROW954", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "PROW954"
prow_maf10_dp20_maxmiss100_ld99_bird[grep("PROW981", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "PROW981"
prow_maf10_dp20_maxmiss100_ld99_bird[grep("PROW984", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "PROW984"

# number of mites
prow_maf10_dp20_maxmiss100_ld99_nummites <- rep(NA, length(prow_maf10_dp20_maxmiss100_ld99_pca$ind))
prow_maf10_dp20_maxmiss100_ld99_nummites[grep("R1", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "R1"
prow_maf10_dp20_maxmiss100_ld99_nummites[grep("R5", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "R5"
prow_maf10_dp20_maxmiss100_ld99_nummites[grep("R20", prow_maf10_dp20_maxmiss100_ld99_pca$ind)] <- "R20"


# With these variables created, we can remake our data.frame like so. Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
prow_maf10_dp20_maxmiss100_ld99_pca <- as_tibble(data.frame(prow_maf10_dp20_maxmiss100_ld99_pca, prow_maf10_dp20_maxmiss100_ld99_bird, prow_maf10_dp20_maxmiss100_ld99_nummites))


# Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).

# first convert to percentage variance explained
prow_maf10_dp20_maxmiss100_ld99_pve <- data.frame(PC = 1:12, pve = prow_maf10_dp20_maxmiss100_ld99_eigenval/sum(prow_maf10_dp20_maxmiss100_ld99_eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(prow_maf10_dp20_maxmiss100_ld99_pve, aes(PC, pve)) + geom_bar(stat = "identity") +
    ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(prow_maf10_dp20_maxmiss100_ld99_pve$pve)


# plot pca by bird/num mites
prow_maf10_dp20_maxmiss100_ld99_pca_plot<-
    ggplot(prow_maf10_dp20_maxmiss100_ld99_pca, aes(PC1, PC2, col = prow_maf10_dp20_maxmiss100_ld99_bird)) +
    geom_point(size = 4, alpha = 0.8, aes(shape=prow_maf10_dp20_maxmiss100_ld99_nummites)) +
    scale_color_manual(values=c("thistle4", "pink3", "darkslategray4", "darkseagreen4")) +
    scale_shape_discrete(labels=c("1", "20", "5")) +
    coord_equal() +
    theme_light() +
    labs(colour = "Bird ID", shape = "Number of Mites") +
    ggtitle(bquote('DP: 20X; LD'~r^2:0.99)) +
    xlab(paste0("PC1 (", signif(prow_maf10_dp20_maxmiss100_ld99_pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(prow_maf10_dp20_maxmiss100_ld99_pve$pve[2], 3), "%)"))

prow_maf10_dp20_maxmiss100_ld99_pca_plot



























#### ++++ ALL PLOTS maf10, variable DPs, variable LD thresholds ####


prow_maf10<-
    ggarrange(
        prow_maf10_dp05_maxmiss100_ld001_pca_plot, prow_maf10_dp05_maxmiss100_ld05_pca_plot, prow_maf10_dp05_maxmiss100_ld99_pca_plot, prow_maf10_dp10_maxmiss100_ld001_pca_plot, prow_maf10_dp10_maxmiss100_ld05_pca_plot, prow_maf10_dp10_maxmiss100_ld99_pca_plot, prow_maf10_dp15_maxmiss100_ld001_pca_plot, prow_maf10_dp15_maxmiss100_ld05_pca_plot, prow_maf10_dp15_maxmiss100_ld99_pca_plot, prow_maf10_dp20_maxmiss100_ld001_pca_plot, prow_maf10_dp20_maxmiss100_ld05_pca_plot, prow_maf10_dp20_maxmiss100_ld99_pca_plot,
        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
        ncol = 3, nrow = 4, common.legend = TRUE)

annotate_figure(prow_maf10, top = text_grob("MAF = 0.10",
                                            color = "black", face = "bold", size = 24))


# 10x14 portrait


