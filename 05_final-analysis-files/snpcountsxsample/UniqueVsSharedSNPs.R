#### Shared SNPs within populations and across populations ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/04_SNP_20221013/snpcountsxsample")

#### LOAD LIBRARIES ####

library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(tidyverse)

# install.packages("ggupset")
library(ggupset)



#### MAF05_DP10_LD001 ####
maf05_dp10_ld001 <- read.table("maf05_dp10_ld001.txt", header = TRUE)
str(maf05_dp10_ld001)
maf05_dp10_ld001$num_mites<-as.factor(maf05_dp10_ld001$num_mites)

maf05_dp10_ld001 <- maf05_dp10_ld001 %>%
    arrange(node.1) %>%
    mutate(nodexsite = factor(nodexsite, unique(nodexsite)))

str(maf05_dp10_ld001)
levels(maf05_dp10_ld001$bird_id)




#### ++ PROW953 ####
maf05_dp10_ld001_PROW953<-subset(maf05_dp10_ld001, bird_id =="PROW953")
levels(maf05_dp10_ld001_PROW953$bird_id)
maf05_dp10_ld001_PROW953$bird_id <- factor(maf05_dp10_ld001_PROW953$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW953$bird_id)

maf05_dp10_ld001_PROW953_list <-
    maf05_dp10_ld001_PROW953 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld001_PROW953_plot <- maf05_dp10_ld001_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW953") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp10_ld001_PROW953_plot
# saved as 5x5 portrait









# make the side plot, if desired.
maf05_dp10_ld001_PROW953_side_plot <- maf05_dp10_ld001_PROW953_list %>%
    select(num_mites) %>%
    unnest(cols = num_mites) %>%
    count(num_mites) %>%
    mutate(num_mites = fct_reorder(as.factor(num_mites), n)) %>%
    ggplot(aes(y = n, x = num_mites)) +
    geom_col() +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    xlab("") + ylab("") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


cowplot::plot_grid(
    cowplot::plot_grid(NULL, maf05_dp10_ld001_PROW953_side_plot + theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1, rel_heights = c(2.5, .5)),
    maf05_dp10_ld001_PROW953_plot, nrow = 1, rel_widths = c(1, 2.5)
)









# Identify which SNP sites are the ones that are shared by all three groups
maf05_dp10_ld001_PROW953_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld001_PROW953_list$num_mites)
print(maf05_dp10_ld001_PROW953_list$nodexsite[maf05_dp10_ld001_PROW953_sharedsnps])
maf05_dp10_ld001_PROW953_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW953_list$nodexsite[maf05_dp10_ld001_PROW953_sharedsnps])
names(maf05_dp10_ld001_PROW953_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW953_sharedsnps_df

View(maf05_dp10_ld001_PROW953_sharedsnps_df)

















#### ++ PROW954 ####
maf05_dp10_ld001_PROW954<-subset(maf05_dp10_ld001, bird_id =="PROW954")
levels(maf05_dp10_ld001_PROW954$bird_id)
maf05_dp10_ld001_PROW954$bird_id <- factor(maf05_dp10_ld001_PROW954$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW954$bird_id)

maf05_dp10_ld001_PROW954_list <-
    maf05_dp10_ld001_PROW954 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld001_PROW954_plot <- maf05_dp10_ld001_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp10_ld001_PROW954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld001_PROW954_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld001_PROW954_list$num_mites)
print(maf05_dp10_ld001_PROW954_list$nodexsite[maf05_dp10_ld001_PROW954_sharedsnps])
maf05_dp10_ld001_PROW954_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW954_list$nodexsite[maf05_dp10_ld001_PROW954_sharedsnps])
names(maf05_dp10_ld001_PROW954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW954_sharedsnps_df

View(maf05_dp10_ld001_PROW954_sharedsnps_df)













#### ++ PROW981 ####
maf05_dp10_ld001_PROW981<-subset(maf05_dp10_ld001, bird_id =="PROW981")
levels(maf05_dp10_ld001_PROW981$bird_id)
maf05_dp10_ld001_PROW981$bird_id <- factor(maf05_dp10_ld001_PROW981$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW981$bird_id)

maf05_dp10_ld001_PROW981_list <-
    maf05_dp10_ld001_PROW981 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld001_PROW981_plot <- maf05_dp10_ld001_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW981") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp10_ld001_PROW981_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld001_PROW981_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld001_PROW981_list$num_mites)
print(maf05_dp10_ld001_PROW981_list$nodexsite[maf05_dp10_ld001_PROW981_sharedsnps])
maf05_dp10_ld001_PROW981_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW981_list$nodexsite[maf05_dp10_ld001_PROW981_sharedsnps])
names(maf05_dp10_ld001_PROW981_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW981_sharedsnps_df

View(maf05_dp10_ld001_PROW981_sharedsnps_df)












#### ++ PROW984 ####
maf05_dp10_ld001_PROW984<-subset(maf05_dp10_ld001, bird_id =="PROW984")
levels(maf05_dp10_ld001_PROW984$bird_id)
maf05_dp10_ld001_PROW984$bird_id <- factor(maf05_dp10_ld001_PROW984$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW984$bird_id)

maf05_dp10_ld001_PROW984_list <-
    maf05_dp10_ld001_PROW984 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld001_PROW984_plot <- maf05_dp10_ld001_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp10_ld001_PROW984_plot
# save as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld001_PROW984_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld001_PROW984_list$num_mites)
print(maf05_dp10_ld001_PROW984_list$nodexsite[maf05_dp10_ld001_PROW984_sharedsnps])
maf05_dp10_ld001_PROW984_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW984_list$nodexsite[maf05_dp10_ld001_PROW984_sharedsnps])
names(maf05_dp10_ld001_PROW984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW984_sharedsnps_df

View(maf05_dp10_ld001_PROW984_sharedsnps_df)






#### ++++ 4 side by sides ####

maf05_dp10_ld001_PROW953_plot_v2 <- maf05_dp10_ld001_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW953_plot_v2



maf05_dp10_ld001_PROW953_plot_v2_cb <- maf05_dp10_ld001_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#90005d") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW953_plot_v2_cb



maf05_dp10_ld001_PROW953_plot_v2_viridis <- maf05_dp10_ld001_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#440154") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 80, by = 20)) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW953_plot_v2_viridis





maf05_dp10_ld001_PROW954_plot_v2 <- maf05_dp10_ld001_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW954_plot_v2



maf05_dp10_ld001_PROW954_plot_v2_cb <- maf05_dp10_ld001_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#b085b9") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW954_plot_v2_cb


maf05_dp10_ld001_PROW954_plot_v2_viridis <- maf05_dp10_ld001_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#414487") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW954_plot_v2_viridis




maf05_dp10_ld001_PROW981_plot_v2 <- maf05_dp10_ld001_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW981_plot_v2


maf05_dp10_ld001_PROW981_plot_v2_cb <- maf05_dp10_ld001_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#174417") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW981_plot_v2_cb


maf05_dp10_ld001_PROW981_plot_v2_viridis <- maf05_dp10_ld001_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#22a884") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 80, by = 20)) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW981_plot_v2_viridis




maf05_dp10_ld001_PROW984_plot_v2 <- maf05_dp10_ld001_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW984_plot_v2

maf05_dp10_ld001_PROW984_plot_v2_cb <- maf05_dp10_ld001_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#7cbf7c") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW984_plot_v2_cb


maf05_dp10_ld001_PROW984_plot_v2_viridis <- maf05_dp10_ld001_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#bddf26") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 80, by = 20)) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld001_PROW984_plot_v2_viridis




maf05_dp10_ld001_4sidebyside_v2<-
    ggarrange(
        maf05_dp10_ld001_PROW953_plot_v2, maf05_dp10_ld001_PROW954_plot_v2, maf05_dp10_ld001_PROW981_plot_v2, maf05_dp10_ld001_PROW984_plot_v2,
        labels = c("A", "B", "C", "D"),
        ncol = 2, nrow = 2)

# 10x10 landscape


maf05_dp10_ld001_4sidebyside_v2_cb<-
    ggarrange(
        maf05_dp10_ld001_PROW953_plot_v2_cb, maf05_dp10_ld001_PROW954_plot_v2_cb, maf05_dp10_ld001_PROW981_plot_v2_cb, maf05_dp10_ld001_PROW984_plot_v2_cb,
        labels = c("A", "B", "C", "D"),
        ncol = 2, nrow = 2)

# 10x10 landscape



maf05_dp10_ld001_4sidebyside_v2_viridis<-
    ggarrange(
        maf05_dp10_ld001_PROW953_plot_v2_viridis, maf05_dp10_ld001_PROW954_plot_v2_viridis, maf05_dp10_ld001_PROW981_plot_v2_viridis, maf05_dp10_ld001_PROW984_plot_v2_viridis,
        labels = c("A", "B", "C", "D"),
        font.label = list(size = 24),
        hjust = -0.1,
        ncol = 2, nrow = 2)

# 10x10 landscape





#### ++ PROW953 and 954 ####
maf05_dp10_ld001_PROW953_954<-subset(maf05_dp10_ld001, maf05_dp10_ld001$bird_id %in% c("PROW953","PROW954"))
levels(maf05_dp10_ld001_PROW953_954$bird_id)
maf05_dp10_ld001_PROW953_954$bird_id <- factor(maf05_dp10_ld001_PROW953_954$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW953_954$bird_id)

maf05_dp10_ld001_PROW953_954_list <-
    maf05_dp10_ld001_PROW953_954 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld001_PROW953_954_plot <- maf05_dp10_ld001_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    xlab("Sampe_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW953_954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld001_PROW953_954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp10_ld001_PROW953_954_sharedsnps<-grep("1:6", maf05_dp10_ld001_PROW953_954_list$sample_id)
print(maf05_dp10_ld001_PROW953_954_list$nodexsite[maf05_dp10_ld001_PROW953_954_sharedsnps])
maf05_dp10_ld001_PROW953_954_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW953_954_list$nodexsite[maf05_dp10_ld001_PROW953_954_sharedsnps])
names(maf05_dp10_ld001_PROW953_954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW953_954_sharedsnps_df

View(maf05_dp10_ld001_PROW953_954_sharedsnps_df)









maf05_dp10_ld001_PROW953_954_plot_v2 <- maf05_dp10_ld001_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    ylab("Number of SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953 and PROW954 (Illinois)") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16)) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld001_PROW953_954_plot_v2
# saved as 5x5














#### ++ PROW981 and 984 ####
maf05_dp10_ld001_PROW981_984<-subset(maf05_dp10_ld001, maf05_dp10_ld001$bird_id %in% c("PROW981","PROW984"))
levels(maf05_dp10_ld001_PROW981_984$bird_id)
maf05_dp10_ld001_PROW981_984$bird_id <- factor(maf05_dp10_ld001_PROW981_984$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld001_PROW981_984$bird_id)

maf05_dp10_ld001_PROW981_984_list <-
    maf05_dp10_ld001_PROW981_984 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld001_PROW981_984_plot <- maf05_dp10_ld001_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    xlab("Sampe_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld001_PROW981_984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))


maf05_dp10_ld001_PROW981_984_plot
# 5 x 9 landscpae



# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp10_ld001_PROW981_984_sharedsnps<-grep("7:12", maf05_dp10_ld001_PROW981_984_list$sample_id)
print(maf05_dp10_ld001_PROW981_984_list$nodexsite[maf05_dp10_ld001_PROW981_984_sharedsnps])
maf05_dp10_ld001_PROW981_984_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_PROW981_984_list$nodexsite[maf05_dp10_ld001_PROW981_984_sharedsnps])
names(maf05_dp10_ld001_PROW981_984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld001_PROW981_984_sharedsnps_df

View(maf05_dp10_ld001_PROW981_984_sharedsnps_df)







maf05_dp10_ld001_PROW981_984_plot_v2 <- maf05_dp10_ld001_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    ylab("Number of SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW981 and PROW984 (Louisiana)") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16)) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))


maf05_dp10_ld001_PROW981_984_plot_v2
# 5 x 9 landscpae














#### ++ ALL FOUR SAMPLES ####
maf05_dp10_ld001_all_list <-
    maf05_dp10_ld001 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld001_all_plot <- maf05_dp10_ld001_all_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "sienna3") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 1000) +
    ggtitle("maf05_dp10_ld001_all_samples") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld001_all_plot
# 6 x 14 landscape







# Identify which SNP sites are the ones that are shared by all six samples
# add the \\b thing on either side of the grep to make sure you don't also grep out the "11:12" matches since it contains "1:12" in the search/grep term!
maf05_dp10_ld001_all_sharedsnps<-grep("\\b1:12\\b", maf05_dp10_ld001_all_list$sample_id)
maf05_dp10_ld001_all_sharedsnps_df<-as.data.frame(maf05_dp10_ld001_all_list$nodexsite[maf05_dp10_ld001_all_sharedsnps])
names(maf05_dp10_ld001_all_sharedsnps_df)[1] <- "sharedsnps"
maf05_dp10_ld001_all_sharedsnps_df

View(maf05_dp10_ld001_all_sharedsnps_df)

























































#### MAF05_DP10_LD99 ####
maf05_dp10_ld99 <- read.table("maf05_dp10_ld99.txt", header = TRUE)
str(maf05_dp10_ld99)
maf05_dp10_ld99$num_mites<-as.factor(maf05_dp10_ld99$num_mites)

maf05_dp10_ld99 <- maf05_dp10_ld99 %>%
    arrange(node.1) %>%
    mutate(nodexsite = factor(nodexsite, unique(nodexsite)))

str(maf05_dp10_ld99)
levels(maf05_dp10_ld99$bird_id)




#### ++ PROW953 ####
maf05_dp10_ld99_PROW953<-subset(maf05_dp10_ld99, bird_id =="PROW953")
levels(maf05_dp10_ld99_PROW953$bird_id)
maf05_dp10_ld99_PROW953$bird_id <- factor(maf05_dp10_ld99_PROW953$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW953$bird_id)

maf05_dp10_ld99_PROW953_list <-
    maf05_dp10_ld99_PROW953 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld99_PROW953_plot <- maf05_dp10_ld99_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld99_PROW953") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp10_ld99_PROW953_plot
# saved as 5x5 portrait









# make the side plot, if desired.
maf05_dp10_ld99_PROW953_side_plot <- maf05_dp10_ld99_PROW953_list %>%
    select(num_mites) %>%
    unnest(cols = num_mites) %>%
    count(num_mites) %>%
    mutate(num_mites = fct_reorder(as.factor(num_mites), n)) %>%
    ggplot(aes(y = n, x = num_mites)) +
    geom_col() +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    xlab("") + ylab("") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


cowplot::plot_grid(
    cowplot::plot_grid(NULL, maf05_dp10_ld99_PROW953_side_plot + theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1, rel_heights = c(2.5, .5)),
    maf05_dp10_ld99_PROW953_plot, nrow = 1, rel_widths = c(1, 2.5)
)









# Identify which SNP sites are the ones that are shared by all three groups
maf05_dp10_ld99_PROW953_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld99_PROW953_list$num_mites)
print(maf05_dp10_ld99_PROW953_list$nodexsite[maf05_dp10_ld99_PROW953_sharedsnps])
maf05_dp10_ld99_PROW953_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW953_list$nodexsite[maf05_dp10_ld99_PROW953_sharedsnps])
names(maf05_dp10_ld99_PROW953_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW953_sharedsnps_df

View(maf05_dp10_ld99_PROW953_sharedsnps_df)

















#### ++ PROW954 ####
maf05_dp10_ld99_PROW954<-subset(maf05_dp10_ld99, bird_id =="PROW954")
levels(maf05_dp10_ld99_PROW954$bird_id)
maf05_dp10_ld99_PROW954$bird_id <- factor(maf05_dp10_ld99_PROW954$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW954$bird_id)

maf05_dp10_ld99_PROW954_list <-
    maf05_dp10_ld99_PROW954 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld99_PROW954_plot <- maf05_dp10_ld99_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld99_PROW954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp10_ld99_PROW954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld99_PROW954_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld99_PROW954_list$num_mites)
print(maf05_dp10_ld99_PROW954_list$nodexsite[maf05_dp10_ld99_PROW954_sharedsnps])
maf05_dp10_ld99_PROW954_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW954_list$nodexsite[maf05_dp10_ld99_PROW954_sharedsnps])
names(maf05_dp10_ld99_PROW954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW954_sharedsnps_df

View(maf05_dp10_ld99_PROW954_sharedsnps_df)













#### ++ PROW981 ####
maf05_dp10_ld99_PROW981<-subset(maf05_dp10_ld99, bird_id =="PROW981")
levels(maf05_dp10_ld99_PROW981$bird_id)
maf05_dp10_ld99_PROW981$bird_id <- factor(maf05_dp10_ld99_PROW981$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW981$bird_id)

maf05_dp10_ld99_PROW981_list <-
    maf05_dp10_ld99_PROW981 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld99_PROW981_plot <- maf05_dp10_ld99_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld99_PROW981") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp10_ld99_PROW981_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld99_PROW981_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld99_PROW981_list$num_mites)
print(maf05_dp10_ld99_PROW981_list$nodexsite[maf05_dp10_ld99_PROW981_sharedsnps])
maf05_dp10_ld99_PROW981_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW981_list$nodexsite[maf05_dp10_ld99_PROW981_sharedsnps])
names(maf05_dp10_ld99_PROW981_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW981_sharedsnps_df

View(maf05_dp10_ld99_PROW981_sharedsnps_df)












#### ++ PROW984 ####
maf05_dp10_ld99_PROW984<-subset(maf05_dp10_ld99, bird_id =="PROW984")
levels(maf05_dp10_ld99_PROW984$bird_id)
maf05_dp10_ld99_PROW984$bird_id <- factor(maf05_dp10_ld99_PROW984$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW984$bird_id)

maf05_dp10_ld99_PROW984_list <-
    maf05_dp10_ld99_PROW984 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp10_ld99_PROW984_plot <- maf05_dp10_ld99_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld99_PROW984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp10_ld99_PROW984_plot
# save as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp10_ld99_PROW984_sharedsnps<-grep(c("1, 3, 2"), maf05_dp10_ld99_PROW984_list$num_mites)
print(maf05_dp10_ld99_PROW984_list$nodexsite[maf05_dp10_ld99_PROW984_sharedsnps])
maf05_dp10_ld99_PROW984_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW984_list$nodexsite[maf05_dp10_ld99_PROW984_sharedsnps])
names(maf05_dp10_ld99_PROW984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW984_sharedsnps_df

View(maf05_dp10_ld99_PROW984_sharedsnps_df)



















#### ++++ 4 side by sides ####

maf05_dp10_ld99_PROW953_plot_v2 <- maf05_dp10_ld99_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW953_plot_v2


maf05_dp10_ld99_PROW953_plot_v2_cb <- maf05_dp10_ld99_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#90005d") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW953_plot_v2_cb




maf05_dp10_ld99_PROW953_plot_v2_viridis <- maf05_dp10_ld99_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#440154") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 200, by = 30)) +
    ggtitle("PROW953") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW953_plot_v2_viridis







maf05_dp10_ld99_PROW954_plot_v2 <- maf05_dp10_ld99_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW954_plot_v2


maf05_dp10_ld99_PROW954_plot_v2_cb <- maf05_dp10_ld99_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#b085b9") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW954_plot_v2_cb


maf05_dp10_ld99_PROW954_plot_v2_viridis <- maf05_dp10_ld99_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#414487") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW954") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW954_plot_v2_viridis






maf05_dp10_ld99_PROW981_plot_v2 <- maf05_dp10_ld99_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW981_plot_v2


maf05_dp10_ld99_PROW981_plot_v2_cb <- maf05_dp10_ld99_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#174417") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW981_plot_v2_cb




maf05_dp10_ld99_PROW981_plot_v2_viridis <- maf05_dp10_ld99_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#22a884") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 200, by = 30)) +
    ggtitle("PROW981") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW981_plot_v2_viridis




maf05_dp10_ld99_PROW984_plot_v2 <- maf05_dp10_ld99_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW984_plot_v2

maf05_dp10_ld99_PROW984_plot_v2_cb <- maf05_dp10_ld99_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#7cbf7c") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW984_plot_v2_cb



maf05_dp10_ld99_PROW984_plot_v2_viridis <- maf05_dp10_ld99_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="#bddf26") +
    ylab("Number of SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    scale_y_continuous(breaks = seq(0, 200, by = 30)) +
    ggtitle("PROW984") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16))

maf05_dp10_ld99_PROW984_plot_v2_viridis




maf05_dp10_ld99_4sidebyside_v2<-
    ggarrange(
        maf05_dp10_ld99_PROW953_plot_v2, maf05_dp10_ld99_PROW954_plot_v2, maf05_dp10_ld99_PROW981_plot_v2, maf05_dp10_ld99_PROW984_plot_v2,
        labels = c("A", "B", "C", "D"),
        ncol = 2, nrow = 2)

# 10x10 landscape




maf05_dp10_ld99_4sidebyside_v2_cb<-
    ggarrange(
        maf05_dp10_ld99_PROW953_plot_v2_cb, maf05_dp10_ld99_PROW954_plot_v2_cb, maf05_dp10_ld99_PROW981_plot_v2_cb, maf05_dp10_ld99_PROW984_plot_v2_cb,
        labels = c("A", "B", "C", "D"),
        ncol = 2, nrow = 2)

# 10x10 landscape







maf05_dp10_ld99_4sidebyside_v2_viridis<-
    ggarrange(
        maf05_dp10_ld99_PROW953_plot_v2_viridis, maf05_dp10_ld99_PROW954_plot_v2_viridis, maf05_dp10_ld99_PROW981_plot_v2_viridis, maf05_dp10_ld99_PROW984_plot_v2_viridis,
        labels = c("A", "B", "C", "D"),
        font.label = list(size = 24),
        hjust = -0.1,
        ncol = 2, nrow = 2)

# 10x10 landscape





















#### ++ PROW953 and 954 ####
maf05_dp10_ld99_PROW953_954<-subset(maf05_dp10_ld99, maf05_dp10_ld99$bird_id %in% c("PROW953","PROW954"))
levels(maf05_dp10_ld99_PROW953_954$bird_id)
maf05_dp10_ld99_PROW953_954$bird_id <- factor(maf05_dp10_ld99_PROW953_954$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW953_954$bird_id)

maf05_dp10_ld99_PROW953_954_list <-
    maf05_dp10_ld99_PROW953_954 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld99_PROW953_954_plot <- maf05_dp10_ld99_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp10_ld99_PROW953_954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld99_PROW953_954_plot
# saved as 5x7

# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp10_ld99_PROW953_954_sharedsnps<-grep("1:6", maf05_dp10_ld99_PROW953_954_list$sample_id)
print(maf05_dp10_ld99_PROW953_954_list$nodexsite[maf05_dp10_ld99_PROW953_954_sharedsnps])
maf05_dp10_ld99_PROW953_954_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW953_954_list$nodexsite[maf05_dp10_ld99_PROW953_954_sharedsnps])
names(maf05_dp10_ld99_PROW953_954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW953_954_sharedsnps_df

View(maf05_dp10_ld99_PROW953_954_sharedsnps_df)





maf05_dp10_ld99_PROW953_954_plot_v2 <- maf05_dp10_ld99_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    ylab("Number of SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("PROW953 and PROW954 (Illinois)") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16)) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld99_PROW953_954_plot_v2
# saved as 5x7






#### ++ PROW981 and 984 ####
maf05_dp10_ld99_PROW981_984<-subset(maf05_dp10_ld99, maf05_dp10_ld99$bird_id %in% c("PROW981","PROW984"))
levels(maf05_dp10_ld99_PROW981_984$bird_id)
maf05_dp10_ld99_PROW981_984$bird_id <- factor(maf05_dp10_ld99_PROW981_984$bird_id) # factor it again to remove other samples
levels(maf05_dp10_ld99_PROW981_984$bird_id)

maf05_dp10_ld99_PROW981_984_list <-
    maf05_dp10_ld99_PROW981_984 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld99_PROW981_984_plot <- maf05_dp10_ld99_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 100) +
    ggtitle("maf05_dp10_ld99_PROW981_984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))


maf05_dp10_ld99_PROW981_984_plot
# 5 x 10 landscpae



# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp10_ld99_PROW981_984_sharedsnps<-grep("7:12", maf05_dp10_ld99_PROW981_984_list$sample_id)
print(maf05_dp10_ld99_PROW981_984_list$nodexsite[maf05_dp10_ld99_PROW981_984_sharedsnps])
maf05_dp10_ld99_PROW981_984_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_PROW981_984_list$nodexsite[maf05_dp10_ld99_PROW981_984_sharedsnps])
names(maf05_dp10_ld99_PROW981_984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp10_ld99_PROW981_984_sharedsnps_df

View(maf05_dp10_ld99_PROW981_984_sharedsnps_df)






maf05_dp10_ld99_PROW981_984_plot_v2 <- maf05_dp10_ld99_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    ylab("Number of SNPs") +
    scale_x_upset(n_intersections = 100) +
    ggtitle("PROW981 and PROW984 (Louisiana)") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size=16)) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))


maf05_dp10_ld99_PROW981_984_plot_v2
# 5 x 10 landscpae











#### ++ ALL FOUR SAMPLES ####
maf05_dp10_ld99_all_list <-
    maf05_dp10_ld99 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp10_ld99_all_plot <- maf05_dp10_ld99_all_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "sienna3") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 1000) +
    ggtitle("maf05_dp10_ld99_all_samples") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp10_ld99_all_plot
# 6 x 24 landscape







# Identify which SNP sites are the ones that are shared by all six samples
# add the \\b thing on either side of the grep to make sure you don't also grep out the "11:12" matches since it contains "1:12" in the search/grep term!
maf05_dp10_ld99_all_sharedsnps<-grep("\\b1:12\\b", maf05_dp10_ld99_all_list$sample_id)
maf05_dp10_ld99_all_sharedsnps_df<-as.data.frame(maf05_dp10_ld99_all_list$nodexsite[maf05_dp10_ld99_all_sharedsnps])
names(maf05_dp10_ld99_all_sharedsnps_df)[1] <- "sharedsnps"
maf05_dp10_ld99_all_sharedsnps_df

View(maf05_dp10_ld99_all_sharedsnps_df)











































#### MAF05_DP15_LD001 ####
maf05_dp15_ld001 <- read.table("maf05_dp15_ld001.txt", header = TRUE)
str(maf05_dp15_ld001)
maf05_dp15_ld001$num_mites<-as.factor(maf05_dp15_ld001$num_mites)

maf05_dp15_ld001 <- maf05_dp15_ld001 %>%
    arrange(node.1) %>%
    mutate(nodexsite = factor(nodexsite, unique(nodexsite)))

str(maf05_dp15_ld001)
levels(maf05_dp15_ld001$bird_id)




#### ++ PROW953 ####
maf05_dp15_ld001_PROW953<-subset(maf05_dp15_ld001, bird_id =="PROW953")
levels(maf05_dp15_ld001_PROW953$bird_id)
maf05_dp15_ld001_PROW953$bird_id <- factor(maf05_dp15_ld001_PROW953$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW953$bird_id)

maf05_dp15_ld001_PROW953_list <-
    maf05_dp15_ld001_PROW953 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld001_PROW953_plot <- maf05_dp15_ld001_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW953") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp15_ld001_PROW953_plot
# saved as 5x5 portrait









# make the side plot, if desired.
maf05_dp15_ld001_PROW953_side_plot <- maf05_dp15_ld001_PROW953_list %>%
    select(num_mites) %>%
    unnest(cols = num_mites) %>%
    count(num_mites) %>%
    mutate(num_mites = fct_reorder(as.factor(num_mites), n)) %>%
    ggplot(aes(y = n, x = num_mites)) +
    geom_col() +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    xlab("") + ylab("") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


cowplot::plot_grid(
    cowplot::plot_grid(NULL, maf05_dp15_ld001_PROW953_side_plot + theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1, rel_heights = c(2.5, .5)),
    maf05_dp15_ld001_PROW953_plot, nrow = 1, rel_widths = c(1, 2.5)
)









# Identify which SNP sites are the ones that are shared by all three groups
maf05_dp15_ld001_PROW953_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld001_PROW953_list$num_mites)
print(maf05_dp15_ld001_PROW953_list$nodexsite[maf05_dp15_ld001_PROW953_sharedsnps])
maf05_dp15_ld001_PROW953_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW953_list$nodexsite[maf05_dp15_ld001_PROW953_sharedsnps])
names(maf05_dp15_ld001_PROW953_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW953_sharedsnps_df

View(maf05_dp15_ld001_PROW953_sharedsnps_df)

















#### ++ PROW954 ####
maf05_dp15_ld001_PROW954<-subset(maf05_dp15_ld001, bird_id =="PROW954")
levels(maf05_dp15_ld001_PROW954$bird_id)
maf05_dp15_ld001_PROW954$bird_id <- factor(maf05_dp15_ld001_PROW954$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW954$bird_id)

maf05_dp15_ld001_PROW954_list <-
    maf05_dp15_ld001_PROW954 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld001_PROW954_plot <- maf05_dp15_ld001_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp15_ld001_PROW954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld001_PROW954_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld001_PROW954_list$num_mites)
print(maf05_dp15_ld001_PROW954_list$nodexsite[maf05_dp15_ld001_PROW954_sharedsnps])
maf05_dp15_ld001_PROW954_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW954_list$nodexsite[maf05_dp15_ld001_PROW954_sharedsnps])
names(maf05_dp15_ld001_PROW954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW954_sharedsnps_df

View(maf05_dp15_ld001_PROW954_sharedsnps_df)













#### ++ PROW981 ####
maf05_dp15_ld001_PROW981<-subset(maf05_dp15_ld001, bird_id =="PROW981")
levels(maf05_dp15_ld001_PROW981$bird_id)
maf05_dp15_ld001_PROW981$bird_id <- factor(maf05_dp15_ld001_PROW981$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW981$bird_id)

maf05_dp15_ld001_PROW981_list <-
    maf05_dp15_ld001_PROW981 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld001_PROW981_plot <- maf05_dp15_ld001_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW981") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp15_ld001_PROW981_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld001_PROW981_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld001_PROW981_list$num_mites)
print(maf05_dp15_ld001_PROW981_list$nodexsite[maf05_dp15_ld001_PROW981_sharedsnps])
maf05_dp15_ld001_PROW981_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW981_list$nodexsite[maf05_dp15_ld001_PROW981_sharedsnps])
names(maf05_dp15_ld001_PROW981_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW981_sharedsnps_df

View(maf05_dp15_ld001_PROW981_sharedsnps_df)












#### ++ PROW984 ####
maf05_dp15_ld001_PROW984<-subset(maf05_dp15_ld001, bird_id =="PROW984")
levels(maf05_dp15_ld001_PROW984$bird_id)
maf05_dp15_ld001_PROW984$bird_id <- factor(maf05_dp15_ld001_PROW984$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW984$bird_id)

maf05_dp15_ld001_PROW984_list <-
    maf05_dp15_ld001_PROW984 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld001_PROW984_plot <- maf05_dp15_ld001_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp15_ld001_PROW984_plot
# save as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld001_PROW984_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld001_PROW984_list$num_mites)
print(maf05_dp15_ld001_PROW984_list$nodexsite[maf05_dp15_ld001_PROW984_sharedsnps])
maf05_dp15_ld001_PROW984_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW984_list$nodexsite[maf05_dp15_ld001_PROW984_sharedsnps])
names(maf05_dp15_ld001_PROW984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW984_sharedsnps_df

View(maf05_dp15_ld001_PROW984_sharedsnps_df)
















#### ++ PROW953 and 954 ####
maf05_dp15_ld001_PROW953_954<-subset(maf05_dp15_ld001, maf05_dp15_ld001$bird_id %in% c("PROW953","PROW954"))
levels(maf05_dp15_ld001_PROW953_954$bird_id)
maf05_dp15_ld001_PROW953_954$bird_id <- factor(maf05_dp15_ld001_PROW953_954$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW953_954$bird_id)

maf05_dp15_ld001_PROW953_954_list <-
    maf05_dp15_ld001_PROW953_954 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld001_PROW953_954_plot <- maf05_dp15_ld001_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW953_954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-20))

maf05_dp15_ld001_PROW953_954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp15_ld001_PROW953_954_sharedsnps<-grep("1:6", maf05_dp15_ld001_PROW953_954_list$sample_id)
print(maf05_dp15_ld001_PROW953_954_list$nodexsite[maf05_dp15_ld001_PROW953_954_sharedsnps])
maf05_dp15_ld001_PROW953_954_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW953_954_list$nodexsite[maf05_dp15_ld001_PROW953_954_sharedsnps])
names(maf05_dp15_ld001_PROW953_954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW953_954_sharedsnps_df

View(maf05_dp15_ld001_PROW953_954_sharedsnps_df)








#### ++ PROW981 and 984 ####
maf05_dp15_ld001_PROW981_984<-subset(maf05_dp15_ld001, maf05_dp15_ld001$bird_id %in% c("PROW981","PROW984"))
levels(maf05_dp15_ld001_PROW981_984$bird_id)
maf05_dp15_ld001_PROW981_984$bird_id <- factor(maf05_dp15_ld001_PROW981_984$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld001_PROW981_984$bird_id)

maf05_dp15_ld001_PROW981_984_list <-
    maf05_dp15_ld001_PROW981_984 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld001_PROW981_984_plot <- maf05_dp15_ld001_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld001_PROW981_984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-20))


maf05_dp15_ld001_PROW981_984_plot
# 5 x 9 landscpae



# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp15_ld001_PROW981_984_sharedsnps<-grep("7:12", maf05_dp15_ld001_PROW981_984_list$sample_id)
print(maf05_dp15_ld001_PROW981_984_list$nodexsite[maf05_dp15_ld001_PROW981_984_sharedsnps])
maf05_dp15_ld001_PROW981_984_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_PROW981_984_list$nodexsite[maf05_dp15_ld001_PROW981_984_sharedsnps])
names(maf05_dp15_ld001_PROW981_984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld001_PROW981_984_sharedsnps_df

View(maf05_dp15_ld001_PROW981_984_sharedsnps_df)











#### ++ ALL FOUR SAMPLES ####
maf05_dp15_ld001_all_list <-
    maf05_dp15_ld001 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld001_all_plot <- maf05_dp15_ld001_all_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "sienna3") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 1000) +
    ggtitle("maf05_dp15_ld001_all_samples") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-20))

maf05_dp15_ld001_all_plot
# 6 x 14 landscape







# Identify which SNP sites are the ones that are shared by all six samples
# add the \\b thing on either side of the grep to make sure you don't also grep out the "11:12" matches since it contains "1:12" in the search/grep term!
maf05_dp15_ld001_all_sharedsnps<-grep("\\b1:12\\b", maf05_dp15_ld001_all_list$sample_id)
maf05_dp15_ld001_all_sharedsnps_df<-as.data.frame(maf05_dp15_ld001_all_list$nodexsite[maf05_dp15_ld001_all_sharedsnps])
names(maf05_dp15_ld001_all_sharedsnps_df)[1] <- "sharedsnps"
maf05_dp15_ld001_all_sharedsnps_df

View(maf05_dp15_ld001_all_sharedsnps_df)

























































#### MAF05_DP15_LD99 ####
maf05_dp15_ld99 <- read.table("maf05_dp15_ld99.txt", header = TRUE)
str(maf05_dp15_ld99)
maf05_dp15_ld99$num_mites<-as.factor(maf05_dp15_ld99$num_mites)

maf05_dp15_ld99 <- maf05_dp15_ld99 %>%
    arrange(node.1) %>%
    mutate(nodexsite = factor(nodexsite, unique(nodexsite)))

str(maf05_dp15_ld99)
levels(maf05_dp15_ld99$bird_id)




#### ++ PROW953 ####
maf05_dp15_ld99_PROW953<-subset(maf05_dp15_ld99, bird_id =="PROW953")
levels(maf05_dp15_ld99_PROW953$bird_id)
maf05_dp15_ld99_PROW953$bird_id <- factor(maf05_dp15_ld99_PROW953$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW953$bird_id)

maf05_dp15_ld99_PROW953_list <-
    maf05_dp15_ld99_PROW953 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld99_PROW953_plot <- maf05_dp15_ld99_PROW953_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="thistle4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld99_PROW953") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp15_ld99_PROW953_plot
# saved as 5x5 portrait









# make the side plot, if desired.
maf05_dp15_ld99_PROW953_side_plot <- maf05_dp15_ld99_PROW953_list %>%
    select(num_mites) %>%
    unnest(cols = num_mites) %>%
    count(num_mites) %>%
    mutate(num_mites = fct_reorder(as.factor(num_mites), n)) %>%
    ggplot(aes(y = n, x = num_mites)) +
    geom_col() +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    xlab("") + ylab("") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


cowplot::plot_grid(
    cowplot::plot_grid(NULL, maf05_dp15_ld99_PROW953_side_plot + theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1, rel_heights = c(2.5, .5)),
    maf05_dp15_ld99_PROW953_plot, nrow = 1, rel_widths = c(1, 2.5)
)









# Identify which SNP sites are the ones that are shared by all three groups
maf05_dp15_ld99_PROW953_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld99_PROW953_list$num_mites)
print(maf05_dp15_ld99_PROW953_list$nodexsite[maf05_dp15_ld99_PROW953_sharedsnps])
maf05_dp15_ld99_PROW953_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW953_list$nodexsite[maf05_dp15_ld99_PROW953_sharedsnps])
names(maf05_dp15_ld99_PROW953_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW953_sharedsnps_df

View(maf05_dp15_ld99_PROW953_sharedsnps_df)

















#### ++ PROW954 ####
maf05_dp15_ld99_PROW954<-subset(maf05_dp15_ld99, bird_id =="PROW954")
levels(maf05_dp15_ld99_PROW954$bird_id)
maf05_dp15_ld99_PROW954$bird_id <- factor(maf05_dp15_ld99_PROW954$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW954$bird_id)

maf05_dp15_ld99_PROW954_list <-
    maf05_dp15_ld99_PROW954 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld99_PROW954_plot <- maf05_dp15_ld99_PROW954_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="pink3") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld99_PROW954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp15_ld99_PROW954_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld99_PROW954_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld99_PROW954_list$num_mites)
print(maf05_dp15_ld99_PROW954_list$nodexsite[maf05_dp15_ld99_PROW954_sharedsnps])
maf05_dp15_ld99_PROW954_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW954_list$nodexsite[maf05_dp15_ld99_PROW954_sharedsnps])
names(maf05_dp15_ld99_PROW954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW954_sharedsnps_df

View(maf05_dp15_ld99_PROW954_sharedsnps_df)













#### ++ PROW981 ####
maf05_dp15_ld99_PROW981<-subset(maf05_dp15_ld99, bird_id =="PROW981")
levels(maf05_dp15_ld99_PROW981$bird_id)
maf05_dp15_ld99_PROW981$bird_id <- factor(maf05_dp15_ld99_PROW981$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW981$bird_id)

maf05_dp15_ld99_PROW981_list <-
    maf05_dp15_ld99_PROW981 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld99_PROW981_plot <- maf05_dp15_ld99_PROW981_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkslategray4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld99_PROW981") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))



maf05_dp15_ld99_PROW981_plot
# saved as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld99_PROW981_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld99_PROW981_list$num_mites)
print(maf05_dp15_ld99_PROW981_list$nodexsite[maf05_dp15_ld99_PROW981_sharedsnps])
maf05_dp15_ld99_PROW981_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW981_list$nodexsite[maf05_dp15_ld99_PROW981_sharedsnps])
names(maf05_dp15_ld99_PROW981_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW981_sharedsnps_df

View(maf05_dp15_ld99_PROW981_sharedsnps_df)












#### ++ PROW984 ####
maf05_dp15_ld99_PROW984<-subset(maf05_dp15_ld99, bird_id =="PROW984")
levels(maf05_dp15_ld99_PROW984$bird_id)
maf05_dp15_ld99_PROW984$bird_id <- factor(maf05_dp15_ld99_PROW984$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW984$bird_id)

maf05_dp15_ld99_PROW984_list <-
    maf05_dp15_ld99_PROW984 %>%
    group_by(nodexsite) %>%
    summarise(num_mites = list(num_mites))

maf05_dp15_ld99_PROW984_plot <- maf05_dp15_ld99_PROW984_list %>%
    ggplot(aes(x=num_mites)) +
    geom_bar(fill="darkseagreen4") +
    xlab("Number of mites") +
    ylab("Number of shared SNPs") +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9))+
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld99_PROW984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold"))

maf05_dp15_ld99_PROW984_plot
# save as 5x5

# Identify which SNP sites are the ones that are shared by all three groups

maf05_dp15_ld99_PROW984_sharedsnps<-grep(c("1, 3, 2"), maf05_dp15_ld99_PROW984_list$num_mites)
print(maf05_dp15_ld99_PROW984_list$nodexsite[maf05_dp15_ld99_PROW984_sharedsnps])
maf05_dp15_ld99_PROW984_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW984_list$nodexsite[maf05_dp15_ld99_PROW984_sharedsnps])
names(maf05_dp15_ld99_PROW984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW984_sharedsnps_df

View(maf05_dp15_ld99_PROW984_sharedsnps_df)
















#### ++ PROW953 and 954 ####
maf05_dp15_ld99_PROW953_954<-subset(maf05_dp15_ld99, maf05_dp15_ld99$bird_id %in% c("PROW953","PROW954"))
levels(maf05_dp15_ld99_PROW953_954$bird_id)
maf05_dp15_ld99_PROW953_954$bird_id <- factor(maf05_dp15_ld99_PROW953_954$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW953_954$bird_id)

maf05_dp15_ld99_PROW953_954_list <-
    maf05_dp15_ld99_PROW953_954 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld99_PROW953_954_plot <- maf05_dp15_ld99_PROW953_954_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "darkorchid4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 50) +
    ggtitle("maf05_dp15_ld99_PROW953_954") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp15_ld99_PROW953_954_plot
# saved as 5x7

# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp15_ld99_PROW953_954_sharedsnps<-grep("1:6", maf05_dp15_ld99_PROW953_954_list$sample_id)
print(maf05_dp15_ld99_PROW953_954_list$nodexsite[maf05_dp15_ld99_PROW953_954_sharedsnps])
maf05_dp15_ld99_PROW953_954_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW953_954_list$nodexsite[maf05_dp15_ld99_PROW953_954_sharedsnps])
names(maf05_dp15_ld99_PROW953_954_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW953_954_sharedsnps_df

View(maf05_dp15_ld99_PROW953_954_sharedsnps_df)








#### ++ PROW981 and 984 ####
maf05_dp15_ld99_PROW981_984<-subset(maf05_dp15_ld99, maf05_dp15_ld99$bird_id %in% c("PROW981","PROW984"))
levels(maf05_dp15_ld99_PROW981_984$bird_id)
maf05_dp15_ld99_PROW981_984$bird_id <- factor(maf05_dp15_ld99_PROW981_984$bird_id) # factor it again to remove other samples
levels(maf05_dp15_ld99_PROW981_984$bird_id)

maf05_dp15_ld99_PROW981_984_list <-
    maf05_dp15_ld99_PROW981_984 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld99_PROW981_984_plot <- maf05_dp15_ld99_PROW981_984_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "turquoise4") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 100) +
    ggtitle("maf05_dp15_ld99_PROW981_984") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))


maf05_dp15_ld99_PROW981_984_plot
# 5 x 10 landscpae



# Identify which SNP sites are the ones that are shared by all six samples
maf05_dp15_ld99_PROW981_984_sharedsnps<-grep("7:12", maf05_dp15_ld99_PROW981_984_list$sample_id)
print(maf05_dp15_ld99_PROW981_984_list$nodexsite[maf05_dp15_ld99_PROW981_984_sharedsnps])
maf05_dp15_ld99_PROW981_984_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_PROW981_984_list$nodexsite[maf05_dp15_ld99_PROW981_984_sharedsnps])
names(maf05_dp15_ld99_PROW981_984_sharedsnps_df)[1] <- "sharedsnps"

maf05_dp15_ld99_PROW981_984_sharedsnps_df

View(maf05_dp15_ld99_PROW981_984_sharedsnps_df)











#### ++ ALL FOUR SAMPLES ####
maf05_dp15_ld99_all_list <-
    maf05_dp15_ld99 %>%
    group_by(nodexsite) %>%
    summarise(sample_id = list(sample_id))

maf05_dp15_ld99_all_plot <- maf05_dp15_ld99_all_list %>%
    ggplot(aes(x=sample_id)) +
    geom_bar(fill = "sienna3") +
    xlab("Sample_Number of mites") +
    ylab("Number of shared SNPs") +
    scale_x_upset(n_intersections = 1000) +
    ggtitle("maf05_dp15_ld99_all_samples") +
    theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14, face = "bold")) +
    geom_text(aes(label = ..count..), stat = "count", vjust=-0.3, position = position_dodge(width=0.9)) +
    theme_combmatrix(combmatrix.label.make_space = TRUE) +
    theme(axis.title.y = element_text(vjust=-15))

maf05_dp15_ld99_all_plot
# 6 x 24 landscape







# Identify which SNP sites are the ones that are shared by all six samples
# add the \\b thing on either side of the grep to make sure you don't also grep out the "11:12" matches since it contains "1:12" in the search/grep term!
maf05_dp15_ld99_all_sharedsnps<-grep("\\b1:12\\b", maf05_dp15_ld99_all_list$sample_id)
maf05_dp15_ld99_all_sharedsnps_df<-as.data.frame(maf05_dp15_ld99_all_list$nodexsite[maf05_dp15_ld99_all_sharedsnps])
names(maf05_dp15_ld99_all_sharedsnps_df)[1] <- "sharedsnps"
maf05_dp15_ld99_all_sharedsnps_df

View(maf05_dp15_ld99_all_sharedsnps_df)
