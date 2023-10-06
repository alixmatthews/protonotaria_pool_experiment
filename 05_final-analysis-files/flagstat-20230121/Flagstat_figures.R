#### Number of SNPs versus LD threshold ####

setwd("~/Desktop/Mites/Genomics/Projects/20210816_projects/20210816_exp/05_Analyses/Flagstat_figures/")


#### LOAD LIBRARIES ####



library(ggplot2)
library(scales)
library(ggpubr)
library(dplyr)
library(ggthemes)


# 1. Using "total" (this is wrong) ----
## 1.1. load data
flagstat <- read.csv("./20230113/Flagstat_figures.csv", header=TRUE, sep = ",")
str(flagstat)
flagstat$Num_Mites <- as.factor(flagstat$Num_Mites)

## 1.2 Take a look at everything at once, separated by reference genome
ggplot(flagstat, aes(x = Sample_ID, y = Num_Reads, fill = Category)) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Reference) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 1.3 subset on category2 - remove total_bac

flagstat2<-subset(flagstat, Category2 !="total_bac")
str(flagstat2)
flagstat2$Category2<-as.character(flagstat2$Category2)
flagstat2<-subset(flagstat2, Category2 !="total_bac")
str(flagstat2)
flagstat2$Category2<-as.factor(flagstat2$Category2)
str(flagstat2)


flagstat3<-subset(flagstat, Reference == "Mite")
flagstat3$Reference<-as.character(flagstat3$Reference)
flagstat3<-subset(flagstat3, Reference == "Mite")
flagstat3$Reference<-as.factor(flagstat3$Reference)
str(flagstat3)




ggplot(flagstat, aes(x = Sample_ID, y = Num_Reads, fill = Category2)) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Reference) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(flagstat, aes(x = Sample, y = Num_Reads, fill = Category2)) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Take a look at mapped reads with only using total number of reads from mite genome (slightly different than when using bacteria genomes)
ggplot(flagstat2, aes(x = Sample_ID, y = Num_Reads, fill = Category2)) +
    geom_bar(stat = "identity", position="dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(flagstat2, aes(x = Sample, y = Num_Reads, fill = Category2)) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### now let's try to get the categories in a nice readable order
category_order <- c("total", "mapped", "total_bac", "mapped_bac")


ggplot(flagstat, aes(x = Sample, y = Num_Reads, fill = factor(Category2, category_order))) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Yes! Now let's make it pretty

facet_label_num_mites <- as_labeller(
    c(`1` = "1 mite", `5` = "5 mites",`20` = "20 mites"))


ggplot(flagstat, aes(x = Sample, y = Num_Reads, fill = factor(Category2, category_order))) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites, labeller = facet_label_num_mites) +
    scale_y_continuous(labels = scales::comma) +
    ylab(label = "Number of Reads") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    scale_fill_manual(name = "Category", values=c("slategray4", "slategray2", "bisque4", "cornsilk3"), labels=c('Total (Mite Reference)', 'Mapped to Mite Reference', "Total (Bacteria Reference)", "Mapped to Bacteria Reference"))

# saved pdf landscape 7x14

#### Here it is without the total bacteria category

#### now let's try to get the categories in a nice readable order
category_order_wo_total_bac <- c("total", "mapped", "mapped_bac")

ggplot(flagstat2, aes(x = Sample, y = Num_Reads, fill = factor(Category2, category_order_wo_total_bac))) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites, labeller = facet_label_num_mites) +
    scale_y_continuous(labels = scales::comma) +
    ylab(label = "Number of Reads") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    scale_fill_manual(name = "Category", values=c("slategray4", "slategray2", "cornsilk3"), labels=c('Total', 'Mapped to Mite Reference', "Mapped to Bacteria Reference"))

# saved pdf landscape 7x14





category_order_miterefonly <- c("total", "mapped")

ggplot(flagstat3, aes(x = Sample, y = Num_Reads, fill = factor(Category2, category_order_miterefonly))) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites, labeller = facet_label_num_mites) +
    scale_y_continuous(labels = scales::comma) +
    ylab(label = "Number of Reads") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    scale_fill_manual(name = "Category", values=c("slategray4", "slategray2"), labels=c('Total', 'Mapped to Mite Reference'))

# saved pdf landscape 7x14


























# 2. The right figures ----
## 2.1. Load data -----
flagstat <- read.csv("./20230121/flagstat_results_20230121.csv", header=TRUE, sep = ",")
str(flagstat)
flagstat$Num_Mites <- as.factor(flagstat$Num_Mites)

## 2.2. Set up figures -----
facet_label_num_mites <- as_labeller(
    c(`1` = "1 mite", `5` = "5 mites",`20` = "20 mites"))

category_order <- c("primary", "mapped_mite", "mapped_bac")

ggplot(flagstat, aes(x = Sample, y = Num_Reads, fill = factor(Category_specific, category_order))) +
    geom_bar(stat = "identity", position="dodge") +
    facet_grid(.~ Num_Mites, labeller = facet_label_num_mites) +
    scale_y_continuous(labels = scales::comma) +
    ylab(label = "Number of Reads") +
    theme_pubclean() +
    theme(axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 20)) +
    scale_fill_manual(name = "Category", values=c("slategray4", "slategray2", "cornsilk3"), labels=c('Total', 'Mapped (Mite)', "Mapped (Bacteria)")) +
    theme(legend.position="bottom",
          legend.text = element_text(size=20),
          legend.title=element_blank())

# saved pdf landscape 7x14
