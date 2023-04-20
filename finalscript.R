# This is a script to satisfy the final project requirement in ENTM 642, Analysis of ecological data
# copy/paste into a new R script
# install required packages if needed
# download the ‘fishpertrawl.csv’ data to your working directory

# clear R environment
rm(list=ls())

# set working directory
setwd("write_path_to_data_here")

# load packages
library(vegan)
library(dendextend)
library(tidyverse)
library(ape)
library(cluster)

# load data
fishpertrawl <- read_csv("fishpertrawl.csv")

# create function to use later
`%notin%` <- Negate(`%in%`)

# keep only species present in the catch for 15 or more years
species <- 
  fishpertrawl %>% 
  group_by(SYear, Species) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Species) %>%
  summarise(count = n()) %>%
  filter(count >= 15) %>%
  select(Species)

common.fish <- 
  fishpertrawl %>%
  right_join(species)

# transform for clustering
fish <- 
  common.fish %>% 
  filter(Common_name %notin% NA) %>% 
  group_by(SYear, Common_name) %>% 
  summarise(biomass = sum(BioPerTrawl)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "SYear", values_from = biomass)

# NAs are instances where no fish of that species were caught in a year, convert to 0s
fish <- 
  fish %>%
  replace(is.na(.), 0)

# select only species data for clustering
fishy <- 
  fish %>%
  as.data.frame()

community <-
  fishy[-1]

row.names(community) <-
  fishy$Common_name

#  Create bray-curtis distance association matrix
com.dist <- vegdist(community, method="bray")	

# Create cluster using Ward's D2
com.clust <- hclust(com.dist, method="ward.D2")

# Create cophenetic matrix
coph <- cophenetic(com.clust)		
# cophenetic correlation
coph.corr <- mantel(coph, com.dist)	
coph.corr
# Mantel's r: 0.6857 significance: 0.001

# Calculate where to prune tree
optclus <- sapply(2:12, function(x) summary(silhouette(cutree(com.clust, k = x), com.dist))$avg.width)
optclus # inspect results
# 3 groups

# convert to dendrogram for plotting
com.den <- as.dendrogram(com.clust)
com.den <- set(com.den, "labels_cex", 1.5)

# plot dendrogram
# tiff("dendrogram.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
par(cex=0.5, col="black", mar=c(6,2,2,10)) 					# # par changes aesthetic parameters
plot(com.den, horiz=TRUE, main="Saginaw Bay Fish Communities", xlab="Distance in Species Space", cex.lab=1.5, cex.axis=1.5, cex.main=2)
rect.dendrogram(com.den, 3, horiz=TRUE, border="red", lty=5, lwd=1)
# dev.off()

# conduct principal coordinates analysis on bray-curtis matrix
hey <- pcoa(com.dist)
# plot pcoa
# tiff("mtxpcoa.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
biplot(hey, main='PCOA of Saginaw Bay Species')
points(hey$vectors[,1],hey$vectors[,2], 
       col=ifelse(hey$vectors[,1] < 0, 'red', 
            ifelse(hey$vectors[,1] > 0 & hey$vectors[,2] > 0, "#CFB53B", "#56B4E9")), cex=2, 
       pch = ifelse(hey$vectors[,1] < 0, 16, 
                    ifelse(hey$vectors[,1] > 0 & hey$vectors[,2] > 0, 15, 17)))	
# dev.off()
# axis 1 - 25.8%, axis 2 - 14.0%, axis 3 - 10.1%

# common species ----------------------------------------------------------
forplot <- 
  common.fish %>% 
  group_by(Common_name) %>% 
  summarise(avgbiomass = sum(BioPerTrawl)/n()) %>% 
  arrange(avgbiomass)

ggplot(forplot) +
  aes(x =fct_rev(fct_reorder(Common_name, avgbiomass)), y = avgbiomass) +
  geom_col()+
  xlab("") +
  ylab(expression(paste("Average biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Species biomass per tow minute") +
  scale_y_continuous(expand =c(0,0.01)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "averagebiomass.jpg", device='jpg', dpi=700)  
