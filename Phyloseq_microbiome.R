getwd()
##Load packages -----

library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)
library(dplyr)
 
##Read in files -----
getwd()
abund_table <- read.csv("data/All_Good_P2_C03.csv", row.names = 1, check.names = F)

#Transpose
abund_table <- t(abund_table)

#Read in taxonomy and metadata
OTU_taxonomy <- read.csv("data/All_Good_P2_C03_Taxonomy.csv", row.names = 1, check.names = F)
meta_table <- read.csv("data/ENV_pitlatrine.csv", row.names = 1, check.names = F)

#Set group info
grouping_info <- data.frame(row.names = rownames(meta_table), t(as.data.frame(strsplit(rownames(meta_table), "_"))))

?strsplit

colnames(grouping_info) <- c("Country", "Latrine", "Depth") # Add column names to grouping_info, and create data.frame by combining meta_table and grouping_info 
meta_table <- data.frame(meta_table, grouping_info)

#Filter - if row names (samples) are not present in meta_table
abund_table1 <- abund_table[rownames(abund_table) %in% rownames(meta_table), ]

#Read in phylogenetic tree
OTU_tree <- read_tree("data/All_Good_P2_C03.tre")

head(abund_table)
str(abund_table)
rowMeans(abund_table1) 
colMeans(abund_table1)

head(OTU_tree)

print(max(abund_table1$C0))
#C0 is a Bacteria from Firmicutes, Clostridia, Clostridiales, Clostridiaceae_1. The abundance mean is low, 0.02469136. Min is 0. Max is ?. 

##Phyloseq -----
#Create phyloseq object
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)

#Richness -- Calculate microbial diversity in the dataset
richness<-estimate_richness(physeq, split = TRUE, measures = NULL)

#Visualize and save plot as p 
p <- plot_richness(physeq, x = "Country", color = "Depth",
                   measures = c("Observed", "Shannon",
                                "InvSimpson"))
p2 <- p+theme_bw()  

#The alpha diversity as shown by the Shannon index is higher in latrine pits in Vietnam than in Tanzania.
  
#Calculate and visualize community composition
distances_bray<-distance(physeq, method="bray")
distances_unifrac<-distance(physeq, method="unifrac")
  
str(distances_bray) #Check structure of tables
str(distances_unifrac)


ord <- ordinate(physeq, method="PCoA", 
                distance = "unifrac", 
                weighted=TRUE)
p <- plot_ordination(physeq, ord, 
                     color="Country",
                     title="Phyloseq's Weighted Unifrac")
p <- p + geom_point(size=5) + theme_bw()



