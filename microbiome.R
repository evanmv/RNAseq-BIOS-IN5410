##Set directory and load packages -----
setwd("/projects/ec34/biosin5410/microbiome")
getwd()
library(tidyverse)
install.packages()
library(RColorBrewer)

##Read in data -----
data <- read.csv("data/town_data.csv")

townA <- data$townA
townB <- data$townB

##Sort and summarize -----
countsA <- table(townA)
countsB <- table(townB)

countsA <- sort(countsA, decreasing = T)
countsB <- sort(countsB, decreasing = T)

countsA.df <- as_tibble(countsA) %>%
  `colnames<-`(c("Breeds", "n"))

countsB.df <- as_tibble(countsB) %>%
  `colnames<-`(c("Breeds", "n"))

barplot(countsA.df$n, names.arg = countsA.df$Breeds, col = colors(distinct=T), args.legend = ncol(6) )
?barplot
barplot(countsB)

## Alpha-diversity -----

# Observed species
obs_sp <- function(counts) {
  obs_ind <- length(counts)
  return(obs_ind)
}

# Shannon diversity index
shannon <- function(counts) {
  total_count <- sum(counts)
  prop <- counts / total_count
  shan_ind <- -sum(prop * log(prop))
  return(shan_ind)
}

# Inversed Simpson index
simpson <- function(counts) {
  total_count <- sum(counts)
  prop <- counts / total_count
  simp_ind <- 1 / sum(prop^2)
  return(simp_ind)
}

#Run functions on Counts A and B
obs_indA <- obs_sp(countsA) # Output is number of observed species. 49 for each group.
shan_indA <- shannon(countsA) # Output is diversity based on Shannon index, equal in both groups. 
simp_indA <- simpson(countsA) # Output is diversity based on Inversed Simpson index, higher in group A. 

obs_indB <- obs_sp(countsB)
shan_indB <- shannon(countsB)
simp_indB <- simpson(countsB)

#Simpson index varies more because it puts heavier weighting on dominant breeds, whereas Shannons treats rare and abundant breeds equally. 
#Shannon uses the log of proportion value and Simpson squares the proportion value
#Diversity is higher in townA

## Beta-diversity -----

#Jaccard only focuses on how many breeds from A are also present in B
jaccard <- function(townA, townB) {
  intersection <- length(intersect(townA, townB))
  combined_community <- unique(c(townA, townB))
  jac <- 1 - intersection / length(combined_community)
  return(jac)
}

#Bray-Curtis also takes into account the numbers of each breed
bray_curtis <- function(countsA, countsB) {
  countsA <- data.frame(countsA)
  countsB <- data.frame(countsB)
  counts <- merge(countsA, countsB, by.x = 'townA', by.y='townB', all = TRUE)
  colnames(counts) <- c('Dog breed','townA','townB')
  counts$Diff <- abs(counts$townA - counts$townB)
  counts$Total <- counts$townA + counts$townB
  bray <- sum(counts$Diff) / sum(counts$Total)
  
  return(bray)
}

jac_ind <- jaccard(townA, townB)
bray_ind <- bray_curtis(countsA, countsB)

#Idicates that all of the same breeds are present in both towns (jaccard = 0), but there are differences in quantities (bray-curtis = 0.282)

## Rarefaction curves ----- 

subsample_sets <- function(sample, step, lbl) {
  div_est <- data.frame(Depth = integer(),
                        Observed_species = numeric(),
                        Shannon = numeric(),
                        Simpson = numeric(),
                        stringsAsFactors = FALSE)
  
  depths <- seq(step, length(sample), by = step)
  
  for (dep in depths) {
    subsample <- sample(sample, dep)
    sub_counts <- table(subsample)
    
    div_est <- rbind(div_est, c(dep, obs_sp(sub_counts), shannon(sub_counts), simpson(sub_counts)))
  }
  
  colnames(div_est) <- c('Depth', 'Observed_species', 'Shannon', 'Simpson')
  
#Plot rarefaction curves
  par(mfrow=c(3, 1))
  plot(div_est$Depth, div_est$Observed_species, type='o', col='lightcoral', xlab='Depth', ylab='Observed species')
  title(paste0('Observed species; ', lbl))
  
  plot(div_est$Depth, div_est$Shannon, type='o', col='lightcoral', xlab='Depth', ylab='Shannon index')
  title(paste0('Shannon index; ', lbl))
  
  plot(div_est$Depth, div_est$Simpson, type='o', col='lightcoral', xlab='Depth', ylab='Inversed Simpson index')
  title(paste0('Inversed Simpson index; ', lbl))
  
  return(div_est)
}

#Run subsample on Towns A and B
step <- 10 # Define minimal number of reads (steps)
div_estA <- subsample_sets(townA, step, 'townA') #Curves look good, reaching plateau at subsample size of < 200 (observed species)

div_estB <- subsample_sets(townB, step, 'townB') #Town B reaches plateau at about the same time (observed species) ~ 200.

# Not a big difference, but B has a large number of individuals within a few breeds, which would increase the required sampling depth

#The Shannon index curves are identical. But the Inversed Simpson curve reaches the peak very early in TownB (~ 100), though it isn't a smooth curve. This is due to the heavy weighting on overrepresented breeds.

