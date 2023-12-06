##Load packages -----
library(tidyverse)
library(ggplot2)
library(plotly)
##Read in data -----
dataset.df <- read_csv("/Users/evanvallenari/RNAseq-BIOS-IN5410/miRNAseq/lc_dataset.csv", delim= ";")
?read_csv
dataset.df <- read.csv("/Users/evanvallenari/RNAseq-BIOS-IN5410/miRNAseq/lc_dataset.csv", sep = ";")
sigmRNA.df <- read.csv("/Users/evanvallenari/RNAseq-BIOS-IN5410/miRNAseq/sig_mRNA.df", sep = " ")
miRNAcounts.df <- read.csv("/Users/evanvallenari/RNAseq-BIOS-IN5410/miRNAseq/lc_mirna_counts.csv", sep = ";")

row.names(sigmRNA.df) <- miRNAcounts.df$ID

##Volcano plot -----

vplot <- ggplot(sigmRNA.df) +
  aes(x=sigmRNA.df$log2FoldChange, y=sigmRNA.df$log10pvalue, text = paste("Transcript ID:", row.names(sigmRNA.df))) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  labs(title="Volcano plot",
       subtitle = "miRNA expression in Cancer vs. Control") +
       xlab("Log2 Fold Change") +
       ylab("-Log10 p-value") +
  theme_bw()

ggplotly(vplot)
ggsave("vplot_gg2.png")
