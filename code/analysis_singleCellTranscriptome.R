# analysis of anti-PD1 immunotherapy RNA-seq data ()

# Scott Dos Santos
# Last edited: 2025-09-03

#################################### setup ####################################

# # install packages if needed
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ALDEx2")
# 
# install.packages("dplyr")
#
# install.packages("ggplot2")
#
# install.packages("patchwork")

library(ALDEx2)
library(dplyr)
library(ggplot2)
library(patchwork)

# read in feature table and metadata for the single cell transcriptome dataset
# from this study's repository
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/singleCell_counts.txt"
sc.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/singleCell_metadata.txt"
sc.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)

# indicate whether to generate files from scratch or load .Rda from github
# (default, set to FALSE)
scratch <- FALSE

# Note: This single-cell dataset is a subset of the original feature table
#       from 'Zheng_cytotoxic_t.txt' and 'Zheng_memory_t.txt' used in Skinnider
#       et al. 2019 whereby only the cells in the upper quartile of read counts
#       were retained (2,551 from memory T cells, 2,544 from cytoxic T cells).
#       Likewise, only reads with >300 reads across all ~10,000 cells were 
#       retained. The final dataset was produced by randomly selecting 1,000 
#       cells from each lineage for a final dataset of 1,508 genes across
#       2,000 cells. For detailed code producing this feature table, see 
#       'https://github.com/ggloor/datasets/blob/main/make-single-cell.R'.

######################### demonstrating scale problem #########################
