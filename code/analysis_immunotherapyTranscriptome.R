# analysis of anti-PD1 immunotherapy RNA-seq data (Riaz et al. 2017; also used
# in Li et al. 2022)

# Scott Dos Santos
# Last edited: 2025-09-03

#################################### setup ####################################

# # install packages if needed
#
# install.packages("seqgendiff")
#
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ALDEx2")
# 
# install.packages("dplyr")
#
# install.packages("ggplot2")

library(seqgendiff)
library(ALDEx2)
library(dplyr)
library(ggplot2)

# read in feature table and metadata for the immunotherapy transcriptome dataset
# from this study's repository
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/immuno_counts.txt"
imm.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/immuno_metadata.txt"
imm.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)

# indicate whether to generate files from scratch or load .Rda from github
# (default, set to FALSE)
scratch <- FALSE

######################### choosing scale values to use #########################

# in a study of cells exposed (or not exposed) to the PD1 cell-cycle checkpoint
# inhibitor, nivolumab, DESeq2 and edgeR return sets of differentially expressed
# genes that barely overlap, suggesting high levels of false-positives. The same
# dataset was used by Li et al. to demonstrate this very point.

# Gloor et al. used the same data to demonstrate that tools like DESeq2 have a
# very high false discovery rates. They used an approach that permuted the group
# variable and injected an artificial signal into 5% of the genes to simulate 
# 'true positive' findings in a differential expression analysis of this data.
# They showed that even a small amount of scale uncertainty can remove many of 
# these false positives, albeit at the cost of sensitivity. To see how different
# amounts of scale uncertainty can affect a differential abundance analysis, one
# can look at the relationship between the differences in abundance of features
# within and between groups, and their corresponding BH-corrected P values

# permute the immunotherapy datasetusing the same approach as Gloor et al. and
# run ALDEx2 on the simulated data using four different amounts of scale 
# uncertainty: scale-naïve, and gamma values of 0.1, 0.2 and 0.5 (load from the
# GH repo if not generating from scratch)

if(scratch == TRUE){
  
  # use binomial thinning to permute data and ensure that 5% of these genes are
  # significantly different (amount of difference between groups added to the
  # actual read counts is drawn from a log normal distribution)
  
  
  
  # scale naïve
  set.seed(2025)
  imm.0.clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                         denom = "all", gamma = 1e-3, verbose = TRUE)
  
  imm.0.clr.e <- aldex.effect(clr = imm.0.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.0.clr.t <- aldex.ttest(clr = imm.0.clr, verbose = TRUE)
  
  imm.0.clr.all <- cbind(imm.0.clr.e, imm.0.clr.t)
  save(imm.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale0.Rda")
  
  # gamma = 0.1
  set.seed(2025)
  imm.1.clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.1, verbose = TRUE)
  
  imm.1.clr.e <- aldex.effect(clr = imm.1.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.1.clr.t <- aldex.ttest(clr = imm.1.clr, verbose = TRUE)
  
  imm.1.clr.all <- cbind(imm.1.clr.e, imm.1.clr.t)
  save(imm.1.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale1.Rda")
  
  # gamma = 0.2
  set.seed(2025)
  imm.2.clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.2, verbose = TRUE)
  
  imm.2.clr.e <- aldex.effect(clr = imm.2.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.2.clr.t <- aldex.ttest(clr = imm.2.clr, verbose = TRUE)
  
  imm.2.clr.all <- cbind(imm.2.clr.e, imm.2.clr.t)
  save(imm.2.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale2.Rda")
  
  # gamma = 0.5
  set.seed(2025)
  imm.5.clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.5, verbose = TRUE)
  
  imm.5.clr.e <- aldex.effect(clr = imm.5.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.5.clr.t <- aldex.ttest(clr = imm.5.clr, verbose = TRUE)
  
  imm.5.clr.all <- cbind(imm.5.clr.e, imm.5.clr.t)
  save(imm.5.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale5.Rda")
  
} else{
  
  for(i in c(0,1,2,5)){
    load(paste0("~/Documents/GitHub/currprotSDS/data/immunotherapy_scale",i,".Rda")) 
  }
  
}

# identify the genes that are significant only in the scale naïve analysis





