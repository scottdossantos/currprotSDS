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
  # actual read counts is drawn from a normal distribution)
  set.seed(2025)
  imm.thin <- thin_2group(mat = as.matrix(imm.gene), prop_null = 0.95,
                          signal_fun = stats::rnorm, alpha = 0,
                          signal_params = list(mean = 0, sd = 2))
  
  # extract group permutation and altered count table
  imm.meta$group.thin <- as.vector(imm.thin$designmat)
  
  imm.gene.thin <- imm.thin$mat
  
  # ALDEx2: scale naïve
  set.seed(2025)
  imm.0.clr <- aldex.clr(reads = imm.gene.thin, conds = imm.meta$group.thin, mc.samples = 128,
                         denom = "all", gamma = 1e-3, verbose = TRUE)
  
  imm.0.clr.e <- aldex.effect(clr = imm.0.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.0.clr.t <- aldex.ttest(clr = imm.0.clr, verbose = TRUE)
  
  imm.0.clr.all <- cbind(imm.0.clr.e, imm.0.clr.t)
  save(imm.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale0.Rda")
  
  # ALDEx2: gamma = 0.1
  set.seed(2025)
  imm.1.clr <- aldex.clr(reads = imm.gene.thin, conds = imm.meta$group.thin, mc.samples = 128,
                         denom = "all", gamma = 0.1, verbose = TRUE)
  
  imm.1.clr.e <- aldex.effect(clr = imm.1.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.1.clr.t <- aldex.ttest(clr = imm.1.clr, verbose = TRUE)
  
  imm.1.clr.all <- cbind(imm.1.clr.e, imm.1.clr.t)
  save(imm.1.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale1.Rda")
  
  # ALDEx2: gamma = 0.2
  set.seed(2025)
  imm.2.clr <- aldex.clr(reads = imm.gene.thin, conds = imm.meta$group.thin, mc.samples = 128,
                         denom = "all", gamma = 0.2, verbose = TRUE)
  
  imm.2.clr.e <- aldex.effect(clr = imm.2.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  imm.2.clr.t <- aldex.ttest(clr = imm.2.clr, verbose = TRUE)
  
  imm.2.clr.all <- cbind(imm.2.clr.e, imm.2.clr.t)
  save(imm.2.clr.all, file = "~/Documents/GitHub/currprotSDS/data/immunotherapy_scale2.Rda")
  
  # ALDEx2: gamma = 0.5
  set.seed(2025)
  imm.5.clr <- aldex.clr(reads = imm.gene.thin, conds = imm.meta$group.thin, mc.samples = 128,
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
genes.ns <- imm.0.clr.all %>% 
  filter(we.eBH >=0.05) %>% 
  mutate(scale = "ns")

genes.s.0 <- imm.0.clr.all %>% 
  filter(we.eBH <0.05,
         imm.1.clr.all$we.eBH >=0.05,
         imm.2.clr.all$we.eBH >=0.05,
         imm.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s0")

# identify genes significant only when scale is 0.1, 0.2 or 0.5
genes.s.1 <- imm.0.clr.all %>% 
  filter(imm.1.clr.all$we.eBH <0.05,
         imm.2.clr.all$we.eBH >=0.05,
         imm.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s1")

genes.s.2 <- imm.0.clr.all %>% 
  filter(imm.2.clr.all$we.eBH <0.05,
         imm.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s2")

genes.s.5 <- imm.0.clr.all %>% 
  filter(imm.5.clr.all$we.eBH <0.05) %>% 
  mutate(scale = "s5")

# bind rows together for plotting, then calculate log10 Qval
genes.plot <- rbind(genes.ns, genes.s.0, genes.s.1, genes.s.2, genes.s.5) %>% 
  mutate(title = "Nivolumab exposed vs. unexposed cells: ALDEx2", qval = -log10(we.eBH))

# some q values are either Inf or incredibly small (>265), so change these to
# an appropriate 'max' value to display on plot (next largest value is 32, so
# change these to 35)
genes.plot$qval <- case_when(genes.plot$qval == Inf ~ 70, .default = genes.plot$qval)

# colours/shapes/alpha vals/point sizes for manual legend
leg.col <- c("ns" = "grey50", "s0" = "gold2", "s1" = "cyan2", "s2" = "red3", "s5" = "magenta2")
leg.shape <- c("ns" = 19, "s0" = 21, "s1" = 21, "s2" = 21, "s5" = 21)
leg.alpha <- c("ns" = 0.15, "s0" = 1, "s1" = 1, "s2" = 1, "s5" = 1)
leg.size <- c("ns" = 1, "s0" = 1.5, "s1" = 1.5, "s2" = 1.5, "s5" = 1.5)

# sort plotting data by diff.btw and restrict to qval <20 to isolate desired red
# point on graph (18th point from right)
highlight <- genes.plot %>% 
  arrange(desc(diff.btw)) %>% 
  filter(qval <20) %>% 
  select(c(diff.btw,qval))

highlight[18,] # coordinates are x = 3.510897 / y = 2.416733

# make volcano plot
# png("~/Documents/GitHub/currprotSDS/figs/immuno_volcano.png",
#     units = "in", height = 4, width = 4.5, res = 600)

ggplot(data = genes.plot, aes(x = diff.btw, y = qval))+
  geom_point(aes(fill = scale, size = scale, shape = scale, alpha = scale), colour = "black", stroke = 0.25)+
  annotate("point", x = 3.510897, y = 2.416733, shape = 21, fill = rgb(0,0,0,0), size = 3, colour = "blue")+
  scale_fill_manual(name = "Genes", values = leg.col, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not significant", "\u03b3 = 0", "\u03b3 = 0.1", "\u03b3 = 0.2", "\u03b3 = 0.5"),)+
  scale_shape_manual(name = "Genes", values = leg.shape, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                     labels = c("Not significant", "\u03b3 = 0", "\u03b3 = 0.1", "\u03b3 = 0.2", "\u03b3 = 0.5"),)+
  scale_size_manual(name = "Genes", values = leg.size, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not significant", "\u03b3 = 0", "\u03b3 = 0.1", "\u03b3 = 0.2", "\u03b3 = 0.5"),)+
  scale_alpha_manual(name = "Genes", values = leg.alpha, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                     labels = c("Not significant", "\u03b3 = 0", "\u03b3 = 0.1", "\u03b3 = 0.2", "\u03b3 = 0.5"),)+
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = "blue")+
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "blue")+
  labs(x = expression("Log"[2]*" difference between groups"), y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  theme(legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.1,"cm"), 
        legend.margin = margin(0,0,0,0,"cm"), legend.position = "top",
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))

# dev.off()
