# analysis of single cell RNA-seq data ()

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

library(ALDEx2)
library(dplyr)
library(ggplot2)

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

######################### choosing scale values to use #########################

# having shown now that scale can be used to dramatically reduce the number of
# false-positive findings, how does one choose an appropriate scale value? By
# changing the amount of scale uncertainty incorporated into the ALDEx2 model,
# one will also change the set of differentially abundant genes returned and we
# will illustrate this below.

# # run ALDEx2 on the single-cell dataset using four different amounts of scale
# uncertainty: scale-naïve, and gamma values of 0.1, 0.2 and 0.5 (load from 
# GH repo if not generating from scratch)

if(scratch == TRUE){
  
  # scale naïve
  set.seed(2025)
  sc.0.clr <- aldex.clr(reads = sc.gene, conds = sc.meta$group, mc.samples = 128,
                         denom = "all", gamma = 1e-3, verbose = TRUE)
  
  sc.0.clr.e <- aldex.effect(clr = sc.0.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  sc.0.clr.t <- aldex.ttest(clr = sc.0.clr, verbose = TRUE)
  
  sc.0.clr.all <- cbind(sc.0.clr.e, sc.0.clr.t)
  save(sc.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/singleCell_scale0.Rda")
  
  # gamma = 0.1
  set.seed(2025)
  sc.1.clr <- aldex.clr(reads = sc.gene, conds = sc.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.1, verbose = TRUE)
  
  sc.1.clr.e <- aldex.effect(clr = sc.1.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  sc.1.clr.t <- aldex.ttest(clr = sc.1.clr, verbose = TRUE)
  
  sc.1.clr.all <- cbind(sc.1.clr.e, sc.1.clr.t)
  save(sc.1.clr.all, file = "~/Documents/GitHub/currprotSDS/data/singleCell_scale1.Rda")
  
  # gamma = 0.2
  set.seed(2025)
  sc.2.clr <- aldex.clr(reads = sc.gene, conds = sc.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.2, verbose = TRUE)
  
  sc.2.clr.e <- aldex.effect(clr = sc.2.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  sc.2.clr.t <- aldex.ttest(clr = sc.2.clr, verbose = TRUE)
  
  sc.2.clr.all <- cbind(sc.2.clr.e, sc.2.clr.t)
  save(sc.2.clr.all, file = "~/Documents/GitHub/currprotSDS/data/singleCell_scale2.Rda")
  
  # gamma = 0.5
  set.seed(2025)
  sc.5.clr <- aldex.clr(reads = sc.gene, conds = sc.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.5, verbose = TRUE)
  
  sc.5.clr.e <- aldex.effect(clr = sc.5.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  sc.5.clr.t <- aldex.ttest(clr = sc.5.clr, verbose = TRUE)
  
  sc.5.clr.all <- cbind(sc.5.clr.e, sc.5.clr.t)
  save(sc.5.clr.all, file = "~/Documents/GitHub/currprotSDS/data/singleCell_scale5.Rda")
  
} else{
  
  for(i in c(0,1,2,5)){
    load(paste0("~/Documents/GitHub/currprotSDS/data/singleCell_scale",i,".Rda")) 
  }
  
}

# identify the genes that are not significant, and genes significant only in the
# scale naïve analysis
genes.ns <- sc.0.clr.all %>% 
  filter(we.eBH >=0.05) %>% 
  mutate(scale = "ns")

genes.s.0 <- sc.0.clr.all %>% 
  filter(we.eBH <0.05,
         sc.1.clr.all$we.eBH >=0.05,
         sc.2.clr.all$we.eBH >=0.05,
         sc.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s0")

# identify genes significant only when scale is 0.1, 0.2 or 0.5
genes.s.1 <- sc.0.clr.all %>% 
  filter(sc.1.clr.all$we.eBH <0.05,
         sc.2.clr.all$we.eBH >=0.05,
         sc.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s1")

genes.s.2 <- sc.0.clr.all %>% 
  filter(sc.2.clr.all$we.eBH <0.05,
         sc.5.clr.all$we.eBH >=0.05) %>% 
  mutate(scale = "s2")

genes.s.5 <- sc.0.clr.all %>% 
  filter(sc.5.clr.all$we.eBH <0.05) %>% 
  mutate(scale = "s5")

# bind rows together for plotting, then calculate log10 Qval
genes.plot <- rbind(genes.ns, genes.s.0, genes.s.1, genes.s.2, genes.s.5) %>% 
  mutate(title = "Memory T cells vs. Cytotoxic T cells: ALDEx2", qval = -log10(we.eBH))

# convert scale column to factor and reorder levels for correct order
lvls <- levels(factor(genes.plot$scale))

genes.plot$scale <- factor(genes.plot$scale, 
                           levels = c(lvls[1], lvls[5], lvls[2], lvls[3], lvls[4]))

# some q values are either Inf or incredibly small (>265), so change these to
# an appropriate 'max' value to display on plot (next largest value is 32, so
# change these to 35)
genes.plot$qval <- case_when(genes.plot$qval == Inf ~ 35,
                             genes.plot$qval >200 ~ 35,
                             .default = genes.plot$qval)

# colours/shapes/alpha vals/point sizes for manual legend
leg.col <- c("ns" = "grey50", "s0" = "gold2", "s1" = "cyan2", "s2" = "red3", "s5" = "magenta2")
leg.shape <- c("ns" = 19, "s0" = 21, "s1" = 21, "s2" = 21, "s5" = 21)
leg.alpha <- c("ns" = 0.15, "s0" = 1, "s1" = 1, "s2" = 1, "s5" = 1)
leg.size <- c("ns" = 1, "s0" = 1.5, "s1" = 1.5, "s2" = 1.5, "s5" = 1.5)

# make volcano plot
# png("~/Documents/GitHub/currprotSDS/figs/sc_volcano.png",
#     units = "in", height = 4, width = 4.5, res = 600)

ggplot(data = genes.plot, aes(x = diff.btw, y = qval))+
  geom_point(aes(fill = scale, size = scale, shape = scale, alpha = scale), colour = "black", stroke = 0.25)+
  scale_fill_manual(name = "Genes", values = leg.col, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not signif.", "Signif. (\u03b3 = 0)", "Signif. (\u03b3 = 0.1)", "Signif. (\u03b3 = 0.2)", "Signif. (\u03b3 = 0.5)"),)+
  scale_shape_manual(name = "Genes", values = leg.shape, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not signif.", "Signif. (\u03b3 = 0)", "Signif. (\u03b3 = 0.1)", "Signif. (\u03b3 = 0.2)", "Signif. (\u03b3 = 0.5)"),)+
  scale_size_manual(name = "Genes", values = leg.size, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not signif.", "Signif. (\u03b3 = 0)", "Signif. (\u03b3 = 0.1)", "Signif. (\u03b3 = 0.2)", "Signif. (\u03b3 = 0.5)"),)+
  scale_alpha_manual(name = "Genes", values = leg.alpha, breaks = c("ns", "s0", "s1", "s2", "s5"), 
                    labels = c("Not signif.", "Signif. (\u03b3 = 0)", "Signif. (\u03b3 = 0.1)", "Signif. (\u03b3 = 0.2)", "Signif. (\u03b3 = 0.5)"),)+
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
