# analysis of yeast transcriptome data (Gierlinski et al. 2015)

# Scott Dos Santos
# Last edited: 2025-08-27

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
#install.packages("patchwork")

library(ALDEx2)
library(dplyr)
library(ggplot2)
library(patchwork)

# read in feature table and metadata for the yeast transcriptome dataset from
# this study's repository
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_counts.txt"
yst.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_metadata.txt"
yst.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)

# indicate whether to generate files from scratch or load .Rda from github
# (default, set to FALSE)
scratch <- FALSE

######################### demonstrating scale problem #########################

# in a highly-replicated yeast transcriptome dataset (knockout SNF2 vs. WT),
# 65-80 % of the genes are called as differential by commonly used tools for
# differential expression analysis of RNA-seq data (e.g. DESeq2, edgeR, ALDEx2).
# The authors of this benchmarking paper suggested using a dual-cutoff approach
# to identify genes which are likely to be truly different. Adding even a modest
# amount of scale uncertainty can solve this problem without the need for any
# sort of dual cutoff approach

# run ALDEx2 with virtually no scale uncertainty included, and again with a
# gamma value of 0.5 standard deviations

if(scratch == TRUE){
  
  set.seed(2025)
  yst.0.clr <- aldex.clr(reads = yst.gene, conds = yst.meta$group, mc.samples = 128,
                           denom = "all", gamma = 1e-3, verbose = TRUE)
  
  yst.0.clr.e <- aldex.effect(clr = yst.0.clr, verbose = TRUE,
                                include.sample.summary = TRUE)
  
  yst.0.clr.t <- aldex.ttest(clr = yst.0.clr, verbose = TRUE)
  
  yst.0.clr.all <- cbind(yst.0.clr.e, yst.0.clr.t)
  save(yst.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/yeast_scale0.Rda")
  
  
  set.seed(2025)
  yst.5.clr <- aldex.clr(reads = yst.gene, conds = yst.meta$group, mc.samples = 128,
                         denom = "all", gamma = 0.5, verbose = TRUE)
  
  yst.5.clr.e <- aldex.effect(clr = yst.5.clr, verbose = TRUE,
                              include.sample.summary = TRUE)
  
  yst.5.clr.t <- aldex.ttest(clr = yst.5.clr, verbose = TRUE)
  
  yst.5.clr.all <- cbind(yst.5.clr.e, yst.5.clr.t)
  save(yst.5.clr.all, file = "~/Documents/GitHub/currprotSDS/data/yeast_scale5.Rda")
  
} else{
  
  url.scale0 <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_scale0.Rda"
  load(url(url.scale0))
  
  url.scale5 <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_scale5.Rda"
  load(url(url.scale5))
  
}

# convert BH-corrected P values into -log10()
yst.0.clr.all$qval <- -log10(yst.0.clr.all$we.eBH)
yst.5.clr.all$qval <- -log10(yst.5.clr.all$we.eBH)

# subset the output from the scale-naïve model for transcripts which are: 1) not 
# significantly different in the scale-naïve analysis; 2) significantly 
# different in the scale-naïve analysis only and; 3) significantly different
# in both analyses
subset.0.ns.0 <- yst.0.clr.all %>% 
  filter(we.eBH >=0.05) %>% 
  mutate(df = "Not significant")

subset.0.s.0 <- yst.0.clr.all %>% 
  filter(we.eBH <0.05, yst.5.clr.all$we.eBH >=0.05) %>% 
  mutate(df = "Significant (\u03b3 = 0)")

subset.0.s.5 <- yst.0.clr.all %>% 
  filter(we.eBH <0.05, yst.5.clr.all$we.eBH <0.05) %>% 
  mutate(df = "Significant (\u03b3 = 0.5)")

# isolate row for SNF2 (YOR290C; i.e. the actual gene that was knocked out; 
# see: https://www.sciencedirect.com/science/article/pii/S088875430900158X )
subset.snf2.0 <- subset.0.s.5["YOR290C",]

# bind rows for plotting
subset.all.0 <- rbind(subset.0.ns.0, subset.0.s.0, subset.0.s.5) %>% 
  mutate(title = "WT vs. \u0394SNF2 yeast transcriptome: scale-naïve ALDEx2")

# several -log10(Q) values are Inf as their corresponding p values are 0, so
# convert these to '80' so they are shown legibly on the plot
subset.all.0$qval <- case_when(subset.all.0$qval == Inf ~ 80, .default = subset.all.0$qval)
subset.snf2.0$qval <- 80

# reorder 'df' factor so 0.5 is after 0
subset.all.0$df <- factor(subset.all.0$df,
                        levels = c("Not significant","Significant (\u03b3 = 0)", "Significant (\u03b3 = 0.5)"))

# plot unscaled yeast data as volcano plot (note: annotate throws a warning but
# okay to ignore: https://stackoverflow.com/questions/77219398/ )
# png("~/Documents/GitHub/currprotSDS/figs/yst_volcano_scale0.png",
#     units = "in", height = 4, width = 5, res = 600)

p1 <- ggplot(data = subset.all.0, aes(x = diff.btw, y = qval))+
  geom_vline(xintercept = 1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_vline(xintercept = -1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_point(aes(colour = df), size = 1.5, alpha = 0.5)+
  geom_point(data = subset.snf2.0, shape = 21, colour = "blue", size = 3.5, stroke = 1)+
  annotate("text", x = 7.5, y = 76, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  scale_x_continuous(limits = c(-3.5, 8), breaks = seq(-4,8, 2), expand = c(0.01,0.01))+
  scale_y_continuous(limits = c(0,82), expand = c(0.02,0.05))+
  scale_colour_manual(name = "Transcripts", values = c("black", "gold2", "red3"))+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

p1

# dev.off()

# 4,169/5,891 genes (~71 %) are called as significantly different between WT
# and SNF2-KO cells, which originally necessitated a dual cutoff of P value
# and log2 fold change to reduce this to a more 'realistic' or manageable 
# number of genes which could then be analysed downstream. This fold-change
# threshold is completely arbitrary and adding even a small amount of scale 
# removes the need for such a threshold

####################### yeast transcriptome scale model #######################

# subset the output from the scale model for transcripts which are: 1) not 
# significantly different; 2) significantly different in the scale-naïve 
# analysis only and; 3) significantly different in both analyses
subset.5.ns.0 <- yst.5.clr.all %>% 
  filter(we.eBH >=0.05) %>% 
  mutate(df = "Not significant")

# subset.5.s.0 <- yst.5.clr.all %>% 
#   filter(we.eBH >=0.05, yst.0.clr.all$we.eBH >0.05) %>% 
#   mutate(df = "Significant (\u03b3 = 0)")

subset.5.s.5 <- yst.5.clr.all %>% 
  filter(we.eBH <0.05, yst.0.clr.all$we.eBH <0.05) %>% 
  mutate(df = "Significant (\u03b3 = 0.5)")

# isolate row for SNF2 (YOR290C; i.e. the actual gene that was knocked out)
subset.snf2.5 <- subset.5.s.5["YOR290C",]

# bind rows for plotting
subset.all.5 <- rbind(subset.5.ns.0, subset.5.s.5) %>% 
  mutate(title = "WT vs. \u0394SNF2 yeast transcriptome: ALDEx2 (\u03b3 = 0.5)")

# several -log10(Q) values are Inf as their corresponding p values are 0, so
# convert these to '80' so they are shown legibly on the plot
subset.all.5$qval <- case_when(subset.all.5$qval == Inf ~ 55, .default = subset.all.5$qval)
subset.snf2.5$qval <- 55

# reorder 'df' factor so 0.5 is after 0
subset.all.5$df <- factor(subset.all.5$df,
                          levels = c("Not significant", "Significant (\u03b3 = 0.5)"))

# plot data
# png("~/Documents/GitHub/currprotSDS/figs/yst_volcano_scale5.png",
#     units = "in", height = 4, width = 5, res = 600)

p2 <- ggplot(data = subset.all.5, aes(x = diff.btw, y = qval))+
  geom_vline(xintercept = 1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_vline(xintercept = -1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_point(aes(colour = df), size = 1.5, alpha = 0.5)+
  geom_point(data = subset.snf2.5, shape = 21, colour = "blue", size = 3.5, stroke = 1)+
  annotate("text", x = 7.5, y = 52, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  scale_x_continuous(limits = c(-3.5, 8), breaks = seq(-4,8, 2), expand = c(0.01,0.01))+
  scale_y_continuous(limits = c(0,57), expand = c(0.02,0.05))+
  scale_colour_manual(name = "Transcripts", values = c("black", "red3"))+
  labs(x = "Log difference between groups", y = expression("-Log"[10]*" adjusted P-value"))+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

p2

# dev.off()

# one can see that adding a modest amount of scale to the model reduces the
# unfeasibly large number of significantly different genes to less than 200. 
# Equally, there is no need any longer for an additional fold-change threshold
# as almost all of these 175 differentially expressed genes have a difference
# between groups greater than the previous arbitrary threshold

# plot of both side by side
p2.edit <- ggplot(data = subset.all.5, aes(x = diff.btw, y = qval))+
  geom_vline(xintercept = 1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_vline(xintercept = -1.4, linetype = 2, linewidth = 0.5, colour = "blue")+
  geom_point(aes(colour = df), size = 1.5, alpha = 0.5)+
  geom_point(data = subset.snf2.5, shape = 21, colour = "blue", size = 3.5, stroke = 1)+
  annotate("text", x = 7.5, y = 52, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  scale_x_continuous(limits = c(-3.5, 8), breaks = seq(-4,8, 2), expand = c(0.01,0.01))+
  scale_y_continuous(limits = c(0,57), expand = c(0.02,0.05))+
  scale_colour_manual(name = "Transcripts", values = c("black", "red3"))+
  labs(x = "Log difference between groups", y = "")+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

# png("~/Documents/GitHub/currprotSDS/figs/yst_volcano_scaleBoth.png",
#     units = "in", height = 4, width = 10, res = 600)

p1 | p2.edit

# dev.off()

################################# effect plots #################################

# make effect plots for the same data subsets

# scale = 0
# png("~/Documents/GitHub/currprotSDS/figs/yst_effect_scale0.png",
#     units = "in", height = 4, width = 5, res = 600)

p3 <- ggplot(data = subset.all.0, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(colour = df), alpha = 0.5, size = 1.5)+
  geom_point(data = subset.snf2.0, col = "blue", shape = 21, size = 3.5, stroke = 1)+
  annotate("text", x = 1.18, y = 8, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_colour_manual(name = "Transcripts", values = c("black", "gold2", "red3"))+
  scale_x_continuous(limits = c(0,4.5), expand = c(0.0001,0.0025))+
  scale_y_continuous(limits = c(-4,8.75), expand = c(0,0))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

p3

# dev.off()

# scale = 0.5
# png("~/Documents/GitHub/currprotSDS/figs/yst_effect_scale5.png",
#     units = "in", height = 4, width = 5, res = 600)

p4 <- ggplot(data = subset.all.5, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(colour = df), alpha = 0.5, size = 1.5)+
  geom_point(data = subset.snf2.5, col = "blue", shape = 21, size = 3.5, stroke = 1)+
  annotate("text", x = 1.39, y = 8, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_colour_manual(name = "Transcripts", values = c("black", "red3"))+
  scale_x_continuous(limits = c(0,4.5), expand = c(0.0001,0.0025))+
  scale_y_continuous(limits = c(-4,8.75), expand = c(0,0))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

p4

# dev.off()

# both side-by-side
# png("~/Documents/GitHub/currprotSDS/figs/yst_effect_scaleBoth.png",
#     units = "in", height = 4, width = 10, res = 600)

p4.edit <- ggplot(data = subset.all.5, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(colour = df), alpha = 0.5, size = 1.5)+
  geom_point(data = subset.snf2.5, col = "blue", shape = 21, size = 3.5, stroke = 1)+
  annotate("text", x = 1.39, y = 8, colour = "blue", label = expression(italic("SNF2")), size = 3)+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_colour_manual(name = "Transcripts", values = c("black", "red3"))+
  scale_x_continuous(limits = c(0,4.5), expand = c(0.0001,0.0025))+
  scale_y_continuous(limits = c(-4,8.75), expand = c(0,0))+
  labs(x = "Log difference within groups", y = "")+
  theme_bw()+
  facet_wrap(~title)+
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        legend.position = "top", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7), strip.text = element_text(size = 9, face = "bold"),
        legend.margin = margin(0,0,0,0,"cm"), legend.box.margin = margin(0,0,0,0,"cm"),
        legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.2,"cm"))

p3 | p4.edit

# dev.off()
