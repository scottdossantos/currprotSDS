# analysis of yeast transcriptome data (Gierlinski et al. 2015)

# Scott Dos Santos
# Last edited: 2025-09-11

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
# gamma value of 0.5 standard deviations. Also run a sensitivity analysis on
# the scale-naive data

if(scratch == TRUE){
  
  set.seed(2023)
  yst.0.clr <- aldex.clr(reads = yst.gene, conds = yst.meta$group, mc.samples = 128,
                           denom = "all", gamma = 1e-3, verbose = TRUE)
  
  yst.0.clr.e <- aldex.effect(clr = yst.0.clr, verbose = TRUE,
                                include.sample.summary = TRUE)
  
  yst.0.clr.t <- aldex.ttest(clr = yst.0.clr, verbose = TRUE)
  
  yst.0.clr.all <- cbind(yst.0.clr.e, yst.0.clr.t)
  # save(yst.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/yeast_scale0.Rda")
  
  
  set.seed(2025)
  yst.sens <- aldex.senAnalysis(yst.0.clr, test = "t", effect = TRUE, verbose = TRUE,
                                gamma = c(1e-3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
  
  # save(yst.sens, file = "~/Documents/GitHub/currprotSDS/data/yeast_sensitivity.Rda")
  
  
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
  
  url.sens <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_sensitivity.Rda"
  load(url(url.sens))
  
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
  mutate(title = "WT vs. \u0394Snf2 yeast transcriptome: scale-naïve ALDEx2")

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
  annotate("text", x = 7.5, y = 76, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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
  mutate(title = "WT vs. \u0394Snf2 yeast transcriptome: ALDEx2 (\u03b3 = 0.5)")

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
  annotate("text", x = 7.5, y = 52, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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
  annotate("text", x = 7.5, y = 52, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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
  annotate("text", x = 1.18, y = 8, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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
  annotate("text", x = 1.39, y = 8, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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
  annotate("text", x = 1.39, y = 8, colour = "blue", label = expression(italic("Snf2")), size = 3)+
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

################################# sensitivity analysis #################################

# ALDEx2 has a built in function for assessing the effect of gamma values on the
# number of significant transcripts. ALDEx2 is run several times with a user-
# specified selection of gamma values and the effect sizes for each feature are
# plotted against the gamma values as individual lines. Lines are plotted in
# black up until the BH-corrected P-value for a given feature drops below the
# user-defined significance threshold, whereupon the line turns grey.

plotGamma(yst.sens, test = "t", thresh = 0.05, blackWhite = T)

# these data can be extracted from the sensitivity analysis output and plotted
# with ggplot2 for more control over the appearance

# obtain data in tidy format
yst.gamma <- do.call(rbind, yst.sens)
yst.gamma$gamma <- rep(c(1e-3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), each = 5891)
yst.gamma <- yst.gamma %>% 
  select(c(diff.btw, diff.win, we.eBH, effect, gamma))

# create gene label
yst.gamma$gene <- rownames(yst.gamma)
yst.gamma$gene <- gsub("gamma_0\\..*\\.", "", yst.gamma$gene)
yst.gamma$gene <- gsub("gamma_1\\.", "", yst.gamma$gene)
length(levels(factor(yst.gamma$gene))) == nrow(yst.gene) # TRUE

# want to make a data frame contaning start and end values for line segments
# showing each gene's effect sizes between two gamma values, from 1e0-3 to 0.1,
# all the way to 0.75 to 1.0 need data in a specific format for that
yst.sens.plot <- data.frame(gam.start = rep(c(1e-3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), each = 5891),
                            gam.end = rep(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), each = 5891),
                            gene = rep(unique(yst.gamma$gene), 10),
                            eff.start = yst.gamma$effect[1:(nrow(yst.gene)*10)],
                            eff.end = yst.gamma$effect[5892:(nrow(yst.gene)*11)],
                            pvl.start = yst.gamma$we.eBH[1:(nrow(yst.gene)*10)],
                            pvl.end = yst.gamma$we.eBH[5892:(nrow(yst.gene)*11)])

# set colour value of all line segments based on pvalue at starting gamma
yst.sens.plot$col <- case_when(yst.sens.plot$pvl.start >=0.05 ~ "grey", .default = "black")

# arrange by gene and gamma start values
yst.sens.plot <- yst.sens.plot %>% 
  arrange(gene, gam.start)

# filter for non-significant and significant segments
yst.sens.ns <- yst.sens.plot %>% 
  filter(col == "grey")

yst.sens.s <- yst.sens.plot %>% 
  filter(col == "black")

# add title
yst.sens.ns$title <- "WT vs. \u0394Snf2 yeast transcriptome: sensitivity analysis"

# make colour vector for manual legend
sens.cols <- c("ns" = "grey", "s" = "black")

# count number of significant genes at each gamma value and store in dataframe
sig.genes <- vector()

for(i in levels(factor(yst.gamma$gamma))){
  sig.genes[i] <- length(which(yst.gamma$gamma == i & yst.gamma$we.eBH <0.05))
}

sig.genes <- as.data.frame(sig.genes)
sig.genes$gamma <- as.numeric(rownames(sig.genes))

# make df for blue lines on x axis
lines <- data.frame(xstart = seq(0,1,0.1), xend = seq(0,1,0.1),
                    ystart = rep(-17,11), yend = rep(-15.5,11))

# plot sensitivity analysis results
# png("~/Documents/GitHub/currprotSDS/figs/yst_sensitivity.png",
#     units = "in", height = 4, width = 5, res = 600)

p1 <- ggplot(data = yst.gamma, aes(x = gamma, y = effect))+
  geom_segment(data = yst.sens.ns, linewidth = 0.1, aes(x = gam.start, xend = gam.end, y = eff.start, yend = eff.end, colour = "ns",))+
  geom_segment(data = yst.sens.s, linewidth = 0.1, aes(x = gam.start, xend = gam.end, y = eff.start, yend = eff.end, colour = "s"))+
  geom_segment(data = lines, aes(x = xend, y = ystart, yend = yend), colour = "blue")+
  geom_text(data = sig.genes, aes(x = gamma, y = -14.5, label = sig.genes), colour = "blue", size = 2.75)+
  scale_colour_manual(name = "Transcripts", values = sens.cols, labels = c("Non-significant", "Significant"), breaks = c("ns", "s"))+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "blue")+
  scale_y_continuous(limits = c(-17, 17), expand = c(0,0))+
  xlab("Gamma")+ ylab("Effect size")+
  theme_bw()+
  facet_wrap(~title)+
  theme(legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.5,"cm"), 
        legend.margin = margin(0,0,0,0,"cm"), legend.position = "top",
        legend.title = element_text(size = 0, face = "bold"), 
        legend.text = element_text(size = 8), legend.key.spacing.x = unit(0.75, "cm"),
        axis.title = element_text(size = 9), axis.text = element_text(size = 8),
        strip.text = element_text(size = 8, face = "bold"))

p1

# dev.off()

# finally, want to show the minimum log-fold change necessary for a significant
# result vs. gamma
lfc <- vector()
for(i in levels(factor(yst.gamma$gamma))){
  tmpdf <- yst.gamma[yst.gamma$gamma == i,]
  lfc[i] <- min(abs(tmpdf$diff.btw[tmpdf$we.eBH <0.05]))
}

# linear regression on gamma vs. minimum fold change
lfc <- data.frame(gamma = as.numeric(names(lfc)), fc = lfc)
lfc$title <- "Minimum log difference required for P <0.05"

# slope is 2.3 (we have seen that the slope is close to 2 for 10 different 
# datasets, over 100 iterations, including bulk RNA-seq, metatranscriptomic data
# and single-cell datasets)
lm(fc~gamma, data = lfc)

# plot data w/ regression line, plus inferred gamma value from minimum diff
# between groups of 0.5
# png("~/Documents/GitHub/currprotSDS/figs/yst_gammaMinDiff.png",
#     units = "in", height = 4, width = 5, res = 600)

p2 <- ggplot(data = lfc, aes(x = gamma, y = fc))+
  geom_point(colour = "black", fill = "royalblue", shape = 21, stroke = 0.3, size = 2.5)+
  xlab("Scale uncertainty (\u03b3)") + ylab("Minimum difference")+
  geom_smooth(method = "lm", linewidth = 0.5, colour = "red")+
  theme_bw()+
  facet_wrap(~title)+
  theme(strip.text = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 8), axis.title = element_text(size = 9))

p2

# dev.off()


# png("~/Documents/GitHub/currprotSDS/figs/yst_sensGammaMin.png",
#     units = "in", height = 4, width = 8, res = 600)

p1 | p2

# dev.off()
