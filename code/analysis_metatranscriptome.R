# analysis of vaginal metatranscriptome data (Dos Santos et al. 2024)

# Scott Dos Santos
# Last edited: 2025-08-26

#################################### setup ####################################

# # install packages if needed
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ALDEx2")
# 
# install.packages('devtools')
# devtools::install_github('scottdossantos/CoDaSeq/CoDaSeq')
# 
# install.packages("dplyr")
#
# install.packages("ggplot2")
#
# install.packages("ggnewscale")

library(ALDEx2)
library(CoDaSeq)
library(dplyr)
library(ggplot2)
library(ggnewscale)

# read in feature table, metadata, KO term -> pathway lookup table and virgo
# gene -> KO term lookup table from this study's repository
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_counts.txt"
mts.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_metadata.txt"
mts.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)

url.func <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_KOlookup.txt"
mts.func <- read.table(url.func, sep = "\t", header = T, quote = "", row.names = 1)

url.virg <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/virgo_kegg_ortholog.txt"
mts.virg <- read.table(url.virg, sep = "\t", header = F, quote = "", row.names = 1)
# NOTE: this is identical to the file '8.A.kegg.ortholog.txt' provided in the
#       VIRGO github repo

# indicate whether to generate files from scratch or load .Rda from github
# (default, set to FALSE)
scratch <- FALSE

############################## aggregate KO terms ##############################

# currently, feature table contains counts for 20,109 genes across 42 samples
# but for functional analyses, we want to aggregate by KEGG orthology (KO) term

# make vector containing KO assignments of each gene in the feature table
# (find index of gene IDs in the gene -> KO lookup provided by virgo)
ko.virgo <- mts.virg[[1]][match(rownames(mts.gene), rownames(mts.virg))]
names(ko.virgo) <- rownames(mts.gene)

# aggregate by counts and clean up dataframe
mts.ko <- aggregate(mts.gene, by = list(ko.virgo), FUN = sum)
rownames(mts.ko) <- mts.ko$Group.1
mts.ko$Group.1 <- NULL

# remove suspect eukaryotic KO terms (identified in our previous analyses)
eukaryotic <- which(grepl(paste("K03260","K06237","K00599", sep = "|"), rownames(mts.ko)))
mts.ko <- mts.ko[-eukaryotic,]

# 1,658 KO terms across 42 samples in same order as lookup
all(rownames(mts.ko) == rownames(mts.func))

######################### demonstrating scale problem #########################

# problem with vaginal metatranscriptome data: marked asymmetry in data in terms
# of both gene content and scale. Causes issues when normalising as described in
# article introduction and illustrated below

# run ALDEx2 using a scale value of virtually zero (or load in combined object
# from repository if not running)
if(scratch == TRUE){
  
  set.seed(2025)
  scale.0.clr <- aldex.clr(reads = mts.ko, conds = mts.meta$group, mc.samples = 128,
                           denom = "all", gamma = 1e-3, verbose = TRUE)
  
  scale.0.clr.e <- aldex.effect(clr = scale.0.clr, verbose = TRUE,
                                include.sample.summary = TRUE)
  
  scale.0.clr.t <- aldex.ttest(clr = scale.0.clr, verbose = TRUE)
  
  scale.0.clr.all <- cbind(scale.0.clr.e, scale.0.clr.t)
  # save(scale.0.clr.all, file = "~/Documents/GitHub/currprotSDS/data/metatranscriptome_scale0.Rda")
  
} else{
  
  url.scale0 <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_scale0.Rda"
  load(url(url.scale0))
  
}
  
# pull the indices of housekeeping functions in the KO -> pathway lookup (using
# ribosomal, glycolysis and tRNA biosynthesis genes as HK functions as we expect
# them to be invariant between groups)
hk.func <- which(mts.func$pathway %in% c("Ribosome","Glycolysis / Gluconeogenesis","Aminoacyl-tRNA biosynthesis"))

# subset CLR data by significance and housekeeping functions
scale0.s <- scale.0.clr.all %>% 
  filter(we.eBH <0.05)

scale0.ns <- scale.0.clr.all %>% 
  filter(we.eBH >= 0.05) %>% 
  mutate(title = "Metatranscriptome: ALDEx2 dispersion vs. difference (\u03b3 = 0)")

scale0.hk <- scale.0.clr.all[hk.func,]

# colours for manual legend
leg.col <- c("ns" = "grey50", "s" = "red3", "hk" = "gold2")
leg.shape <- c("ns" = 19, "s" = 21, "hk" = 21)
leg.alpha <- c("ns" = 0.3, "s" = 1, "hk" = 1)
leg.size <- c("ns" = 0.9, "s" = 1.25, "hk" = 1.25)

# plot scale = 0
# png("~/Documents/GitHub/currprotSDS/figs/mts_effect_scale0.png",
#     units = "in", height = 3, width = 5, res = 600)

ggplot(data = scale0.ns, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(fill = "ns", shape = "ns", alpha = "ns", size = "ns"))+
  geom_point(data = scale0.s, aes(fill = "s", shape = "s", alpha = "s", size = "s"), stroke = 0.15)+
  geom_point(data = scale0.hk, aes(fill = "hk", shape = "hk", alpha = "hk", size = "hk"), stroke = 0.15)+
  scale_fill_manual(name = "Functions", values = leg.col, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_shape_manual(name = "Functions", values = leg.shape, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_alpha_manual(name = "Functions", values = leg.alpha, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_size_manual(name = "Functions", values = leg.size, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_x_continuous(limits = c(0,9.5), expand = c(0.00001,0.01))+
  scale_y_continuous(limits = c(-9.25,12.5), expand = c(0.00001,0.01))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  theme(legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.1,"cm"), 
        legend.margin = margin(0,0,0,0,"cm"), legend.position = "top",
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))

# dev.off()

# housekeeping functions are very off-centre (blue dashed line) and many are 
# called (erroneously!) as significantly different. We can address this problem
# using a scale uncertainty model in ALDEx2

############################# scale model: simple #############################

# we will illustrate how using simple and full scale models affect the above
# plot in this metatranscriptome dataset

# simple scale model: adds 0.5 standard deviation to the scale model across all
# samples
if(scratch == TRUE){
  
  set.seed(2025)
  scale.5.clr <- aldex.clr(reads = mts.ko, conds = mts.meta$group, mc.samples = 128,
                           denom = "all", gamma = 0.5, verbose = TRUE)
  
  scale.5.clr.e <- aldex.effect(clr = scale.5.clr, verbose = TRUE,
                                include.sample.summary = TRUE)
  
  scale.5.clr.t <- aldex.ttest(clr = scale.5.clr, verbose = TRUE)
  
  scale.5.clr.all <- cbind(scale.5.clr.e, scale.5.clr.t)
  # save(scale.5.clr.all, file = "~/Documents/GitHub/currprotSDS/data/metatranscriptome_scale5.Rda")
  
} else{
  
  url.scale5 <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_scale5.Rda"
  load(url(url.scale5))
  
}

# subset CLR data by significance and housekeeping functions
scale5.s <- scale.5.clr.all %>% 
  filter(we.eBH <0.05)

scale5.ns <- scale.5.clr.all %>% 
  filter(we.eBH >= 0.05) %>% 
  mutate(title = "Metatranscriptome: ALDEx2 dispersion vs. difference (\u03b3 = 0.5)")

scale5.hk <- scale.5.clr.all[hk.func,]

# plot scale = 0.5
# png("~/Documents/GitHub/currprotSDS/figs/mts_effect_scale5.png",
#     units = "in", height = 3, width = 5, res = 600)

ggplot(data = scale5.ns, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(fill = "ns", shape = "ns", alpha = "ns", size = "ns"))+
  geom_point(data = scale5.s, aes(fill = "s", shape = "s", alpha = "s", size = "s"), stroke = 0.15)+
  geom_point(data = scale5.hk, aes(fill = "hk", shape = "hk", alpha = "hk", size = "hk"), stroke = 0.15)+
  scale_fill_manual(name = "Functions", values = leg.col, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_shape_manual(name = "Functions", values = leg.shape, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_alpha_manual(name = "Functions", values = leg.alpha, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_size_manual(name = "Functions", values = leg.size, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_x_continuous(limits = c(0,9.5), expand = c(0.00001,0.01))+
  scale_y_continuous(limits = c(-9.25,12.5), expand = c(0.00001,0.01))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  theme(legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.1,"cm"), 
        legend.margin = margin(0,0,0,0,"cm"), legend.position = "top",
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))

# dev.off()

# this looks largely the same as the previous model (with virtually no scale
# uncertainty). Perhaps a few features along the line of equivalence between
# diff within and diff between are no longer 'significant', but the problematic
# asymmetry is still present

############################## scale model: full ##############################

# we can specify a difference in the amount of scale uncertainty per-group in
# ALDEx2. While the simple scale model increased the dispersion and pushed a
# few HK functions towards the centre, it did not address the issue at large;
# therefore, we have to use a full scale model to really rectify the issue

# full scale model: add a 15% difference in scale between groups
if(scratch == TRUE){
  
  # make matrix of scale values to pass to ALDEx2; absolute values given
  # to mu parameter do not matter, but ratio of the values DOES (i.e. c(1,5) is
  # would give the same result on the MW plot as c(5,25), using the same seed)
  set.seed(2025)
  scale.f <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15),
                                   conditions = mts.meta$group, mc.samples = 128)
  
  scale.f.clr <- aldex.clr(reads = mts.ko, conds = mts.meta$group, mc.samples = 128,
                           denom = "all", gamma = scale.f, verbose = TRUE)
  
  scale.f.clr.e <- aldex.effect(clr = scale.f.clr, verbose = TRUE,
                                include.sample.summary = TRUE)
  
  scale.f.clr.t <- aldex.ttest(clr = scale.f.clr, verbose = TRUE)
  
  scale.f.clr.all <- cbind(scale.f.clr.e, scale.f.clr.t)
  save(scale.f.clr.all, file = "~/Documents/GitHub/currprotSDS/data/metatranscriptome_scaleFull.Rda")
  
} else{
  
  url.scaleFull <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/metatranscriptome_scaleFull.Rda"
  load(url(url.scaleFull))
  
}

# subset CLR data by significance and housekeeping functions
scalef.s <- scale.f.clr.all %>% 
  filter(we.eBH <0.05)

scalef.ns <- scale.f.clr.all %>% 
  filter(we.eBH >= 0.05) %>% 
  mutate(title = "Metatranscriptome: ALDEx2 dispersion vs. difference (\u03b3 = 0.5, \u03bc = 15 %)")

scalef.hk <- scale.f.clr.all[hk.func,]

# plot scale = 0.5, difference = 15 %
# png("~/Documents/GitHub/currprotSDS/figs/mts_effect_scaleFull.png",
#     units = "in", height = 3, width = 5, res = 600)

ggplot(data = scalef.ns, aes(x = diff.win, y = diff.btw))+
  geom_point(aes(fill = "ns", shape = "ns", alpha = "ns", size = "ns"))+
  geom_point(data = scalef.s, aes(fill = "s", shape = "s", alpha = "s", size = "s"), stroke = 0.15)+
  geom_point(data = scalef.hk, aes(fill = "hk", shape = "hk", alpha = "hk", size = "hk"), stroke = 0.15)+
  scale_fill_manual(name = "Functions", values = leg.col, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_shape_manual(name = "Functions", values = leg.shape, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_alpha_manual(name = "Functions", values = leg.alpha, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  scale_size_manual(name = "Functions", values = leg.size, labels = c("Non-significant","Significant","Housekeeping"), breaks = c("ns","s", "hk"))+
  geom_abline(intercept = 0, slope = 1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = -1, colour = "grey30",linewidth = 0.5, linetype = 2)+
  geom_abline(intercept = 0, slope = 0, colour = "blue",linewidth = 0.5, linetype = 2)+
  scale_x_continuous(limits = c(0,9.5), expand = c(0.00001,0.01))+
  scale_y_continuous(limits = c(-9.25,12.5), expand = c(0.00001,0.01))+
  labs(x = "Log difference within groups", y = "Log difference between groups")+
  theme_bw()+
  facet_wrap(~title)+
  theme(legend.box.spacing = unit(0,"cm"), legend.key.width = unit(0.1,"cm"), 
        legend.margin = margin(0,0,0,0,"cm"), legend.position = "top",
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))

# dev.off()

# adding a scale difference of just 15% between groups fixes the asymmetry 
# entirely! HK genes are now centred around the line of no-difference. As a 
# result, many functions previously deemed to be significantly upregulated in
# the healthy group are no longer so; likewise, many genes previously deemed to
# be not significantly different are now significantly more expressed in the BV
# group. For datasets like this, where there is a marked disparity in both gene
# content and scale between conditions, a full-scale model is appropriate

############################## visualisation: PCA ##############################

# multiple different dimensionality reduction methods are in widespread use for 
# visualising beta-diversity among microbiomes. Singular value decomposition of
# is implemented in R via the 'prcomp()' function. Using the scaled log-ratio
# count data as input allows us to create a plot of the major principal
# components (i.e. a compositional biplot) capturing the most variation in the
# dataset

# the log-ratio data are found in the output of 'aldex.effect()' but are only
# produced when the 'include.sample.summary' argument is set to 'TRUE'. Values
# for each feature in a given sample are found in 'rab.sample.' columns. Extract
# these data
scale.f.lr <- scale.f.clr.all[, grep("rab.sample.", colnames(scale.f.clr.all))]
colnames(scale.f.lr) <- gsub("rab\\.sample\\.", "", colnames(scale.f.lr))

# principal component analysis (transposed input as samples expected in rows)
pca.f <- prcomp(t(scale.f.lr))

# basic biplot can be achieved with biplot()
# png("~/Documents/GitHub/currprotSDS/figs/mts_biplot_basic.png",
#     units = "in", height = 6, width = 6, res = 600)

biplot(pca.f, cex = c(0.55, 0.25), col = c("blue", rgb(0,0,0,0.1)))

# dev.off()

# samples separate out very well based on healthy vs. BV. Can use the CoDaSeq
# package to make a prettier plot; codaSeq.PCAplot() requires lists contaning
# row/column indices of features/samples (respectively) in the data that you
# want to show and/or colour
ind.grp <- list(Healthy = which(mts.meta$group == "Healthy"),
                BV = which(mts.meta$group == "BV"))

# highlight several pathways from our previous analysis, as well as HK functions
pathways <- c("Aminoacyl-tRNA biosynthesis","Bacterial chemotaxis","Butanoate metabolism",
              "CAMP resistance","Exopolysaccharide biosynthesis","Flagellar assembly",
              "Porphyrin metabolism","Ribosome","Starch and sucrose metabolism",
              "Two-component system")

ind.load <- list()
for(i in pathways){
  ind.load[[i]] <- which(mts.func$pathway == i)
}

# add all other pathways as 'Other'
ind.load[["Other"]] <- setdiff(1:nrow(scale.f.lr), unlist(ind.load))

# colour vector for groups and pathways
cols.group <- c("dodgerblue2", "orange2")
cols.load <- c("skyblue1","red3","gold2","chocolate4","greenyellow","purple3",
               "olivedrab","black","maroon","cyan3", rgb(0,0,0,0.075))

# plot nicer looking biplot with CoDaSeq
# png("~/Documents/GitHub/currprotSDS/figs/mts_biplot_CoDaSeq.png",
#     units = "in", height = 6, width = 10, res = 600)

codaSeq.PCAplot(pca.f, plot.groups = T, plot.loadings = T, plot.density = "groups",
                PC = c(1,2), grp = ind.grp, grp.col = cols.group, grp.cex = 0.75,
                load.grp = ind.load, load.col = cols.load, load.sym = 19, load.cex = 0.6,
                plot.legend = "loadings", leg.position = "bottomright", leg.columns = 2, leg.cex = 0.625,
                title = "Vaginal metatranscriptome: PCA plot (\u03b3 = 0.5, \u03bc = 15 %)")

# dev.off()


# alternatively, one can extract the raw plotting data for the feature loadings 
# using codaSeq.PCAvalues(), and extract the raw sample loadings from the prcomp
# object, the then use your own desired method for plotting
load.f.samp <- data.frame(pca.f$x)
load.f.feat <- codaSeq.PCAvalues(pca.f)

# add column for setting pathway
for(i in 1:length(ind.load)){
  load.f.feat[ind.load[[i]], "path"] <- names(ind.load)[i]
}

load.f.feat$col <- case_when(load.f.feat$path == "Other" ~ rgb(0,0,0,0.05), .default = "black")

# reorder levels of pathway factor to ensure 'Other' is plotted last & colours
# correspond with correct pathway
lvs <- levels(factor(load.f.feat$path))
load.f.feat$path <- factor(load.f.feat$path, levels = c(lvs[1:6], lvs[8:11], lvs[7]))

# add columns for setting group colours and strip title
for(i in 1:length(ind.grp)){
  load.f.samp[ind.grp[[i]], "group"] <- names(ind.grp)[i]
}

load.f.samp$title <- "Vaginal metatranscriptome: PCA plot (\u03b3 = 0.5, \u03bc = 15 %)"

# make vector for ensuring order of grouping colours and changing legend symbol
grp_labels<- c(Healthy = "h", BV = "v")

# plot with ggplot 
# png("~/Documents/GitHub/currprotSDS/figs/mts_biplot_ggplot.png",
#     units = "in", height = 6, width = 10, res = 600)

ggplot(data = load.f.samp, aes(x = PC1, y = PC2))+
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5, linetype = 2)+
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5, linetype = 2)+
  geom_text(aes(colour = group), label = rownames(load.f.samp))+
  scale_colour_manual(name = "Group", values = cols.group, breaks = names(grp_labels), 
                      guide = guide_legend(override.aes = list(label = grp_labels)))+
  new_scale_colour()+
  geom_point(data = load.f.feat, aes(fill = path, colour = col), stroke = 0.25, shape = 21, size = 2)+
  scale_fill_manual(name = "Pathway", values = cols.load,)+
  scale_colour_manual(name = "Pathway", values = c(rgb(0,0,0,0.05), "black"))+
  xlab("PC1: 49.4 % variance explained")+
  ylab("PC2: 11.0 % variance explained")+
  guides(colour = "none")+
  theme_bw()+
  facet_wrap(~title)+
  theme(strip.text = element_text(size = 9, face = "bold"),
        panel.grid.major = element_blank())

# dev.off()
