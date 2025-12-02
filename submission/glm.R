library(ALDEx2)

# pull KO-aggregated count table and metadata for vaginal metatranscriptomes
# used in Dos Santos et al. (2024)
url.df <- "https://github.com/scottdossantos/dossantos2024study/raw/refs/heads/main/Rdata/ko.both.Rda"
tf.df <- tempfile(fileext = ".Rda")
download.file(url = url.df, destfile = tf.df, mode = "wb")
load(tf.df)
unlink(tf.df)

url.md <- "https://github.com/scottdossantos/dossantos2024study/raw/refs/heads/main/Rdata/hm.metadata.Rda"
tf.md <- tempfile(fileext = ".Rda")
download.file(url = url.md, destfile = tf.md, mode = "wb")
load(tf.md)
unlink(tf.md)

# delete temp files
rm(list = c(ls(pattern = "url."),ls(pattern = "tf.")))

# remove two weird samples from the read count table and delete suspect reads
# of eukaryotic origin as in original paper
ko.both <- ko.both[,-c(9,22)]

euk <- paste("K03364","K13963","K01173","K12373", "K14327","K00863","K00599",
             "K13993","K00811","K03260","K00985", sep = "|")

ko.both <- ko.both[-which(grepl(euk, rownames(ko.both))),]

# set seed and generate matrix of (reproducible) scale values as in original
# paper
set.seed(2023)
mu.hbv <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15), mc.samples = 128,
                                conditions = hm.metadata$`BV status`)

# make model matrix for a GLM including BV status (healthy vs. BV) and dataset
# (London or Europe datasets; this should have no impact on differential
# expression)
mm <- model.matrix(~ `BV status` + Dataset, hm.metadata)

# transform data w/ different amounts of scale added to each group (15 %
# difference between groups on average)
set.seed(2023)
x <- aldex.clr(reads = ko.both, conds = mm, mc.samples = 128,
               denom = "all", verbose = TRUE, gamma = mu.hbv)

# run the GLM function (will take 1-2 minutes to run)
x.glm <- aldex.glm(clr = x, verbose = T, fdr.method = "BH")

# calculate effect sizes from the GLM
x.glm.e <- aldex.glm.effect(clr = x, verbose = T, include.sample.summary = T)

# extract effect sizes and CLR abundances for both variables
effect.bv <- x.glm.e[[1]]
effect.dataset <- x.glm.e[[2]]

# analysis shows only a single significantly different gene expressed at a
# higher level in samples from the London dataset compared to Europe Dataset,
# which is expected because dataset shouldn't really influence this at all!
length(which(x.glm$`DatasetLondon:pval.padj` <0.05))

# look at this gene in particular: present in all London samples at low counts
# (~20-30 reads per sample) but absent from all Europe samples. This is very
# clearly an artifact of sequencing depth and not truly differential expression
ko.both[which(x.glm$`DatasetLondon:pval.padj` <0.05),]

# many differentially expressed genes based on BV status (note: original paper
# restricted analysis of differentially expressed genes to those with absolute
# effect size >1)
length(which(x.glm$`\`BV status\`Healthy:pval.padj` <0.05)) # 911
