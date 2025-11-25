# time taken to run ALDEx2 modular functions on all 3 datasets

# Scott Dos Santos
# Last edited: 2025-09-24

#################################### setup ####################################

# # install packages if needed
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ALDEx2")

library(ALDEx2)

# /!\ only load data for one dataset at a time or else  /!\
# /!\ times will be inflated by unnecessary RAM usage   /!\

# immunotherapy
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/immuno_counts.txt"
imm.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/immuno_metadata.txt"
imm.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)


# yeast
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_counts.txt"
yst.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/yeast_metadata.txt"
yst.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)


# metatranscriptome
url.gene <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/mts_counts.txt"
mts.gene <- read.table(url.gene, sep = "\t", header = T, quote = "", row.names = 1)

url.meta <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/mts_metadata.txt"
mts.meta <- read.table(url.meta, sep = "\t", header = T, quote = "", row.names = 1)

url.func <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/mts_KOlookup.txt"
mts.func <- read.table(url.func, sep = "\t", header = T, quote = "", row.names = 1)

url.virg <- "https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/mts_virgo_KOs.txt"
mts.virg <- read.table(url.virg, sep = "\t", header = F, quote = "", row.names = 1)

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

######################## timing: immunotherapy dataset ########################

# run ALDEx2 clr, effect and ttest functions on immunotherapy dataset, checking
# the time before and after running, then calculate time in seconds. Run the
# three command a total of 5 times and store times in a dataframe

times.imm <- data.frame(iteration = 1:5, start = rep(NA, 5),
                        finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  times.imm[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE, verbose = TRUE)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  times.imm[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  times.imm[i,4] <- times.imm[i,3] - times.imm[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(times.imm, file = "~/Documents/GitHub/currprotSDS/data/times.imm.Rda")

# same but for multithreaded aldex.clr
mc.times.imm <- data.frame(iteration = 1:5, start = rep(NA, 5),
                           finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  mc.times.imm[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = imm.gene, conds = imm.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE, useMC = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE,
                         verbose = TRUE,useMC = T)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  mc.times.imm[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  mc.times.imm[i,4] <- mc.times.imm[i,3] - mc.times.imm[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(mc.times.imm, file = "~/Documents/GitHub/currprotSDS/data/mc.times.imm.Rda")

##################### timing: yeast transcriptome dataset #####################

# run ALDEx2 clr, effect and ttest functions on yeast transcriptome dataset, 
# checking the time before and after running, then calculate time in seconds.
# Run the three command a total of 5 times and store times in a dataframe
times.yst <- data.frame(iteration = 1:5, start = rep(NA, 5),
                        finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  times.yst[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = yst.gene, conds = yst.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE, verbose = TRUE)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  times.yst[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  times.yst[i,4] <- times.yst[i,3] - times.yst[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(times.yst, file = "~/Documents/GitHub/currprotSDS/data/times.yst.Rda")


# now do the same for multicore
mc.times.yst <- data.frame(iteration = 1:5, start = rep(NA, 5),
                           finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  mc.times.yst[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = yst.gene, conds = yst.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE, useMC = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE,
                         verbose = TRUE, useMC = TRUE)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  mc.times.yst[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  mc.times.yst[i,4] <- mc.times.yst[i,3] - mc.times.yst[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(mc.times.yst, file = "~/Documents/GitHub/currprotSDS/data/mc.times.yst.Rda")

###################### timing: metatranscriptome dataset ######################

# run ALDEx2 clr, effect and ttest functions on vaginal metatranscriptome
# dataset checking the time before and after running, then calculate time in 
# seconds. Run the three command a total of 5 times and store times in a 
# dataframe
times.mts <- data.frame(iteration = 1:5, start = rep(NA, 5),
                        finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  times.mts[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = mts.ko, conds = mts.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE, verbose = TRUE)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  times.mts[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  times.mts[i,4] <- times.mts[i,3] - times.mts[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(times.mts, file = "~/Documents/GitHub/currprotSDS/data/times.mts.Rda")

# do the same but for multicore processing
mc.times.mts <- data.frame(iteration = 1:5, start = rep(NA, 5),
                           finish. = rep(NA, 5), total = rep(NA, 5))

for(i in 1:5){
  
  # returns number of seconds past the Unix epoch of Jan 1st, 1970, 00:00:00 UTC
  mc.times.mts[i,2] <- Sys.time()
  
  set.seed(2025)
  clr <- aldex.clr(reads = mts.ko, conds = mts.meta$group, mc.samples = 128,
                   denom = "all", gamma = 0.5, verbose = TRUE, useMC = TRUE)
  
  effect <- aldex.effect(clr = clr, include.sample.summary = TRUE,
                         verbose = TRUE, useMC = TRUE)
  
  ttest <- aldex.ttest(clr = clr, verbose = TRUE)
  
  mc.times.mts[i,3] <- Sys.time()
  
  # calculates time taken in number of seconds
  mc.times.mts[i,4] <- mc.times.mts[i,3] - mc.times.mts[i,2]
  
  rm(clr, effect, ttest)
  
  gc()
  
}

# save(mc.times.mts, file = "~/Documents/GitHub/currprotSDS/data/mc.times.mts.Rda")

####################### timing: calculating mean & s.d.  #######################

# load all .Rda objects containing time information from GitHub
mc.t <- paste0("https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/",
               "mc.times.", c("imm", "mts", "yst"),".Rda")

mc.f <- paste0("https://github.com/scottdossantos/currprotSDS/raw/refs/heads/main/data/",
               "times.", c("imm", "mts", "yst"),".Rda")

for(i in c(mc.t, mc.f)){
  load(url(i))
  rm(i)
}

# calculate means and standard deviations across 5 instances of ALDEx2
for(i in ls(pattern = "times")){
  message(paste0("Mean time taken (", i, ") \n",
                 mean(get(i)$total), "\n",
                 "Std dev. time taken (", i, ") \n", sd(get(i)$total), "\n"))
}
