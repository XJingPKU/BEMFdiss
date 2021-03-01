###########################################################
# beta-diversity and beta-multifunctionality
# 
# Data are from Nutrient Network
# 
# First created by Xin Jing 9/5/2016
# Last modified by Xin Jing 8/28/2020
#  
# Contact: Xin Jing <jingxin0123@gmail.com>
###########################################################

rm(list = ls())

# load libraries
library(tidyverse)
library(TeachingDemos)

char2seed("betadiversity")

# load data
nut.meta <- read.csv("./data/data_processing/NutNet_metadata_controls_updated_noNN10.csv")
nut.plant <- data.table::fread("./data/data_processing/NutNet_species_cover_data.csv", na.strings = 'NULL')
nut.bac <- read.csv("./data/data_processing/NutNet_16S_otu_table_wTax.csv")
nut.fun <- read.csv("./data/data_processing/NutNet_ITS_otu_table_unite75_wTax.csv")
nut.env <- data.table::fread("./data/data_processing/NutNet_plot_summary_data.csv", na.strings = "NULL")
nut.NP <- read.csv("./data/data_processing/Xin-NutNet-bulk-nutrients-biomass-28Apr2017.csv")


###########################################################
# Part 1: data cleaning 
###########################################################

# create a unique id
nut.meta$id <- paste(nut.meta$site_code, nut.meta$plot, sep = "_")
nut.meta.id <- data.frame(id = nut.meta$id, sampleID = nut.meta$sampleID)
nut.env <- data.frame(nut.env)
nut.env$env.id <- paste(nut.env$site_code, nut.env$plot, sep = "_")
# plant
nut.plant <- data.table::dcast(nut.plant, site_code + plot ~ Taxon, value.var = 'max_cover')
nut.plant <- data.frame(nut.plant)
nut.plant[is.na(nut.plant)] <- 0
nut.plant$id <- paste(nut.plant$site_code, nut.plant$plot, sep = "_")
nut.plant <- nut.plant[nut.plant$id %in% nut.env$env.id, ]
# bacteria 
rownames(nut.bac) <- nut.bac$OUT_ID
nut.bac <- nut.bac[names(nut.bac) != "OUT_ID" & names(nut.bac) != "taxonomy"]
nut.bac <- t(nut.bac)
nut.bac <- data.frame(nut.bac)
nut.bac$sampleID <- rownames(nut.bac)
nut.bac <- merge(nut.bac, nut.meta.id, by = "sampleID")
nut.bac <- nut.bac[nut.bac$id %in% nut.env$env.id, ]
# fungi
rownames(nut.fun) <- nut.fun$OUT_ID
nut.fun <- nut.fun[names(nut.fun) != "OUT_ID" & names(nut.fun) != "taxonomy"]
nut.fun <- t(nut.fun)
nut.fun <- data.frame(nut.fun)
nut.fun$sampleID <- rownames(nut.fun)
nut.fun <- merge(nut.fun, nut.meta.id, by = "sampleID")
nut.fun <- nut.fun[nut.fun$id %in% nut.env$env.id, ]
# aboveground N and P data
nut.NP <- nut.NP %>% 
  filter(trt == "Control" & live == 1) %>% 
  mutate(id = paste(site_code, plot, sep = "_")) %>% 
  mutate(id = factor(id)) %>% 
  dplyr::select(id, mass, pct_N, pct_P) %>% 
  mutate(pct_N = as.numeric(as.character(pct_N)),
         pct_P = as.numeric(as.character(pct_P))) %>%  
  group_by(id) %>% 
  mutate(PTN = sum(pct_N * mass) / sum(mass, na.rm = T),
         PTP = sum(pct_P * mass) / sum(mass, na.rm = T)) %>% 
  filter(!duplicated(id)) %>% 
  dplyr::select(id, PTN, PTP)
nut.NP <- as.data.frame(nut.NP)  # 61 obs. of 3 variables
# re-clean all the data
# nut.NP table has small sample size
setdiff(nut.NP$id, nut.env$env.id)
setdiff(nut.env$env.id, nut.NP$id)
id2 <- nut.NP$id[which(nut.NP$id %in% nut.env$env.id)]  # 52 samples
# fungal species table has the lowest samples
id2 <- nut.fun$id[which(nut.fun$id %in% id2)]  # 47 samples
# env
nut.env <- nut.env[nut.env$env.id %in% id2, ]
# plant
nut.plant <- nut.plant[nut.plant$id %in% id2, ]
rownames(nut.plant) <- nut.plant$id
nut.plant <- nut.plant %>% 
  dplyr::select(-id, -site_code, -plot)
nut.plant <- nut.plant[rowSums(nut.plant) > 0, ]
nut.plant <- nut.plant[, colSums(nut.plant) > 0]
names(nut.plant) <- sub("\\.", "_", names(nut.plant))  # replace the first . to _, use sub function
names(nut.plant) <- tolower(names(nut.plant))  # to the lower case
names(nut.plant) <- gsub(pattern = "\\b([a-z])", replacement = "\\U\\1", 
                         names(nut.plant), perl = TRUE)  # capitalize the first letter
# bacteria
nut.bac <- nut.bac[nut.bac$id %in% id2, ]
rownames(nut.bac) <- nut.bac$id
nut.bac <- nut.bac[names(nut.bac) != "id" & names(nut.bac) != "sampleID"]
nut.bac <- nut.bac[rowSums(nut.bac) > 0, ]
nut.bac <- nut.bac[, colSums(nut.bac) > 0]
# fungi
nut.fun <- nut.fun[nut.fun$id %in% id2, ]
rownames(nut.fun) <- nut.fun$id
nut.fun <- nut.fun[names(nut.fun) != "id" & names(nut.fun) != "sampleID"]
nut.fun <- nut.fun[rowSums(nut.fun) > 0, ]
nut.fun <- nut.fun[, colSums(nut.fun) > 0]
# nut.NP
nut.NP <- nut.NP[nut.NP$id %in% id2, ]
nut.env <- merge(nut.env, nut.NP, by.x = "env.id", by.y = "id")
nut.env$id <- nut.env$env.id
# save data for later use
nut.emf <- nut.env %>% 
  rename(AGB = aboveground_live_mass, PTN = PTN, PTP = PTP,
         BGB = belowground_mass, STN = pct_N, avaP = ppm_P) %>% 
  dplyr::select(id, AGB, PTN, PTP,
                BGB, STN, avaP)
nut.abio <- nut.env %>% 
  dplyr::select(id, latitude, longitude, MAP, MAT, pH)

# sort samples to be the same order
nut.bac <- nut.bac[sort(rownames(nut.bac)), ]
nut.fun <- nut.fun[sort(rownames(nut.fun)), ]
nut.plant <- nut.plant[sort(rownames(nut.plant)), ]

# inspect the order of samples
nut.abio$id
nut.emf$id
rownames(nut.bac)
rownames(nut.fun)
rownames(nut.plant)

# remove unused data and package
rm(list = c("id2", "nut.meta", "nut.meta.id", "nut.NP", "nut.env"))

###########################################################
# Part 2: generate data for later use ---
###########################################################
library(fields)  # Calls: rdist.earth
library(ecodist)  # Calls: distance
source("./R/beta.div.comp.R")

# 2.1 calculate distance matrices
# varScale function
# range-relevant standardization
varScale <- function(x) {
  x = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  return(x)
}
# Geographic distance
geo <- data.frame(nut.abio$longitude, nut.abio$latitude)
nut.geo <- as.dist(rdist.earth(geo, geo, miles = FALSE))
# Environmental distance
env <- apply(nut.abio[4:6], 2, varScale)
env <- data.frame(env)
nut.clim <- distance(env[1:2], method = "euclidean")
nut.pH <- distance(env[3], method = "euclidean")
nut.env <- distance(env, method = "euclidean")
# Ecosystem functions
nut.emf$BGB[is.na(nut.emf$BGB)] <- 545 # missing value was replaced by the mean 
emf <- apply(nut.emf[-1], 2, varScale)
emf <- data.frame(emf)
nut.emf <- distance(emf, method = "euclidean")
AGB <- distance(emf["AGB"], method = "euclidean")
PTN <- distance(emf["PTN"], method = "euclidean")
PTP <- distance(emf["PTP"], method = "euclidean")
BGB <- distance(emf["BGB"], method = "euclidean")
STN <- distance(emf["STN"], method = "euclidean")
avaP <- distance(emf["avaP"], method = "euclidean")

# Community dissimilarity
# calculate dissimilarity summary table using vegetarian package
plant.q1 <- 1 - vegetarian::sim.table(nut.plant, q = 1, half = FALSE)  # Horn
plant.q2 <- 1 - vegetarian::sim.table(nut.plant, q = 2, half = FALSE)  # Morisita-Horn
bac.q1 <- 1 - vegetarian::sim.table(nut.bac, q = 1, half = FALSE)  # Horn
bac.q2 <- 1 - vegetarian::sim.table(nut.bac, q = 2, half = FALSE)  # Morisita-Horn
fun.q1 <- 1 - vegetarian::sim.table(nut.fun, q = 1, half = FALSE)  # Horn
fun.q2 <- 1 - vegetarian::sim.table(nut.fun, q = 2, half = FALSE)  # Morisita-Horn
# calculate dissimilarity summary table using beta.div.comp function
nut.plant <- ifelse(nut.plant > 0, 1, 0)
nut.plant <- data.frame(nut.plant)
nut.plant.beta <- beta.div.comp(nut.plant, coef = "S", quant = FALSE)
nut.bac <- ifelse(nut.bac > 0, 1, 0)
nut.bac <- data.frame(nut.bac)
nut.bac.beta <- beta.div.comp(nut.bac, coef = "S", quant = FALSE)
nut.fun <- ifelse(nut.fun > 0, 1, 0)
nut.fun <- data.frame(nut.fun)
nut.fun.beta <- beta.div.comp(nut.fun, coef = "S", quant = FALSE)
# Combine data by columns
df.nut <- cbind(data.frame(cbind(nut.geo, nut.env, nut.clim, nut.pH,
                                 nut.emf, AGB, PTN, PTP, BGB, STN, avaP)),
                cbind(nut.plant.beta$repl, nut.plant.beta$rich, nut.plant.beta$D),
                cbind(nut.bac.beta$repl, nut.bac.beta$rich, nut.bac.beta$D),
                cbind(nut.fun.beta$repl, nut.fun.beta$rich, nut.fun.beta$D))
names(df.nut) <- c("geo", "env", "clim", "pH",
                   "EMF", "AGB", "PTN", "PTP", "BGB", "STN", "avaP",
                   "plant.repl", "plant.rich", "plant.sor",
                   "bac.repl", "bac.rich", "bac.sor",
                   "fun.repl", "fun.rich", "fun.sor")
m <- nut.plant.beta$repl
m <- as.matrix(m)
m2 <- data.frame(row = rownames(m)[row(m)[lower.tri(m)]],
                 col = colnames(m)[col(m)[lower.tri(m)]])
df.nut <- cbind(m2, df.nut)
df.nut$plant.q1 <- plant.q1[lower.tri(plant.q1)]
df.nut$plant.q2 <- plant.q2[lower.tri(plant.q2)]
df.nut$bac.q1 <- bac.q1[lower.tri(bac.q1)]
df.nut$bac.q2 <- bac.q2[lower.tri(bac.q2)]
df.nut$fun.q1 <- fun.q1[lower.tri(fun.q1)]
df.nut$fun.q2 <- fun.q2[lower.tri(fun.q2)]

# 2.2 MC permutation tests for ecosystem functions
nSim <- 1000
obs <- as.vector(distance(emf, method = "euclidean"))
sims <- matrix(NA, length(obs), nSim)
for (k in seq_len(nSim)) {
  indd <- sample(1:dim(emf)[1], replace = FALSE)
  temp <- emf[indd, ]
  sims[, k] <- as.vector(distance(temp, method = "euclidean"))
}

# calculate means and sds of the simulated data
mu <- data.frame(mu = apply(sims, 1, mean))
sd <- data.frame(sd = apply(sims, 1, sd))

# calculate standardized effect size (SES)
SES <- (obs - mu) / sd
names(SES) <- "SES"
SES$SES.sig <- ifelse(SES$SES > 1.96, "pos", ifelse(SES$SES < -1.96, "neg", "neutral"))
df.nut <- cbind(df.nut, SES)

# 2.3 sensitive analysis
nut.vars4func <- names(emf)
nut.res <- NULL
for (i in seq_along(nut.vars4func)) {
  N <- dim(combn(nut.vars4func, i))[2]
  for (j in 1:N) {
    dat <- emf[, combn(nut.vars4func, i)[, j]]
    emf.dist <- distance(dat, method = "euclidean")
    df.dist <- data.frame(sor = df.nut$plant.sor, emf = as.vector(emf.dist))
    res <- ecodist::mantel(df.dist$emf ~ df.dist$sor, nperm = 1, mrank = TRUE)
    res <- data.frame(t(res))
    print(paste(i, "functions out of", length(nut.vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    nut.res <- rbind(nut.res, res)
  }
}

nut.res %>% 
  data.frame() %>% 
  ggplot(aes(n.funcs, mantelr)) +
  geom_point(position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  labs(x = "Number of functions",
       y = "Mantel coefficient (r)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("./outs/NutNet_no_of_functions.pdf")

write.csv(df.nut, "./data/data_processing/nutnet_distance_matrix.csv")


###########################################################
#                End of the Script                        #
###########################################################
