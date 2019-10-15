# Monte Carlo permutation tests for beta-multifunctionality
# XJ
# 1/10/2019

rm(list = ls())

# load library
library(TeachingDemos)
library(ecodist)
library(ggplot2)
# library(tidyverse)

char2seed("beta-functionality")

# load data ---
source("./R/BEMFdiss_load_source_data.R")

# FUNCTION: varScale, variables to be scaled from 0 to 1
varScale <- function(x) {
  x = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


###########################################################
# Monte Carlo permutaion tests
# Shuffle ecosystem functions column-wise

# by nutrient network ---
# clean data
d.raw <- nut.emf[, -which(names(nut.emf) %in% c("site_code", "plot"))]
d.std <- data.frame(apply(d.raw, 2, varScale))

# calculate Euclidean distance
raw.obs <- as.vector(distance(d.raw, method = "euclidean"))
std.obs <- as.vector(distance(d.std, method = "euclidean"))

# monte carlo permutations
nSim <- 1000
raw.sims <- matrix(NA, length(raw.obs), nSim)
std.sims <- matrix(NA, length(std.obs), nSim)
for (k in seq_len(nSim)) {
  indd <- sample(1:dim(d.raw)[1], replace = FALSE)
  raw.temp <- d.raw[indd, ]
  std.temp <- d.std[indd, ]
  raw.sims[, k] <- as.vector(distance(raw.temp, method = "euclidean"))
  std.sims[, k] <- as.vector(distance(std.temp, method = "euclidean"))
}

# calculate means and sds of the simulated data
raw.mu <- data.frame(mu = apply(raw.sims, 1, mean))
raw.sd <- data.frame(sd = apply(raw.sims, 1, sd))
std.mu <- data.frame(mu = apply(std.sims, 1, mean))
std.sd <- data.frame(sd = apply(std.sims, 1, sd))

# calculate SES
SES.raw <- (raw.obs - raw.mu) / raw.sd
SES.std <- (std.obs - std.mu) / std.sd
df <- cbind(SES.raw, SES.std)
names(df) <- c("SES.raw", "SES.std")
df.nut <- df
ggplot(df.nut, aes(SES.raw, SES.std)) +
  geom_point(shape = 1) +
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())

# by Tibetan grasslands ---
# clean data
d.raw <- tibet.emf[, -which(names(tibet.emf) %in% c("site", "plot"))]
d.std <- data.frame(apply(d.raw, 2, varScale))

# calculate Euclidean distance
raw.obs <- as.vector(distance(d.raw, method = "euclidean"))
std.obs <- as.vector(distance(d.std, method = "euclidean"))

# monte carlo permutations
nSim <- 1000
raw.sims <- matrix(NA, length(raw.obs), nSim)
std.sims <- matrix(NA, length(std.obs), nSim)
for (k in seq_len(nSim)) {
  indd <- sample(1:dim(d.raw)[1], replace = FALSE)
  raw.temp <- d.raw[indd, ]
  std.temp <- d.std[indd, ]
  raw.sims[, k] <- as.vector(distance(raw.temp, method = "euclidean"))
  std.sims[, k] <- as.vector(distance(std.temp, method = "euclidean"))
}

# calculate means and sds of the simulated data
raw.mu <- data.frame(mu = apply(raw.sims, 1, mean))
raw.sd <- data.frame(sd = apply(raw.sims, 1, sd))
std.mu <- data.frame(mu = apply(std.sims, 1, mean))
std.sd <- data.frame(sd = apply(std.sims, 1, sd))

# calculate SES
SES.raw <- (raw.obs - raw.mu) / raw.sd
SES.std <- (std.obs - std.mu) / std.sd
df <- cbind(SES.raw, SES.std)
names(df) <- c("SES.raw", "SES.std")
df.tibet <- df
ggplot(df.tibet, aes(SES.raw, SES.std)) +
  geom_point(shape = 1) +
  # geom_smooth(method = "lm") +
  geom_abline(slope = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())
# # save the results
# write.csv(df.nut, "./outs/Euclidean_null_estimate_nutnet.csv")
# write.csv(df.tibet, "./outs/Euclidean_null_estimate_tibet.csv")

# # reload the data
df.nut <- read.csv("./outs/Euclidean_null_estimate_nutnet.csv")
df.tibet <- read.csv("./outs/Euclidean_null_estimate_tibet.csv")


# load beta-diversity data
nut.beta <- read.csv("./outs/nut_beta.csv")
tibet.beta <- read.csv("./outs/tibet_beta.csv")

df.nut.long <- nut.beta %>% 
  filter(q.order == 0) %>% 
  mutate(SES.raw = df.nut$SES.raw,
         SES.std = df.nut$SES.std) %>% 
  dplyr::select(plant, bac, fun, emf, SES.raw, SES.std) %>% 
  reshape2::melt(variable.name = "diversity", 
                 value.name = "beta.div",
                 id.vars = c("emf", "SES.raw", "SES.std")) %>% 
  reshape2::melt(id.vars = c("emf", "diversity", "beta.div"))

df.tibet.long <- tibet.beta %>% 
  filter(q.order == 0) %>% 
  mutate(SES.raw = df.tibet$SES.raw,
         SES.std = df.tibet$SES.std) %>% 
  dplyr::select(plant, bac, fun, emf, SES.raw, SES.std) %>% 
  reshape2::melt(variable.name = "diversity", 
                 value.name = "beta.div",
                 id.vars = c("emf", "SES.raw", "SES.std")) %>% 
  reshape2::melt(id.vars = c("emf", "diversity", "beta.div"))

df.nut.long$group <- rep(0, dim(df.nut.long)[1])
df.tibet.long$group <- rep(1, dim(df.tibet.long)[1])

# generate a dataframe for labels
df.labs <- expand.grid(diversity = c("Plant", "Bacteria", "Fungi"),
                       group = c("Nutrient Network", "Tibetan grasslands"))
df.labs$x <- 0.35
df.labs$y <- 1.85
df.labs$labs <- c("a) r = 0.36, P < 0.001",
                   "b) r = 0.06, P = 0.416",
                   "c) r = 0.05, P = 0.366",
                   "d) r = 0.18, P < 0.001",
                   "e) r = 0.29, P < 0.001",
                   "f) r = 0.15, P = 0.002")

# plot the beta-diversity and beta-EMF relationships
p.nuttibet <- rbind(df.nut.long, df.tibet.long) %>% 
  mutate(diversity = factor(diversity, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  mutate(group = factor(group, levels = c("0", "1"),
                        labels = c("Nutrient Network", "Tibetan grasslands"))) %>%
  mutate(fsig = ifelse(value < -1.96, -1, ifelse(value >= -1.96 & value <= 1.96, 0, 1))) %>% 
  mutate(fsig = factor(fsig, levels = c("-1", "0", "1"))) %>% 
  mutate(flty = factor(paste(group, diversity))) %>% 
  filter(variable == "SES.std") %>% 
  ggplot(aes(beta.div, emf)) +
  geom_point(aes(color = fsig), 
             shape = 1) +
  geom_smooth(method = "lm", aes(lty = flty),
              lwd = 0.3, se = FALSE) +
  geom_vline(xintercept = 0.5, lty = 2, color = "black", alpha = 0.5) +
  scale_linetype_manual(values = c(2, 2, 1, 1, 1, 1)) +
  geom_text(data = df.labs, aes(x = x, y = y),
             label = df.labs$labs) +
  facet_grid(group ~ diversity) +
  lims(y = c(0, 2.0)) +
  scale_color_manual(values = c("#bc80bd", "gray70", "#e41a1c")) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0), 
                     labels = c("0.0", "0.5", "1.0")) +
  labs(x = "Community beta-diversity",
       y = "Ecosystem beta-multifunctionality") +
  expand_limits(x = 0, y = 0) +
  theme_bw(base_size = 14.5) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 15),
        panel.grid = element_blank())
# print(p.nuttibet)
ggsave("./outs/NutTibet_beta_EMFdiv_std.pdf", width = 8.5, height = 4.8)

# p.nuttibet <- rbind(df.nut.long, df.tibet.long) %>% 
#   mutate(diversity = factor(diversity, levels = c("plant", "bac", "fun"),
#                             labels = c("Plant", "Bacteria", "Fungi"))) %>% 
#   mutate(group = factor(group, levels = c("0", "1"),
#                         labels = c("Nutrient Network", "Tibetan grasslands"))) %>%
#   mutate(fsig = ifelse(value < -1.96, -1, ifelse(value >= -1.96 & value <= 1.96, 0, 1))) %>% 
#   mutate(fsig = factor(fsig, levels = c("-1", "0", "1"))) %>% 
#   mutate(flty = factor(paste(group, diversity))) %>% 
#   filter(variable == "SES.raw") %>% 
#   ggplot(aes(beta.div, emf)) +
#   geom_point(aes(color = fsig, alpha = 0.5), 
#              shape = 1) +
#   geom_smooth(method = "lm", aes(lty = flty),
#               lwd = 0.3, se = FALSE) +
#   scale_linetype_manual(values = c(2, 2, 1, 1, 1, 1)) +
#   facet_grid(group ~ diversity) +
#   scale_color_manual(values = c("gray70", "#e41a1c")) +
#   labs(x = "Community beta-diversity",
#        y = "Ecosystem beta-multifunctionality") +
#   expand_limits(x = 0, y = 0) +
#   theme_bw(base_size = 14.5) +
#   theme(legend.position = 'none',
#         strip.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 30))
# # print(p.nuttibet)
# ggsave("./outs/NutTibet_beta_EMFdiv_raw.pdf", width = 9, height = 6)

###########################################################
# Are those pairs that are more similar than expected closer to one another
# geographically? Envrionmentally?

df.nut.long2 <- nut.beta %>% 
  filter(q.order == 0) %>% 
  mutate(SES.std = df.nut$SES.std) %>% 
  dplyr::select(plant, bac, fun, emf, geo, env, SES.std) %>% 
  reshape2::melt(variable.name = "geo.env", 
                 value.name = "GEvalue",
                 id.vars = c("plant", "bac", "fun", "emf", "SES.std"))

df.tibet.long2 <- tibet.beta %>% 
  filter(q.order == 0) %>% 
  mutate(SES.std = df.tibet$SES.std) %>% 
  dplyr::select(plant, bac, fun, emf, geo, env, SES.std) %>% 
  reshape2::melt(variable.name = "geo.env", 
                 value.name = "GEvalue",
                 id.vars = c("plant", "bac", "fun", "emf", "SES.std"))

df.nut.long2$group <- rep(0, dim(df.nut.long2)[1])
df.tibet.long2$group <- rep(1, dim(df.tibet.long2)[1])

df.labs <- expand.grid(geo.env = c("Geographic distance", "Environmental distance"),
                       group = c("Nutrient Network", "Tibetan grasslands"))
df.labs$x <- c(17907 * 0.32, 1.50 * 0.3, 1250 * 0.32, 1.50 *0.3)
df.labs$y <- 1.85
df.labs$labs <- c("a) r = 0.21, P = 0.002",
                  "b) r = 0.18, P = 0.007",
                  "c) r = 0.16, P < 0.001",
                  "d) r = 0.23, P < 0.001")

rbind(df.nut.long2, df.tibet.long2) %>% 
  mutate(group = factor(group, levels = c("0", "1"),
                        labels = c("Nutrient Network", "Tibetan grasslands"))) %>%
  mutate(geo.env = factor(geo.env, levels = c("geo", "env"),
                          labels = c("Geographic distance", "Environmental distance"))) %>% 
  mutate(fsig = ifelse(SES.std < -1.96, -1, 
                       ifelse(SES.std >= -1.96 & SES.std <= 1.96, 0, 1))) %>% 
  mutate(fsig = factor(fsig, levels = c("-1", "0", "1"))) %>% 
  ggplot(aes(GEvalue, emf)) +
  geom_point(aes(color = fsig), 
             shape = 1) +
  geom_smooth(method = "lm",
              lwd = 0.3, se = FALSE) +
  facet_wrap(group ~ geo.env, scales = "free_x") +
  lims(y = c(0, 2)) +
  geom_text(data = df.labs, aes(x = x, y = y),
            label = df.labs$labs, size = 5.0) +
  scale_color_manual(values = c("#bc80bd", "gray70", "#e41a1c")) +
  labs(x = "",
       y = "Ecosystem beta-multifunctionality") +
  expand_limits(x = 0, y = 0) +
  theme_bw(base_size = 14.5) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16))
ggsave("./outs/NutTibet_beta_EMF_geoenv.pdf", width = 7.6, height = 7.6)
