###########################################################
# On controls of Ecosystem beta-mutlifunctionality (EMF)
# 
# Leading authors: Xin Jing, Aimee Classen, Nathan Sanders
# Data are from Nutrient Network and Tibetan grasslands
# 
# First created by Xin Jing on September 5, 2016
# Last modified by Xin Jing
#  
# Contact: Xin Jing <Xin.Jing@uvm.edu>
###########################################################

rm(list = ls())

# load data ---
source("./R/BEMFdiss_load_source_data.R")

# load library ---
library(ggplot2)
library(fields)  # Calls: rdist.earth
library(ecodist)  # Calls: distance, mantel
library(TeachingDemos)  # Calls: char2seed
library(gridExtra)  # Calls: grid.arrange
library(vegetarian)  # Calls: sim.table
library(picante)
# library(adespatial)  # Calls: beta.div.comp
source("./R/beta.div.comp.R")

# read tree
nut.tree.sp <- phytools::read.newick("./data/data_processing/NutNet_SpeciesTree_Zanne_Lengths.txt")
nut.tree.sp <- ape::collapse.singles(nut.tree.sp)
# plot(nut.tree.sp, type = "fan", cex = 0.4)
tibet.tree.sp <- phytools::read.newick("./data/data_processing/TibetanPlateau_SpeciesTree_Lengths_Zanne.txt")
tibet.tree.sp <- ape::collapse.singles(tibet.tree.sp)
# plot(tibet.tree.sp, type = "fan", cex = 0.4)

char2seed("beta-multifunctionality")

###########################################################
# Map figure ---

# generate a dataframe containing coordinates
df.map <- rbind(cbind(site = as.character(nut.abio$site_code),
                      latitude = nut.abio$latitude, 
                      longitude = nut.abio$longitude,
                      projects = rep("NutNet", dim(nut.abio)[1])),
                cbind(site = as.character(tibet.abio$site),
                      latitude = tibet.abio$latitude, 
                      longitude = tibet.abio$longitude,
                      projects = rep("Tibet", dim(tibet.abio)[1]))) %>% 
  data.frame() %>% 
  mutate(projects = factor(projects),
         latitude = as.numeric(as.character(latitude)),
         longitude = as.numeric(as.character(longitude)))
summary(df.map)
# keep unique rows
df.map <- df.map %>% 
  distinct()
summary(df.map)
# ggmap
df.map %>% 
  ggplot(aes(x = longitude, y = latitude)) +
  borders(database = "world", 
          colour = "gray80", fill = NA, lwd = 0.02) +
  geom_point(aes(color = projects), size = 0.8, alpha = 0.6) +
  labs(x = expression("Longitude ("*degree*")"),
       y = expression("Latitude ("*degree*")")) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.15),
        legend.title = element_blank())
# ggsave("./outs/NutTibet_map_figure.pdf", width = 7.5, height = 4.5)

###########################################################
# inspect spatial patterns of ecosystem beta-multifunctionality ---

# FUNCTION: varScale, variables to be scaled from 0 to 1
varScale <- function(x) {
  x = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
# FUNCTION: calDist, calculate matrics for geographic distance, environmental distance and beta-EMF
calDist <- function(geo, env, emf) {
  # geo is a matrix with longitude in the first column and latitude in the second column 
  geo.dist <- as.dist(rdist.earth(geo, geo, miles = FALSE))
  # abiotic factors, e.g., climate and soil pH
  env <- data.frame(apply(env, 2, varScale))
  env.dist <- distance(env, method = "euclidean")
  # ecosystem functions
  emf <- data.frame(apply(emf, 2, varScale))
  emf.dist <- distance(emf, method = "euclidean")
  dist.df <- cbind(geo.dist, env.dist, emf.dist)
  dist.df <- data.frame(dist.df)
  return(dist.df)
}

#-----------------------------
# Relationships between geographic distance and beta-EMF, and environmental distance and beta-EMF
# by NutNet
nut.emf.sp <- calDist(geo = data.frame(nut.abio$longitude, nut.abio$latitude),
                      env = data.frame(nut.abio$MAP, nut.abio$MAT, nut.abio$pH),
                      emf = nut.emf[-c(1, 2)])
mt.geo <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(nut.emf.sp$geo.dist) - min(nut.emf.sp$geo.dist)) * 0.25, 
          (max(nut.emf.sp$env.dist) - min(nut.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[1:2], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")

p1 <- nut.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist)) +
  geom_point(alpha = 0.6, shape = 1, aes(color = variable)) +
  geom_smooth(method = "lm", size = 0.2) +
  facet_grid(~ variable, scales = "free_x") +
  labs(x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), color = "black", label = mt$labs) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# by Tibet
tibet.emf.sp <- calDist(geo = data.frame(tibet.abio$longitude, tibet.abio$latitude),
                        env = data.frame(tibet.abio$MAP, tibet.abio$MAT, tibet.abio$pH),
                        emf = tibet.emf[-c(1, 2)])
mt.geo <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(tibet.emf.sp$geo.dist) - min(tibet.emf.sp$geo.dist)) * 0.25, 
          (max(tibet.emf.sp$env.dist) - min(tibet.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[3:4], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")
mt$labs <- gsub("P = 0", "P < 0.001", mt$labs)

p2 <- tibet.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist)) +
  geom_point(alpha = 0.6, shape = 1, aes(color = variable)) +
  geom_smooth(method = "lm", size = 0.2) +
  facet_grid(~ variable, scales = "free_x") +
  labs(x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
pdf("./outs/NutTibet_GeoEnv_EMF.pdf", width = 8.5, height = 7.5)
grid.arrange(p1, p2, ncol = 1)
dev.off()

#-----------------------------
# Sensitivity analysis to test whether the number of function matters
# by Nutrient network
nut.sa <- nut.emf.sp %>% 
  reshape2::melt(., id.vars = "emf.dist")
nut.emf.std <- data.frame(apply(nut.emf[-c(1, 2)], 2, varScale))
nut.vars4func <- names(nut.emf.std)
nut.res <- NULL
for (i in seq_along(nut.vars4func)) {
  N <- dim(combn(nut.vars4func, i))[2]
  for (j in 1:N) {
    dat <- nut.emf.std[, combn(nut.vars4func, i)[, j]]
    emf.dist <- distance(dat, method = "euclidean")
    df.dist <- data.frame(nut.sa[-1], emf = rep(emf.dist, 2))
    res <- ddply(df.dist, .(variable), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(nut.vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    nut.res <- rbind(nut.res, res)
  }
}
names(nut.res)[6:7] <- c("llim", "ulim") 
nut.res$variable <- factor(nut.res$variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))
# add labs
nut.labs <- ddply(nut.res, .(variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1],
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
nut.labs$labs <- paste(letters[1:2], 
                       ") Slope = ", round(nut.labs$slope, 3), 
                       ", R^2 = ", round(nut.labs$r.squared, 2),
                       ", P = ", round(nut.labs$P, 3), sep = "")
nut.labs$labs[1] <- gsub("P = 0", "P < 0.001", nut.labs$labs[1])
nut.labs$x <- c(6 * 0.5, 6 * 0.5)
nut.labs$y <- c(0.5 * 0.66, 0.5 * 0.66)
p.nut.sa <- ggplot(nut.res, aes(x = n.funcs, y = mantelr, color = variable,
                                lty = variable)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Nutrient network", 
       x = "Number of functions", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(c(-0.1, 0.35)) +
  geom_text(data = nut.labs, aes(x = x, y = y), label = nut.labs$labs) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# by Tibetan grasslands
tibet.sa <- tibet.emf.sp %>% 
  reshape2::melt(., id.vars = "emf.dist")
tibet.emf.std <- data.frame(apply(tibet.emf[-c(1, 2)], 2, varScale))
tibet.vars4func <- names(tibet.emf.std)
tibet.res <- NULL
for (i in seq_along(tibet.vars4func)) {
  N <- dim(combn(tibet.vars4func, i))[2]
  for (j in 1:N) {
    dat <- tibet.emf.std[, combn(tibet.vars4func, i)[, j]]
    emf.dist <- distance(dat, method = "euclidean")
    df.dist <- data.frame(tibet.sa[-1], emf = rep(emf.dist, 2))
    res <- ddply(df.dist, .(variable), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(tibet.vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    tibet.res <- rbind(tibet.res, res)
  }
}
names(tibet.res)[6:7] <- c("llim", "ulim") 
tibet.res$variable <- factor(tibet.res$variable, levels = c("geo.dist", "env.dist"),
                             labels = c("Geographic distance", "Environmental distance"))
# add labs
tibet.labs <- ddply(tibet.res, .(variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1],
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
tibet.labs$labs <- paste(letters[3:4],
                         ") Slope = ", round(tibet.labs$slope, 3),
                         ", R^2 = ", round(tibet.labs$r.squared, 2),
                         ", P = ", round(tibet.labs$P, 3), sep = "")
tibet.labs$labs <- gsub("P = 0", "P < 0.001", tibet.labs$labs)
tibet.labs$x <- c(6 * 0.5, 6 * 0.5)
tibet.labs$y <- c(0.5 * 0.66, 0.5 * 0.66)
p.tibet.sa <- ggplot(tibet.res, aes(x = n.funcs, y = mantelr, color = variable)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Tibetan grasslands", 
       x = "Number of functions", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(c(-0.1, 0.35)) +
  geom_text(data = tibet.labs, aes(x = x, y = y), label = tibet.labs$labs) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# pdf("./outs/NutTibet_GeoEnv_Sensitivity_analysis.pdf", width = 8, height = 8)
grid.arrange(p.nut.sa, p.tibet.sa)
# dev.off()

#-----------------------------
# Sensitivity analysis on the environmental variables
# by Nutrient network
# use different combinations of the environmental variables
nut.emf.env <- calDist(geo = data.frame(nut.abio$longitude, nut.abio$latitude),
                       env = data.frame(nut.abio$MAP, nut.abio$MAT, nut.abio$pH),
                       emf = nut.emf[-c(1, 2)])
ecodist::mantel(nut.emf.env$env ~ nut.emf.env$emf, nperm = 10000)
# by Tibetan grasslands
tibet.emf.env <- calDist(geo = data.frame(tibet.abio$longitude, tibet.abio$latitude),
                         env = data.frame(tibet.abio$MAP, tibet.abio$MAT, tibet.abio$pH),
                         emf = tibet.emf[-c(1, 2)])
ecodist::mantel(tibet.emf.env$env ~ tibet.emf.env$emf, nperm = 10000)

#-----------------------------
# Geo and env distance with above- and belowground beta-EMF
# by Aboveground ecosysem beta-EMF
nut.emf.sp <- calDist(geo = data.frame(nut.abio$longitude, nut.abio$latitude),
                      env = data.frame(nut.abio$MAP, nut.abio$MAT, nut.abio$pH),
                      emf = nut.emf[, c("AGB", "PTN", "PTP")])
mt.geo <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(nut.emf.sp$geo.dist) - min(nut.emf.sp$geo.dist)) * 0.25, 
          (max(nut.emf.sp$env.dist) - min(nut.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[1:2], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")

p1 <- nut.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist, color = variable)) +
  geom_point(alpha = 0.5, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  facet_grid(~ variable, scales = "free_x") +
  labs(title = "Nutrient network", 
       x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# by Tibet
tibet.emf.sp <- calDist(geo = data.frame(tibet.abio$longitude, tibet.abio$latitude),
                        env = data.frame(tibet.abio$MAP, tibet.abio$MAT, tibet.abio$pH),
                        emf = tibet.emf[, c("AGB", "PTN", "PTP")])
mt.geo <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(tibet.emf.sp$geo.dist) - min(tibet.emf.sp$geo.dist)) * 0.25, 
          (max(tibet.emf.sp$env.dist) - min(tibet.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[3:4], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")
mt$labs <- gsub("P = 0", "P < 0.001", mt$labs)

p2 <- tibet.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist, color = variable)) +
  geom_point(alpha = 0.2, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  facet_grid(~ variable, scales = "free_x") +
  labs(title = "Tibetan grasslands", 
       x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# pdf("./outs/NutTibet_GeoEnv_EMF_aboveground.pdf", width = 8, height = 8)
grid.arrange(p1, p2, ncol = 1)
# dev.off()

# by Belowground ecosysem beta-EMF
nut.emf.sp <- calDist(geo = data.frame(nut.abio$longitude, nut.abio$latitude),
                      env = data.frame(nut.abio$MAP, nut.abio$MAT, nut.abio$pH),
                      emf = nut.emf[, c("BGB", "STN", "STP")])
mt.geo <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(nut.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(nut.emf.sp$geo.dist) - min(nut.emf.sp$geo.dist)) * 0.25, 
          (max(nut.emf.sp$env.dist) - min(nut.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[1:2], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")

p1 <- nut.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist, color = variable)) +
  geom_point(alpha = 0.5, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  facet_grid(~ variable, scales = "free_x") +
  labs(title = "Nutrient network", 
       x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# by Tibet
tibet.emf.sp <- calDist(geo = data.frame(tibet.abio$longitude, tibet.abio$latitude),
                        env = data.frame(tibet.abio$MAP, tibet.abio$MAT, tibet.abio$pH),
                        emf = tibet.emf[, c("BGB", "STN", "STP")])
mt.geo <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ geo.dist, nperm = 10000))
mt.env <- with(tibet.emf.sp, ecodist::mantel(emf.dist ~ env.dist, nperm = 10000))
mt <- data.frame(rbind(mt.geo, mt.env))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, levels = c("mt.geo", "mt.env"),
                      labels = c("Geographic distance", "Environmental distance"))
mt$x <- c((max(tibet.emf.sp$geo.dist) - min(tibet.emf.sp$geo.dist)) * 0.25, 
          (max(tibet.emf.sp$env.dist) - min(tibet.emf.sp$env.dist)) * 0.25)
mt$y <- c(1.85, 1.85)
mt$labs <- paste(letters[3:4], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")
mt$labs <- gsub("P = 0", "P < 0.001", mt$labs)

p2 <- tibet.emf.sp %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("geo.dist", "env.dist"),
                           labels = c("Geographic distance", "Environmental distance"))) %>% 
  ggplot(aes(x = value, y = emf.dist, color = variable)) +
  geom_point(alpha = 0.2, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  facet_grid(~ variable, scales = "free_x") +
  labs(title = "Tibetan grasslands", 
       x = "", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# pdf("./outs/NutTibet_GeoEnv_EMF_belowground.pdf", width = 8, height = 8)
grid.arrange(p1, p2, ncol = 1)
# dev.off()

###########################################################
# inspect the relationship between geographic and environmental distance ---
with(nut.emf.sp, ecodist::mantel(env.dist ~ geo.dist, nperm = 10000))
p3 <- nut.emf.sp %>% 
  ggplot(aes(x = geo.dist, y = env.dist)) +
  geom_point(alpha = 0.1, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  labs(title = "Nutrient network",
       x = "Geographic distance", y = "Environmental distance") +
  lims(y = c(0, 1.6)) +
  annotate("text", x = (max(nut.emf.sp$geo.dist) - min(nut.emf.sp$geo.dist))*0.25, 
           y = 1.5,
           label = "a) r = 0.24, P < 0.001") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
with(tibet.emf.sp, ecodist::mantel(env.dist ~ geo.dist, nperm = 10000))
p4 <- tibet.emf.sp %>% 
  ggplot(aes(x = geo.dist, y = env.dist)) +
  geom_point(alpha = 0.1, shape = 1) +
  geom_smooth(method = "lm", size = 0.2, color = "black") +
  labs(title = "Tibetan grasslands",
       x = "Geographic distance", y = "Environmental distance") +
  lims(y = c(0, 1.6)) +
  annotate("text", x = (max(tibet.emf.sp$geo.dist) - min(tibet.emf.sp$geo.dist))*0.25, 
           y = 1.5,
           label = "b) r = 0.20, P < 0.001") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# pdf("./outs/NutTibet_Geo_Env.pdf", width = 4.0, height = 6)
grid.arrange(p3, p4, ncol = 1)
# dev.off()
rm(list = c("p1", "p2", "p3", "p4", "p.nut.sa", "p.tibet.sa"))

###########################################################
# inspect relationship between beta-diversity and beta-multifunctionality ---

# calculate community dissimilarity according to Jost 2007
# Partitioning diversity into independent alpha and beta components
# I skipped this step and saved the data for later use
# plant.q0 <- 1 - sim.table(nut.plant, q = 0, half = FALSE)  # Sorensen
# plant.q1 <- 1 - sim.table(nut.plant, q = 1, half = FALSE)  # Horn
# plant.q2 <- 1 - sim.table(nut.plant, q = 2, half = FALSE)  # Morisita-Horn
# bac.q0 <- 1 - sim.table(nut.bac, q = 0, half = FALSE)
# bac.q1 <- 1 - sim.table(nut.bac, q = 1, half = FALSE)
# bac.q2 <- 1 - sim.table(nut.bac, q = 2, half = FALSE)
# fun.q0 <- 1 - sim.table(nut.fun, q = 0, half = FALSE)
# fun.q1 <- 1 - sim.table(nut.fun, q = 1, half = FALSE)
# fun.q2 <- 1 - sim.table(nut.fun, q = 2, half = FALSE)

# distCombine <- function() {
#   plant.q0 <- plant.q0[lower.tri(plant.q0)]
#   plant.q1 <- plant.q1[lower.tri(plant.q1)]
#   plant.q2 <- plant.q2[lower.tri(plant.q2)]
#   bac.q0 <- bac.q0[lower.tri(bac.q0)]
#   bac.q1 <- bac.q1[lower.tri(bac.q1)]
#   bac.q2 <- bac.q2[lower.tri(bac.q2)]
#   fun.q0 <- fun.q0[lower.tri(fun.q0)]
#   fun.q1 <- fun.q1[lower.tri(fun.q1)]
#   fun.q2 <- fun.q2[lower.tri(fun.q2)]
#   
#   plant.q0 <- data.frame(plant.q0)
#   plant.q1 <- data.frame(plant.q1)
#   plant.q2 <- data.frame(plant.q2)
#   bac.q0 <- data.frame(bac.q0)
#   bac.q1 <- data.frame(bac.q1)
#   bac.q2 <- data.frame(bac.q2)
#   fun.q0 <- data.frame(fun.q0)
#   fun.q1 <- data.frame(fun.q1)
#   fun.q2 <- data.frame(fun.q2)
#   
#   dist.df0 <- data.frame(cbind(plant.q0, bac.q0, fun.q0))
#   dist.df1 <- data.frame(cbind(plant.q1, bac.q1, fun.q1))
#   dist.df2 <- data.frame(cbind(plant.q2, bac.q2, fun.q2))
#   
#   names(dist.df0) <- c("plant", "bac", "fun")
#   names(dist.df1) <- c("plant", "bac", "fun")
#   names(dist.df2) <- c("plant", "bac", "fun")
#   
#   dist.df <- rbind(dist.df0, dist.df1, dist.df2)
#   dist.df <- data.frame(dist.df)
#   return(dist.df)
# }

# nut.beta <- distCombine()
# nut.beta$q.order <- rep(c(0, 1, 2), each = dim(nut.beta)[1] / 3)
# nut.beta$emf <- rep(nut.emf.sp$emf.dist, 3)
# nut.beta$geo <- rep(nut.emf.sp$geo.dist, 3)
# nut.beta$env <- rep(nut.emf.sp$env.dist, 3)
# write.csv(nut.beta, "./outs/nut_beta.csv")

# plant.q0 <- 1 - sim.table(tibet.plant, q = 0, half = FALSE)  # Sorensen
# plant.q1 <- 1 - sim.table(tibet.plant, q = 1, half = FALSE)  # Horn
# plant.q2 <- 1 - sim.table(tibet.plant, q = 2, half = FALSE)  # Morisita-Horn
# bac.q0 <- 1 - sim.table(tibet.bac, q = 0, half = FALSE)
# bac.q1 <- 1 - sim.table(tibet.bac, q = 1, half = FALSE)
# bac.q2 <- 1 - sim.table(tibet.bac, q = 2, half = FALSE)
# fun.q0 <- 1 - sim.table(tibet.fun, q = 0, half = FALSE)
# fun.q1 <- 1 - sim.table(tibet.fun, q = 1, half = FALSE)
# fun.q2 <- 1 - sim.table(tibet.fun, q = 2, half = FALSE)

# tibet.beta <- distCombine()
# tibet.beta$q.order <- rep(c(0, 1, 2), each = dim(tibet.beta)[1] / 3)
# tibet.beta$emf <- rep(tibet.emf.sp$emf.dist, 3)
# tibet.beta$geo <- rep(tibet.emf.sp$geo.dist, 3)
# tibet.beta$env <- rep(tibet.emf.sp$env.dist, 3)
# 
# write.csv(tibet.beta, "./outs/tibet_beta.csv")

# # load data
# nut.beta <- read.csv("./outs/nut_beta.csv")
# tibet.beta <- read.csv("./outs/tibet_beta.csv")

# # calculate beta-phylogenetic diversity
# nut.clean <- match.phylo.comm(nut.tree.sp, nut.plant)
# tibet.clean <- match.phylo.comm(tibet.tree.sp, tibet.plant)
# nut.mntd.awf <- comdistnt(nut.clean$comm, cophenetic(nut.clean$phy), abundance.weighted = FALSE)
# nut.mntd.awt <- comdistnt(nut.clean$comm, cophenetic(nut.clean$phy), abundance.weighted = TRUE)
# tibet.mntd.awf <- comdistnt(tibet.clean$comm, cophenetic(tibet.clean$phy), abundance.weighted = FALSE)
# tibet.mntd.awt <- comdistnt(tibet.clean$comm, cophenetic(tibet.clean$phy), abundance.weighted = TRUE)
# combine the data
# nut.beta$mntd.awf <- rep(nut.mntd.awf, 3)  # abundance weighted = FALSE
# nut.beta$mntd.awt <- rep(nut.mntd.awt, 3)   # abundance weighted = TRUE
# tibet.beta$mntd.awf <- rep(tibet.mntd.awf, 3)  # abundance weighted = FALSE
# tibet.beta$mntd.awt <- rep(tibet.mntd.awt, 3)   # abundance weighted = TRUE
# save the data 
# write.csv(nut.beta, "./outs/nut_beta.csv")
# write.csv(tibet.beta, "./outs/tibet_beta.csv")

# load the saved data
nut.beta <- read.csv("./outs/nut_beta.csv")
tibet.beta <- read.csv("./outs/tibet_beta.csv")
nut.beta <- nut.beta[-1]
tibet.beta <- tibet.beta[-1]

# inspect the distribution of plant beta-diversity
nut.beta.dist <- nut.beta %>% 
  filter(q.order == 0)
tibet.beta.dist <- tibet.beta %>% 
  filter(q.order == 0)

median(nut.beta.dist$plant)
median(tibet.beta.dist$plant)

pdf("./outs/taxonomic_beta_diversity_distribution.pdf", 
    height = 3.5, width = 5.5)
par(mfrow = c(1, 2))
hist(nut.beta.dist$plant, main = "Nutrient Network",
     xlab = "Taxonomic beta-diversity")
hist(tibet.beta.dist$plant, main = "Tibetan grasslands",
     xlab = "Taxonomic beta-diversity")
dev.off()
par(mfrow = c(1, 1))

# inspect the distribution of bacterial beta-diversity
median(nut.beta.dist$bac)
median(tibet.beta.dist$bac)

pdf("./outs/bacteria_taxonomic_beta_diversity_distribution.pdf", 
    height = 3.5, width = 5.5)
par(mfrow = c(1, 2))
hist(nut.beta.dist$bac, main = "Nutrient Network",
     xlab = "Taxonomic beta-diversity")
hist(tibet.beta.dist$bac, main = "Tibetan grasslands",
     xlab = "Taxonomic beta-diversity")
dev.off()
par(mfrow = c(1, 1))

# mantel test for beta-diversity and beta-EMF
nut.mt <- nut.beta %>% 
  select(-mntd.awf, -mntd.awt) %>% 
  reshape2::melt(., id.vars = c("X", "emf", "q.order", "geo", "env"),
                 variable.name = "beta.var",
                 value.name = "beta.val") %>%
  reshape2::melt(., id.vars = c("X", "q.order", "beta.var", "beta.val")) %>% 
  ddply(.(variable, beta.var, q.order), function(x) {
    mt <- ecodist::mantel(x$value ~ x$beta.val, nperm = 10000)
  })
names(nut.mt)[8:9] <- c("llim", "ulim")
nut.mt$variable <- factor(nut.mt$variable, levels = c("emf", "geo", "env"),
                          labels = c("Beta-multifunctionality", "Geographic distance", 
                                     "Environmental distance"))
nut.mt$beta.var <- factor(nut.mt$beta.var, levels = c("plant", "bac", "fun"),
                          labels = c("Plant", "Bacteria", "Fungi"))

tibet.mt <- tibet.beta %>% 
  select(-mntd.awf, -mntd.awt) %>% 
  reshape2::melt(., id.vars = c("X", "emf", "q.order", "geo", "env"),
                 variable.name = "beta.var",
                 value.name = "beta.val") %>%
  reshape2::melt(., id.vars = c("X", "q.order", "beta.var", "beta.val")) %>% 
  ddply(.(variable, beta.var, q.order), function(x) {
    mt <- ecodist::mantel(x$value ~ x$beta.val, nperm = 10000)
  })
names(tibet.mt)[8:9] <- c("llim", "ulim")
tibet.mt$variable <- factor(tibet.mt$variable, levels = c("emf", "geo", "env"),
                            labels = c("Beta-multifunctionality", "Geographic distance", 
                                       "Environmental distance"))
tibet.mt$beta.var <- factor(tibet.mt$beta.var, levels = c("plant", "bac", "fun"),
                            labels = c("Plant", "Bacteria", "Fungi"))

p5 <- nut.mt %>% 
  mutate(variable = factor(variable,
                           levels = c("Geographic distance",
                                      "Environmental distance",
                                      "Beta-multifunctionality"))) %>% 
  ggplot(aes(x = beta.var, y = mantelr, color = factor(q.order))) +
  geom_point(position = position_dodge(0.5), size = 2.5) +
  geom_errorbar(aes(ymin = llim, ymax = ulim),
                width = 0.1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Nutrient network", x = "", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(-0.03, 0.6) +
  theme_bw() +
  scale_color_discrete(name = "q order") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.25, 0.65))
p6 <- tibet.mt %>% 
  mutate(variable = factor(variable,
                           levels = c("Geographic distance",
                                      "Environmental distance",
                                      "Beta-multifunctionality"))) %>% 
  ggplot(aes(x = beta.var, y = mantelr, color = factor(q.order))) +
  geom_point(position = position_dodge(0.5), size = 2.5) +
  geom_errorbar(aes(ymin = llim, ymax = ulim),
                width = 0.1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Tibetan grasslands", x = "", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(c(-0.03, 0.6)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

# pdf("./outs/NutTibet_beta_EMFGeoEnv.pdf", width = 7.5, height = 5.5)
grid.arrange(p5, p6, ncol = 1)
# dev.off()

#----------------------------
# Sensitivity analysis: number of functions
nut.df <- nut.beta %>% 
  dplyr::select(plant, bac, fun, q.order) %>% 
  reshape2::melt(., id.vars = c("q.order"))
emf <- data.frame(apply(nut.emf[-c(1, 2)], 2, varScale))
vars4func <- names(emf)
res.all <- NULL
for (i in 1:6) {
  N <- dim(combn(vars4func, i))[2]
  for (j in 1:N) {
    dat <- emf[, combn(vars4func, i)[, j]]
    emf.dist <- distance(dat, method = "euclidean")
    emf.dist <- rep(emf.dist, 9)
    dist.df <- data.frame(nut.df, emf.dist)
    res <- ddply(dist.df, .(variable, q.order), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    res.all <- rbind(res.all, res)
  }
}
names(res.all)[7:8] <- c("llim", "ulim")
res.all$fq <- factor(res.all$q.order)
res.all$variable <- factor(res.all$variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))
df.labs <- ddply(res.all, .(q.order, variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1], 
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
df.labs$labs <- paste(letters[1:9], ") Slope = ", round(df.labs$slope, 3),
                      ", R^2 = ", round(df.labs$r.squared, 2),
                      ", P = ", round(df.labs$P, 3), sep = "")
df.labs$labs[c(1, 4, 7)] <- gsub("P = 0", "P < 0.001", df.labs$labs[c(1, 4, 7)])
df.labs$x <- rep(3.5, 9)
df.labs$y <- rep(0.45, 9)
df.labs$fq <- factor(df.labs$q.order)

ggplot(res.all, aes(x = n.funcs, y = mantelr, color = fq)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(x = "Number of functions", y = "Mantel coefficient (r)") +
  facet_grid(fq ~ variable) +
  lims(y = c(-0.1, 0.5)) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# ggsave("./outs/Nut_Sensitivity_analysis.pdf", width = 10, height = 6)

tibet.df <- tibet.beta %>% 
  dplyr::select(plant, bac, fun, q.order) %>% 
  reshape2::melt(., id.vars = c("q.order"))
emf <- data.frame(apply(tibet.emf[-c(1, 2)], 2, varScale))
vars4func <- names(emf)
res.all <- NULL
for (i in 1:6) {
  N <- dim(combn(vars4func, i))[2]
  for (j in 1:N) {
    dat <- emf[, combn(vars4func, i)[, j]]
    emf.dist <- distance(dat, method = "euclidean")
    emf.dist <- rep(emf.dist, 9)
    dist.df <- data.frame(tibet.df, emf.dist)
    res <- ddply(dist.df, .(variable, q.order), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    res.all <- rbind(res.all, res)
  }
}
names(res.all)[7:8] <- c("llim", "ulim")
res.all$fq <- factor(res.all$q.order)
res.all$variable <- factor(res.all$variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))
df.labs <- ddply(res.all, .(q.order, variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1],
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
df.labs$labs <- paste(letters[1:9], ") Slope = ", round(df.labs$slope, 3), 
                      ", R^2 = ", round(df.labs$r.squared, 2),
                      ", P = ", round(df.labs$P, 3), sep = "")
df.labs$labs[-c(7, 9)] <- gsub("P = 0", "P < 0.001", df.labs$labs[-c(7, 9)])
df.labs$x <- rep(3.5, 9)
df.labs$y <- rep(0.42, 9)
df.labs$fq <- factor(df.labs$q.order)

ggplot(res.all, aes(x = n.funcs, y = mantelr, color = fq)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", size = 0.2, se = FALSE, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(x = "Number of functions", y = "Mantel coefficient (r)") +
  ylim(c(-0.05, 0.45)) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  facet_grid(fq ~ variable) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# ggsave("./outs/Tibet_Sensitivity_analysis.pdf", width = 10, height = 6)

#----------------------------
# Sensitivity analysis: environmental variables
# by Nutrient network
nut.emf.env <- rep(distance(nut.abio[, c("MAP", "MAT", "pH")], method = "euclidean"), 9)
df.nut <- nut.beta %>% 
  select(plant, bac, fun, q.order) %>% 
  reshape2::melt(id.vars = "q.order")
df.nut$env <- nut.emf.env
ddply(df.nut, .(variable, q.order), function(x) {
  ecodist::mantel(x$value ~ x$env, nperm = 10000)
})

# by Tibetan grasslands
tibet.emf.env <- rep(distance(tibet.abio[, c("MAP", "MAT", "pH")], method = "euclidean"), 9)
df.tibet <- tibet.beta %>% 
  select(plant, bac, fun, q.order) %>% 
  reshape2::melt(id.vars = "q.order")
df.tibet$env <- tibet.emf.env
ddply(df.tibet, .(variable, q.order), function(x) {
  ecodist::mantel(x$value ~ x$env, nperm = 10000)
})

#
# by Nutrient network
df.nut <- nut.beta %>% 
  select(plant, bac, fun, q.order, emf) %>% 
  filter(q.order == 0) %>% 
  reshape2::melt(id.vars = c("q.order", "emf")) %>% 
  mutate(variable = factor(variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi")))
(range(df.nut$emf)[2] - range(df.nut$emf)[1])/2
df.labs <- nut.mt %>% 
  filter(variable == "Beta-multifunctionality" & q.order == 0)
df.labs$labs <- paste(letters[1:3], ") r = ", round(df.labs$mantelr, 2),
                      ", P = ", round(df.labs$pval3, 3), sep = "")  
df.labs$labs[1] <- gsub("P = 0", "P < 0.001", df.labs$labs[1])
df.labs$x <- rep(.23, 3)
df.labs$y <- rep((range(df.nut$emf)[2] - range(df.nut$emf)[1])*0.99, 3)
df.labs <- df.labs[-1]
names(df.labs)[1] <- "variable"

p1 <- ggplot(df.nut, aes(x = value, y = emf)) +
  geom_point(aes(color = variable, alpha = 0.6)) +
  geom_smooth(method = "lm", aes(lty = variable), lwd = 0.35) +
  scale_linetype_manual(values = c(1, 2, 2)) +
  geom_vline(xintercept = 0.5, color = "gray50", lty = 2) +
  geom_hline(yintercept = 0.80, color = "gray50", lty = 2) +
  facet_grid(~ variable) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  labs(x = "Community beta-diversity",
       y = "Ecosystem beta-multifunctionality") +
  expand_limits(x = 0, y = 0) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

# by Tibetan grassland
df.tibet <- tibet.beta %>% 
  select(plant, bac, fun, q.order, emf) %>% 
  filter(q.order == 0) %>% 
  reshape2::melt(id.vars = c("q.order", "emf")) %>% 
  mutate(variable = factor(variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi")))
(range(df.tibet$emf)[2] - range(df.tibet$emf)[1])/2
df.labs <- tibet.mt %>% 
  filter(variable == "Beta-multifunctionality" & q.order == 0)
df.labs$labs <- paste(letters[4:6], ") r = ", round(df.labs$mantelr, 2),
                      ", P = ", round(df.labs$pval3, 3), sep = "")  
df.labs$labs[1:2] <- gsub("P = 0", "P < 0.001", df.labs$labs[1:2])
df.labs$x <- rep(.23, 3)
df.labs$y <- rep((range(df.tibet$emf)[2] - range(df.tibet$emf)[1])*0.99, 3)
df.labs <- df.labs[-1]
names(df.labs)[1] <- "variable"

p2 <- ggplot(df.tibet, aes(x = value, y = emf)) +
  geom_point(aes(color = variable, alpha = 0.6)) +
  geom_smooth(method = "lm", aes(lty = variable), lwd = 0.35) +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_vline(xintercept = 0.5, color = "gray50", lty = 2) +
  geom_hline(yintercept = 0.82, color = "gray50", lty = 2) +
  facet_grid(~ variable) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  labs(x = "Community beta-diversity",
       y = "Ecosystem beta-multifunctionality") +
  expand_limits(x = 0, y = 0) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

pdf("./outs/NutTibet_beta_EMFdiv.pdf", width = 10.5, height = 7.5)
grid.arrange(p1, p2, ncol = 1)
dev.off()

###########################################################
# Plant beta-phylogenetic diversity and ecosystem multifunctionality
# caculate plant genus and family table for Nutrient network
nut.species.list <- read.csv("./data/data_processing/NutNet_SpeciesList.csv")
nut.plant.phy <- data.frame(t(nut.plant))
nut.plant.phy$species <- rownames(nut.plant.phy)
nut.plant.phy <- merge(nut.plant.phy, nut.species.list, by = "species") 
# by genus
nut.plant.genus <- nut.plant.phy %>% 
  select(-family, -species) %>% 
  ddply(., .(genus), colwise(sum)) %>% 
  data.frame()
rownames(nut.plant.genus) <- nut.plant.genus$genus 
nut.plant.genus <- nut.plant.genus %>% 
  select(-genus) %>% 
  t() %>% 
  data.frame()
nut.plant.genus <- data.frame(ifelse(nut.plant.genus > 0, 1, 0))
# by family
nut.plant.family <- nut.plant.phy %>% 
  select(-genus, -species) %>% 
  ddply(., .(family), colwise(sum)) %>% 
  data.frame()
rownames(nut.plant.family) <- nut.plant.family$family 
nut.plant.family <- nut.plant.family %>% 
  select(-family) %>% 
  t() %>% 
  data.frame()
nut.plant.family <- data.frame(ifelse(nut.plant.family > 0, 1, 0))
# caculate plant genus and family table for Tibetan grasslands
tibet.species.list <- read.csv("./data/data_processing/Tibet_plant_species_taxa_table.csv")
tibet.plant.phy <- data.frame(t(tibet.plant))
tibet.plant.phy$species <- rownames(tibet.plant.phy)
tibet.plant.phy <- merge(tibet.plant.phy, tibet.species.list, by = "species") 
# by genus
tibet.plant.genus <- tibet.plant.phy %>% 
  select(-family, -species) %>% 
  ddply(., .(genus), colwise(sum)) %>% 
  data.frame()
rownames(tibet.plant.genus) <- tibet.plant.genus$genus 
tibet.plant.genus <- tibet.plant.genus %>% 
  select(-genus) %>% 
  t() %>% 
  data.frame()
tibet.plant.genus <- data.frame(ifelse(tibet.plant.genus > 0, 1, 0))
# by family
tibet.plant.family <- tibet.plant.phy %>% 
  select(-genus, -species) %>% 
  ddply(., .(family), colwise(sum)) %>% 
  data.frame()
rownames(tibet.plant.family) <- tibet.plant.family$family 
tibet.plant.family <- tibet.plant.family %>% 
  select(-family) %>% 
  t() %>% 
  data.frame()
tibet.plant.family <- data.frame(ifelse(tibet.plant.family > 0, 1, 0))

# caculate plant genus and family tree for Nutrient network
# genus tree
nut.tips <- nut.tree.sp$tip.label
nut.genus <- unique(sapply(strsplit(nut.tips, "_"), function(x) x[1]))
ii <- sapply(nut.genus, function(x, y) grep(x, y)[1], y = nut.tips)
nut.tree.genus <- drop.tip(nut.tree.sp, setdiff(nut.tree.sp$tip.label, nut.tips[ii]))
nut.tree.genus$tip.label <- sapply(strsplit(nut.tree.genus$tip.label, "_"), function(x) x[1])
# family tree
nut.tree.family <- nut.tree.genus
nut.family.list <- nut.species.list[, -3]
nut.family.list <- nut.family.list[which(!duplicated(nut.family.list[, 2])), ]
nut.family.list[, 1] <- paste(nut.family.list$family, nut.family.list$genus, sep = "_")
nut.family.list <- nut.family.list[which(nut.family.list$genus %in% nut.tree.family$tip.label), ]
for (i in 1:length(nut.tree.family$tip.label)){
  nut.tree.family$tip.label[i] <- nut.family.list[, 1][which(nut.family.list[, 2] == nut.tree.family$tip.label[i])]
}
nut.tips <- nut.tree.family$tip.label
nut.family <- unique(sapply(strsplit(nut.tips, "_"), function(x) x[1]))
ii <- sapply(nut.family, function(x, y) grep(x, y)[1], y = nut.tips)
nut.tree.family <- drop.tip(nut.tree.family, setdiff(nut.tree.family$tip.label, nut.tips[ii]))
nut.tree.family$tip.label <- sapply(strsplit(nut.tree.family$tip.label, "_"), function(x) x[1])
# caculate plant genus and family tree for Tibetan grasslands
# genus tree
tibet.tips <- tibet.tree.sp$tip.label
tibet.genus <- unique(sapply(strsplit(tibet.tips, "_"), function(x) x[1]))
ii <- sapply(tibet.genus, function(x, y) grep(x, y)[1], y = tibet.tips)
tibet.tree.genus <- drop.tip(tibet.tree.sp, setdiff(tibet.tree.sp$tip.label, tibet.tips[ii]))
tibet.tree.genus$tip.label <- sapply(strsplit(tibet.tree.genus$tip.label, "_"), function(x) x[1])
# family tree
tibet.tree.family <- tibet.tree.genus
tibet.family.list <- tibet.species.list[, -3]
tibet.family.list <- tibet.family.list[which(!duplicated(tibet.family.list[, 2])), ]
tibet.family.list[, 1] <- paste(tibet.family.list$family, tibet.family.list$genus, sep = "_")
tibet.family.list <- tibet.family.list[which(tibet.family.list$genus %in% tibet.tree.family$tip.label), ]
for (i in 1:length(tibet.tree.family$tip.label)){
  tibet.tree.family$tip.label[i] <- tibet.family.list[, 1][which(tibet.family.list[, 2] == tibet.tree.family$tip.label[i])]
}
tibet.tips <- tibet.tree.family$tip.label
tibet.family <- unique(sapply(strsplit(tibet.tips, "_"), function(x) x[1]))
ii <- sapply(tibet.family, function(x, y) grep(x, y)[1], y = tibet.tips)
tibet.tree.family <- drop.tip(tibet.tree.family, setdiff(tibet.tree.family$tip.label, tibet.tips[ii]))
tibet.tree.family$tip.label <- sapply(strsplit(tibet.tree.family$tip.label, "_"), function(x) x[1])

# # calculate beta-phylogenetic diversity at species level
# nut.clean.sp <- match.phylo.comm(nut.tree.sp, nut.plant)
# tibet.clean.sp <- match.phylo.comm(tibet.tree.sp, tibet.plant)
# nut.mpd.sp <- comdist(nut.clean.sp$comm, cophenetic(nut.clean.sp$phy), abundance.weighted = FALSE)
# tibet.mpd.sp <- comdist(tibet.clean.sp$comm, cophenetic(tibet.clean.sp$phy), abundance.weighted = FALSE)
# # calculate beta-phylogenetic diversity at genus and family level
# nut.clean.genus <- match.phylo.comm(nut.tree.genus, nut.plant.genus)
# nut.clean.family <- match.phylo.comm(nut.tree.family, nut.plant.family)
# tibet.clean.genus <- match.phylo.comm(tibet.tree.genus, tibet.plant.genus)
# tibet.clean.family <- match.phylo.comm(tibet.tree.family, tibet.plant.family)
# nut.mpd.genus <- comdist(nut.clean.genus$comm, cophenetic(nut.clean.genus$phy), abundance.weighted = FALSE)
# tibet.mpd.genus <- comdist(tibet.clean.genus$comm, cophenetic(tibet.clean.genus$phy), abundance.weighted = FALSE)
# nut.mpd.family <- comdist(nut.clean.family$comm, cophenetic(nut.clean.family$phy), abundance.weighted = FALSE)
# tibet.mpd.family <- comdist(tibet.clean.family$comm, cophenetic(tibet.clean.family$phy), abundance.weighted = FALSE)

# combine the beta-phylogenetic diversity into a distance data frame
# nut.beta.phy <- nut.beta %>% 
#   filter(q.order == 0) %>% 
#   select(emf, geo, env)
# nut.beta.phy <- data.frame(cbind(nut.beta.phy, 
#                                  as.vector(nut.mpd.sp), 
#                                  as.vector(nut.mpd.genus), 
#                                  as.vector(nut.mpd.family)))
# names(nut.beta.phy)[4:6] <- c("Species", "Genus", "Family")
# write.csv(nut.beta.phy, "./outs/nut_beta_phy.csv")
# 
# tibet.beta.phy <- tibet.beta %>% 
#   filter(q.order == 0) %>% 
#   select(emf, geo, env)
# tibet.beta.phy <- data.frame(cbind(tibet.beta.phy, 
#                                  as.vector(tibet.mpd.sp), 
#                                  as.vector(tibet.mpd.genus), 
#                                  as.vector(tibet.mpd.family)))
# names(tibet.beta.phy)[4:6] <- c("Species", "Genus", "Family")
# write.csv(tibet.beta.phy, "./outs/tibet_beta_phy.csv")

# inspect the distribution of phylogenetic beta-diversity
nut.beta.phy <- read.csv("./outs/nut_beta_phy.csv")
tibet.beta.phy <- read.csv("./outs/tibet_beta_phy.csv")
nut.beta.phy <- nut.beta.phy[-1]
tibet.beta.phy <- tibet.beta.phy[-1]

pdf("./outs/phylogenetic_beta_diversity_distribution.pdf", 
    height = 3.5, width = 5.5)
par(mfrow = c(1, 2))
hist(nut.beta.phy$Species, main = "Nutrient Network",
     xlab = "Phylogenetic beta-diversity")
hist(tibet.beta.phy$Species, main = "Tibetan grasslands",
     xlab = "Phylogenetic beta-diversity")
dev.off()
par(mfrow = c(1, 1))

# mantel test for beta-diversity and beta-EMF
nut.mt <- nut.beta.phy %>% 
  reshape2::melt(., id.vars = c("emf", "geo", "env"),
                 variable.name = "beta.var",
                 value.name = "beta.val") %>%
  reshape2::melt(., id.vars = c("beta.var", "beta.val")) %>% 
  ddply(.(variable, beta.var), function(x) {
    mt <- ecodist::mantel(x$value ~ x$beta.val, nperm = 10000)
  })
names(nut.mt)[7:8] <- c("llim", "ulim")
nut.mt$variable <- factor(nut.mt$variable, levels = c("geo", "env", "emf"),
                          labels = c("Geographic distance", 
                                     "Environmental distance",
                                     "Beta-multifunctionality"))
nut.mt$beta.var <- factor(nut.mt$beta.var, levels = c("Species", "Genus", "Family"))

tibet.mt <- tibet.beta.phy %>% 
  reshape2::melt(., id.vars = c("emf", "geo", "env"),
                 variable.name = "beta.var",
                 value.name = "beta.val") %>%
  reshape2::melt(., id.vars = c("beta.var", "beta.val")) %>% 
  ddply(.(variable, beta.var), function(x) {
    mt <- ecodist::mantel(x$value ~ x$beta.val, nperm = 10000)
  })
names(tibet.mt)[7:8] <- c("llim", "ulim")
tibet.mt$variable <- factor(tibet.mt$variable, levels = c("geo", "env", "emf"),
                            labels = c("Geographic distance", 
                                       "Environmental distance",
                                       "Beta-multifunctionality"))
tibet.mt$beta.var <- factor(tibet.mt$beta.var, levels = c("Species", "Genus", "Family"))

p5 <- nut.mt %>% 
  ggplot(aes(x = beta.var, y = mantelr)) +
  geom_point(position = position_dodge(0.5), size = 2.5) +
  geom_errorbar(aes(ymin = llim, ymax = ulim),
                width = 0.1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Nutrient network", x = "", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(-0.03, 0.45) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.25, 0.65))
p6 <- tibet.mt %>% 
  ggplot(aes(x = beta.var, y = mantelr)) +
  geom_point(position = position_dodge(0.5), size = 2.5) +
  geom_errorbar(aes(ymin = llim, ymax = ulim),
                width = 0.1, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Tibetan grasslands", x = "", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  ylim(c(-0.03, 0.45)) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

# pdf("./outs/NutTibet_beta_EMFGeoEnv_phy.pdf", width = 9, height = 7)
grid.arrange(p5, p6, ncol = 1)
# dev.off()

#----------------------------
# Sensitivity analysis: number of functions
nut.df.phy <- nut.beta.phy %>% 
  dplyr::select(Species, Genus, Family) %>% 
  reshape2::melt()
emf <- data.frame(apply(nut.emf[-c(1, 2)], 2, varScale))
vars4func <- names(emf)
res.all <- NULL
for (i in 1:6) {
  N <- dim(combn(vars4func, i))[2]
  for (j in 1:N) {
    dat <- emf[, combn(vars4func, i)[, j]]
    emf.dist <- as.vector(distance(dat, method = "euclidean"))
    emf.dist <- rep(emf.dist, 3)
    dist.df <- data.frame(nut.df.phy, emf.dist)
    res <- ddply(dist.df, .(variable), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    res.all <- rbind(res.all, res)
  }
}
names(res.all)[6:7] <- c("llim", "ulim")
res.all$variable <- factor(res.all$variable, levels = c("Species", "Genus", "Family"))
df.labs <- ddply(res.all, .(variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1],
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
df.labs$labs <- paste(letters[1:3], 
                      ") Slope = ", round(df.labs$slope, 3),
                      ", R^2 = ", round(df.labs$r.squared, 2),
                      ", P = ", round(df.labs$P, 3), sep = "")
df.labs$x <- rep(3.5, 3)
df.labs$y <- rep(0.425, 3)

nut.phy <- ggplot(res.all, aes(x = n.funcs, y = mantelr)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Nutrient network", 
       x = "Number of functions", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  lims(y = c(-0.05, 0.45)) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# ggsave("./outs/Nut_Sensitivity_analysis_phy.pdf", width = 8, height = 6)

tibet.df.phy <- tibet.beta.phy %>% 
  dplyr::select(Species, Genus, Family) %>% 
  reshape2::melt()
emf <- data.frame(apply(tibet.emf[-c(1, 2)], 2, varScale))
vars4func <- names(emf)
res.all <- NULL
for (i in 1:6) {
  N <- dim(combn(vars4func, i))[2]
  for (j in 1:N) {
    dat <- emf[, combn(vars4func, i)[, j]]
    emf.dist <- as.vector(distance(dat, method = "euclidean"))
    emf.dist <- rep(emf.dist, 3)
    dist.df <- data.frame(tibet.df.phy, emf.dist)
    res <- ddply(dist.df, .(variable), function(x) {
      mt <- ecodist::mantel(x$value ~ x$emf, nperm = 1)
    })
    print(paste(i, "functions out of", length(vars4func), "completed", sep = " "))
    res$n.funcs <- i
    res$n.combn <- j
    res.all <- rbind(res.all, res)
  }
}
names(res.all)[6:7] <- c("llim", "ulim")
res.all$variable <- factor(res.all$variable, levels = c("Species", "Genus", "Family"))
df.labs <- ddply(res.all, .(variable), function(x) {
  mod <- lm(mantelr ~ n.funcs, data = x)
  cbind(slope = summary(mod)$coefficients[2, 1],
        r.squared = summary(mod)$r.squared,
        P = summary(mod)$coefficients[2, 4])
})
df.labs$labs <- paste(letters[4:6], ") Slope = ", round(df.labs$slope, 3),
                      ", R^2 = ", round(df.labs$r.squared, 2),
                      ", P = ", round(df.labs$P, 3), sep = "")
df.labs$labs <- gsub("P = 0", "P < 0.001", df.labs$labs)
df.labs$x <- rep(3.5, 3)
df.labs$y <- rep(0.425, 3)

tibet.phy <- ggplot(res.all, aes(x = n.funcs, y = mantelr)) +
  geom_point(shape = 1, position = position_jitter(0.09), size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.2, color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "gray70") +
  labs(title = "Tibetan grasslands", 
       x = "Number of functions", y = "Mantel coefficient (r)") +
  facet_grid(~ variable) +
  lims(y = c(-0.05, 0.45)) +
  geom_text(data = df.labs, aes(x = x, y = y), 
            label = df.labs$labs, inherit.aes = FALSE) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
# ggsave("./outs/Tibet_Sensitivity_analysis.pdf", width = 8, height = 6)
# pdf("./outs/NutTibet_Sensitivity_analysis_phy.pdf", width = 9.5, height = 7)
grid.arrange(nut.phy, tibet.phy, ncol = 1)
# dev.off()

###########################################################
# Sensitivity analyses
# Nutrient network ---
# keep only samples from the North America
nut.emf2 <- nut.emf[grep(".us", nut.emf$site_code), ]
nut.abio2 <- nut.abio[grep(".us", nut.abio$site_code), ]
nut.plant2 <- nut.plant[grep(".us", rownames(nut.plant)), ]
nut.bac2 <- nut.bac[grep(".us", rownames(nut.bac)), ]
nut.fun2 <- nut.fun[grep(".us", rownames(nut.fun)), ]
nut.plant2 <- nut.plant2[rowSums(nut.plant2) > 0, ]
nut.plant2 <- nut.plant2[, colSums(nut.plant2) > 0]
nut.bac2 <- nut.bac2[rowSums(nut.bac2) > 0, ]
nut.bac2 <- nut.bac2[, colSums(nut.bac2) > 0]
nut.fun2 <- nut.fun2[rowSums(nut.fun2) > 0, ]
nut.fun2 <- nut.fun2[, colSums(nut.fun2) > 0]
# calculate distance matrix
nut.dist <- calDist(geo = data.frame(nut.abio2$longitude, nut.abio2$latitude),
                    env = data.frame(nut.abio2$MAP, nut.abio2$MAT, nut.abio2$pH),
                    emf = nut.emf2[-c(1, 2)])
# nut.plant.temp <- 1 - sim.table(nut.plant2, q = 0, half = FALSE) 
# nut.dist$plant <- nut.plant.temp[lower.tri(nut.plant.temp)]
# nut.bac.temp <- 1 - sim.table(nut.bac2, q = 0, half = FALSE) 
# nut.dist$bac <- nut.bac.temp[lower.tri(nut.bac.temp)]
# nut.fun.temp <- 1 - sim.table(nut.fun2, q = 0, half = FALSE) 
# nut.dist$fun <- nut.fun.temp[lower.tri(nut.fun.temp)]
# write.csv(nut.dist, "./outs/nut_dist_Sensitivity_analysis.csv")
nut.dist <- read.csv("./outs/nut_dist_Sensitivity_analysis.csv")

df.nut.dist <- nut.dist %>% 
  dplyr::select(emf.dist, plant, fun, bac) %>% 
  reshape2::melt(id.vars = "emf.dist") %>% 
  mutate(variable = factor(variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))) 
mt <- ddply(df.nut.dist, .(variable), function(x) {
  ecodist::mantel(x$value ~ x$emf.dist, nperm = 10000)
})
mt$x <- c(0.27, 0.27, 0.27)
mt$y <- c(1.95, 1.95, 1.95)
mt$labs <- paste(letters[1:3], ") r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")
mt$labs[1] <- gsub("P = 0", "P < 0.001", mt$labs[1])

ggplot(df.nut.dist, aes(x = value, y = emf.dist)) +
  geom_point(alpha = 0.8, aes(color = variable)) +
  geom_smooth(method = "lm", size = 0.2, aes(lty = variable)) +
  scale_linetype_manual(values = c(1, 2, 2)) +
  geom_vline(xintercept = 0.5, color = "gray50", lty = 2) +
  geom_hline(yintercept = 0.88, color = "gray50", lty = 2) +
  facet_grid(~ variable, scales = "free_x") +
  expand_limits(x = 0, y = 0) +
  labs(x = "Community beta-diversity", y = "Ecosystem beta-multifunctionality") +
  lims(y = c(0, 2)) +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  theme_bw(base_size = 14.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# ggsave("./outs/Nut_Sensitivity_Analysis_North_America.pdf", width = 9, height = 4)

# Tibetan grasslands ---
calDistTibet <- function(emf, plant, bacteria, fungi, N) {
  plant$site <- emf$site
  bacteria$site <- emf$site
  fungi$site <- emf$site
  ids <- sample(unique(emf$site), N)
  # EMF
  emf <- emf[which(emf$site %in% ids), ]
  emf <- emf[, -c(1:2)]
  emf <- data.frame(apply(emf, 2, varScale))
  emf.dist <- distance(emf, method = "euclidean")
  # plant community
  plant <- plant[which(plant$site %in% ids), ]
  plant <- plant[, names(plant) != "site"]
  # plant <- ifelse(plant > 0, 1, 0)
  plant <- plant[rowSums(plant) > 0, ]
  plant <- plant[, colSums(plant) > 0]
  plant.dist <- distance(plant, method = "sorensen")
  # bacteria community
  bacteria <- bacteria[which(bacteria$site %in% ids), ]
  bacteria <- bacteria[, names(bacteria) != "site"]
  # bacteria <- ifelse(bacteria > 0, 1, 0)
  bacteria <- bacteria[rowSums(bacteria) > 0, ]
  bacteria <- bacteria[, colSums(bacteria) > 0]
  bacteria.dist <- distance(bacteria, method = "sorensen")
  # fungi community
  fungi <- fungi[which(fungi$site %in% ids), ]
  fungi <- fungi[, names(fungi) != "site"]
  # fungi <- ifelse(fungi > 0, 1, 0)
  fungi <- fungi[rowSums(fungi) > 0, ]
  fungi <- fungi[, colSums(fungi) > 0]
  fungi.dist <- distance(fungi, method = "sorensen")
  # combine the distance matrix
  dist.df <- cbind(plant.dist, bacteria.dist, fungi.dist, emf.dist)
  dist.df <- data.frame(dist.df)
  dist.df <- reshape2::melt(dist.df, id.vars = "emf.dist")
  mantel.test <- ddply(dist.df, .(variable), function(x) {
    ecodist::mantel(x$value ~ x$emf.dist, nperm = 1)
  })
  return(mantel.test)
}

tibet.plant.PA <- data.frame(ifelse(tibet.plant > 0, 1, 0))
tibet.bac.PA <- data.frame(ifelse(tibet.bac > 0, 1, 0))
tibet.fun.PA <- data.frame(ifelse(tibet.fun > 0, 1, 0))

# calDistTibet(emf = tibet.emf, 
#              plant = tibet.plant.PA,
#              bacteria = tibet.bac.PA,
#              fungi = tibet.fun.PA, N = 18)
# nSim <- 1000
# Xsim <- list()
# for (i in seq_len(nSim)) {
#   Xsim[[i]] <- calDistTibet(emf = tibet.emf, 
#                             plant = tibet.plant.PA,
#                             bacteria = tibet.bac.PA,
#                             fungi = tibet.fun.PA, N = 18)
#   print(paste(i, "loops out of", nSim, "completed", sep = " "))
# }
# Xsim2 <- do.call(rbind, Xsim)
# write.csv(Xsim2, "./outs/Xsim2.csv")
Xsim <- read.csv("./outs/Xsim2.csv")

Xobs <- tibet.beta %>%
  filter(q.order == 0) %>% 
  select(emf, plant, bac, fun) %>% 
  reshape2::melt(id.vars = "emf") %>% 
  ddply(., .(variable), function(x) {
    ecodist::mantel(x$value ~ x$emf, nperm = 1)
  })
Xsim$variable <- factor(Xsim$variable,
                        levels = c("plant.dist", "bacteria.dist", "fungi.dist"),
                        labels = c("Plant", "Bacteria", "Fungi"))
Xobs$variable <- factor(Xobs$variable,
                        levels = c("plant", "bac", "fun"),
                        labels = c("Plant", "Bacteria", "Fungi"))
Xobs$nut <- c(0.36, 0.06, 0.05)
ggplot(Xsim, aes(x = mantelr, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.8) +
  geom_vline(data = Xobs, aes(xintercept = mantelr),
             lty = 1, color = "black") +
  geom_vline(data = Xobs, aes(xintercept = nut),
             lty = 2, color = "gray50") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Mantel coefficient (r)",
       y = "Count") +
  facet_grid(~ variable) +
  theme_bw(base_size = 16.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
ggsave("./outs/Tibet_SES.pdf", width = 9.5, height = 5.5)
# calculate standardized index
Xsim.sd <- ddply(Xsim, .(variable), function(x) {
  data.frame(mean = mean(x$mantelr),
             sd = sd(x$mantelr))
})
(Xobs$mantelr - Xsim.sd$mean)/Xsim.sd$sd
# anova and TukeyHSD analysis
summary(aov(mantelr ~ variable, data = Xsim))
TukeyHSD(aov(mantelr ~ variable, data = Xsim))

##########################################################
# Structural equation model: distance based-piecewiseSEM
# scale the distance variables
nut.beta.sem1 <- nut.beta %>% 
  filter(q.order == 0)
nut.beta.sem1[, c("plant", "bac", "fun", "emf", "geo", "env")] <- apply(nut.beta.sem1[, c("plant", "bac", "fun", "emf", "geo", "env")], 2, varScale)
tibet.beta.sem1 <- tibet.beta %>% 
  filter(q.order == 0)
tibet.beta.sem1[, c("plant", "bac", "fun", "emf", "geo", "env")] <- apply(tibet.beta.sem1[, c("plant", "bac", "fun", "emf", "geo", "env")], 2, varScale)
# significance tests using MRM in ecodist package
MRM(plant ~ env, nperm = 10000, data = nut.beta.sem1)
MRM(bac ~ env + plant, nperm = 10000, data = nut.beta.sem1)
MRM(fun ~ env + plant, nperm = 10000, data = nut.beta.sem1)
MRM(emf ~ env + plant + bac + fun, nperm = 10000, data = nut.beta.sem1)

MRM(plant ~ env, nperm = 10000, data = tibet.beta.sem1)
MRM(bac ~ env + plant, nperm = 10000, data = tibet.beta.sem1)
MRM(fun ~ env + plant, nperm = 10000, data = tibet.beta.sem1)
MRM(emf ~ env + plant + bac + fun, nperm = 10000, data = tibet.beta.sem1)

#----------------------------
# correct for geographic distance
nut.beta.sem2 <- nut.beta.sem1
m1 <- lm(plant ~ geo, data = nut.beta.sem2)
m2 <- lm(bac ~ geo, data = nut.beta.sem2)
m3 <- lm(fun ~ geo, data = nut.beta.sem2)
m4 <- lm(env ~ geo, data = nut.beta.sem2)
m5 <- lm(emf ~ geo, data = nut.beta.sem2)
nut.beta.sem2$plant <- varScale(residuals(m1))
nut.beta.sem2$bac <- varScale(residuals(m2))
nut.beta.sem2$fun <- varScale(residuals(m3))
nut.beta.sem2$env <- varScale(residuals(m4))
nut.beta.sem2$emf <- varScale(residuals(m5))

MRM(plant ~ env, nperm = 10000, data = nut.beta.sem2)
MRM(bac ~ env + plant, nperm = 10000, data = nut.beta.sem2)
MRM(fun ~ env + plant, nperm = 10000, data = nut.beta.sem2)
MRM(emf ~ env + plant + bac + fun, nperm = 10000, data = nut.beta.sem2)

tibet.beta.sem2 <- tibet.beta.sem1
m1 <- lm(plant ~ geo, data = tibet.beta.sem2)
m2 <- lm(bac ~ geo, data = tibet.beta.sem2)
m3 <- lm(fun ~ geo, data = tibet.beta.sem2)
m4 <- lm(env ~ geo, data = tibet.beta.sem2)
m5 <- lm(emf ~ geo, data = tibet.beta.sem2)
tibet.beta.sem2$plant <- varScale(residuals(m1))
tibet.beta.sem2$bac <- varScale(residuals(m2))
tibet.beta.sem2$fun <- varScale(residuals(m3))
tibet.beta.sem2$env <- varScale(residuals(m4))
tibet.beta.sem2$emf <- varScale(residuals(m5))

MRM(plant ~ env, nperm = 10000, data = tibet.beta.sem2)
MRM(bac ~ env + plant, nperm = 10000, data = tibet.beta.sem2)
MRM(fun ~ env + plant, nperm = 10000, data = tibet.beta.sem2)
MRM(emf ~ env + plant + bac + fun, nperm = 10000, data = tibet.beta.sem2)

#----------------------------
# include richness difference to the path analysis
# by Nutrient Network
nut.plant.rich <- as.vector(beta.div.comp(nut.plant, coef = "S", quant = FALSE)$rich)
nut.bac.rich <- as.vector(beta.div.comp(nut.bac, coef = "S", quant = FALSE)$rich)
nut.fun.rich <- as.vector(beta.div.comp(nut.fun, coef = "S", quant = FALSE)$rich)

nut.beta.sem3 <- nut.beta.sem1
nut.beta.sem3 <- cbind(nut.beta.sem3, nut.plant.rich, nut.bac.rich, nut.fun.rich)
nut.beta.sem3[, c("nut.plant.rich", "nut.bac.rich", "nut.fun.rich")] <- apply(nut.beta.sem3[, c("nut.plant.rich", "nut.bac.rich", "nut.fun.rich")], 2, varScale)

# by tibetan grasslands
tibet.plant.rich <- as.vector(beta.div.comp(tibet.plant, coef = "S", quant = FALSE)$rich)
tibet.bac.rich <- as.vector(beta.div.comp(tibet.bac, coef = "S", quant = FALSE)$rich)
tibet.fun.rich <- as.vector(beta.div.comp(tibet.fun, coef = "S", quant = FALSE)$rich)

tibet.beta.sem3 <- tibet.beta.sem1
tibet.beta.sem3 <- cbind(tibet.beta.sem3, tibet.plant.rich, tibet.bac.rich, tibet.fun.rich)
tibet.beta.sem3[, c("tibet.plant.rich", "tibet.bac.rich", "tibet.fun.rich")] <- apply(tibet.beta.sem3[, c("tibet.plant.rich", "tibet.bac.rich", "tibet.fun.rich")], 2, varScale)

MRM(nut.plant.rich ~ env, nperm = 10000, data = nut.beta.sem3)
MRM(nut.bac.rich ~ env + nut.plant.rich, nperm = 10000, data = nut.beta.sem3)
MRM(nut.fun.rich ~ env + nut.plant.rich, nperm = 10000, data = nut.beta.sem3)
MRM(emf ~ env + nut.plant.rich + nut.bac.rich + nut.fun.rich, nperm = 10000, data = nut.beta.sem3)

MRM(tibet.plant.rich ~ env, nperm = 10000, data = tibet.beta.sem3)
MRM(tibet.bac.rich ~ env + tibet.plant.rich, nperm = 10000, data = tibet.beta.sem3)
MRM(tibet.fun.rich ~ env + tibet.plant.rich, nperm = 10000, data = tibet.beta.sem3)
MRM(emf ~ env + tibet.plant.rich + tibet.bac.rich + tibet.fun.rich, nperm = 10000, data = tibet.beta.sem3)

#----------------------------
# include replacement to the path analysis
# by Nutrient Network
nut.plant.repl <- as.vector(beta.div.comp(nut.plant, coef = "S", quant = FALSE)$repl)
nut.bac.repl <- as.vector(beta.div.comp(nut.bac, coef = "S", quant = FALSE)$repl)
nut.fun.repl <- as.vector(beta.div.comp(nut.fun, coef = "S", quant = FALSE)$repl)

nut.beta.sem4 <- nut.beta.sem1
nut.beta.sem4 <- cbind(nut.beta.sem4, nut.plant.repl, nut.bac.repl, nut.fun.repl)
nut.beta.sem4[, c("nut.plant.repl", "nut.bac.repl", "nut.fun.repl")] <- apply(nut.beta.sem4[, c("nut.plant.repl", "nut.bac.repl", "nut.fun.repl")], 2, varScale)

# by tibetan grasslands
tibet.plant.repl <- as.vector(beta.div.comp(tibet.plant, coef = "S", quant = FALSE)$repl)
tibet.bac.repl <- as.vector(beta.div.comp(tibet.bac, coef = "S", quant = FALSE)$repl)
tibet.fun.repl <- as.vector(beta.div.comp(tibet.fun, coef = "S", quant = FALSE)$repl)

tibet.beta.sem4 <- tibet.beta.sem1
tibet.beta.sem4 <- cbind(tibet.beta.sem4, tibet.plant.repl, tibet.bac.repl, tibet.fun.repl)
tibet.beta.sem4[, c("tibet.plant.repl", "tibet.bac.repl", "tibet.fun.repl")] <- apply(tibet.beta.sem4[, c("tibet.plant.repl", "tibet.bac.repl", "tibet.fun.repl")], 2, varScale)

MRM(nut.plant.repl ~ env, nperm = 10000, data = nut.beta.sem4)
MRM(nut.bac.repl ~ env + nut.plant.repl, nperm = 10000, data = nut.beta.sem4)
MRM(nut.fun.repl ~ env + nut.plant.repl, nperm = 10000, data = nut.beta.sem4)
MRM(emf ~ env + nut.plant.repl + nut.bac.repl + nut.fun.repl, nperm = 10000, data = nut.beta.sem4)

MRM(tibet.plant.repl ~ env, nperm = 10000, data = tibet.beta.sem4)
MRM(tibet.bac.repl ~ env + tibet.plant.repl, nperm = 10000, data = tibet.beta.sem4)
MRM(tibet.fun.repl ~ env + tibet.plant.repl, nperm = 10000, data = tibet.beta.sem4)
MRM(emf ~ env + tibet.plant.repl + tibet.bac.repl + tibet.fun.repl, nperm = 10000, data = tibet.beta.sem4)

# save the data for Glasso analysis
nut.beta2 <- nut.beta %>% 
  filter(q.order == 0) %>% 
  select(plant, bac, fun, emf, geo, env)
nut.beta2 <- cbind(nut.beta2, nut.plant.rich, nut.bac.rich, nut.fun.rich,
                   nut.plant.repl, nut.bac.repl, nut.fun.repl)
write.csv(nut.beta2, "./data/data_processing/nut.beta.comp.csv")

tibet.beta2 <- tibet.beta %>% 
  filter(q.order == 0) %>% 
  select(plant, bac, fun, emf, geo, env)
tibet.beta2 <- cbind(tibet.beta2, tibet.plant.rich, tibet.bac.rich, tibet.fun.rich,
                   tibet.plant.repl, tibet.bac.repl, tibet.fun.repl)
write.csv(tibet.beta2, "./data/data_processing/tibet.beta.comp.csv")

###########################################################
# multiple threshold approach
thresholds <- 99

# create a list of arrays
thresh.array <- NULL
for (i in seq_len(thresholds)) {
  temp <- nut.emf.std
  j <- i/100
  temp[temp >= j] <- 1
  temp[temp < j] <- 0
  thresh.array[[i]] <- temp
}
nut.beta.thresh <- nut.beta %>%
  filter(q.order == 0) %>% 
  select(plant, bac, fun, env)
nut.null <- matrix(NA, nrow = 99, ncol = 8)
for (i in seq_len(thresholds)) {
  emf.temp <- thresh.array[[i]]
  emf.temp <- distance(emf.temp, method = "euclidean")
  dat.temp <- cbind(nut.beta.thresh, as.vector(emf.temp))
  dat.temp <- data.frame(apply(dat.temp, 2, varScale))
  names(dat.temp) <- c("plant", "bac", "fun", "env", "emf")
  mrm.temp <- MRM(emf ~ env + plant + bac + fun, data = dat.temp, nperm = 10000)
  nut.null[i, ] <- cbind(t(mrm.temp$coef[-1, 1]), t(mrm.temp$coef[-1, 2]))
  print(paste(i, "% completed", sep = " "))
}
nut.null <- data.frame(nut.null)
nut.null$thresh <- (seq_len(thresholds)) / 100
names(nut.null)[1:8] <- c("env", "plant", "bac", "fun", 
                          "env.sig", "plant.sig", "bac.sig", "fun.sig")
nut.null.coefs <- nut.null %>% 
  select(plant, bac, fun, thresh) %>% 
  reshape2::melt(id.var = "thresh")
nut.null.sig <- nut.null %>%
  select(plant.sig, bac.sig, fun.sig, thresh) %>% 
  reshape2::melt(id.var = "thresh")
nut.null.df <- data.frame(cbind(nut.null.coefs, nut.null.sig))
nut.null.df$sig <- NA
nut.null.df$sig[nut.null.df$value < 0 & nut.null.df$value.1 <= 0.05] <- "sig.neg"
nut.null.df$sig[nut.null.df$value > 0 & nut.null.df$value.1 <= 0.05] <- "sig.pos"
nut.null.df$sig[is.na(nut.null.df$sig)] <- "neutral"

nut.null.p <- nut.null.df %>% 
  mutate(sig = factor(sig, levels = c("sig.pos", "neutral", "sig.neg"))) %>%
  mutate(variable = factor(variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  ggplot(aes(x = thresh, y = value, color = sig)) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2, col = "gray70") +
  facet_grid(~ variable) +
  scale_color_manual(values = c("#e41a1c", "gray80", "#377eb8")) +
  labs(title = "Nutrient network", x = "Thresholds", 
       y = "Slope estimate") +
  ylim(c(-0.50, 0.70)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')

# create a list of arrays
thresh.array <- NULL
for (i in seq_len(thresholds)) {
  temp <- tibet.emf.std
  j <- i/100
  temp[temp >= j] <- 1
  temp[temp < j] <- 0
  thresh.array[[i]] <- temp
}
tibet.beta.thresh <- tibet.beta %>%
  filter(q.order == 0) %>% 
  select(plant, bac, fun, env)
tibet.null <- matrix(NA, nrow = 99, ncol = 8)
for (i in seq_len(thresholds)) {
  emf.temp <- thresh.array[[i]]
  emf.temp <- distance(emf.temp, method = "euclidean")
  dat.temp <- cbind(tibet.beta.thresh, as.vector(emf.temp))
  dat.temp <- data.frame(apply(dat.temp, 2, varScale))
  names(dat.temp) <- c("plant", "bac", "fun", "env", "emf")
  mrm.temp <- MRM(emf ~ env + plant + bac + fun, data = dat.temp, nperm = 10000)
  tibet.null[i, ] <- cbind(t(mrm.temp$coef[-1, 1]), t(mrm.temp$coef[-1, 2]))
  print(paste(i, "% completed", sep = " "))
}
tibet.null <- data.frame(tibet.null)
tibet.null$thresh <- (seq_len(thresholds)) / 100
names(tibet.null)[1:8] <- c("env", "plant", "bac", "fun", 
                            "env.sig", "plant.sig", "bac.sig", "fun.sig")
tibet.null.coefs <- tibet.null %>% 
  select(plant, bac, fun, thresh) %>% 
  reshape2::melt(id.var = "thresh")
tibet.null.sig <- tibet.null %>%
  select(plant.sig, bac.sig, fun.sig, thresh) %>% 
  reshape2::melt(id.var = "thresh")
tibet.null.df <- data.frame(cbind(tibet.null.coefs, tibet.null.sig))
tibet.null.df$sig <- NA
tibet.null.df$sig[tibet.null.df$value < 0 & tibet.null.df$value.1 <= 0.05] <- "sig.neg"
tibet.null.df$sig[tibet.null.df$value > 0 & tibet.null.df$value.1 <= 0.05] <- "sig.pos"
tibet.null.df$sig[is.na(tibet.null.df$sig)] <- "neutral"

tibet.null.p <- tibet.null.df %>% 
  mutate(sig = factor(sig, levels = c("sig.pos", "neutral", "sig.neg"))) %>%
  mutate(variable = factor(variable, levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  ggplot(aes(x = thresh, y = value, color = sig)) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2, col = "gray70") +
  facet_grid(~ variable) +
  scale_color_manual(values = c("#e41a1c", "gray80", "#377eb8")) +
  labs(title = "Tibetan grasslands", x = "Thresholds", 
       y = "Slope estimate") +
  ylim(c(-0.50, 0.70)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
pdf("./outs/NutTibet_thresh.pdf", width = 9, height = 6)
grid.arrange(nut.null.p, tibet.null.p)
dev.off()

###########################################################
#                    End of Script                        #
###########################################################