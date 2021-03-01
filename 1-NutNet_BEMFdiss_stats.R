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
library(cowplot)

char2seed("multifunctionality")

# load data
df.nut <- read.csv("./data/data_processing/nutnet_distance_matrix.csv")


# location map
# load site coordinate data
nutnet.dat <- read.table("./data/data_processing/nutnet_site_locs.txt")
names(nutnet.dat) <- c("regional", "site", "plot", "longitude", "latitude")
sites <- levels(factor(gsub("_\\d*", "", df.nut$col)))
nutnet.dat <- nutnet.dat[nutnet.dat$site %in% sites, ]
nutnet.dat <- nutnet.dat %>% droplevels() %>% 
  select(site, longitude, latitude) %>% 
  distinct()

# Map figure
mp <- NULL
mapWorld <- borders("world", colour = "gray80", fill = "gray100", lwd = 0.2) # create a layer of borders
mp <- ggplot() + mapWorld
# site.labs <- unique(env.df$site)
x <- nutnet.dat$longitude
y <- nutnet.dat$latitude

#Now Layer the sites on top
mp <- mp +
  geom_point(aes(x = x, y = y), 
             size = 2.5,
             shape = 1) +
  xlab(expression("Longitude ("*degree*")")) + 
  ylab(expression("Lattitude ("*degree*")")) +
  # xlim(c(-200, 200)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))
print(mp)
ggsave("./outs/NutNet_site_maps.pdf", width = 7.5, height = 4.5)


###########################################################
# Part 1: Mantel analysis
###########################################################

# 1.1 sensitivity analysis for beta-diversity index
names(df.nut)
id.vars <- c("plant.sor", "plant.repl", "plant.rich", "plant.q1", "plant.q2",
             "bac.sor", "bac.repl", "bac.rich", "bac.q1", "bac.q2",
             "fun.sor", "fun.repl", "fun.rich", "fun.q1", "fun.q2")
res.beta <- NULL
for (i in id.vars) {
  temp <- df.nut[c("EMF", i)]
  names(temp) <- c("y", "x")
  res <- ecodist::mantel(temp$y ~ temp$x, nperm = 9999, mrank = TRUE)
  res <- data.frame(t(res))
  res$id <- i
  res.beta <- rbind(res.beta, res)
}
res.beta

# 1.2 partial Mantel analysis
names(df.nut)

# Geo + Clim + pH --> beta-diversity
partialMantel1 <- function(data) {
  temp <- data
  names(temp) <- c("y", "geo", "clim", "pH")
  geo1 <- ecodist::mantel(y ~ geo, data = temp, nperm = 9999, mrank = TRUE)
  geo2 <- ecodist::mantel(y ~ geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  geo3 <- ecodist::mantel(y ~ geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo4 <- ecodist::mantel(y ~ geo + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim1 <- ecodist::mantel(y ~ clim, data = temp, nperm = 9999, mrank = TRUE)
  clim2 <- ecodist::mantel(y ~ clim + geo, data = temp, nperm = 9999, mrank = TRUE)
  clim3 <- ecodist::mantel(y ~ clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim4 <- ecodist::mantel(y ~ clim + geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  pH1 <- ecodist::mantel(y ~ pH, data = temp, nperm = 9999, mrank = TRUE)
  pH2 <- ecodist::mantel(y ~ pH + geo, data = temp, nperm = 9999, mrank = TRUE)
  pH3 <- ecodist::mantel(y ~ pH + clim, data = temp, nperm = 9999, mrank = TRUE)
  pH4 <- ecodist::mantel(y ~ pH + geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  res <- rbind(geo1, geo2, geo3, geo4, clim1, clim2, clim3, clim4,
               pH1, pH2, pH3, pH4)
  res <- data.frame(res)
  res$id <- c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
              "pH1", "pH2", "pH3", "pH4")
  return(res)
}

plant.part <- partialMantel1(df.nut[c("plant.sor", "geo", "clim", "pH")])
bac.part <- partialMantel1(df.nut[c("bac.sor", "geo", "clim", "pH")])
fun.part <- partialMantel1(df.nut[c("fun.sor", "geo", "clim", "pH")])

p1 <- rbind(plant.part, bac.part, fun.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("Plant", "Bacteria", "Fungi"), each = 12)) %>% 
  mutate(variable = factor(variable,
                           levels = c("Plant", "Bacteria", "Fungi")),
         id = factor(id,
                     levels = c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
                                "pH1", "pH2", "pH3", "pH4"),
                     labels = c("Geo", "Geo | Clim", "Geo | pH", "Geo | Clim+pH",
                                "Clim", "Clim | Geo", "Clim | pH", "Clim | Geo+pH",
                                "pH", "pH | Geo", "pH | Clim", "pH | Geo+Clim"))) %>% 
  filter(!id %in% c("Geo", "Clim", "pH")) %>% 
  ggplot(aes(forcats::fct_rev(id), mantelr)) +
  geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.), 
                color = "black", alpha = 0.66,
                width = 0.15) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  geom_vline(xintercept = 3.5, color = "gray") +
  geom_vline(xintercept = 6.5, color = "gray") +
  geom_point(size = 3) +
  facet_grid(~ variable) +
  labs(x = NULL, y = "Mantel coefficient (r)") +
  coord_flip() +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
p2 <- df.nut %>% 
  select(plant.sor, bac.sor, fun.sor, geo, clim, pH) %>% 
  pivot_longer(geo:pH, names_to = "abiotic.variable", values_to = "abiotic.value") %>% 
  pivot_longer(plant.sor:fun.sor, names_to = "biotic.variable", values_to = "biotic.value") %>% 
  mutate(biotic.variable = factor(biotic.variable,
                           levels = c("plant.sor", "bac.sor", "fun.sor"),
                           labels = c("Plant", "Bacteria", "Fungi")),
         abiotic.variable = factor(abiotic.variable,
                                   levels = c("geo", "clim", "pH"),
                                   labels = c("Geographic distance", 
                                              "Climatic distance",
                                              "pH distance"))) %>% 
  ggplot(aes(abiotic.value, biotic.value)) +
  geom_point(size = 1.5, color = "grey") +
  geom_smooth(method = "lm", se = FALSE, size = 0.6) +
  facet_grid(biotic.variable ~ abiotic.variable, scales = "free_x") +
  labs(y = expression(paste(beta, "-diversity")),
       x = "Abiotic factors") +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
plot_grid(p2, p1, ncol = 1, labels = c("a)", "b)"), rel_heights = c(1/2, 1/2))
ggsave("./outs/NutNet_beta-diversity.pdf", width = 6.5, height = 9.5)

# beta-diversity --> beta-multifunctionality
# EMF vs. Sorensen
partialMantel2 <- function(data) {
  temp <- data
  names(temp) <- c("y", "x", "geo", "clim", "pH")
  nul <- ecodist::mantel(y ~ x, data = temp, nperm = 9999, mrank = TRUE)
  geo <- ecodist::mantel(y ~ x + geo, data = temp, nperm = 9999, mrank = TRUE)
  clim <- ecodist::mantel(y ~ x + clim, data = temp, nperm = 9999, mrank = TRUE)
  pH <- ecodist::mantel(y ~ x + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo.clim <- ecodist::mantel(y ~ x + geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  geo.pH <- ecodist::mantel(y ~ x + geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim.pH <- ecodist::mantel(y ~ x + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo.clim.pH <- ecodist::mantel(y ~ x + geo + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  res <- rbind(nul, geo, clim, pH, geo.clim, geo.pH, clim.pH, geo.clim.pH)
  res <- data.frame(res)
  res$id <- c("nul", "geo", "clim", "pH", 
              "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH")
  return(res)
}

plant.part <- partialMantel2(df.nut[c("EMF", "plant.sor", "geo", "clim", "pH")])
bac.part <- partialMantel2(df.nut[c("EMF", "bac.sor", "geo", "clim", "pH")])
fun.part <- partialMantel2(df.nut[c("EMF", "fun.sor", "geo", "clim", "pH")])

p3 <- rbind(plant.part, bac.part, fun.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("Plant", "Bacteria", "Fungi"), each = 8)) %>% 
  mutate(variable = factor(variable,
                           levels = c("Plant", "Bacteria", "Fungi")),
         id = factor(id,
                     levels = c("nul", "geo", "clim", "pH", 
                                "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH"),
                     labels = c("X", "X | Geo", "X | Clim", "X | pH", 
                                "X | Geo+Clim", "X | Geo+pH", "X | Clim+pH", "X | Geo+Clim+pH"))) %>% 
  filter(!id %in% c("X")) %>% 
  ggplot(aes(forcats::fct_rev(id), mantelr)) +
  geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.), 
                color = "black", alpha = 0.66,
                width = 0.15) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  geom_point(size = 3) +
  facet_grid(~ variable) +
  labs(x = NULL, y = "Mantel coefficient (r)") +
  coord_flip() +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
p4 <- df.nut %>% 
  select(EMF, SES.sig, plant.sor, bac.sor, fun.sor) %>% 
  pivot_longer(plant.sor:fun.sor, names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable,
                           levels = c("plant.sor", "bac.sor", "fun.sor"),
                           labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  ggplot(aes(value, EMF)) +
  geom_point(aes(color = SES.sig), size = 1.5) +
  scale_color_manual(values = c("#2b83ba", "gray", "#d7191c")) +
  geom_smooth(method = "lm", se = FALSE, size = 0.6) +
  facet_grid(~ variable) +
  labs(x = expression(paste(beta, "-diversity")),
       y = expression(paste(beta, "-multifunctionality"))) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# plot_grid(p4, p3, ncol = 1, labels = c("a)", "b)"), rel_heights = c(2/5, 3/5))
# ggsave("./outs/NutNet_beta-BEMF.pdf", width = 7.5, height = 7.5)


# EMF vs. abiotic factors
mt.geo <- with(df.nut, ecodist::mantel(EMF ~ geo, nperm = 9999, mrank = TRUE))
mt.clim <- with(df.nut, ecodist::mantel(EMF ~ clim, nperm = 9999, mrank = TRUE))
mt.pH <- with(df.nut, ecodist::mantel(EMF ~ pH, nperm = 9999, mrank = TRUE))
mt <- data.frame(rbind(mt.geo, mt.clim, mt.pH))
mt$variable <- rownames(mt)
mt$variable <- factor(mt$variable, 
                      levels = c("mt.geo", "mt.clim", "mt.pH"),
                      labels = c("Geographic distance", 
                                 "Climatic distance", 
                                 "pH distance"))
mt$x <- c((max(df.nut$geo) - min(df.nut$geo)) * 0.35, 
          (max(df.nut$env) - min(df.nut$env)) * 0.35,
          (max(df.nut$pH) - min(df.nut$pH)) * 0.35)
mt$y <- c(1.85, 1.85, 1.85)
mt$labs <- paste("r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")

p5 <- df.nut %>% 
  select(EMF, SES.sig, geo, clim, pH) %>% 
  pivot_longer(geo:pH, names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable,
                           levels = c("geo", "clim", "pH"),
                           labels = c("Geographic distance", 
                                      "Climatic distance", 
                                      "pH distance"))) %>% 
  ggplot(aes(value, EMF)) +
  geom_point(aes(color = SES.sig), size = 1.5) +
  scale_color_manual(values = c("#2b83ba", "gray", "#d7191c")) +
  geom_smooth(method = "lm", se = FALSE, size = 0.6, aes(lty = variable)) +
  scale_linetype_manual(values = c(1, 1, 1)) +
  facet_grid(~ variable, scales = "free_x") +
  geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
  labs(x = "Abiotic factors",
       y = expression(paste(beta, "-multifunctionality"))) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none')
# ggsave("./outs/NutNet_beta-abio-EMF.pdf", width = 7.0, height = 4.0)
plot_grid(p5, p4, p3, ncol = 1, labels = c("a)", "b)", "c)"), 
          rel_heights = c(2/6, 2/6, 4/6))
ggsave("./outs/NutNet_beta-BEMF.pdf", width = 6.5, height = 9.5)


# AGB vs. Sorensen (replace AGB with other functions, PTN, PTP, BGB, STN, avaP)
var4functions <- c("AGB", "PTN", "PTP", "BGB", "STN", "avaP")
res <- NULL
for (i in var4functions) {
  plant.part <- partialMantel2(df.nut[c(i, "plant.sor", "geo", "clim", "pH")])
  bac.part <- partialMantel2(df.nut[c(i, "bac.sor", "geo", "clim", "pH")])
  fun.part <- partialMantel2(df.nut[c(i, "fun.sor", "geo", "clim", "pH")])
  print(i)
  res <- rbind(res, cbind(rbind(plant.part, bac.part, fun.part), i))
}
res$variable <- rep(c("Plant", "Bacteria", "Fungi"), each = 8)
write.csv(res, "./outs/NutNet_single_functions.csv")

# # EMF vs. Sorensen in USA
# df.usa <- df.nut[grep(".us", df.nut$row), ]
# df.usa <- df.usa[grep(".us", df.usa$col), ]
# 
# plant.part <- partialMantel2(df.usa[c("EMF", "plant.sor", "geo", "clim", "pH")])
# bac.part <- partialMantel2(df.usa[c("EMF", "bac.sor", "geo", "clim", "pH")])
# fun.part <- partialMantel2(df.usa[c("EMF", "fun.sor", "geo", "clim", "pH")])
# 
# p1.usa <- rbind(plant.part, bac.part, fun.part) %>%
#   data.frame() %>%
#   mutate(variable = rep(c("Plant", "Bacteria", "Fungi"), each = 8)) %>%
#   mutate(variable = factor(variable,
#                            levels = c("Plant", "Bacteria", "Fungi")),
#          id = factor(id,
#                      levels = c("nul", "geo", "clim", "pH",
#                                 "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH"),
#                      labels = c("X", "X | Geo", "X | Clim", "X | pH",
#                                 "X | Geo+Clim", "X | Geo+pH", "X | Clim+pH", "X | Geo+Clim+pH"))) %>%
#   ggplot(aes(forcats::fct_rev(id), mantelr)) +
#   geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.),
#                 color = "black", alpha = 0.66,
#                 width = 0.15) +
#   geom_hline(yintercept = 0, color = "gray", lty = 2) +
#   geom_point(size = 3) +
#   facet_grid(~ variable) +
#   labs(x = NULL, y = "Mantel coefficient (r)") +
#   coord_flip() +
#   theme_bw(base_size = 13.5) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank())
# p2.usa <- df.usa %>%
#   select(EMF, SES.sig, plant.sor, bac.sor, fun.sor) %>%
#   pivot_longer(plant.sor:fun.sor, names_to = "variable", values_to = "value") %>%
#   mutate(variable = factor(variable,
#                            levels = c("plant.sor", "bac.sor", "fun.sor"),
#                            labels = c("Plant", "Bacteria", "Fungi"))) %>%
#   ggplot(aes(value, EMF)) +
#   geom_point(aes(color = SES.sig), size = 1.5) +
#   scale_color_manual(values = c("#2b83ba", "gray", "#d7191c")) +
#   geom_smooth(method = "lm", se = FALSE, size = 0.6) +
#   facet_grid(~ variable) +
#   labs(x = expression(paste(beta, "-diversity")),
#        y = expression(paste(beta, "-multifunctionality"))) +
#   theme_bw(base_size = 13.5) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = 'none')
# plot_grid(p2, p1, ncol = 1, labels = c("a)", "b)"), rel_heights = c(2/5, 3/5))
# ggsave("./outs/NutNet_beta-BEMF_USA.pdf", width = 7.5, height = 7.5)
# 
# mt.geo <- with(df.usa, ecodist::mantel(EMF ~ geo, nperm = 9999, mrank = TRUE))
# mt.clim <- with(df.usa, ecodist::mantel(EMF ~ clim, nperm = 9999, mrank = TRUE))
# mt.pH <- with(df.usa, ecodist::mantel(EMF ~ pH, nperm = 9999, mrank = TRUE))
# mt <- data.frame(rbind(mt.geo, mt.clim, mt.pH))
# mt$variable <- rownames(mt)
# mt$variable <- factor(mt$variable, 
#                       levels = c("mt.geo", "mt.clim", "mt.pH"),
#                       labels = c("Geographic distance", 
#                                  "Climatic distance", 
#                                  "pH distance"))
# mt$x <- c((max(df.usa$geo) - min(df.usa$geo)) * 0.35, 
#           (max(df.usa$env) - min(df.usa$env)) * 0.35,
#           (max(df.usa$pH) - min(df.usa$pH)) * 0.35)
# mt$y <- c(1.85, 1.85, 1.85)
# mt$labs <- paste("r = ", round(mt$mantelr, 2), ", P = ", round(mt$pval3, 3), sep = "")
# 
# p3.usa <- df.usa %>% 
#   select(EMF, SES.sig, geo, clim, pH) %>% 
#   pivot_longer(geo:pH, names_to = "variable", values_to = "value") %>% 
#   mutate(variable = factor(variable,
#                            levels = c("geo", "clim", "pH"),
#                            labels = c("Geographic distance", 
#                                       "Climatic distance", 
#                                       "pH distance"))) %>% 
#   ggplot(aes(value, EMF)) +
#   geom_point(aes(color = SES.sig), size = 1.5) +
#   scale_color_manual(values = c("#2b83ba", "gray", "#d7191c")) +
#   geom_smooth(method = "lm", se = FALSE, size = 0.6, aes(lty = variable)) +
#   scale_linetype_manual(values = c(1, 1, 2)) +
#   facet_grid(~ variable, scales = "free_x") +
#   geom_text(data = mt, aes(x = x, y = y), label = mt$labs) +
#   labs(x = "Abiotic factors",
#        y = expression(paste(beta, "-multifunctionality"))) +
#   theme_bw(base_size = 13.5) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = 'none')
# plot_grid(p3.usa, p2.usa, p1.usa, ncol = 1, labels = c("a)", "b)", "c)"), 
#           rel_heights = c(2/6, 2/6, 4/6))
# ggsave("./outs/NutNet_beta-BEMF_USA.pdf", width = 6.5, height = 9.5)


###########################################################
# Part 2: Structural equation modelling
###########################################################

library(lavaan)
library(lavaanPlot)

sem.data <- df.nut %>% 
  select(EMF, geo, clim, plant.sor, bac.sor, fun.sor) %>% 
  scale() %>% 
  data.frame()

mod_formula_A <- '
# regression
EMF ~ plant.sor + fun.sor
plant.sor ~ geo + clim
bac.sor ~ geo + clim + plant.sor
fun.sor ~ geo + clim + plant.sor
# residual correlation
geo ~~ clim
bac.sor ~~ fun.sor
'
fit <- sem(mod_formula_A, data = sem.data, se = "boot", bootstrap = 9999)
summary(fit, standardized = T, rsq = T, fit.measures = F)
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
lavaanPlot(model = fit, 
           node_options = list(shape = "box", fontname = "Arial"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE, stars = TRUE)

mod_formula_B <- '
# regression
EMF ~ plant.sor + fun.sor
plant.sor ~ clim + fun.sor
bac.sor ~ geo + clim
fun.sor ~ geo + clim
# residual correlation
geo ~~ clim
bac.sor ~~ fun.sor
EMF 
'
fit <- sem(mod_formula_B, data = sem.data, se = "boot", bootstrap = 9999)
summary(fit, standardized = T, rsq = T, fit.measures = F)
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
lavaanPlot(model = fit, 
           node_options = list(shape = "box", fontname = "Arial"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE, stars = TRUE)


###########################################################
#                End of the Script                        #
###########################################################