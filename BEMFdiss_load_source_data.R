###########################################################
# On controls of Ecosystem beta-mutlifunctionality
# 
# Leading authors: Xin Jing, Aimee Classen, Nathan Sanders
# Data are from Nutrient Network and Tibetan grasslands
# 
# First created by Xin Jing on September 5, 2016
# Last modified by Xin Jing on April 25, 2018
#  
# Contact: Xin Jing <Xin.Jing@uvm.edu>
###########################################################

rm(list = ls())

# load library ---
library(raster)  # Calls: raster, extract etc.
library(data.table)  # Calls: fread, dcast
library(plyr)  # Calls: ddply
library(dplyr)  # Calls: select

# data loading and cleaning ---

# NutNet data
# load data
nut.meta <- read.csv("./data/data_processing/NutNet_metadata_controls_updated_noNN10.csv")
nut.plant <- fread("./data/data_processing/NutNet_species_cover_data.csv", na.strings = 'NULL')
nut.bac <- read.csv("./data/data_processing/NutNet_16S_otu_table_wTax.csv")
nut.fun <- read.csv("./data/data_processing/NutNet_ITS_otu_table_unite75_wTax.csv")
nut.env <- fread("./data/data_processing/NutNet_plot_summary_data.csv", na.strings = "NULL")
nut.NP <- read.csv("./data/data_processing/Xin-NutNet-bulk-nutrients-biomass-28Apr2017.csv")
# clean data
# create a unique id
nut.meta$id <- paste(nut.meta$site_code, nut.meta$plot, sep = "_")
nut.meta.id <- data.frame(id = nut.meta$id, sampleID = nut.meta$sampleID)
nut.env <- data.frame(nut.env)
nut.env$env.id <- paste(nut.env$site_code, nut.env$plot, sep = "_")
# plant
nut.plant <- dcast(nut.plant, site_code + plot ~ Taxon, value.var = 'max_cover')
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
# nut.NP table has the lowest samples
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
# save data for later use
nut.emf <- nut.env %>% 
  rename(AGB = aboveground_live_mass, PTN = PTN, PTP = PTP,
         BGB = belowground_mass, STN = pct_N, STP = ppm_P) %>% 
  mutate(STP = STP / 10000) %>% 
  dplyr::select(site_code, plot, AGB, PTN, PTP, BGB, STN, STP)
nut.emf$BGB[is.na(nut.emf$BGB)] <- mean(nut.emf$BGB, na.rm = TRUE)
nut.abio <- nut.env %>% 
  dplyr::select(site_code, plot, latitude, longitude, MAP, MAT, pH)
# remove unused data and package
rm(list = c("id2", "nut.meta", "nut.meta.id", "nut.NP", "nut.env"))
detach(package:data.table)

# Tibet data
# load data
tibet.plant <- read.csv("./data/data_processing/tibet_plant_species_table_cover.csv", header = TRUE, row.names = 1)
tibet.bac <- read.csv("./data/data_processing/tibet_bacteria_otu_table_raw.csv")
tibet.fun <- read.csv("./data/data_processing/tibet_fungi_otu_table_raw.csv")
tibet.env <- read.csv("./data/data_processing/tibet_env_table.csv")
# clean data
# plant
tibet.plant <- tibet.plant[, -which(names(tibet.plant) %in% paste("unknown", 1:6, sep = ""))]
tibet.plant <- tibet.plant[rowSums(tibet.plant) > 0, ]
tibet.plant <- tibet.plant[, colSums(tibet.plant) > 0]
names(tibet.plant) <- sub("\\.", "_", names(tibet.plant))  # replace '.' by '_'
names(tibet.plant)[which(names(tibet.plant) %in% "Stipa_capillate")] <- "Stipa_capillata"
names(tibet.plant)[which(names(tibet.plant) %in% "Stellera_chamsejasme")] <- "Stellera_chamaejasme"
names(tibet.plant)[which(names(tibet.plant) %in% "Helictotrichon_schellianu")] <- "Helictotrichon_schellianum"
names(tibet.plant)[which(names(tibet.plant) %in% "Pedicularis_cheilanthifol")] <- "Pedicularis_cheilanthifolia"

# bacteria
rownames(tibet.bac) <- tibet.bac$OUT_ID
tibet.bac <- tibet.bac[names(tibet.bac) != "OUT_ID" & names(tibet.bac) != "taxonomy"]
tibet.bac <- t(tibet.bac)
tibet.bac <- data.frame(tibet.bac)
tibet.bac <- tibet.bac[rowSums(tibet.bac) > 0, ]
tibet.bac <- tibet.bac[, colSums(tibet.bac) > 0]
# fungi
rownames(tibet.fun) <- tibet.fun$OUT_ID
tibet.fun <- tibet.fun[names(tibet.fun) != "OUT_ID" & names(tibet.fun) != "taxonomy"]
tibet.fun <- t(tibet.fun)
tibet.fun <- data.frame(tibet.fun)
rownames(tibet.fun) <- rownames(tibet.bac)
tibet.fun <- tibet.fun[rowSums(tibet.fun) > 0, ]
tibet.fun <- tibet.fun[, colSums(tibet.fun) > 0]
# load and extract climate data
MAT <- raster("./data/data_processing/bio_1/w001001.adf")
MAP <- raster("./data/data_processing/bio_12/w001001.adf")
xy.locs <- tibet.env %>% dplyr :: select(longitude, latitude)
coordinates(xy.locs) <- c("longitude", "latitude")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
projection(MAT) <- crs.geo
projection(MAP) <- crs.geo
tibet.env$MAT <- extract(MAT, xy.locs)/10
tibet.env$MAP <- extract(MAP, xy.locs)

# save data for later use
tibet.emf <- tibet.env %>% 
  dplyr::select(site, plot, AGB, PTN, PTP, BGB, STN, STP) %>% 
  mutate(BGB = 0.56 * BGB)
tibet.abio <- tibet.env %>% 
  dplyr::select(site, plot, latitude, longitude, MAP, MAT, pH)
# remove unused data
rm(list = c("crs.geo", "MAT", "MAP", "tibet.env", "xy.locs"))
detach(package:raster)

###########################################################
#                End of the Script                        #
###########################################################