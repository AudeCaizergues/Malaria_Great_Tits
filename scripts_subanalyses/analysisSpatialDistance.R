#~~~~~~~~~~~~~~~
# Linking patterns of diversity to spatial organisation
#~~~~~~~~~~~~~~~

# What does this script do?

# This script computes several analyses/data processing/figures for investigating the pattern of parasitology with Malaria germs in great tits of the long-term study based at the CEFE, Montpellier.
# # In particular you will find three main parts:
# 1) The plot of the distribution of sampling sites. Two key things are missing: the exact nest locations, as well as one site ("La rouviere"), a little urbanised study site excentred from the main city (Montpellier). Inserts must be added (one for la rouviere, one for the city) on the map of France to localise this site.
# 2) Analysing whether some strains are actually more linked to urbanised or less urbanised environments
# 3) Analysing whether there is a "contagious" pattern, at the site or nest level: closer sites/nests are more likely (or not!) to host similar strains.

# Libraries ---------------------------------------------------------------

library(readr)
library(dplyr)
library(sf)
library(vegan)
library(ape)

# Data --------------------------------------------------------------------

#Remove inexploitable and negatif
leucoMalaria_df <- leucoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")

PlasmoMalaria_df <- PlasmoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")

#Some summary data
souche_df <- leucoMalaria_df %>% 
  dplyr::select(station, 
                strainID)

urbanisationLevel_df <- leucoMalaria_df %>% 
  select(station, average_urb_level) %>% 
  unique()

#Sample size:
Nleuco <- leucoMalaria_df %>% 
  group_by(station) %>% 
  summarise(
    N_indiv = length(unique(CODE.MIVEGEC))
  )
NPlasmo <- PlasmoMalaria_df %>% 
  group_by(station) %>% 
  summarise(
    N_indiv = length(unique(Code.MIVEGEC))
  )

# Strain overlap

#Count data
subdata1 <-
  table(leucoMalaria_df$strainID, leucoMalaria_df$station) %>% as.data.frame.matrix()
subdata2 <-
  table(PlasmoMalaria_df$strainID, PlasmoMalaria_df$station) %>% as.data.frame.matrix()

data <- rbind(subdata1, subdata2)
data[is.na(data)] <- 0
data <- data %>% as.data.frame()
# library(tidyverse)
data2 <- data
data3 <- t(data2)
data4 <- as.data.frame(data3)

#Get only binary data: if strain is present or not
data_binary <- data4
data_binary[data_binary > 0] <- 1

#Prop of individuals infected
subdata1_prop <- t(apply(subdata1, 1, function(x){x/Nleuco$N_indiv}))#Are naturally in the alphabetical order for both
subdata2_prop <- t(apply(subdata2, 1, function(x){x/NPlasmo$N_indiv}))

data_prop <- rbind(subdata1_prop, subdata2_prop)
data_prop[is.na(data_prop)] <- 0
data_prop <- data_prop %>% as.data.frame()
# library(tidyverse)
data2_prop <- data_prop
data3_prop <- t(data2_prop)
data4_prop <- as.data.frame(data3_prop)

data_percent <- data4_prop

brayBinary <- vegdist(data_binary, method = "bray", binary = TRUE) %>% as.matrix()
brayPercent <- vegdist(data_percent, method = "bray", binary = TRUE) %>% as.matrix()
diag(brayBinary) <- NA
diag(brayPercent) <- NA

strainNest_leuco <- cbind(leucoMalaria_df$strainID, leucoMalaria_df$stanic)
strainNest_plasmo <- cbind(PlasmoMalaria_df$strainID, PlasmoMalaria_df$stanic)
colnames(strainNest_leuco) <- c("strain", "stanic")
colnames(strainNest_plasmo) <- c("strain", "stanic")

strainNest <- rbind(strainNest_leuco, strainNest_plasmo) %>% as.data.frame()
countStrainNest <- table(strainNest$strain, strainNest$stanic)
countStrainNest[countStrainNest > 0] <- 1

## Import polygon data ------------------------------------------------------
polygonsToAdd_sf <- st_read("../data/nestboxes_data/nestboxes_zones.kml")
polygonsToAdd_sf$Name <- tolower(polygonsToAdd_sf$Name)
polygonsToAdd_sf$description <- 1:nrow(polygonsToAdd_sf)

# Add value urbanisation
polygonsToAdd_sf <- left_join(polygonsToAdd_sf,
                              urbanisationLevel_df,
                              by = c("Name" = "station"))

# Calculate the centroid for labelling
centroid_sf <- st_centroid(polygonsToAdd_sf)

## Strain overlap and distance between sites: do they relate? -------------

matrixOverlapstrain <- brayBinary
matrixOverlapstrainPercent <- brayPercent

# Matrix of distance between sites
matrixOfDistances <- st_distance(polygonsToAdd_sf %>% 
                                   arrange(Name) %>% 
                                   sf::st_transform(paste("+proj=utm +zone=",31," +ellps=WGS84", " +datum=WGS84", " +units=m", " +towgs84:0,0,0", sep=""))
)#I reproject the data and calculate the distances
diag(matrixOfDistances) <- NA
matrixOfDistances <- matrix(matrixOfDistances, nrow = length(unique(leucoMalaria_df$station)))
colnames(matrixOfDistances) <- (polygonsToAdd_sf %>% arrange(Name))$Name
rownames(matrixOfDistances) <- colnames(matrixOfDistances) 

# Compare correlation of matrices with mantel test

mantelTestResultsSiteStrainId <- mantel.test(matrixOfDistances, matrixOverlapstrain, nperm = 999)
mantelTestResultsSiteStrainId

mantelTestResultsSitePercent <- mantel.test(matrixOfDistances, matrixOverlapstrainPercent, nperm = 999)
mantelTestResultsSitePercent

## Strain overlap and distance between nests: do they relate? -------------

#Import nest data
nestLocation_sf <- read_csv("../data/nestboxes_data/Montpellier_NB.csv") %>%
  mutate(
    NestboxID = gsub(" ", "_", NestboxID),
    NestboxID = gsub("_NB", "", NestboxID),
    NestboxID = gsub("_", ".", NestboxID)
  ) %>% 
  dplyr::filter(!is.na(Longitude) & !is.na(Longitude)) %>% 
  dplyr::select(Longitude, Latitude, NestboxID) %>% 
  unique() %>% 
  arrange(NestboxID) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"))
st_crs(nestLocation_sf) <- st_crs(polygonsToAdd_sf)

#Add missing nests
nestLocation_missing_sf <- read_csv("../data/nestboxes_data/netsboxes_extra.csv") %>%
  dplyr::rename(
    NestboxID = "nichoir"
  ) %>% 
  dplyr::filter(!is.na(Longitude) & !is.na(Longitude)) %>% 
  dplyr::select(Longitude, Latitude, NestboxID) %>% 
  unique() %>% 
  arrange(NestboxID) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"))
st_crs(nestLocation_missing_sf) <- st_crs(polygonsToAdd_sf)

toAdd <- nestLocation_sf[nestLocation_sf$NestboxID == "fac.22",]
toAdd$NestboxID <- "fac.23" #From Aude, same loc.

nestLocation_sf <- rbind(
  nestLocation_sf,
  nestLocation_missing_sf,
  toAdd
)
## First summarise strains for each nest

countStrainNest <- countStrainNest[, which(colnames(countStrainNest) %in% nestLocation_sf$NestboxID)]
brayBinary_nest <- vegdist(t(countStrainNest), method = "bray", binary = TRUE) %>% as.matrix()
matrixNestOverlapstrain_nest <- brayBinary_nest

## Spatial
colnames(countStrainNest)[which(!(colnames(countStrainNest) %in% nestLocation_sf$NestboxID))]

matrixNestOfDistances <- st_distance(nestLocation_sf %>% 
                                       arrange(NestboxID) %>% 
                                       filter(NestboxID %in% colnames(countStrainNest)) %>% 
                                       sf::st_transform(paste("+proj=utm +zone=",31," +ellps=WGS84", " +datum=WGS84", " +units=m", " +towgs84:0,0,0", sep=""))
)#I reproject the data and calculate the distances
diag(matrixNestOfDistances) <- NA
matrixNestOfDistances <- matrix(matrixNestOfDistances, nrow = ncol(countStrainNest))

colnames(matrixNestOfDistances) <- (nestLocation_sf %>% arrange(NestboxID) %>% filter(NestboxID %in% colnames(countStrainNest)))$NestboxID 
rownames(matrixNestOfDistances) <- colnames(matrixNestOfDistances) 

# Compare correlation of matrices with mantel test

mantelTestResultsNest <- mantel.test(matrixNestOfDistances, matrixNestOverlapstrain_nest, nperm = 999)
mantelTestResultsNest

# Strain overlap and level of urbanisation: do they related? --------------

#Site level
## Calculate distance matrix urbanisation
urbanisation_stations <- leucoMalaria_df %>% 
  dplyr::select(station, average_urb_level) %>% 
  arrange(station) %>% 
  unique() %>% 
  mutate(
    average_urb_level = gsub(",", ".", average_urb_level)
  )
site_names <- urbanisation_stations[,1]
urbanisation_stations <- as.numeric(as.character(urbanisation_stations$average_urb_level)) %>% as.matrix()
rownames(urbanisation_stations) <- as.character(site_names$station) 
matrixOfDistancesUrbanisationSite <- vegdist(urbanisation_stations, method = "manhattan", binary = FALSE) %>% as.matrix()
diag(matrixOfDistancesUrbanisationSite) <- NA
##Run test correlation
mantelTestResultsSiteStrainId_urbanisation <- mantel.test(matrixOfDistancesUrbanisationSite, matrixOverlapstrain, nperm = 999)
mantelTestResultsSiteStrainId_urbanisation
mantelTestResultsSitePercent_urbanisation <- mantel.test(matrixOfDistancesUrbanisationSite, matrixOverlapstrainPercent, nperm = 999)
mantelTestResultsSitePercent_urbanisation

#Nest level
## Calculate distance matrix urbanisation
urbanisation_stanic <- PlasmoMalaria_df %>% #Plasmo ahs the largest stanic nb
  dplyr::select(stanic, urb_level) %>% 
  arrange(stanic) %>% 
  unique() %>% 
  mutate(
    urb_level = gsub(",", ".", urb_level)
  )
site_names <- urbanisation_stanic[,1]
urbanisation_stanic <- as.numeric(as.character(urbanisation_stanic$urb_level)) %>% as.matrix()
rownames(urbanisation_stanic) <- as.character(site_names$stanic)
urbanisation_stanic <- urbanisation_stanic[which(rownames(urbanisation_stanic) %in% nestLocation_sf$NestboxID),] %>%  as.matrix()
ncol(urbanisation_stanic)

matrixOfDistancesUrbanisationNest <- vegdist(urbanisation_stanic, method = "manhattan", binary = FALSE) %>% as.matrix()
diag(matrixOfDistancesUrbanisationNest) <- NA

##Run test correlation
mantelTestResultsNest_urbanisation <- mantel.test(matrixOfDistancesUrbanisationNest, matrixNestOverlapstrain_nest, nperm = 999)
mantelTestResultsNest_urbanisation