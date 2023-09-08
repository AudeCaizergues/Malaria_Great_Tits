#~~~~~~~~~~~
# DIVERSITY & ABUNDANCE MALARIA URBAN
#~~~~~~~~~~~

# What does this script do?
# This script calculates the different indices of parasitic diversity

# LIBRARIES ---------------------------------------------------------------
library(vegan)
library(BiodiversityR)
library(SpadeR)
library(ggplot2)
library(ggrepel)
library(reshape)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DATA  -------------------------------------------------------------------
# setwd("~/Dropbox/PostDoc_Malaria/Malaria_urbanization/1_data")
# data<-read.csv("1_data/table_diversity.csv",sep=";", row.names=1)#Not needed anymore

#Remove inexploitable and negatif
leucoMalaria_df <- leucoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")

PlasmoMalaria_df <- PlasmoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")

#Sample size:
Nleuco <- leucoMalaria_df %>% 
  group_by(station) %>% 
  summarise(
    N_indiv = length(unique(X))
  )
NPlasmo <- PlasmoMalaria_df %>% 
  group_by(station) %>% 
  summarise(
    N_indiv = length(unique(code))
  )

#Compute diversity tables

#Count data
subdata1 <-
  table(leucoMalaria_df$strainID, leucoMalaria_df$station) %>% as.data.frame.matrix()
subdata2 <-
  table(PlasmoMalaria_df$strainID, PlasmoMalaria_df$station) %>% as.data.frame.matrix()

data <- rbind(subdata1, subdata2)
data[is.na(data)] <- 0
data <- data %>% as.data.frame()
# library(tidyverse)
data2 <-
  cbind(
    data$rou,
    data$zoo,
    data$cef,
    data$gram,
    data$bot,
    data$fac,
    data$font,
    data$mos,
    data$mas
  )
colnames(data2) <-
  c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")
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
data2_prop <-
  cbind(
    data_prop$rou,
    data_prop$zoo,
    data_prop$cef,
    data_prop$gram,
    data_prop$bot,
    data_prop$fac,
    data_prop$font,
    data_prop$mos,
    data_prop$mas
  )
colnames(data2_prop) <-
  c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")
data3_prop <- t(data2_prop)
data4_prop <- as.data.frame(data3_prop)

data_prevalence <- data4_prop

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# VARIOUS DIVERSITY INDEXES -----------------------------------------------

#richness per site
div <- diversityresult(
  data3,
  y = NULL,
  factor = NULL,
  level = NULL,
  index = c("richness"),
  method = c("each site"),
  sortit = FALSE,
  digits = 3
)

div2 <- diversityresult(
  data3,
  y = NULL,
  factor = NULL,
  level = NULL,
  index = c("abundance"),
  method = c("each site"),
  sortit = FALSE,
  digits = 3
)

div3 <- diversityresult(
  data3,
  y = NULL,
  factor = NULL,
  level = NULL,
  index = c("Shannon"),
  method = c("each site"),
  sortit = FALSE,
  digits = 3
)

div4 <- diversityresult(
  data3,
  y = NULL,
  factor = NULL,
  level = NULL,
  index = c("inverseSimpson"),
  method = c("each site"),
  sortit = FALSE,
  digits = 3
)

# div5 <- diversityresult(
#   data2,
#   y = NULL,
#   factor = NULL,
#   level = NULL,
#   index = c("Logalpha"),
#   #Error here in previous script, was written logalpha
#   method = c("each site"),
#   sortit = FALSE,
#   digits = 3
# )

diversityResults <-
  cbind(div,
        div2$abundance,
        div3$Shannon,
        div4$inverseSimpson)
rownames(div)

### Expected diversity
expectedDiversity1 <-
  diversityresult(
    data3,
    y = NULL,
    factor = NULL,
    level = NULL,
    index = c("jack1"),
    method = c("pooled"),
    sortit = FALSE,
    digits = 8
  )

expectedDiversity2 <-
  diversityresult(
    data3,
    y = NULL,
    factor = NULL,
    level = NULL,
    index = c("chao"),
    method = c("pooled"),
    sortit = FALSE,
    digits = 8
  )

expectedDiversity3 <-
  diversityresult(
    data3,
    y = NULL,
    factor = NULL,
    level = NULL,
    index = c("boot"),
    method = c("pooled"),
    sortit = FALSE,
    digits = 8
  )
expectedDiversity <-
  c(expectedDiversity1, expectedDiversity2, expectedDiversity3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RANK ABUNDANCE CURVES ---------------------------------------------------
count_leuco_rur <- as.data.frame(data2[, c(8)])
count_leuco_urb <- data2[, c(1, 2, 3, 4, 5, 6, 7, 9)]
data_env$site <- as.factor(data_env$site)
data_env$habitat <- as.factor(data_env$habitat)

data_env <- data_env %>% arrange(site)
data2 <- data2 %>% as.data.frame() %>% dplyr::select(order(colnames(.)))

RA.data <- rankabuncomp(
  t(data2) %>%  as.data.frame(),
  y = data_env,
  factor = 'site',
  return.data = TRUE,
  legend = FALSE
)#Put legend = FALSE because with TRUE takes for ever, why?

BioR.theme <- theme(
 panel.background = element_blank(),
 panel.border = element_blank(),
 #panel.grid = element_blank(),
 axis.line = element_line("gray25"),
 #text = element_text(size = 12, family="Arial"),
 axis.text = element_text(size = 10, colour = "gray25"),
 axis.title = element_text(size = 14, colour = "gray25"),
 legend.title = element_text(size = 14),
 legend.text = element_text(size = 14),
 legend.key = element_blank()
)

# BRAY-CURTIS INDEX ------------------------------------------------------

# Bray-Curtis difference matrix: binary
brayBinary <-
  vegdist(data_binary, method = "bray", binary = TRUE) %>% as.matrix()
brayBinary[lower.tri(brayBinary,diag=TRUE)] <- NA

brayBinary <- melt(brayBinary, na.rm = FALSE) %>%
  filter(!is.na(value)) %>%
  dplyr::rename(A = "X1",
                B = "X2")
brayBinary

# Bray-Curtis difference matrix: prevalence
brayPercent <-
  vegdist(data_prevalence, method = "bray", binary = FALSE) %>% as.matrix()
brayPercent[upper.tri(brayPercent,diag=TRUE)] <- NA

brayPercent <- melt(brayPercent, na.rm = FALSE) %>%
  filter(!is.na(value)) %>%
  dplyr::rename(A = "X1",
                B = "X2")
brayPercent


# Save prevalence of each strain ------------------------------------------

#leucoMalaria_df = dataSets_subsampled[[1]][[1]]
#PlasmoMalaria_df = dataSets_subsampled[[2]][[1]]
  
leucoMalaria_df <- leucoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")
PlasmoMalaria_df <- PlasmoMalaria_df %>%
  filter(strainID != "inexploitable" & strainID != "negatif")


leucoMalaria_df2<-leucoMalaria_df %>% 
  group_by(station) %>% 
  mutate(nStrain = length(strainID)) %>% 
  group_by(station,strainID) %>% 
  summarize(prevalenceStrain = length(strainID)/nStrain) %>% 
  mutate(type= "leuco") %>% unique()


PlasmoMalaria_df2<-PlasmoMalaria_df %>% 
  group_by(station) %>% 
  mutate(nStrain = length(strainID)) %>% 
  group_by(station,strainID) %>% 
  summarize(prevalenceStrain = length(strainID)/nStrain) %>% 
  mutate(type= "plasmo") %>% unique()



