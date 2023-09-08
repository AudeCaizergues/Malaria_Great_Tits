## ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Running diversity analysis ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Prerequisite: 
# This script must be run AFTER running 0_randomSampling.R

# Goal of the script:
# This script runs the "script_diversity_malaria.R" script on 1000 subsampled datasets to account for leucocytozoon strain incertitude.
# Estimates diversity metrics for each sampled site, rank abundance curve and disimilarity indices between sites.

rm(list = ls())

# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(viridis)
library(ggtext)
library(RColorBrewer)
library(cowplot)

# Load subsampled data ---------------------------------------------------------------
data_env <- read.csv("../data/data_env_rank_plot.csv", sep = ";")
dataSets_subsampled <- readRDS("../data/dataSets_subsampled.rds")

# Function to run script on subsampled data ---------------------------------------------------------------

runAnalyseDiversitySubsample <- function(leucoMalaria_df,
                                         PlasmoMalaria_df,
                                         data_env) {
  ## Running analysis ------------------------------------------------------
  source("../scripts_subanalyses/script_diversity_malaria.R",
         local = TRUE)
  
  ## Save results ------------------------------------------------------------
  return(list(
    diversityResults,
    expectedDiversity,
    RA.data,
    brayBinary,
    brayPercent,
    leucoMalaria_df2,
    PlasmoMalaria_df
  ))
}

# Run analysis on subsampled data -----------------------------------------
outputSampling_AnalyseDiversity <-
  apply(
    dataSets_subsampled,
    1,
    FUN = function(x) {
      runAnalyseDiversitySubsample(leucoMalaria_df = x[1][[1]],
                                   PlasmoMalaria_df = x[2][[1]],
                                   data_env = data_env)
    }
  )
# leucoMalaria_df = dataSets_subsampled[[1]][[1]]
# PlasmoMalaria_df = dataSets_subsampled[[2]][[1]]


outputSampling_AnalyseDiversity <-
  do.call(rbind, outputSampling_AnalyseDiversity) %>%  as.data.frame()
colnames(outputSampling_AnalyseDiversity) <-
  c("diversity",
    "expectedDiversity",
    "RA.data",
    "brayBinary",
    "brayPercent",
    "leucoMalaria_df2",
    "PlasmoMalaria_df2")

# Summarise results for diversity indices ---------------------------------

mergedDiversity <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$diversity, function(x) {
            x %>%  mutate(site = rownames(.))
          })) %>%  as.data.frame()

mergedDiversity_summarised <- mergedDiversity %>%
  dplyr::rename(
    abundance = "div2$abundance",
    shannon = "div3$Shannon",
    inverseSimpson = "div4$inverseSimpson"
  ) %>%
  pivot_longer(
    !site,
    names_to = "variable",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  group_by(site, variable) %>%
  summarise(
    mean = mean(value),
    margin = qt(0.92,df=999)*sd(value)/sqrt(1000),
    CIlow = mean(value) - margin,
    CIhigh = mean(value) + margin
  )
mergedDiversity_summarised
richness<-subset(mergedDiversity_summarised, mergedDiversity_summarised$variable=="richness")
shannon<-subset(mergedDiversity_summarised, mergedDiversity_summarised$variable=="shannon")
inverseSimpson<-subset(mergedDiversity_summarised, mergedDiversity_summarised$variable=="inverseSimpson")


plot_richness<- ggplot(richness, aes(x = site,y = mean))+
  geom_point(data = richness, aes(x = site, y = mean),pch = 21, fill = "white", size = 3)+
  geom_errorbar(aes(ymin  = CIlow, ymax  = CIhigh), width = 0.2,size  = 0.7) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.title = element_text(face = "bold"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        #legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93"))
plot_richness

plot_shannon<- ggplot(shannon, aes(x = site,y = mean))+
  geom_point(data = shannon, aes(x = site, y = mean),pch = 21, fill = "white", size = 3)+
  geom_errorbar(aes(ymin  = CIlow, ymax  = CIhigh), width = 0.2,size  = 0.7) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.title = element_text(face = "bold"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        #legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93"))
plot_shannon

plot_inverseSimpson<- ggplot(inverseSimpson, aes(x = site,y = mean))+
  geom_point(data = inverseSimpson, aes(x = site, y = mean),pch = 21, fill = "white", size = 3)+
  geom_errorbar(aes(ymin  = CIlow, ymax  = CIhigh), width = 0.2,size  = 0.7) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.title = element_text(face = "bold"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        #legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93"))
plot_inverseSimpson

mergedExpectedDiversity <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$expectedDiversity, function(x) {
            x %>%  as.data.frame()
          })) %>%  as.data.frame()
mergedExpectedDiversity <- mergedExpectedDiversity %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value),
    margin = qt(0.92,df=999)*sd(value)/sqrt(1000),
    CIlow = mean(value) - margin,
    CIhigh = mean(value) + margin
  )
mergedExpectedDiversity

# Plot rank abundance figure ----------------------------------------------
mergedRank <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$RA.data, function(x) {
            x %>%  as.data.frame()
          })) %>%  as.data.frame()
mergedRank <- mergedRank %>%
  group_by(Grouping, rank) %>%
  summarise(
    abundance = median(abundance),
    abundanceLow = quantile(abundance, 0.025),
    abundanceHigh = quantile(abundance, 0.975)
  )

BioR.theme <- theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_line(colour = "grey97"),
    panel.grid.major = element_line(colour = "grey93"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank()
  )

# Summarise distance matrices ---------------------------------------------
mergedbrayBinary <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$brayBinary, function(x) {
            return(x %>%  as.data.frame())
          })) %>%  as.data.frame()
head(mergedbrayBinary)

mergedbrayPercent <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$brayPercent, function(x) {
            return(x %>%  as.data.frame())
          })) %>%  as.data.frame()

mergedbrayBinary_summarised <- mergedbrayBinary %>%
  group_by(A, B) %>%
  summarise(
    distances = median(value),
    distancesLow = quantile(value, 0.025),
    distancesHigh = quantile(value, 0.975),
    type = "brayBinary"
  ) #%>%
# filter(
#   A != B
# )

mergedbrayPercent_summarised <- mergedbrayPercent %>%
  group_by(A, B) %>%
  summarise(
    distances = median(value),
    distancesLow = quantile(value, 0.025),
    distancesHigh = quantile(value, 0.975),
    type = "brayPercent"
  ) #%>%
# filter(
#   A != B
# )

# Plotting: different colour for the two different matrices --------------------
mergedAll <-
  rbind(mergedbrayBinary_summarised, mergedbrayPercent_summarised) %>% unique()

ggplotMatrixDistance <- ggplot(mergedAll,
                               aes(x = A, y = B, fill = distances)) +
  # The first layer, with its own fill scale
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayBinary", NA, distances
  )),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "viridis",
    limits = c(0, 1),
    direction = -1,
    name = "Bottom: Bray-Curtis\n (composition)",
    na.value = NA
  ) +
  new_scale_fill() +
  # Declare new fill scale for the second layer
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayPercent", NA, distances
  )), #~ subset(.x, type == "OverlapStrainPrevalence"),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "magma",
    limits = c(0, 1),
    direction = -1,
    name = "Top: Bray-Curtis\n (prevalence)",
    na.value = NA
  ) +

  # geom_text(aes(label = paste0(
  #   round(distancesLow, digits = 2),
  #   "-",
  #   round(distancesHigh, digits = 2)
  # )), color = "white", size = 4) +
  #scale_fill_viridis_c(option = "magma",limits=c(0, 1),direction=-1, name=("Top\nOverlap of strain presence\nBottom\nOverlap of strain prevalence"))+
  #scale_x_discrete(labels = c("ROU", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  #scale_y_discrete(limits = c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    text=element_text(family="Times New Roman", size=12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_x_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRAM", "BOT", "FAC", "FONT", "MOS", "MAS")) +
  scale_y_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRAM", "BOT", "FAC", "FONT", "MOS", "MAS"))
ggplotMatrixDistance


ggsave(
  "heatmap_bray.png",
  plot = ggplotMatrixDistance,
  device = "png",
  path = "../Figures",
  scale = 1,
  width = 17,
  height = 11,
  units = "cm")



ggplotMatrixDistance2 <- ggplot(mergedAll,
                               aes(x = A, y = B, fill = distances)) +
  # The first layer, with its own fill scale
  # Declare new fill scale for the second layer
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayBinary", NA, distances
  )), #~ subset(.x, type == "OverlapStrainPrevalence"),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "viridis",
    limits = c(0, 1),
    direction = -1,
    name = "Top: Bray-Curtis\n (prevalence)",
    na.value = NA
  ) +
  new_scale_fill() +
  
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayPercent", NA, distances
  )),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "magma",
    limits = c(0, 1),
    direction = -1,
    name = "Bottom: Bray-Curtis\n (composition)",
    na.value = NA
  ) +

  geom_text(aes(label = paste0(
    round(distancesLow, digits = 2),
    "\n",
    round(distancesHigh, digits = 2)
  )), color = "white", size = 3) +
  #scale_fill_viridis_c(option = "magma",limits=c(0, 1),direction=-1, name=("Top\nOverlap of strain presence\nBottom\nOverlap of strain prevalence"))+
  #scale_x_discrete(labels = c("ROU", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  #scale_y_discrete(limits = c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    text=element_text(family="Times New Roman", size=12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_x_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS")) +
  scale_y_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS"))
ggplotMatrixDistance2


ggsave(
  "heatmap_bray_withCI.png",
  plot = ggplotMatrixDistance2,
  device = "png",
  path = "~/Dropbox/PostDoc_Malaria/Malaria_urbanization/Pipeline_Malaria_03032023/3_results/Figures",
  scale = 1,
  width = 17,
  height = 11,
  units = "cm")

ggplotMatrixDistance <- ggplot(mergedAll,
                               aes(x = A, y = B, fill = distances)) +
  # The first layer, with its own fill scale
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayBinary", NA, distances
  )),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "viridis",
    limits = c(0, 1),
    direction = -1,
    name = "Bottom: Bray-Curtis\n (composition)",
    na.value = NA
  ) +
  new_scale_fill() +
  # Declare new fill scale for the second layer
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayPercent", NA, distances
  )), #~ subset(.x, type == "OverlapStrainPrevalence"),
  aes(fill = distances)) +
  scale_fill_viridis(
    option = "magma",
    limits = c(0, 1),
    direction = -1,
    name = "Top: Bray-Curtis\n (prevalence)",
    na.value = NA
  ) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    text=element_text(family="Times New Roman", size=12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_x_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS")) +
  scale_y_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS"))
ggplotMatrixDistance




# # Plot heatmap test colors ----------------------------------------------
#colourBrewer1 <- colorRampPalette(c(brewer.pal(n = 7, name = "Set3")[2], brewer.pal(n = 7, name = "Set3")[3]))(3)#brewer.pal(n = 9, name = "YlOrRd")
#colourBrewer2 <- colorRampPalette(c(brewer.pal(n = 7, name = "Set3")[2], brewer.pal(n = 7, name = "Set3")[4]))(3)#brewer.pal(n = 9, name = "YlOrRd")
colourBrewer1 <- colorRampPalette(c("#CAB2D6","#6A3D9A"))(3)#brewer.pal(n = 9, name = "YlOrRd")
colourBrewer2 <- colorRampPalette(c("#FB9A99", "#E31A1C"))(3)#brewer.pal(n = 9, name = "YlOrRd")

ggplotMatrixDistance3 <- ggplot(mergedAll,
                               aes(x = A, y = B, fill = distances)) +
  # The first layer, with its own fill scale
 
  # Declare new fill scale for the second layer
  geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayPercent", NA, distances
  )), #~ subset(.x, type == "OverlapStrainPrevalence"),
  aes(fill = distances)) +
  scale_fill_gradientn(
    colours = colourBrewer1,
    limits = c(0.2,0.7),
    na.value = NA,
    name = "Top: Bray-Curtis\n (prevalence)"
  ) +
  new_scale_fill() +
  
   geom_tile(data = mergedAll %>% mutate(distances = ifelse(
    type == "brayBinary", NA, distances
  )),
  aes(fill = distances)) +
  scale_fill_gradientn(
    colours = colourBrewer2,
    limits = c(0.2,0.7),
    na.value = NA,
    name = "Bottom: Bray-Curtis\n (composition)"
  ) +

  geom_text(aes(label = paste0(
    round(distancesLow, digits = 2),
    "\n",
    round(distancesHigh, digits = 2)
  )), color = "white", size = 3) +
  #scale_fill_viridis_c(option = "magma",limits=c(0, 1),direction=-1, name=("Top\nOverlap of strain presence\nBottom\nOverlap of strain prevalence"))+
  #scale_x_discrete(labels = c("ROU", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  #scale_y_discrete(limits = c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    text=element_text(family="Times New Roman", size=12),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    panel.background = element_blank(),
    #panel.border = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_x_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS")) +
  scale_y_discrete(expand = c(0, 0),labels = c("ROU", "ZOO", "CEF", "GRA", "BOT", "FAC", "FON", "MOS", "MAS"))
ggplotMatrixDistance3



# Plotting: same colour for the two different matrices --------------------

mergedData <-
  rbind(mergedbrayBinary_summarised, mergedbrayPercent_summarised) %>%
  mutate(
    distances = ifelse(A == B, NA, distances),
    distancesLow = ifelse(A == B, "", round(distancesLow, digits = 2)),
    distancesHigh = ifelse(A == B, "", round(distancesHigh, digits = 2))
  )

# ggplotMatrixDistance <-
#   ggplot(mergedData, aes(x = A, y = B, fill = distances)) +
#   geom_tile() +
#   scale_x_discrete(limits = c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
#   scale_y_discrete(limits = c("rou", "zoo", "cef", "gram", "bot", "fac", "font", "mos", "mas")) +
#   geom_text(aes(label = paste0(distancesLow, "-", distancesHigh)), color = "white", size = 4) +
#   scale_fill_viridis_c(
#     option = "magma",
#     limits = c(0, ),
#     direction = -1,
#     name = "Top: Bray-Curtis (composition)\nBottom: Bray-Curtis (prevalence)",
#     na.value = "white"
#   ) +
#   xlab("") +
#   ylab("") +
#   theme_bw() +
#   theme(
#     axis.title = element_text(face = "bold"),
#     axis.text = element_text(face = "bold", size = 10),
#     legend.position = "right",
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 14),
#     legend.key = element_blank(),
#     panel.background = element_blank(),
#     #panel.border = element_blank(),
#     panel.grid = element_blank(),
#   ) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0))
# ggplotMatrixDistance3


# Souche prevalence per siter for barplot ---------------------------------
outputSampling_leucoStrain <-
  do.call(rbind, outputSampling_AnalyseDiversity$leucoMalaria_df2) %>% group_by(station, strainID) %>% 
  summarize(prevalenceStrainMoy = sum(prevalenceStrain)/1000)
outputSampling_leucoStrain


strainImportance <- outputSampling_leucoStrain %>% group_by( strainID)  %>%   
  summarize(prevalenceStrainMoy_tot = sum(prevalenceStrainMoy)/length(strainID)) %>% 
  arrange(., prevalenceStrainMoy_tot)
strainImportance


outputSampling_leucoStrain$condition <- outputSampling_leucoStrain$strainID!="PARUS74" & outputSampling_leucoStrain$strainID!="PARUS4" & outputSampling_leucoStrain$strainID!="PARUS72" & outputSampling_leucoStrain$strainID!="PARUS77" & outputSampling_leucoStrain$strainID!="SYBOR22" & outputSampling_leucoStrain$strainID!="PARUS18" & outputSampling_leucoStrain$strainID!="PARUS73" & outputSampling_leucoStrain$strainID!="PYCGPI03" & outputSampling_leucoStrain$strainID!="PARUS16" & outputSampling_leucoStrain$strainID!="PARUS70"

outputSampling_leucoStrain$strainOK = ifelse (outputSampling_leucoStrain$condition, "other", outputSampling_leucoStrain$strainID )

# Palette tests
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n = 12, name = 'Set3')
palette1<- c("#FFFFB3" ,"#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#FFED6F","#BC80BD" ,"#CCEBC5","#D9D9D9" ,"#8DD3C7")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n = 12, name = 'Paired')
palette2<- c("#1F78B4" ,"#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928")


outputSampling_leucoStrain$station <- factor(outputSampling_leucoStrain$station,levels = c("rou", "zoo", "cef","gram", "bot","fac","font","mos","mas"))
outputSampling_leucoStrain$strainOK<-factor(outputSampling_leucoStrain$strainOK, levels = c("PARUS74","PARUS4","PARUS72","PARUS77","SYBOR22","PARUS18","PARUS73","PYCGPI03","PARUS16","PARUS70","other"))

outputSampling_leucoStrain <- outputSampling_leucoStrain %>% group_by(station, strainOK) %>%  summarize(prevalenceStrainMoy_tot = sum(prevalenceStrainMoy))

leuco_strain<- ggplot(outputSampling_leucoStrain, aes(x = station, y=prevalenceStrainMoy_tot, fill = strainOK)) + 
  geom_col(width=0.85,color="black")+
  scale_fill_manual(values = palette1) +
  scale_y_continuous(labels = c("0%","25%","50%","75%","100%")) +
  guides(
    fill = guide_legend(title = "Leucocytozoon \n      lineages"),
    nrow = 5,
    byrow = TRUE
  ) +
  ylab("Relative abundance \n of lineages") +
  xlab("") +
  ggtitle("Leucocytozoon") +
  theme_bw() +
  scale_x_discrete(labels= c("ROU", "ZOO", "CEF","GRA", "BOT","FAC","FON","MOS","MAS"))+
  xlab("Site")+
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.minor = element_line(colour = "grey97"),
    panel.grid.major = element_line(colour = "grey93"),
    text = element_text(family = "Times New Roman")
  )
leuco_strain

#PlasmoMalaria_pos<-subset(PlasmoMalaria, PlasmoMalaria$souche!="negatif")

### renaming "inexploitable" into "unknown"
PlasmoMalaria <- read.csv("../data/data_malaria_adultes_ok_barplot.csv",sep=";")
PlasmoMalaria_pos<-subset(PlasmoMalaria, PlasmoMalaria$souche!="negatif")

PlasmoMalaria_pos$condition2 <- PlasmoMalaria_pos$souche=="inexploitable" 
PlasmoMalaria_pos$souche = ifelse (PlasmoMalaria_pos$condition2, "unknown", PlasmoMalaria_pos$souche )
PlasmoMalaria_pos$station <- factor(PlasmoMalaria_pos$station,levels = c("rou", "zoo", "cef","gram", "bot","fac","font","mos","mas"))
PlasmoMalaria_pos$souche<-factor(PlasmoMalaria_pos$souche, levels=c("SGS1","YWT4","GRW11","DELURB4","AFR065","unknown"))
PlasmoMalaria_pos_ok<-PlasmoMalaria_pos[,c("station","souche")]

Plasmo_strain<- ggplot(na.omit(PlasmoMalaria_pos_ok), aes(x = station, fill = souche)) + 
  geom_bar(position="fill", width=0.85,color="black")+
  scale_fill_manual(values=palette1)+
  scale_y_continuous(labels = scales::percent)+
  guides(fill=guide_legend(title="Plasmodium \n   lineages",ncol=1))+
  ylab("Relative abundance \n of lineages")+
  ggtitle("Plasmodium")+
  xlab("Site")+
  scale_x_discrete(labels= c("ROU", "ZOO", "CEF","GRA", "BOT","FAC","FON","MOS","MAS"))+
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93"),
        text = element_text(family = "Times New Roman")) 
Plasmo_strain

pannel<-plot_grid(Plasmo_strain,leuco_strain,labels = "AUTO", ncol = 1, rel_heights = c(1,1.125),align = "v")
pannel

setwd(
  "../Figures"
)
library(svglite)
save_plot(
  "barplot_souche.png",
  pannel,
  base_height = 6,
  base_asp = 1.5,
  base_width = 7
)

# Summarise rank abundance ------------------------------------------------
head(outputSampling_AnalyseDiversity[[3]])

#Averaging
mergedDiversity <-
  do.call(rbind,
          lapply(outputSampling_AnalyseDiversity$diversity, function(x) {
            x %>%  mutate(site = rownames(.))
          })) %>%  as.data.frame()

mergedRankAbundance <- do.call(rbind,
                               lapply(outputSampling_AnalyseDiversity$RA.data, function(x) {
                                  x %>% 
                                   as.data.frame() %>% 
                                   dplyr::select(Grouping, rank, abundance)
                               }
                               )
                        )

mergedRankAbundanceSummarised <- mergedRankAbundance %>% 
  group_by(Grouping, rank) %>% 
  summarise(
    n = length(rank),
    meanAbundance = mean(abundance),
    lowerAbundance = meanAbundance + qt(0.05/2, df = n-1)*sd(abundance)/sqrt(n-1),
    upperAbundance = meanAbundance - qt(0.05/2, df = n-1)*sd(abundance)/sqrt(n-1)
  )

#Plotting
mergedRankAbundanceSummarised

plotRankAbundance <- ggplot(mergedRankAbundanceSummarised, aes(x = rank, y = meanAbundance, colour = Grouping)) +
  #CI
  geom_ribbon(aes(x = rank, ymin = lowerAbundance, ymax = upperAbundance), alpha = 0.2, colour = NA) + 
  #Mean
  geom_line(position=position_dodge(width=0.4)) +
  geom_point(position=position_dodge(width=0.4)) +
  scale_color_discrete(name = "Site", labels = c("BOT", "CEF", "FAC", "FONT", "GRAM", "MAS", "MOS","ROU","ZOO"))+
  xlab("Lineage rank") +
  ylab("Lineage abundance") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        #legend.title = element_text(face = "bold"),
        #legend.position = "right",
        #legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93"),
        text = element_text(family = "Times New Roman", size =15)
        ) 
plotRankAbundance
