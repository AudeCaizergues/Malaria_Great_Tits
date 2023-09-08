#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~        ANALYSES MALARIA STRAINS ~ JURBANISATION LEVEL     #~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#What does this script?
#This script aims to assess whether some strains specific to urbanised vs less urbanised areas? 

# Methods

# To investigate whether one strain occurred more frequently than randomly expected in urban vs rural environment, 
# we compared the proportion of sites at which the strain was present that were urban to the proportion of urban sites overall using a binomial test ("binon.test" function). 
# We computed the test only for strains for which the type II error was below 0.20.

# Libraries ---------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(cowplot)

# Analyses ----------------------------------------------------------------

## LEUCOCYTOZOON  ----------------------------------------------------------

strainSummary_df <- leucoMalaria_df %>% 
  mutate(
    numberNestTot = length(unique(stanic)),
    numberNestUrb = length(unique(stanic[habitat == "urb"]))
  ) %>% 
  select(stanic, 
         strainID, 
         habitat, 
         numberNestTot,
         numberNestUrb
  ) %>% 
  group_by(strainID) %>% 
  unique() %>% 
  summarise(
    numberNestTot = numberNestTot,
    numberNestUrb = numberNestUrb,
    numberNestUrbWithStrain = length(unique(stanic[habitat == "urb"])),
    numberNestWithStrain = length(unique(stanic)),
    criticalAreaBinomTestLow = min(qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot), numberNestWithStrain - qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot)),
    criticalAreaBinomTestHigh = max(qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot), numberNestWithStrain - qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot)),
    typeII = 1 - (pbinom(criticalAreaBinomTestHigh, numberNestWithStrain, numberNestUrb/numberNestTot) - pbinom(criticalAreaBinomTestLow, numberNestWithStrain, numberNestUrb/numberNestTot))
  ) %>% 
  ungroup() %>% 
  unique() 

### Run binomial test for strains with sufficient sampling  -----------------

#### Calculate the minimum observation to have to see difference  ------------

strainSummary_df_rdc <- strainSummary_df %>% 
  filter(
    typeII >= 0.8 & # Filter when power is insufficient
      strainID != "inexploitable" &
      strainID != "" &
      strainID != "negatif") %>% 
  filter( #Remove if 100% power because of low sample
    typeII != 1 | numberNestWithStrain >= 10
  )

strainSummary_df_rdc <- as.data.frame(strainSummary_df_rdc)

#### Compute binomial test only when sufficient power  -----------------------

resultsEstimatedParametersBinomialTest_summary <- apply(strainSummary_df_rdc, 1, function(x){
  testBinom80 <- binom.test(x = as.numeric(x[4]), 
                            n = as.numeric(x[5]), 
                            p = as.numeric(x[3])/as.numeric(x[2]),
                            conf.level = 0.80
  )
  testBinom95 <- binom.test(x = as.numeric(x[4]), 
                            n = as.numeric(x[5]), 
                            p = as.numeric(x[3])/as.numeric(x[2]),
                            conf.level = 0.95
  )
  
  output <- c(
    x[1],
    testBinom80$estimate,
    testBinom80$'conf.int',
    testBinom80$'p.value',
    testBinom95$'conf.int',
    testBinom95$'p.value'
  )
  return(output)
}
)
resultsEstimatedParametersBinomialTest_summary <- as.data.frame(t(resultsEstimatedParametersBinomialTest_summary))
colnames(resultsEstimatedParametersBinomialTest_summary) <-
  c("Strain", "Estimate", "CILow80", "CIHigh80", "p80", "CILow95", "CIHigh95", "p95")

resultsEstimatedParametersBinomialTest_summary[, 2:ncol(resultsEstimatedParametersBinomialTest_summary)] <- apply(resultsEstimatedParametersBinomialTest_summary[, 2:ncol(resultsEstimatedParametersBinomialTest_summary)], 2, function(x){
  x <- as.numeric(x)
  return(x)
})
resultsEstimatedParametersBinomialTest_summary <-
  resultsEstimatedParametersBinomialTest_summary %>% 
  arrange(Estimate) %>% 
  group_by(Strain) %>% 
  mutate(
    colour80 = if(strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] < CILow80 | strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] > CIHigh80){
      if(Estimate > strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]){
        "darkgrey"
      }else{
        "darkolivegreen3"
      }}else{
        "black"
      },
    colour95 = if(strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] < CILow95 | strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] > CIHigh95){
      if(Estimate > strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]){
        "darkgrey"
      }else{
        "darkolivegreen3"
      }}else{
        "black"
      }
  )

### Reorder to have increased occurrence prob
resultsEstimatedParametersBinomialTest_summary$Strain <- factor(resultsEstimatedParametersBinomialTest_summary$Strain,
                                                                levels=
                                                                  resultsEstimatedParametersBinomialTest_summary$Strain)

### Plot results

# LEGEND:
# The mean probability of occurrence of a given strain in a urban environment is depicted by the black dots.
# The vertical segments depict the 95% (dashed) and 80% (plain) confidence intervals.
# The expected probability (number of urban nests over total number of nests monitored) is highlighted by the horizontal black line. 
# When the confidence interval does not intersect this line, this is highlighted by coloured dots and segments 
# (green if below, which means that rural habitats preponderantly host the strain, or grey, which means that urban habitats 
# preponderantly host the strain).

uniqueValuesColours <- unique(c(resultsEstimatedParametersBinomialTest_summary$colour95, resultsEstimatedParametersBinomialTest_summary$colour80))
uniqueLabelsColours <- ifelse(
  uniqueValuesColours == "darkolivegreen3", 
  "Less urban",
  "Unspecific"
)
uniqueLabelsColours <- ifelse(
  uniqueValuesColours == "darkgrey", 
  "More urban",
  uniqueLabelsColours
)

# plotBinomialTest <- ggplot(resultsEstimatedParametersBinomialTest_summary,                ### The data frame to use.
#                            aes(x = Strain,
#                                y = Estimate)) +
#   # CI
#   geom_errorbar(aes(ymin  = CILow95,
#                     ymax  = CIHigh95, 
#                     color = colour95),
#                 width = 0.2,
#                 size  = 0.35) +
#   geom_errorbar(aes(ymin  = CILow80,
#                     ymax  = CIHigh80, 
#                     color = colour80),
#                 width = 0,
#                 size  = 1) +
#   # Expected probability
#   geom_hline(yintercept = strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2], linetype = "dashed") +
#   # Estimate
#   geom_point(data = resultsEstimatedParametersBinomialTest_summary, aes(x = Strain, y = Estimate, color = colour80),
#              pch = 21,
#              fill = "white",
#              size = 2.5) +
#   scale_colour_manual(breaks = uniqueValuesColours,
#                       values = uniqueValuesColours,
#                       labels = uniqueLabelsColours,
#                       name = "Habitat specificity") +
#   ggtitle("Leucocytozoon") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         legend.title = element_text(face = "bold"),
#         legend.position = "bottom",
#         legend.box="vertical",
#         panel.grid.minor=element_line(colour="grey97"),
#         panel.grid.major=element_line(colour="grey93")) +
#   ylim(0,1) +
#   ylab("Occurrence probability") +
#   xlab("Lineage")
# plotBinomialTest

resultsEstimatedParametersBinomialTest_summary_LEUCO <- resultsEstimatedParametersBinomialTest_summary
valueExpectedLeuco <- strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]

## PLASMODIUM  -------------------------------------------------------------

strainSummary_df <- PlasmoMalaria_df %>% 
  mutate(
    numberNestTot = length(unique(stanic)),
    numberNestUrb = length(unique(stanic[habitat == "urb"]))
  ) %>% 
  select(stanic, 
         strainID, 
         habitat, 
         numberNestTot,
         numberNestUrb
  ) %>% 
  group_by(strainID) %>% 
  unique() %>% 
  summarise(
    numberNestTot = numberNestTot,
    numberNestUrb = numberNestUrb,
    numberNestUrbWithStrain = length(unique(stanic[habitat == "urb"])),
    numberNestWithStrain = length(unique(stanic)),
    criticalAreaBinomTestLow = min(qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot), numberNestWithStrain - qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot)),
    criticalAreaBinomTestHigh = max(qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot), numberNestWithStrain - qbinom(.025, numberNestWithStrain, numberNestUrb/numberNestTot)),
    typeII = 1 - (pbinom(criticalAreaBinomTestHigh, numberNestWithStrain, numberNestUrb/numberNestTot) - pbinom(criticalAreaBinomTestLow, numberNestWithStrain, numberNestUrb/numberNestTot))
  ) %>% 
  ungroup() %>% 
  unique() 

### Run binomial test for strains with sufficient sampling  -----------------

#### Calculate the minimum observation to have to see difference  ------------

strainSummary_df_rdc <- strainSummary_df %>% 
  filter(
    typeII >= 0.8 & # Filter when power is insufficient
      strainID != "inexploitable" &
      strainID != "negatif") %>% 
  filter( #Remove if 100% power because of low sample
    typeII != 1 | numberNestWithStrain >= 10
  )
strainSummary_df_rdc <- as.data.frame(strainSummary_df_rdc)

#### Compute binomial test only when sufficient power  -----------------------

resultsEstimatedParametersBinomialTest_summary <- apply(strainSummary_df_rdc, 1, function(x){
  testBinom80 <- binom.test(x = as.numeric(x[4]), 
                            n = as.numeric(x[5]), 
                            p = as.numeric(x[3])/as.numeric(x[2]),
                            conf.level = 0.80
  )
  testBinom95 <- binom.test(x = as.numeric(x[4]), 
                            n = as.numeric(x[5]), 
                            p = as.numeric(x[3])/as.numeric(x[2]),
                            conf.level = 0.95
  )
  
  output <- c(
    x[1],
    testBinom80$estimate,
    testBinom80$'conf.int',
    testBinom80$'p.value',
    testBinom95$'conf.int',
    testBinom95$'p.value'
  )
  return(output)
}
)
resultsEstimatedParametersBinomialTest_summary <- as.data.frame(t(resultsEstimatedParametersBinomialTest_summary))
colnames(resultsEstimatedParametersBinomialTest_summary) <-
  c("Strain", "Estimate", "CILow80", "CIHigh80", "p80", "CILow95", "CIHigh95", "p95")

resultsEstimatedParametersBinomialTest_summary[, 2:ncol(resultsEstimatedParametersBinomialTest_summary)] <- apply(resultsEstimatedParametersBinomialTest_summary[, 2:ncol(resultsEstimatedParametersBinomialTest_summary)], 2, function(x){
  x <- as.numeric(x)
  return(x)
})
resultsEstimatedParametersBinomialTest_summary <-
  resultsEstimatedParametersBinomialTest_summary %>% 
  arrange(Estimate) %>% 
  group_by(Strain) %>% 
  mutate(
    colour80 = if(strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] < CILow80 | strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] > CIHigh80){
      if(Estimate > strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]){
        "darkgrey"
      }else{
        "darkolivegreen3"
      }}else{
        "black"
      },
    colour95 = if(strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] < CILow95 | strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2] > CIHigh95){
      if(Estimate > strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]){
        "darkgrey"
      }else{
        "darkolivegreen3"
      }}else{
        "black"
      }
  )

## Reorder to have increased occurrence prob -----------------
resultsEstimatedParametersBinomialTest_summary$Strain <- factor(resultsEstimatedParametersBinomialTest_summary$Strain,
                                                                levels=
                                                                  resultsEstimatedParametersBinomialTest_summary$Strain)

resultsEstimatedParametersBinomialTest_summary_PLASMO <- resultsEstimatedParametersBinomialTest_summary
valueExpectedPlasmo <- strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2]

# Plot results ------------------------------------------------------------

# LEGEND:
# The mean probability of occurrence of a given strain in a urban environment is depicted by the black dots.
# The vertical segments depict the 95% (dashed) and 80% (plain) confidence intervals.
# The expected probability (number of urban nests over total number of nests monitored) is highlighted by the horizontal black line. 
# When the confidence interval does not intersect this line, this is highlighted by coloured dots and segments 
# (green if below, which means that rural habitats preponderantly host the strain, or grey, which means that urban habitats 
# preponderantly host the strain).

uniqueValuesColours <- unique(c(resultsEstimatedParametersBinomialTest_summary$colour95, resultsEstimatedParametersBinomialTest_summary$colour80))
uniqueLabelsColours <- ifelse(
  uniqueValuesColours == "darkolivegreen3", 
  "Less urban",
  "Unspecific"
)
uniqueLabelsColours <- ifelse(
  uniqueValuesColours == "darkgrey", 
  "More urban",
  uniqueLabelsColours
)

plotBinomialTest2 <- ggplot(resultsEstimatedParametersBinomialTest_summary,                ### The data frame to use.
                            aes(x = Strain,
                                y = Estimate)) +
  # Estimate
  geom_point(data = resultsEstimatedParametersBinomialTest_summary, aes(x = Strain, y = Estimate, color = colour80),
             pch = 19,
             size = 2.5) +
  # CI
  geom_errorbar(aes(ymin  = CILow95,
                    ymax  = CIHigh95,
                    color = colour95),
                width = 0.2,
                size  = 0.7) +
  # Expected probability
  geom_hline(yintercept = strainSummary_df_rdc[1,3]/strainSummary_df_rdc[1,2], linetype = "dashed") +
  # Estimate
  geom_point(data = resultsEstimatedParametersBinomialTest_summary, aes(x = Strain, y = Estimate, color = colour80),
             pch = 21,
             fill = "white",
             size = 3) +
  scale_colour_manual(breaks = uniqueValuesColours,
                      values = uniqueValuesColours,
                      labels = uniqueLabelsColours,
                      name = "Habitat specificity") +
  ggtitle("Plasmodium") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "bottom",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))+
  ylim(0,1) +
  ylab("Occurrence probability\n in urban areas")+
  xlab("Lineage")
plotBinomialTest2

# plot_final<-plot_grid(plotBinomialTest2,plotBinomialTest, align = "h", rel_widths = c(1.3, 5))

# Save plots --------------------------------------------------------------

# setwd("~/Dropbox/PostDoc_Malaria/Malaria_urbanization/3_results/4_souche_repartition")
# library(svglite)
# save_plot(
#   "plot_souche.svg",
#   plot_final,
#   ncol=2,
#   base_height=5,
#   base_asp = 1.5,
#   base_width=4
# )