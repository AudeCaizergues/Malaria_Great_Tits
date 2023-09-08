## ~~~~~~~~~~~~~~~~~~~~~~~ ##
## Running strain analysis ##
## ~~~~~~~~~~~~~~~~~~~~~~~ ##

# Prerequisite: 
# This script must be run AFTER running 0_randomSampling.R 

# Goal of the script:
# This script runs the "lineage_analysis.R" script on 1000 subsampled datasets to account for leucocytozoon strain incertitude.
# Here we investigate if some malarial lineages are associated with some sites/levels of urbanization.
rm(list = ls())

# Libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Load subsampled data ---------------------------------------------------------------
dataSets_subsampled <- readRDS("../data/dataSets_subsampled.rds")

# Function to run script on subsampled data ---------------------------------------------------------------

runAnalyseSoucheSubsample <- function(
    leucoMalaria_df,
    PlasmoMalaria_df,
    onlyPlotPlasmo = FALSE){

  ## Running analysis ------------------------------------------------------
  source("../scripts_subanalyses/lineage_analysis.R", local=TRUE)
  
  ## Save results ------------------------------------------------------------
  
  if(onlyPlotPlasmo){
    return(plotBinomialTest2)
  }else{
    results <- list(
      leucoMalaria_df,
      valueExpectedLeuco,
      PlasmoMalaria_df,
      valueExpectedPlasmo,
      resultsEstimatedParametersBinomialTest_summary_LEUCO,
      resultsEstimatedParametersBinomialTest_summary_PLASMO
    )
    return(results)
  }
}

# Run analysis on subsampled data -----------------------------------------

outputSampling_AnalyseSouche <- apply(dataSets_subsampled, 1, FUN = function(x){
  runAnalyseSoucheSubsample(
    leucoMalaria_df = x[1][[1]],
    PlasmoMalaria_df = x[2][[1]],
    onlyPlotPlasmo = FALSE)}
  )
outputSampling_AnalyseSouche <- do.call(rbind, outputSampling_AnalyseSouche) %>% as.data.frame()
colnames(outputSampling_AnalyseSouche) <- c("dataLeuco",  "valueExpectedLeuco", "dataPlasmo", "valueExpectedPlasmo", "resultBinomLeuco", "resultBinomPlasmo")

# For LEUCO: From subsample, extract only median result ----------------------------------------

valueExpectedLeuco <- outputSampling_AnalyseSouche$valueExpectedLeuco[1][[1]]
# Get median estimates for all
mergedResults_LEUCO <- do.call(rbind, outputSampling_AnalyseSouche$resultBinomLeuco) %>% 
  as.data.frame()

removeInsufficientlySampledStrains <- table(mergedResults_LEUCO$Strain)

whichToKeepGlobalInsufficientSamples <- names(removeInsufficientlySampledStrains[removeInsufficientlySampledStrains>=50])

medianResults_LEUCO <- mergedResults_LEUCO %>% 
  filter(Strain %in% whichToKeepGlobalInsufficientSamples) %>%
  group_by(Strain) %>% 
  mutate(nobs=length(Strain)) %>% 
  summarise(
    #Calculating the CI of the median estimate based on sampling
    EstimateLow = median(Estimate) - qt(0.95,df=nobs-1)*sd(Estimate)/sqrt(nobs),
    EstimateHigh = median(Estimate) + qt(0.95,df=nobs-1)*sd(Estimate)/sqrt(nobs),
    #Median estimate
    Estimate = median(Estimate)
  ) %>% 
  ungroup() %>% 
  left_join(
    mergedResults_LEUCO %>% unique(),
    by = c("Strain" = "Strain", "Estimate" = "Estimate")
  ) %>% 
  unique()

# Replot median plot
resultsEstimatedParametersBinomialTest_summary <- medianResults_LEUCO

#If multiple samples give the median, all are selected: hence subsample one by taking the one with the least sensitive to error I (i.e., if non  significant results)
countStrain_v <-
  table(resultsEstimatedParametersBinomialTest_summary$Strain)
countStrain_v <- countStrain_v[countStrain_v > 1]
strainWithDuplicates <- names(countStrain_v)

for(i in 1:length(strainWithDuplicates)){
  rowsToFilter <-
    which(resultsEstimatedParametersBinomialTest_summary$Strain ==
            strainWithDuplicates[i])
  resCI80 <-
    resultsEstimatedParametersBinomialTest_summary$colour80[rowsToFilter]
  resCI95 <-
    resultsEstimatedParametersBinomialTest_summary$colour95[rowsToFilter]
  
  if(length(unique(resCI80)) == 1 & length(unique(resCI95)) == 1){
    toKeep <- sample(rowsToFilter, 1)
  } else if(length(unique(resCI80)) == 2){
    if("black" %in% unique(resCI80)){
      whichNonSignificant <- which(resCI80 == "black")
      toKeep <- sample(rowsToFilter[whichNonSignificant], 1)
    }else{
      toKeep <- sample(rowsToFilter, 1)
    }
  }
  #print(toKeep)
  if(length(toKeep) > 1){
    print(paste0("Error in keeping one value among all estimates being
equal to the median for step ", i))
  }else{
    toRemove <- rowsToFilter[rowsToFilter != toKeep]
    resultsEstimatedParametersBinomialTest_summary <-
      resultsEstimatedParametersBinomialTest_summary[-toRemove,]
  }
}

#If for any reason, subsample gave a "even" number of simulations then, no median is found: take the existing value just above

mergedDataBinom <- do.call(rbind, outputSampling_AnalyseSouche$resultBinomLeuco) %>% 
  as.data.frame()

strainWithNoOutput <- resultsEstimatedParametersBinomialTest_summary$Strain[
  which(is.na(resultsEstimatedParametersBinomialTest_summary$colour80))
]
if(length(strainWithNoOutput)>0){
  for(i in 1:length(strainWithNoOutput)){
   transitoryData <- mergedDataBinom[
     mergedDataBinom$Strain == strainWithNoOutput[i],
   ]%>% 
      dplyr::arrange(Estimate)
   toKeep <- which(transitoryData$Estimate >= median(transitoryData$Estimate))[1]
  
   resultsEstimatedParametersBinomialTest_summary[resultsEstimatedParametersBinomialTest_summary$Strain == strainWithNoOutput[1], 4:ncol(resultsEstimatedParametersBinomialTest_summary)] <- transitoryData[toKeep,-1]
  }
}





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

plotBinomialTest <- ggplot(resultsEstimatedParametersBinomialTest_summary %>% arrange(Estimate) %>% mutate(Strain=factor(Strain,unique(Strain))),### The data frame to use.
                           aes(x = Strain,
                               y = Estimate)) +
  # CI of the median mean
  geom_errorbar(aes(ymin  = CILow95,
                    ymax  = CIHigh95, 
                    color = colour95),
                width = 0.2,
                size  = 0.7) +
  # geom_errorbar(aes(ymin  = CILow80,
  #                  ymax  = CIHigh80, 
  #                 color = colour80),
  #            width = 0,
  #           size  = 1) +
  # #CI of the median
  # geom_errorbar(aes(x = as.numeric(as.factor(Strain)) - 0.25,
  #                 ymin  = EstimateLow,
  #                 ymax  = EstimateHigh,
  #                 color = colour80),
  #                 width = 0, size  = 0.35) +
# Expected probability
geom_hline(yintercept = valueExpectedLeuco, linetype = "dashed") +
  # Estimate
  geom_point(data = resultsEstimatedParametersBinomialTest_summary, aes(x = Strain, y = Estimate, color = colour80),
             pch = 21,
             fill = "white",
             size = 3) +
  scale_colour_manual(breaks = uniqueValuesColours,
                      values = uniqueValuesColours,
                      labels = uniqueLabelsColours,
                      name = "Habitat specificity") +
  ggtitle("Leucocytozoon") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        #legend.box="vertical",
        panel.grid.minor=element_line(colour="grey97"),
        panel.grid.major=element_line(colour="grey93")) +
  ylim(0,1) +
  ylab("") +
  xlab("Lineage")
plotBinomialTest
# For PLASMO:  ----------------------------------------

plotBinomialTest2 <- runAnalyseSoucheSubsample(
  leucoMalaria_df = dataSets_subsampled[1,1][[1]],
  PlasmoMalaria_df = dataSets_subsampled[1,2][[1]],
  onlyPlotPlasmo = TRUE)
plotBinomialTest2

#Saving the merge plot
library(ggpubr)
plot_final <-ggarrange(plotBinomialTest2,
                       plotBinomialTest,
                       nrow=1,
                       common.legend = TRUE,
                       legend = "bottom",
                       widths = c(1.3, 5))
plot_final 
  
  
  
  
plot_grid(
    plotBinomialTest2,
    plotBinomialTest,
    align = "h",
    rel_widths = c(1.3, 5),
    
  )
plot_final

setwd(
  "../Figure"
)
library(svglite)
save_plot(
  "plot_strains.png",
  plot_final,
  ncol = 2,
  base_height = 5,
  base_asp = 1.5,
  base_width = 4
)
