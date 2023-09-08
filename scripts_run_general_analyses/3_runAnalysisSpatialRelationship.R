## ~~~~~~~~~~~~~~~~~~~~~~~ ##
## Running spatial analysis ##
## ~~~~~~~~~~~~~~~~~~~~~~~ ##

# Prerequisite: 
# This script must be run AFTER running 0_randomSampling.R 

# Goal of the script:
# This script runs the "analysisSpatialDistance.R" script on 1000 subsampled datasets to account for leucocytozoon strain incertitude.
# Here we investigate if there is spatial correlation in malarial lineages infections.

rm(list = ls())

# Libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

# Load subsampled data ---------------------------------------------------------------


dataSets_subsampled <- readRDS("../data/dataSets_subsampled.rds")

# Function to run script on subsampled data ---------------------------------------------------------------

runSpatialDistanceSubsample <- function(
    leucoMalaria_df,
    PlasmoMalaria_df){
  ### IMPORTANT QUESTION: Here, individuals are considered different between years (even if they are the same). Shouldn't we?
  
  ## Running analysis ------------------------------------------------------
  source("../scripts_subanalyses/analysisSpatialDistance.R", local = TRUE)
  
  ## Save results ------------------------------------------------------------
  
  results <- list(
    mantelTestResultsSiteStrainId,
    mantelTestResultsSitePercent,
    mantelTestResultsNest,
    matrixOverlapstrain,
    matrixOverlapstrainPercent,
    matrixNestOverlapstrain_nest,
    mantelTestResultsSiteStrainId_urbanisation,
    mantelTestResultsSitePercent_urbanisation,
    mantelTestResultsNest_urbanisation,
    matrixOfDistancesUrbanisationSite,
    matrixOfDistancesUrbanisationNest
  )
  return(results)
}

# Run analysis on subsampled data -----------------------------------------

outputSampling_SpatialDistance <- apply(dataSets_subsampled, 1, FUN = function(x){
  runSpatialDistanceSubsample(
    leucoMalaria_df = x[1][[1]],
    PlasmoMalaria_df = x[2][[1]])}
)
outputSampling_SpatialDistance <- do.call(rbind, outputSampling_SpatialDistance) %>% as.data.frame()
colnames(outputSampling_SpatialDistance) <- c("mantelSiteStrain", 
                                              "mantelSitePercent", 
                                              "mantelNest", 
                                              "OverlapStrainSite", 
                                              "OverlapStrainSitePercent", 
                                              "OverlapStrainNest",
                                              "mantelSiteStrain_urbanisation",
                                              "mantelSitePercent_urbanisation",
                                              "mantelNest_urbanisation",
                                              "DissimilaritySite_urbanisation",
                                              "DissimilarityNest_urbanisation")

# Summarise mantel test (i.e. how many had significant p-value?) with FDR---------------------------------------------------

### FDR
resultsMantelSiteStrain <- lapply(outputSampling_SpatialDistance$mantelSiteStrain, function(x){
  x$p 
})
vector_pvalue_corrected <- p.adjust(unlist(resultsMantelSiteStrain), method="fdr")
resultsMantelSiteStrain <- sum(ifelse(vector_pvalue_corrected > 0.05, 0, 1))
resultsMantelSiteStrain
## Distance and strain diversity -------------------------------------------

resultsMantelSiteStrain <- lapply(outputSampling_SpatialDistance$mantelSiteStrain, function(x){
  x$p
})
vector_pvalue_correctedSiteStrain <- p.adjust(unlist(resultsMantelSiteStrain), method="fdr")
resultsMantelSiteStrain <- sum(ifelse(vector_pvalue_correctedSiteStrain > 0.05, 0, 1))
resultsMantelSiteStrain
#0

resultsMantelSitePercent <- lapply(outputSampling_SpatialDistance$mantelSitePercent, function(x){
  x$p 
  
})
vector_pvalue_correctedSitePercent <- p.adjust(unlist(resultsMantelSitePercent), method="fdr")
resultsMantelSitePercent <- sum(ifelse(vector_pvalue_correctedSitePercent > 0.05, 0, 1))
resultsMantelSitePercent
#0

resultsMantelNest <- lapply(outputSampling_SpatialDistance$mantelNest, function(x){
  x$p
})
vector_pvalue_correctedSiteNest <- p.adjust(unlist(resultsMantelNest), method="fdr")
resultsMantelSiteNest <- sum(ifelse(vector_pvalue_correctedSiteNest > 0.05, 0, 1))
resultsMantelSiteNest
#0

## Urbanisation and strain diversity ---------------------------------------

resultsMantelSiteStrain_urbanisation <- lapply(outputSampling_SpatialDistance$mantelSiteStrain_urbanisation, function(x){
  x$p 
})
vector_pvalue_correctedSiteStrain_urbanisation <- p.adjust(unlist(resultsMantelSiteStrain_urbanisation), method="fdr")
resultsMantelSiteStrain_urbanisation <- sum(ifelse(vector_pvalue_correctedSiteStrain_urbanisation > 0.05, 0, 1))
resultsMantelSiteStrain_urbanisation
#0

resultsMantelSitePercent_urbanisation <- lapply(outputSampling_SpatialDistance$mantelSitePercent_urbanisation, function(x){
  x$p
})
vector_pvalue_correctedSitePercent_urbanisation <- p.adjust(unlist(resultsMantelSitePercent_urbanisation), method="fdr")
resultsMantelSitePercent_urbanisation <- sum(ifelse(vector_pvalue_correctedSitePercent_urbanisation > 0.05, 0, 1))
resultsMantelSitePercent_urbanisation
#0

resultsMantelNest_urbanisation <- lapply(outputSampling_SpatialDistance$mantelNest_urbanisation, function(x){
  x$p
})
vector_pvalue_correctedNest_urbanisation <- p.adjust(unlist(resultsMantelNest_urbanisation), method="fdr")
resultsMantelNest_urbanisation <- sum(ifelse(vector_pvalue_correctedNest_urbanisation > 0.05, 0, 1))
resultsMantelNest_urbanisation
#0



# Summarise mantel test (i.e. how many had significant p-value?) RAW PVALUE ---------------------------------------------------

# 
# ## Distance and strain diversity -------------------------------------------
# 
# resultsMantelSiteStrain <- lapply(outputSampling_SpatialDistance$mantelSiteStrain, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelSiteStrain <- sum(unlist(resultsMantelSiteStrain))
# resultsMantelSiteStrain
# #2
# resultsMantelSitePercent <- lapply(outputSampling_SpatialDistance$mantelSitePercent, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelSitePercent <- sum(unlist(resultsMantelSitePercent))
# resultsMantelSitePercent
# #2
# resultsMantelNest <- lapply(outputSampling_SpatialDistance$mantelNest, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelNest <- sum(unlist(resultsMantelNest))
# resultsMantelNest
# #613
# 
# ## Urbanisation and strain diversity ---------------------------------------
# 
# resultsMantelSiteStrain_urbanisation <- lapply(outputSampling_SpatialDistance$mantelSiteStrain_urbanisation, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelSiteStrain_urbanisation <- sum(unlist(resultsMantelSiteStrain_urbanisation))
# resultsMantelSiteStrain_urbanisation
# #0
# resultsMantelSitePercent_urbanisation <- lapply(outputSampling_SpatialDistance$mantelSitePercent_urbanisation, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelSitePercent_urbanisation <- sum(unlist(resultsMantelSitePercent_urbanisation))
# resultsMantelSitePercent_urbanisation
# #0
# 
# resultsMantelNest_urbanisation <- lapply(outputSampling_SpatialDistance$mantelNest_urbanisation, function(x){
#   if(x$p < 0.05){
#     return(1)
#   }else{
#     return(0)
#   }
# })
# resultsMantelNest_urbanisation <- sum(unlist(resultsMantelNest_urbanisation))
# resultsMantelNest_urbanisation
# #0
# 
# 
