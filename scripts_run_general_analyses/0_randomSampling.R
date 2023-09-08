# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)

# Subsampling data (only leuco) --------------------------------------------------------
setwd("~/Dropbox/PostDoc_Malaria/Malaria_urbanization/Pipeline_Malaria_03032023")

#Number of subsampling
NSimulations = 1000
set.seed(42)#Fix seed for reproducibility

#Function to subsample
randomSampling <- function(){#Subsampling only for leuco (modify commented part if plasmodium has to be subsampled too)
  leucoMalaria_df <-
    read.csv("1_data/data_leuco_adultes_tri_ok.csv", sep = ";")
  PlasmoMalaria_df <-
    read.csv("1_data/data_malaria_adultes_ok_binom.csv", sep = ";")
  
  leucoMalaria_df <- leucoMalaria_df %>%
    dplyr::rename(nest = 'ad.nichoir',
                  stanic = 'sta.nic') %>%
    pivot_longer(
      cols = starts_with("souch"),
      names_to = "strainNumber",
      values_to = "strainID",
      values_drop_na = TRUE
    ) %>%
    filter(strainID != "inexploitable" & strainID != "") %>%
    #Subsampling only one strain by individual
    group_by(X) %>% #Group by id
    sample_n(1) %>%
    ungroup()
  
  PlasmoMalaria_df <- PlasmoMalaria_df %>%
    dplyr::rename(nest = 'ad.nichoir',
                  stanic = 'sta.nic') %>%
    pivot_longer(
      cols = starts_with("souch"),
      names_to = "strainNumber",
      values_to = "strainID",
      values_drop_na = TRUE
    ) #%>%
  # #Subsampling only one strain by individual: NOT TO DO WITH PLASMODIUM
  # group_by(code) %>% #Group by id
  # sample_n(1) %>%
  # ungroup()
  
return(
  list(
  leucoMalaria_df,
  PlasmoMalaria_df
  )
)
}

#Subsample data and store those subsamples to run analyses (see other scripts)
dataSets_subsampled <- lapply(1:NSimulations, FUN = function(x)randomSampling())
dataSets_subsampled <- do.call(rbind, dataSets_subsampled) %>% as.data.frame()
colnames(dataSets_subsampled) <- c("dataLeuco", "dataPlasmo")
saveRDS(dataSets_subsampled, "1_data/dataSets_subsampled.rds")
