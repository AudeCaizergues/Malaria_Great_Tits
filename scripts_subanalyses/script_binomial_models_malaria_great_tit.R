#~~~~~~~~~~~
# FACTORS EXPLAINING MALARIA PREVALENCE
#~~~~~~~~~~~

# What does this script do?
# This script computes linear models to assess the factors associated to variation in parasitic prevalence

# Libraries ---------------------------------------------------------------
library(readr)
library(lme4)
library(lmerTest)
library(lmtest)
library(dplyr)
library(ggplot2)
library(librarian)

#~~~~~~~~~~~

# Home made functions -----------------------------------------------------

## To source package efficiently (and install if not installed)
if (!("librarian" %in% installed.packages())) {
  install.packages("librarian")
}
library(librarian)

#Function to check linear model assumptions
diagnostics.plot.dharma <-
  function(mod.res,
           col = grey(level = 0.25, alpha = 0.5),
           breaks.histo = 20,
           quantreg = TRUE) {
    old.par = par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
    par(mar = c(3, 3, 3, 0.5))
    hist(
      residuals(mod.res),
      probability = T,
      xlab = "",
      ylab = "",
      main = "",
      breaks = breaks.histo
    )
    mtext(text = "Histogram of residuals",
          side = 3,
          line = 0)
    x = seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out =
              100)
    lines(x, dnorm(x, mean = 0, sd = sd(residuals(mod.res))))
    
    shelf(DHARMa)
    simulationOutput <-
      simulateResiduals(fittedModel = mod.res, plot = FALSE)
    
    plotQQunif(simulationOutput) # left plot in plot.DHARMa()
    plotResiduals(simulationOutput, quantreg = quantreg)
    #Old way without dharma, from Roger Mundry
    # qqnorm(residuals(mod.res), main="", pch=19)
    # qqline(residuals(mod.res))
    # mtext(text="qq-plot of residuals", side=3, line=0)
    # plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
    # abline(h=0, lty=2)
    # mtext(text="residuals against fitted values", side=3, line=0)
    par(old.par)
  }

#Function to check linear model validity
check_model <- function(#Still in construction not tested for all model types
  model,
  fittedWith = c("lm", "lmer", "glm", "glmer", "glmmTMB")
){
  
  shelf(performance)
  shelf(lme4)
  shelf(influence.ME)
  
  #Convergence
  singularityTest <- check_singularity(model)
  convergenceTest <- check_convergence(model)
  
  if(fittedWith %in% c("lm", "lmer")){
    family = "gaussian"
  }else{
    family = family(model)$family
  }
  #Overdispersion
  if(family[1] == "poisson"){
    overdispersionTest <- check_overdispersion(model)
  }else{
    overdispersionTest <- NA
  }
  
  #Autocorr
  autocorrTest <- check_autocorrelation(model)
  
  #Check correlation var (VIF)
  collinearityTest <- check_collinearity(model)
  
  
  
  #Check dfbetas
  if(fittedWith == "lm" | fittedWith == "glm"){
    dfbetas <- 
      round(cbind(coefficients(model), 
                  coefficients(model)+t(apply(X=dfbeta(model), 
                                              MARGIN=2, FUN=range))), 5)
    #Check outliers
    outlierTest <- check_outliers(model)
  }else if(fittedWith == "lmer" | fittedWith == "glmer"){
    dfbetas <- round(dfbetas(influence.ME::influence(model, obs=TRUE)), 2) 
    dfbetas <- as.data.frame(cbind(
      rownames(summary(model)$coefficients),
      summary(model)$coefficients[,1],
      summary(model)$coefficients[,1] + apply(dfbetas, 2, min),
      summary(modelDistance)$coefficients[,1] + apply(dfbetas, 2, max)
    )
    )
    colnames(dfbetas) <- c("Estimate", "Min", "Max")
    print("Dfbetas are computed at the singular value level. You may want to run it at the random level too with the influence() function of the lme4 package")
    #Check outliers
    outlierTest <- check_outliers(model)
  }else{
    randomEffects <- ifelse(length(ranef(model)$cond) > 0, TRUE, FALSE)
    if(randomEffects){
      print("Dfbetas are computed at the singular value level. You may want to run it at the random level too with the influence() function of the lme4 package")
      dfToTest <- model$frame
      if(family == "gaussian"){
        modelnonglmmTMB <- lmer(formula = formula(model), data = dfToTest)
      }else{
        modelnonglmmTMB <- glmer(formula = formula(model), data = dfToTest, family = family(model))
      }
      #Check outliers
      outlierTest <- check_outliers(modelnonglmmTMB)
      dfbetas <- round(dfbetas(influence.ME::influence(modelnonglmmTMB, data = dfToTest, obs=TRUE)), 2) 
      dfbetas <- as.data.frame(cbind(
        rownames(summary(modelnonglmmTMB)$coefficients),
        summary(modelnonglmmTMB)$coefficients[,1],
        summary(modelnonglmmTMB)$coefficients[,1] + apply(dfbetas, 2, min),
        summary(modelnonglmmTMBDistance)$coefficients[,1] + apply(dfbetas, 2, max)
      )
      )
    }else{
      if(family == "gaussian"){
        modelnonglmmTMB <- lm(formula = formula(model), data = model$frame)
      }else{
        modelnonglmmTMB <- glm(formula = formula(model), data = model$frame, family = family(model))
      }
      #Check outliers
      outlierTest <- check_outliers(modelnonglmmTMB)
      dfbetas <- 
        round(cbind(coefficients(modelnonglmmTMB), 
                    coefficients(modelnonglmmTMB)+t(apply(X=dfbeta(modelnonglmmTMB), 
                                                          MARGIN=2, FUN=range))), 2)
    }
  }
  colnames(dfbetas) <- c("Estimate", "Min", "Max")
  
  #Check assumptions
  dharmaOutput <- diagnostics.plot.dharma(model)
  heteroskTest <- check_heteroscedasticity(model)
  
  output <- list(
    singularityTest,
    convergenceTest,
    overdispersionTest,
    autocorrTest,
    collinearityTest,
    outlierTest,
    dfbetas,
    heteroskTest,
    dharmaOutput
  )
  names(output) <- 
    c(
      "singularityTest",
      "convergenceTest",
      "overdispersionTest",
      "autocorrTest",
      "collinearityTest",
      "outlierTest",
      "dfbetas",
      "heteroskTest",
      "dharmaOutput"
    )
  return(
    output
  )
}

# Import data adults -------------------------------------------------------------

setwd("~/Dropbox/PostDoc_Malaria/Malaria_urbanization/1_data")

# data <- read.csv("1_data/data_malaria_adultes_ok.csv",
#                  header = TRUE,
#                  sep = ";")
data <- read.csv("data_malaria_adultes_ok.csv",
                 header = TRUE,
                 sep = ";")
data$sex <- as.factor(data$sex)
data$plasmo_infection <- as.factor(data$HaemoPlasmo.1)
data$leuco_infection <- as.factor(data$Leuco.1)
data$year <- as.factor(data$year)
data$age_reel <- as.numeric(data$age_reel)

str(data)

#Relabel for good looking output:

data <- data %>%
  mutate(
    habitat = ifelse(habitat == "urb", "Urban", "Rural"),
    sex = ifelse(sex == "1", "Female", "Male")
  ) %>%
  dplyr::rename(
    `Nest-level naturalness` = "PC1",
    `Site-level naturalness` = "PC1_mean_zone",
    `Infection (Plasmodium/Haemoproteus)` = "plasmo_infection",
    `Infection (Leucocytozoon)` = "leuco_infection",
    Station = "station",
    Sex = "sex",
    Year = "year",
    Age = "age_reel",
    Habitat = "habitat"
  )

#Problem of complete separation because insufficient sampling
table(data$Year)

data <- data %>%
  filter(Year != 2017 & Year != 2018)
#~~~~~~~~~~~



# Adult Run separated models ---------------------------------------------------

## MODELS HAEMO/PLASMO ---------------------------------------------------------
### with Habitat ------------------------------------------------------------
#### full model ----------------------

plasmo_hab1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex + Year ,
    family = "binomial",
    data = data
  )
summary(plasmo_hab1)

#Check and plot + table
modelCheck_plasmo_hab1 <- check_model(plasmo_hab1,
                                      fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_hab1

library(sjPlot)
tab_model(
  plasmo_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_plasmo_hab1<-plot_model(
  plasmo_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Year
plasmo_hab2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex,
    family = "binomial",
    data = data
  )
summary(plasmo_hab2)
lrtest(plasmo_hab1, plasmo_hab2)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex + 
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -37.5377                      
# 2   4 -37.9855 -1 0.89564    0.34395

## Sex
plasmo_hab3 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_hab3)
lrtest(plasmo_hab1, plasmo_hab3)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -37.538
# 2   6 -38.065 -1 1.0539     0.3046

## AGE
plasmo_hab4 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_hab4)
lrtest(plasmo_hab1, plasmo_hab4)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -37.538
# 2   6 -37.924 -1 0.7733     0.3792

## Habitat
plasmo_hab5 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_hab5)
lrtest(plasmo_hab1, plasmo_hab5)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -37.538
# 2   6 -37.539 -1 0.0032     0.9549

### with zone urbanization level --------------------------------------------

#### full model ----------------------

plasmo_urbzone1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_urbzone1)

#Check and plot + table
modelCheck_plasmo_urbzone1 <- check_model(plasmo_urbzone1,
                                          fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  plasmo_urbzone1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

plot_plasmo_urbzone1<-plot_model(
  plasmo_urbzone1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )


## TESTING VARIABLE SIGNIFICANCE
# Sex : IS NOT SIGNIFICANT
plasmo_urbzone2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + Age + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_urbzone2)
lrtest(plasmo_urbzone1, plasmo_urbzone2)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Sex + Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -37.3591                      
# 2   4 -37.8922 -1 1.06604    0.30184

## Year : NOT SIGNIFICANT
plasmo_urbzone3 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + Age + Sex,
    family = "binomial",
    data = data
  )
summary(plasmo_urbzone3)
lrtest(plasmo_urbzone1, plasmo_urbzone3)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Sex + Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Sex
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -37.3591                      
# 2   4 -38.0668 -1 1.41526    0.23419

## AGE : NOT SIGNIFICANT
plasmo_urbzone4 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_urbzone4)
lrtest(plasmo_urbzone1, plasmo_urbzone4)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Sex + Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Sex + Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -37.3591                      
# 2   4 -37.7498 -1 0.78139    0.37672

## Station :
plasmo_urbzone5 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_urbzone5)
lrtest(plasmo_urbzone1, plasmo_urbzone5)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` + 
#   Age + Sex + Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -37.3591                      
# 2   4 -37.5393 -1 0.36036    0.54831

### with nestbox urbanization level -----------------------------------------

#### full model ----------------------
plasmo_nest1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_nest1)

#Check and plot + table
modelCheck_plasmo_nest1 <- check_model(plasmo_nest1,
                                       fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  plasmo_nest1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

plot_plasmo_nest1<-plot_model(
  plasmo_nest1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Sex : NOT SIGNIFICANT
plasmo_nest2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_nest2)
lrtest(plasmo_nest1, plasmo_nest2)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -36.071
# 2   6 -36.533 -1 0.9247     0.3363

## Year : NOT SIGNIFICANT
plasmo_nest3 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Sex,
    family = "binomial",
    data = data
  )
summary(plasmo_nest3)
lrtest(plasmo_nest1, plasmo_nest3)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + 
#   Age + Sex + Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + 
#   Age + Sex
# #Df   LogLik Df  Chisq Pr(>Chisq)
# 1   5 -36.0710                     
# 2   4 -37.3816 -1 2.6212    0.10544

## AGE :
plasmo_nest4 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_nest4)
lrtest(plasmo_nest1, plasmo_nest4)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -36.071
# 2   6 -36.433 -1 0.7241     0.3948

## NESTBOX : MARGINALLY SIGNIFICANT
plasmo_nest5 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_nest5)
lrtest(plasmo_nest1, plasmo_nest5)
# Likelihood ratio test
#
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness` + Age + Sex +
#   Year
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ Age + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -36.071
# 2   6 -37.539 -1 2.9367    0.08659 .
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
##~~~~~~~~

## MODELS LEUCO ------------------------------------------------------------

### with Habitat ------------------------------------------------------------

#### full model ----------------------
leuco_hab1 <-
  glm(
    `Infection (Leucocytozoon)` ~ Habitat + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_hab1)


#Check and plot + table
modelCheck_leuco_hab1 <- check_model(leuco_hab1,
                                     fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  #plasmo_hab1,
  leuco_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

plot_leuco_hab1<-plot_model(
  leuco_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Year
leuco_hab2 <-
  glm(`Infection (Leucocytozoon)` ~ Habitat + Age + Sex,
      family = "binomial",
      data = data)
summary(leuco_hab2)
# lrtest(leuco_hab1, leuco_hab2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ Habitat + Age + Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ Habitat + Age + Sex
# #Df   LogLik Df   Chisq Pr(>Chisq)  
# 1   5 -79.7041                        
# 2   4 -82.2417 -1 5.07512   0.024272 *

## Sex
leuco_hab3 <-
  glm(`Infection (Leucocytozoon)` ~ Habitat + Age + Year,
      family = "binomial",
      data = data)
summary(leuco_hab3)
lrtest(leuco_hab1, leuco_hab3)
# Likelihood ratio test
#
# Model 1: `Infection (Leucocytozoon)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Leucocytozoon)` ~ Habitat + Age + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -79.704
# 2   6 -79.849 -1 0.2887      0.591

## AGE
leuco_hab4 <-
  glm(`Infection (Leucocytozoon)` ~ Habitat + Sex + Year,
      family = "binomial",
      data = data)
summary(leuco_hab4)
lrtest(leuco_hab1, leuco_hab4)
# Likelihood ratio test
#
# Model 1: `Infection (Leucocytozoon)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Leucocytozoon)` ~ Habitat + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -79.704
# 2   6 -80.140 -1 0.8713     0.3506

## Habitat
leuco_hab5 <-
  glm(`Infection (Leucocytozoon)` ~ Age + Sex + Year,
      family = "binomial",
      data = data)
summary(leuco_hab5)
lrtest(leuco_hab1, leuco_hab5)
# Likelihood ratio test
#
# Model 1: `Infection (Leucocytozoon)` ~ Habitat + Age + Sex +
#   Year
# Model 2: `Infection (Leucocytozoon)` ~ Age + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -79.704
# 2   6 -80.559 -1 1.7102      0.191



### with zone urbanization level --------------------------------------------

#### full model  ----------------------
leuco_urbzone1 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_urbzone1)

#Check and plot + table
modelCheck_leuco_urbzone1 <- check_model(leuco_urbzone1,
                                         fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  leuco_urbzone1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

plot_leuco_urbzone1<-plot_model(
  leuco_urbzone1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Year
leuco_urbzone2 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + Sex,
    family = "binomial",
    data = data
  )
summary(leuco_urbzone2)
lrtest(leuco_urbzone1, leuco_urbzone2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Sex
# #Df   LogLik Df   Chisq Pr(>Chisq)  
# 1   5 -80.5530                        
# 2   4 -82.0199 -1 2.93379   0.086744 .

## Sex
leuco_urbzone3 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + Year,
    family = "binomial",
    data = data
  )
summary(leuco_urbzone3)
lrtest(leuco_urbzone1, leuco_urbzone3)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -80.5530                      
# 2   4 -80.7339 -1 0.36174    0.54754

## AGE
leuco_urbzone4 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_urbzone4)
lrtest(leuco_urbzone1, leuco_urbzone4)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Sex + 
#   Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -80.5530                      
# 2   4 -80.8265 -1 0.54692    0.45958

## Station URB
leuco_urbzone5 <-
  glm(`Infection (Leucocytozoon)` ~ Age + Sex + Year,
      family = "binomial",
      data = data)
summary(leuco_urbzone5)
lrtest(leuco_urbzone1, leuco_urbzone5)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Site-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ Age + Sex + Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -80.5530                      
# 2   4 -80.5593 -1 0.01246    0.91113

### with nestbox urbanization level  ----------------------------------------

#### full model ----------------------
leuco_nest1 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_nest1)

#Check and plot + table
modelCheck_leuco_nest1 <- check_model(leuco_nest1,
                                      fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  leuco_nest1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

plot_leuco_nest1<-plot_model(
  leuco_nest1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Year
leuco_nest2 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + Sex,
    family = "binomial",
    data = data
  )
summary(leuco_nest2)
lrtest(leuco_nest1, leuco_nest2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + 
#   Sex
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -80.2760                      
# 2   4 -81.3277 -1 2.10332    0.14698

## Sex
leuco_nest3 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + Year,
    family = "binomial",
    data = data
  )
summary(leuco_nest3)
lrtest(leuco_nest1, leuco_nest3)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + 
#   Sex + Year
# Model 2: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + 
#   Year
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   5 -80.2760                      
# 2   4 -80.4601 -1 0.36815    0.54402

## AGE
leuco_nest4 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_nest4)
lrtest(leuco_nest1, leuco_nest4)
# Likelihood ratio test
#
# Model 1: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + Sex +
#   Year
# Model 2: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -80.276
# 2   6 -80.504 -1 0.4552     0.4999

## NESTBOX
leuco_nest5 <-
  glm(`Infection (Leucocytozoon)` ~ Age + Sex + Year,
      family = "binomial",
      data = data)
summary(leuco_nest5)
lrtest(leuco_nest1, leuco_nest5)
# Likelihood ratio test
#
# Model 1: `Infection (Leucocytozoon)` ~ `Nest-level naturalness` + Age + Sex +
#   Year
# Model 2: `Infection (Leucocytozoon)` ~ Age + Sex + Year
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   7 -80.276
# 2   6 -80.559 -1 0.5666     0.4516


# Adult Run full models --------------------------------------------------------------

## MODELS HAEMO/PLASMO -----------------------------------------------------


### FULL MODEL ----------------------

plasmo_hab1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat + `Site-level naturalness` + `Nest-level naturalness` + Sex + Age + Year,
    family = "binomial",
    data = data
  )
summary(plasmo_hab1)

lrtest_plasmo <- drop1(plasmo_hab1, test = "Chisq")
lrtest_plasmo

#Check and plot + table
modelCheck_plasmo_hab1 <- check_model(plasmo_hab1,
                                      fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_hab1

library(sjPlot)
tab_model(
  plasmo_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 3,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plotModel_plasmo <- plot_model(
  plasmo_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )



## MODELS LEUCO ------------------------------------------------------------
### FULL MODEL ----------------------
leuco_hab1 <-
  glm(
    `Infection (Leucocytozoon)` ~ Habitat + `Site-level naturalness` + `Nest-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_hab1)

leuco_hab1 <-
  glm(
    `Infection (Leucocytozoon)` ~ Habitat + `Site-level naturalness` + Age + Sex + Year,
    family = "binomial",
    data = data
  )
summary(leuco_hab1)


lrtest_leuco <- drop1(leuco_hab1, test = "Chisq")
lrtest_leuco

#Check and plot + table
modelCheck_leuco_hab1 <- check_model(leuco_hab1,
                                     fittedWith = "glm")#A changer en fonction de ce qui est utilise

library(sjPlot)
tab_model(
  leuco_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 3,
  show.ci = 0.95,
  transform = NULL
)

tab_model(
  plasmo_hab1,
  leuco_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 3,
  show.ci = 0.95,
  transform = NULL
)

plotModel_leuco <- plot_model(
  leuco_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )




# Import data nestlings-------------------------------------------------------------

setwd("~/Dropbox/PostDoc_Malaria/Malaria_urbanization/Pipeline_Malaria_03032023/1_data/")
data_pouss<-read.csv("data_malaria_poussins_ok.csv",header=TRUE, sep=";")
data_pouss$plasmo_infection<-as.factor(data_pouss$HaemoPlasmo.1)
data_pouss$leuco_infection<-as.factor(data_pouss$Leuco.1)
data_pouss$station<-as.factor(data_pouss$pouss.station)
is.numeric(data_pouss$`Site-level naturalness`)
str(data_pouss)

#Relabel for good looking output:

data_pouss <- data_pouss %>%
  mutate(
    habitat = ifelse(habitat == "urb", "Urban", "Rural")
  ) %>%
  dplyr::rename(
    `Nest-level naturalness` = "PC1",
    `Site-level naturalness` = "PC1_mean_zone",
    `Infection (Plasmodium/Haemoproteus)` = "plasmo_infection",
    `Infection (Leucocytozoon)` = "leuco_infection",
    Station = "station",
    Habitat = "habitat"
  )
str(data_pouss)
#~~~~~~~~~~~


# Nestlings Run models ---------------------------------------------------------

## MODELS HAEMO/PLASMO ---------------------------------------------------------

### with Habitat ----
#### full model ----
plasmo_nest_hab1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ Habitat,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_hab1)

#Check and plot + table
modelCheck_plasmo_nest_hab1 <- check_model(plasmo_nest_hab1,
                                      fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_nest_hab1

library(sjPlot)
tab_model(
  plasmo_nest_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_plasmo_nest_hab1<-plot_model(
  plasmo_nest_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Habitat : SIGNIFICANT
plasmo_nest_hab2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_hab2)
lrtest(plasmo_nest_hab1, plasmo_nest_hab2)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ Habitat
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)   
# 1   2 -24.3303                         
# 2   1 -29.2575 -1 9.85432  0.0016943 **


### with zone urbanization level ----
#### full model ----
plasmo_nest_urbzone1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness` ,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_urbzone1)

#Check and plot + table
modelCheck_plasmo_nest_urbzone1 <- check_model(plasmo_nest_urbzone1,
                                           fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_nest_urbzone1

library(sjPlot)
tab_model(
  plasmo_nest_urbzone1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_plasmo_nest_urbzone1<-plot_model(
  plasmo_nest_urbzone1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Urb : NS
plasmo_nest_urbzone2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_urbzone2)
lrtest(plasmo_nest_urbzone1, plasmo_nest_urbzone2)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Site-level naturalness`
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   2 -28.6646                      
# 2   1 -29.2575 -1 1.18576    0.27619

### with nest box urbanization level ----
#### full model ----
plasmo_nest_urbnest1 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness`  ,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_urbnest1)

#Check and plot + table
modelCheck_plasmo_nest_urbnest1 <- check_model(plasmo_nest_urbnest1,
                                               fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_nest_urbnest1

library(sjPlot)
tab_model(
  plasmo_nest_urbnest1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_plasmo_nest_urbnest1<-plot_model(
  plasmo_nest_urbnest1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Urb : NS
plasmo_nest_urbnest2 <-
  glm(
    `Infection (Plasmodium/Haemoproteus)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(plasmo_nest_urbnest2)
lrtest(plasmo_nest_urbnest1, plasmo_nest_urbnest2)
# Likelihood ratio test
# 
# Model 1: `Infection (Plasmodium/Haemoproteus)` ~ `Nest-level naturalness`
# Model 2: `Infection (Plasmodium/Haemoproteus)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   2 -29.2507                      
# 2   1 -29.2575 -1 0.01349    0.90754


## MODELS LEUCO ---------------------------------------------------------

### with Habitat ----
#### full model ----
leuco_nest_hab1 <-
  glm(
    `Infection (Leucocytozoon)` ~ Habitat,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_hab1)

#Check and plot + table
modelCheck_leuco_nest_hab1 <- check_model(leuco_nest_hab1,
                                           fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_leuco_nest_hab1

library(sjPlot)
tab_model(
  leuco_nest_hab1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_leuco_nest_hab1<-plot_model(
  leuco_nest_hab1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Habitat : NS
leuco_nest_hab2 <-
  glm(
    `Infection (Leucocytozoon)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_hab2)
lrtest(leuco_nest_hab1, leuco_nest_hab2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ Habitat
# Model 2: `Infection (Leucocytozoon)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   2 -23.4064                      
# 2   1 -24.5977 -1 2.38259    0.12269


### with zone urbanization level ----
#### full model ----
leuco_nest_urbzone1 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Site-level naturalness` ,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_urbzone1)

#Check and plot + table
modelCheck_leuco_nest_urbzone1 <- check_model(leuco_nest_urbzone1,
                                               fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_plasmo_nest_urbzone1

library(sjPlot)
tab_model(
  leuco_nest_urbzone1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_leuco_nest_urbzone1<-plot_model(
  leuco_nest_urbzone1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Urb : SIGNIFICANT
leuco_nest_urbzone2 <-
  glm(
    `Infection (Leucocytozoon)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_urbzone2)
lrtest(leuco_nest_urbzone1, leuco_nest_urbzone2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Site-level naturalness`
# Model 2: `Infection (Leucocytozoon)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   2 -23.9524                      
# 2   1 -24.5977 -1 1.29065    0.25593

### with nest box urbanization level ----
#### full model ----
leuco_nest_urbnest1 <-
  glm(
    `Infection (Leucocytozoon)` ~ `Nest-level naturalness`  ,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_urbnest1)

#Check and plot + table
modelCheck_leuco_nest_urbnest1 <- check_model(leuco_nest_urbnest1,
                                               fittedWith = "glm")#A changer en fonction de ce qui est utilise
modelCheck_leuco_nest_urbnest1

library(sjPlot)
tab_model(
  leuco_nest_urbnest1,
  show.se = TRUE,
  show.r2 = FALSE,
  digits = 2,
  show.ci = 0.95,
  transform = NULL
)

library(ggplot2)
plot_leuco_nest_urbnest1<-plot_model(
  leuco_nest_urbnest1,
  colors = "bw",
  vline.color = "black",
  sort.est = TRUE,
  show.intercept = FALSE,
  show.p = TRUE,
  show.values = TRUE,
  transform = NULL,
  value.offset = .15,
  p.style = "asterisk",
  p.threshold = c(0.05, 0.01, 0.001),
  line.size = 0.5
) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_line(colour = "grey93"),
    panel.grid.major = element_line(colour = "grey93"),
    strip.background = element_rect(colour = "white",
                                    fill = "white"),
    strip.text = element_text(face = "bold", size = 14)
  )

## Urb : NS
leuco_nest_urbnest2 <-
  glm(
    `Infection (Leucocytozoon)` ~ 1,
    family = "binomial",
    data = data_pouss
  )
summary(leuco_nest_urbnest2)
lrtest(leuco_nest_urbnest1, leuco_nest_urbnest2)
# Likelihood ratio test
# 
# Model 1: `Infection (Leucocytozoon)` ~ `Nest-level naturalness`
# Model 2: `Infection (Leucocytozoon)` ~ 1
# #Df   LogLik Df   Chisq Pr(>Chisq)
# 1   2 -23.6794                      
# 2   1 -24.5977 -1 1.83671    0.17534



#~~~~~~~~~~~
# Save analysis -----------------------------------------------------------
#save.image("Renvironments/binomial_models.RData")
save.image("/Users/aude/Dropbox/PostDoc_Malaria/Malaria_urbanization/Pipeline_Malaria_03032023/2_analyses/2_binomial_models/binomial_models.RData")


