#~~~~~~~~~~~
# ESTIMATION OF PREVALENCES CONFIDENCE INTERVALS
#~~~~~~~~~~~

# What does this script do?
# It estimates confidence intervals of malaria infection prevalence in urban and forest Great tits

# Libraries ---------------------------------------------------------------
library(prevalence)

## Wilson score interval methods (with or without continuity 
# correction) have been shown to be the most accurate and the most 
# robus
# sources:
# doi:10.1080/09296174.2013.799918
# doi:10.1214/ss/1009213286
# doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E

# Analyses ----------------------------------------------------------------

## Prevalence Haemo/plasmo in adults ---------------------------------------

# ROU
propCI(x = 51, n = 52, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 51 52 0.98076923 wilson  0.95 0.89879489 0.99659719

# URB TOT
propCI(x = 236, n = 244, method = "wilson", level = c(0.95))
# x   n          p method level      lower      upper
# 1 236 244 0.96721311 wilson  0.95 0.93664841 0.98329453

#ZOO
propCI(x = 39, n = 41, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 39 41 0.95121951 wilson  0.95 0.83861021 0.98651905

#CEF
propCI(x = 8, n = 8, method = "wilson", level = c(0.95))
# x n p method level     lower upper
# 1 8 8 1 wilson  0.95 0.7472764     1

#GRAM
propCI(x = 43, n = 45, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 43 45 0.95555556 wilson  0.95 0.85172489 0.98772587

#BOT
propCI(x = 12, n = 12, method = "wilson", level = c(0.95))
# x  n p method level      lower upper
# 1 12 12 1 wilson  0.95 0.81601881     1

#FAC
propCI(x = 19, n =20, method = "wilson", level = c(0.95))
# x   n          p method level      lower      upper
# 1 236 244 0.96721311 wilson  0.95 0.93664841 0.98329453

#FONT
propCI(x = 38, n = 40, method = "wilson", level = c(0.95))
# x  n    p method level      lower      upper
# 1 38 40 0.95 wilson  0.95 0.83496123 0.98617933

#MOS
propCI(x = 37, n = 37, method = "wilson", level = c(0.95))
# x  n p method level      lower upper
# 1 37 37 1 wilson  0.95 0.93185981     1

#MAS
propCI(x = 40, n = 41, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 40 41 0.97560976 wilson  0.95 0.87404938 0.99568147

#TOTAL 
#MAS
propCI(x = 287, n = 296, method = "wilson", level = c(0.95))
# x   n          p method level      lower      upper
# 1 287 296 0.96959459 wilson  0.95 0.94323393 0.98392271

## Prevalence leuco in adults ----------------------------------------------

# ROU
propCI(x = 48, n = 52, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 48 52 0.92307692 wilson  0.95 0.81826438 0.96968065

# URB TOT
propCI(x = 223, n = 244, method = "wilson", level = c(0.95))
# x   n          p method level      lower      upper
# 1 223 244 0.91393443 wilson  0.95 0.87201732 0.94301985

#ZOO
propCI(x = 36, n = 41, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 36 41 0.87804878 wilson  0.95 0.74455788 0.94676664

#CEF
propCI(x = 8, n = 8, method = "wilson", level = c(0.95))
# x n p method level     lower upper
# 1 8 8 1 wilson  0.95 0.7472764     1

#GRAM
propCI(x = 36, n = 45, method = "wilson", level = c(0.95))
# x  n   p method level      lower      upper
# 1 36 45 0.8 wilson  0.95 0.66177031 0.89103873

#BOT
propCI(x = 12, n = 12, method = "wilson", level = c(0.95))
# x  n p method level      lower upper
# 1 12 12 1 wilson  0.95 0.81601881     1

#FAC
propCI(x = 18, n =20, method = "wilson", level = c(0.95))
# x  n   p method level      lower      upper
# 1 18 20 0.9 wilson  0.95 0.69896635 0.97213352

#FONT
propCI(x = 38, n = 40, method = "wilson", level = c(0.95))
# x  n    p method level      lower      upper
# 1 38 40 0.95 wilson  0.95 0.83496123 0.98617933

#MOS
propCI(x = 37, n = 37, method = "wilson", level = c(0.95))
# x  n p method level      lower upper
# 1 37 37 1 wilson  0.95 0.93185981     1

#MAS
propCI(x = 38, n = 41, method = "wilson", level = c(0.95))
# x  n          p method level      lower      upper
# 1 38 41 0.92682927 wilson  0.95 0.80572551 0.97480217

#TOTAL 
#MAS
propCI(x = 271, n = 296, method = "wilson", level = c(0.95))
# x   n          p method level      lower      upper
# 1 271 296 0.91554054 wilson  0.95 0.87829477 0.94213881

## Prevalence Haemo/plasmo in nestlings ------------------------------------

# ROU
propCI(x = 0, n = 36, method = "wilson", level = c(0.95))
# x  n p method level lower       upper
# 1 0 36 0 wilson  0.95     0 0.069900671

# URB TOT
propCI(x = 9, n = 54, method = "wilson", level = c(0.95))
# x  n          p method level       lower      upper
# 1 9 54 0.16666667 wilson  0.95 0.090243916 0.28736514

#ZOO
propCI(x = 2, n = 14, method = "wilson", level = c(0.95))
# x  n          p method level       lower      upper
# 1 2 14 0.14285714 wilson  0.95 0.040093921 0.39941379

#GRAM
propCI(x = 2, n = 12, method = "wilson", level = c(0.95))
#x  n         p method level      lower     upper
#1 2 12 0.1666667 wilson  0.95 0.04696514 0.4480309

#FONT
propCI(x = 2, n = 7, method = "wilson", level = c(0.95))
#x  n         p method level      lower     upper
#1 2 7 0.2857143 wilson  0.95 0.08221892 0.6410655

#MOS
propCI(x = 3, n = 8, method = "wilson", level = c(0.95))
#x n     p method level     lower     upper
#1 3 8 0.375 wilson  0.95 0.1368443 0.6942576

#MAS
propCI(x = 0, n = 12, method = "wilson", level = c(0.95))
#x  n p method level lower     upper
#1 0 12 0 wilson  0.95     0 0.1839812

#TOTAL 

propCI(x = 9, n = 90, method = "wilson", level = c(0.95))
#   x  n   p method level      lower     upper
#1 9 90 0.1 wilson  0.95 0.05350675 0.1792417

## Prevalence Leuco in nestlings -------------------------------------------

# ROU
propCI(x = 1, n = 36, method = "wilson", level = c(0.95))
#  x  n          p method level       lower     upper
#1 1 36 0.02777778 wilson  0.95 0.004920407 0.1416972

# URB TOT
propCI(x = 6, n = 54, method = "wilson", level = c(0.95))
#   x  n         p method level      lower    upper
#1 6 54 0.1111111 wilson  0.95 0.05193022 0.221947

#ZOO
propCI(x = 0, n = 14, method = "wilson", level = c(0.95))
#  x  n p method level lower     upper
#1 0 14 0 wilson  0.95     0 0.1619548

#GRAM
propCI(x = 4, n = 12, method = "wilson", level = c(0.95))
#  x  n         p method level     lower     upper
#1 4 12 0.3333333 wilson  0.95 0.1381201 0.6093779

#FONT
propCI(x = 0, n = 7, method = "wilson", level = c(0.95))
#  x n p method level lower     upper
#1 0 7 0 wilson  0.95     0 0.2787627

#MOS
propCI(x = 0, n = 8, method = "wilson", level = c(0.95))
# x n p method level lower     upper
#1 0 8 0 wilson  0.95     0 0.2527236

#MAS
propCI(x = 2, n = 12, method = "wilson", level = c(0.95))
#  x  n         p method level      lower     upper
#1 2 12 0.1666667 wilson  0.95 0.04696514 0.4480309

#TOTAL 

propCI(x = 7, n = 90, method = "wilson", level = c(0.95))
#   x  n          p method level      lower     upper
#1 7 90 0.07777778 wilson  0.95 0.03818482 0.1519386

## CHI2 malaria prevalence in nestlings ------------------------------------

setwd("C:/Users/CAIZERGUES/Dropbox/PostDoc_Malaria/Malaria_urbanization/2_analyses/1_prevalences")

data_plasmo<-read.csv("preval_plasmo_poussin_pour_chi2.csv",sep=";", row.names=1)
chisq.test(data_plasmo)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  data_plasmo
# X-squared = 4.9434, df = 1, p-value = 0.02619

data_leuco<-read.csv("preval_leuco_poussin_pour_chi2.csv",sep=";", row.names=1)
chisq.test(data_leuco)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  data_leuco
# X-squared = 1.0908, df = 1, p-value = 0.2963

# Save analysis -----------------------------------------------------------
save.image("Renvironments/results_CI_prevalence.RData")