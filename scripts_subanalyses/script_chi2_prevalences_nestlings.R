#################################################################
########## CHI2 MALARIA PREVALENCE IN NESTLINGS  ################
#################################################################


conting_table_plasmo<-matrix(c(0,36,9,45), ncol=2, byrow=TRUE)
chisq.test(conting_table_plasmo)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  conting_table_plasmo
# X-squared = 4.9434, df = 1, p-value = 0.02619

conting_table_leuco<-matrix(c(1,35,6,48), ncol=2, byrow=TRUE)
chisq.test(conting_table_leuco)
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  conting_table_leuco
# X-squared = 1.0908, df = 1, p-value = 0.2963
