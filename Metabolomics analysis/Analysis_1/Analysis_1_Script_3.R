
#### SOP ####
# run pair-wise correlations between omics data layers using the pearson method
# correlation results will be visualized using a sankey diagram


#### Required Libraries ####
library(MetLipNet)
library(Hmisc)
library(rstatix)



##### Data prep #####
meta_data_names <- NULL

dat_ready <- read.csv("2025-09-26_Cor_ready.csv") %>% 
  select(-Layer) %>% 
  pivot_wider(id_cols = Feature, names_from = Patient, values_from = Value_diff) 


##### Correlations #####
# This function calculates the pairwise Pearson correlations for metabolite intensities 
# and provides the correlation estimate and p-value for each pair.

# view the help page
?cal_met_cor

# The returned data frame will be ready for the igraph function.
dat_cor_result <- cal_met_cor(data = dat_ready, 
                              meta_cols = meta_data_names, 
                              identifier_col = "Feature", 
                              method = "pearson",
                              cor_threshold = 0.2,
                              p_threshold = 0.05) 



##### Correlation Data Export #####

# save the correlation results
# this will be used for the Sankey Diagram
write.csv(dat_cor_result, file=paste(Sys.Date(), "Hourigan_0005_dat_cor_result.csv", sep = "_"), row.names = F)




