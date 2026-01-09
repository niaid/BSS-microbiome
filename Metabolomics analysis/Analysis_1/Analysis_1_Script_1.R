
#### SOP ####
# the purpose of this script was to split the patient cohort depending on their responder status to BSS.
# then, compare each time point specific fold changes back to day0 and similarly determine temporal significant differences.
# finally, this script generates bubble plots for the main figure, and matrix bar plots for the supplement.

# original dataframe for responder status was provided with entries in Response detailing Responder, and Non-Responder.
# the collaborators requested that we refer to them as high and low responders moving forward.


#### Required Libraries ####
library(tidyverse)
library(janitor)
library(tibble)

library(matrixBP)
library(deeptime)
library(patchwork)
library(cowplot)


#### Functions ####
# Function to run wilcox test for one metabolite and one time point
# Important to note: by using pivot wider, time points are moved to columns.
# This strategy innately aligns patients across time points.
# The wilcox test relies on vector input. Therefore, patient IDs are not matched.
# Hence, the importance of the alignment via pivot wider.
wilcox_test <- function(met_name, tp) {
  wide <- df %>%
    filter(metabolite == met_name) %>%
    pivot_wider(names_from = time, values_from = Value) %>%
    drop_na(`Day 0`, !!sym(tp)) # remove pairs with NA
  
  if (nrow(wide) >= 3) { # need at least 3 pairs
    wt <- wilcox.test(wide$`Day 0`, wide[[tp]], paired = TRUE, exact = TRUE)
    return(wt$p.value)
  } else {
    return(NA)
  }
}




#### Preliminary High Responder Analysis ####

# this was the preliminary responder analysis. Visualization included bubble plot and volcano plot
# the MBP chunks are further downstream in the analysis


##### High Responder Coercion #####

dat_t <- dat %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_wider <- dat_t %>% 
  pivot_longer(4:ncol(dat_t), names_to = "Metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>%
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  pivot_wider(id_cols = c(Patient, Metabolite), names_from = "Group", values_from = "Value")


# this is the backbone data frame used for statistics and fold change calculations
dat_longer <- dat_t %>% 
  pivot_longer(4:ncol(dat_t), names_to = "Metabolite", values_to = "Value") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  mutate_at(c("Value", "SubjectID"), as.numeric) %>%
  
  
  # merge patient responder status to treatment (data provided by Lauren)
  left_join(dat_responders, join_by(SubjectID)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient) %>% 
  
  # filter data set for responders
  filter(Response == "Responder")
  

  
##### High Responder Stat calculations #####
df <- dat_longer %>% 
  rename(time = Group, metabolite = Metabolite)

# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### High Responder Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the bubble plot and corr network analysis.


# determine which metabolites show a significant difference for at least 1 time point
responder_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()


# save the sig hits to use as filter in a sankey diagram
#write.csv(responder_sig_hits, "Responder_only_wilcox_day2_sig_hits.csv", row.names = F)


# count the number of significant hits
num_sig_hits <- length(responder_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)

# this strategy pairs the fold changes then takes the median value. 
Day0_value <- df %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  
  # calculate the median per metabolite
  group_by(time, metabolite) %>% 
  mutate(log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # keep only the columns required for downstream plotting
  select(metabolite, time, log2FC)


# calculate paired fold changes maintain patient specificity
# the fc_results_pat_spef dataframe will be used in the correlation network script. 
fc_results_pat_spef <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  ungroup() %>% 
  select(Patient, metabolite, time, fold_change, log2FC)


# export here
#write.csv(fc_results_pat_spef, paste(Sys.Date(), "Hourigansk_0005_TCM_SCFA_Fold_Change.csv", sep = "_"), row.names = FALSE)

# combine the fold change and stat results
dat_final_responders <- fc_results %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% responder_sig_hits) %>% 
  
  # remove Day0 values as they all equal zero
  filter(time != "Day 0") %>% 
  
  # merge pathway level information, used during plotting
  left_join(dat_pathways, join_by(metabolite), relationship = "many-to-many") %>% 
  relocate(Pathway, Pathway_adj, .after = metabolite)



##### High Responder Bubble plot prep #####

# keep only compounds that passed the cut off
dat_bubb_pre <- dat_final_responders %>% 
  arrange(Pathway_adj) %>% 
  mutate(log2FC_adj = case_when(log2FC > 2 ~ 2,
                                log2FC < -2 ~ -2, 
                                TRUE~ log2FC)) %>% 
  mutate(time = factor(time, levels = c("Baseline", "Day 2", "Day 8", "Day 14", "Day 28")))

legend <- data.frame(metabolite = c("zz_low", "zz_high"), neg_logp = c(0,3), log2FC_adj = c(-2, 2))


dat_bubb <- bind_rows(dat_bubb_pre, legend)

##### High Responder Make the bubble plot #####

bubble_plot <- ggplot(data=dat_bubb, 
                      mapping=aes(x=time, y=fct_rev(fct_inorder(metabolite)), 
                                  color=log2FC_adj, size=neg_logp)) +
  
  geom_point() +
  
  ggtitle("Metabolites") +
  xlab("Sample") +
  ylab("Metabolites") +
  labs(size="-Log(p)", color="log2FC D0 Paired") +
  
  scale_size(range=c(0,5)) +
  scale_x_discrete(position="top") +
  
  scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint = 0)+
  
  
  cowplot::theme_cowplot() + 
  theme(#legend.position = "none",
    axis.line  = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0, hjust=0, size = 6), axis.text.y = element_text(size=6),
    axis.ticks = element_blank(),
    text = element_text(size=8))

plot(bubble_plot)


ggsave(plot = bubble_plot, filename = paste(Sys.Date(), "lapointpo_TCM_responder_bubble_plot.svg", sep = "_"), width = 3.5, height = 10, units = "in")      



##### Color Palette Generation #####

metab_path <- unique(dat_pathways$Pathway_adj)

# determine the number of distinct colors needed
print( paste("I need", length(unique(metab_path)), "distinct colors"))

n <- length(unique(metab_path))


# use Rcolor Brewer to generate base colors from Dark2 color palette
base_colors <-  RColorBrewer::brewer.pal(8, "Dark2")

# adjust to accomidate the number of colors required for this particular analysis
expanded_colors <- colorRampPalette(base_colors)(n)

# match with grouping variable
pathway_colors <- setNames(expanded_colors, metab_path)

# specify that insignificant points will be grey
pathway_colors["insignificant"] <- "grey"



#### Preliminary Low Responder Stats and Calcs ####

dat_longer <- dat_t %>% 
  pivot_longer(4:ncol(dat_t), names_to = "Metabolite", values_to = "Value") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  mutate_at(c("Value", "SubjectID"), as.numeric) %>%
  
  
  # merge patient responder status to treatment (data provided by Lauren)
  left_join(dat_responders, join_by(SubjectID)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient) %>% 
  
  # filter data set for responders
  filter(Response == "Non-responder")



##### Low Responder Stat calculations #####
df <- dat_longer %>% 
  rename(time = Group, metabolite = Metabolite)

# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)


##### Low Responder Significance filtering #####
# this section utilizes significant hits determined through the Responder analysis. 
# Such that sig hits identified above for responders will be used here to filter the non-responder data.

# first run the significance test for the non-responders
# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)

# this strategy pairs the fold changes then takes the median value. 
Day0_value <- df %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  
  # calculate the median per metabolite
  group_by(time, metabolite) %>% 
  mutate(log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # keep only the columns required for downstream plotting
  select(metabolite, time, log2FC)


# calculate paired fold changes maintain patient specificity
fc_results_pat_spef <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  ungroup() %>% 
  select(Patient, metabolite, time, fold_change, log2FC)



# combine the statistics and fold change results
dat_final_non_responders <- fc_results %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  
  # filter for significant hits that were detected in the High Responder Analysis
  filter(metabolite %in% responder_sig_hits) %>% 
  filter(time != "Day 0") %>% 
  
  # merge pathway level information
  left_join(dat_pathways, join_by(metabolite), relationship = "many-to-many") %>% 
  relocate(Pathway, Pathway_adj, .after = metabolite) 



##### Low Responder Bubble plot prep #####

# keep only compounds that passed the cut off
dat_bubb_pre <- dat_final_non_responders %>% 
  arrange(Pathway_adj) %>% 
  
  mutate(log2FC_adj = case_when(log2FC > 2 ~ 2,
                                log2FC < -2 ~ -2, 
                                TRUE~ log2FC)) %>% 
  
  mutate(time = factor(time, levels = c("Baseline", "Day 2", "Day 8", "Day 14", "Day 28")))


legend <- data.frame(metabolite = c("zz_low", "zz_high"), neg_logp = c(0,3), log2FC_adj = c(-2, 2))


dat_bubb <- bind_rows(dat_bubb_pre, legend)

##### Low Responder Make the bubble plot #####

bubble_plot <- ggplot(data=dat_bubb, 
                      mapping=aes(x=time, y=fct_rev(fct_inorder(metabolite)), 
                                  color=log2FC_adj, size=neg_logp)) +
  
  geom_point() +
  
  ggtitle("Metabolites") +
  xlab("Sample") +
  ylab("Metabolites") +
  labs(size="-Log(p)", color="log2FC D0 Paired") +
  
  scale_size(range=c(0,5)) +
  scale_x_discrete(position="top") +
  
  scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint = 0)+
  
  
  cowplot::theme_cowplot() + 
  theme(#legend.position = "none",
    axis.line  = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0, hjust=0, size = 6), axis.text.y = element_text(size=6),
    axis.ticks = element_blank(),
    text = element_text(size=8))

plot(bubble_plot)


ggsave(plot = bubble_plot, filename = paste(Sys.Date(), "lapointpo_TCM_non_responder_bubble_plot.svg", sep = "_"), width = 3.5, height = 10, units = "in")      

#### Follow-Up MBP Analysis SOP ####
# below chunks are follow up to the preliminary responder status analysis
# the purpose was to essentially convert the bubble plots to MBPs
# break down MBPs by responder status will go into the supplement
# additionally, Bile Acid and Bile Acid Ratio data were added to the MBP


#### Total Cohort MBP Analysis ####
##### Tot Cohort Metab Stats and Calcs #####
# these fold changes will be used in the correlation matrix


dat_longer <- dat_t %>% 
  pivot_longer(4:ncol(dat_t), names_to = "Metabolite", values_to = "Value") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  mutate_at(c("Value", "SubjectID"), as.numeric) %>%
  
  
  # merge patient responder status to treatment (data provided by Lauren)
  left_join(dat_responders, join_by(SubjectID)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient) 




##### Tot Cohort Metab FC Calcs #####
df <- dat_longer %>% 
  rename(time = Group, metabolite = Metabolite)


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- df %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  group_by(time, metabolite) %>% 
  mutate(log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  select(metabolite, time, log2FC)

# calculate paired fold changes maintain patient specificity
fc_results_pat_spef <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  ungroup() %>% 
  select(Patient, metabolite, time, fold_change, log2FC)

# the fc_results_pat_spef dataframe will be used in the correlation network script. 
# export here
#write.csv(fc_results_pat_spef, paste(Sys.Date(), "Hourigansk_0005_TCM_SCFA_Fold_Change.csv", sep = "_"), row.names = FALSE)







##### Total Cohort Metab Calcs for MBP #####
dat_longer <- dat_t %>% 
  pivot_longer(4:ncol(dat_t), names_to = "Metabolite", values_to = "Value") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  mutate_at(c("Value", "SubjectID"), as.numeric) %>%
  
  
  # merge patient responder status to treatment (data provided by Lauren)
  left_join(dat_responders, join_by(SubjectID)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient) 
  

df <- dat_longer %>% 
  rename(time = Group, metabolite = Metabolite)

# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### Tot Cohort Metab Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# determine which metabolites show a significant difference for at least 1 time point
total_cohort_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()


# save the sig hits to use as filter in a sankey diagram
# write.csv(responder_sig_hits, "Responder_only_wilcox_day2_sig_hits.csv", row.names = F)


# count the number of significant hits
num_sig_hits <- length(responder_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)

# this strategy pairs the fold changes then takes the median value. 
Day0_value <- df %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  group_by(time, metabolite) %>% 
  
  mutate(log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  select(metabolite, time, log2FC)

# calculate paired fold changes maintain patient specificity
fc_results_pat_spef <- df %>% 
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>% 
  ungroup() %>% 
  select(Patient, metabolite, time, fold_change, log2FC)

# the fc_results_pat_spef dataframe will be used in the correlation network script. 
# export here
#write.csv(fc_results_pat_spef, paste(Sys.Date(), "Hourigansk_0005_TCM_SCFA_Fold_Change.csv", sep = "_"), row.names = FALSE)

# combine the results
# contains metabolomics data, scfa data for total cohort.
dat_metab_final <- fc_results %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% total_cohort_sig_hits) %>%
  left_join(dat_pathways, join_by(metabolite), relationship = "many-to-many") %>% 
  filter(time != "Day 0") %>% 
  dplyr::rename(Compound = metabolite, median_log2FC = log2FC, Class = Pathway_adj) %>% 
  select(-Pathway) %>% 
  mutate(Order_col = Class)




##### Tot Cohort Bile Acid Data Coercion #####

dat_bile_t <- dat_bile %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bile_merged <- dat_bile_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) 


##### Tot Cohort Bile Acid Data Stats #####

dat_bile_fc_pre <- dat_bile_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bile_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bile_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bile <- dat_bile_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### Total Cohort Bile Calcs for MBP #####
df <- dat_bile_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### Tot Cohort Bile Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# determine which metabolites show a significant difference for at least 1 time point
total_cohort_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()



# count the number of significant hits
num_sig_hits <- length(total_cohort_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
# contains metabolomics data, scfa data for total cohort.
dat_bile_final <- fc_results_bile %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% total_cohort_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
 
  left_join(dat_bla_conj, join_by(Compound)) %>% 

  mutate(Order_col = "zzz_bile",
         Class = paste("zz", Class, sep = "_")) %>% 
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)





##### Tot Cohort BLA Ratio Coercion #####

dat_bla_ratios_t <- dat_bla_ratios %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bla_ratios_merged <- dat_bla_ratios_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) 


##### Tot Cohort BLA Ratio Stats #####

dat_bla_ratios_fc_pre <- dat_bla_ratios_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bla_ratios_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bla_ratios_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bla_ratios <- dat_bla_ratios_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### Total Cohort BLA Ratio Calcs for MBP #####
df <- dat_bla_ratios_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### Tot Cohort BLA Ratio Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# determine which metabolites show a significant difference for at least 1 time point
total_cohort_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()



# count the number of significant hits
num_sig_hits <- length(total_cohort_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
dat_bla_ratios_final <- fc_results_bla_ratios %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% total_cohort_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
  
  left_join(dat_bla_ratios_conj, join_by(Compound)) %>% 
  
  mutate(Order_col = "zzz_bile_ratios",
         Class = paste("zzzz", Class, sep = "_")) %>%  
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)




##### Tot Cohort MBP prep ####

dat_final <- dat_metab_final %>% 
  bind_rows(dat_bile_final, dat_bla_ratios_final)

# define a capping value
cap <- 2

# make final adjustments specific to the matrix barplot
dat_mbp_pre <- dat_final %>% 
  mutate(median_log2FC = case_when(median_log2FC > cap ~ cap,
                                   median_log2FC < -cap ~ -cap,
                                   TRUE~median_log2FC)) %>% 
  mutate(time = factor(time, levels = unique(time)))


end <- ceiling(max(dat_mbp_pre$neg_logp))

# this data frame is used to generate a significance legend in the matrix barplot
legend <- data.frame("time" = rep("Baseline", end),
                    "Compound" = paste0("zz", sprintf("%02d", 1:end)), 
                    "Order_col" = paste0("ZZZZZ", sprintf("%02d", 1:end)), 
                    "Class" = paste0("ZZZZZ", sprintf("%02d", 1:end)),
                    "median_log2FC" = rep(0, end) , 
                    "neg_logp" = c(1: end))


# combine the legend data frame with the mbp data
dat_mbp <- bind_rows(dat_mbp_pre, legend) %>% 
  drop_na(Class) %>% 
  mutate(time = factor(time, levels = c("Baseline", "Day 2", "Day 8", "Day 14", "Day 28"))) %>% 
  arrange(Order_col) 



##### Tot Cohort Make the MBP #####

# generate the plot
my_tot_cohort_plot <- matrix_barplot(data = dat_mbp, aes(x = neg_logp, y = Compound, color = median_log2FC), 
                          sig_lines_at = NA,
                          show_strip.x.text = TRUE, show_strip.y.text = TRUE, 
                          remove.x.labels = TRUE, remove.y.labels = FALSE, switch = FALSE,
                          cols = vars(time), rows = vars(Class), orientation = "portrait") +
  
  scale_color_gradient2(low = 'blue', mid = 'grey75', high = 'red', midpoint = 0, limits = c(-cap, cap))+
  
  # customize the theme
  theme(#legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black",  size = 4),
    strip.background = element_rect(fill = "black", size = 0.5))

# view the plot
plot(my_tot_cohort_plot)

# save the condensed plot for the main figure
ggsave(plot = my_tot_cohort_plot, paste(Sys.Date(), "Hourigansk_0005_Matrix_Barplot_Arrange_Pathway_TOTAL_COHORT.svg", sep = "_"), width = 2.5, height = 6, units = "in")


# save the extended height plot for the supplement
ggsave(plot = my_tot_cohort_plot, paste(Sys.Date(), "Hourigansk_0005_MBP_TOT_COHORT_EXT_Height.svg", sep = "_"), width = 2.5, height = 9, units = "in")




#### High Responder MBP Analysis ####

##### Metab High Responder Finalization #####
# pull results from above responder analysis
# contains metabolomics data, scfa data for total cohort.
dat_metab_final <- dat_final_responders %>% 
  rename(Compound = metabolite, median_log2FC = log2FC, Class = Pathway_adj) %>% 
  select(-Pathway) %>% 
  mutate(Order_col = Class)




##### High Responder Bile Acid Data Coercion #####

dat_bile_t <- dat_bile %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bile_merged <- dat_bile_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) %>% 
  
  # responder analysis portion
  rename(Patient_id = Patient) %>% 
  left_join(dat_responders, join_by(Patient_id)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient_id) %>% 
  filter(Response == "Responder") %>% 
  select(-c(SubjectID, Response, Colonoscopy)) %>% 
  rename(Patient = Patient_id)


##### High Responder Bile Acid Data Stats #####

dat_bile_fc_pre <- dat_bile_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bile_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bile_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bile <- dat_bile_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### High Responder Bile Calcs for MBP #####
df <- dat_bile_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### High Responder Bile Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# determine which metabolites show a significant difference for at least 1 time point
bile_high_responder_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()



# count the number of significant hits
num_sig_hits <- length(bile_high_responder_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
# contains metabolomics data, scfa data for total cohort.
dat_bile_final <- fc_results_bile %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% bile_high_responder_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
  
  left_join(dat_bla_conj, join_by(Compound)) %>% 
  
  mutate(Order_col = "zzz_bile",
         Class = paste("zz", Class, sep = "_")) %>% 
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)






##### High Responder BLA Ratio Coercion #####

dat_bla_ratios_t <- dat_bla_ratios %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bla_ratios_merged <- dat_bla_ratios_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) %>% 
  
  # responder analysis portion
  rename(Patient_id = Patient) %>% 
  left_join(dat_responders, join_by(Patient_id)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient_id) %>% 
  filter(Response == "Responder") %>% 
  select(-c(SubjectID, Response, Colonoscopy)) %>% 
  rename(Patient = Patient_id)


##### High Responder BLA Ratio Stats #####

dat_bla_ratios_fc_pre <- dat_bla_ratios_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bla_ratios_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bla_ratios_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bla_ratios <- dat_bla_ratios_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### High Responder BLA Ratio Calcs for MBP #####
df <- dat_bla_ratios_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### High Responder BLA Ratio Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# determine which metabolites show a significant difference for at least 1 time point
bla_ratio_high_responder_sig_hits <- combinations %>% 
  filter(timepoint == "Day 0" | timepoint == "Day 2") %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH")) %>% 
  filter(p_adj < 0.1) %>% 
  ungroup() %>% 
  distinct(metabolite) %>% 
  pull()



# count the number of significant hits
num_sig_hits <- length(bla_ratio_high_responder_sig_hits)  

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
dat_bla_ratios_final <- fc_results_bla_ratios %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% bla_ratio_high_responder_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
  
  left_join(dat_bla_ratios_conj, join_by(Compound)) %>% 
  
  mutate(Order_col = "zzz_bile_ratios",
         Class = paste("zzzz", Class, sep = "_")) %>%  
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)






##### High Responder MBP prep ####

# merge the high responder data sets
dat_final <- dat_metab_final %>% 
  bind_rows(dat_bile_final, dat_bla_ratios_final)


# define a capping value
cap <- 2

# make final adjustments specific to the matrix barplot
dat_mbp_pre <- dat_final %>% 
  mutate(median_log2FC = case_when(median_log2FC > cap ~ cap,
                                   median_log2FC < -cap ~ -cap,
                                   TRUE~median_log2FC)) %>% 
  mutate(time = factor(time, levels = unique(time)))

# set manually to match the total cohort MBP
end <- 6

# this data frame is used to generate a significance legend in the matrix barplot
legend <- data.frame("time" = rep("Baseline", end),
                    "Compound" = paste0("zz", sprintf("%02d", 1:end)), 
                    "Order_col" = paste0("ZZZZZ", sprintf("%02d", 1:end)), 
                    "Class" = paste0("ZZZZZ", sprintf("%02d", 1:end)),
                    "median_log2FC" = rep(0, end) , 
                    "neg_logp" = c(1: end))


# combine the legend data frame with the real data
dat_mbp <- bind_rows(dat_mbp_pre, legend) %>% 
  drop_na(Class) %>% 
  mutate(time = factor(time, levels = c("Baseline", "Day 2", "Day 8", "Day 14", "Day 28"))) %>% 
  arrange(Order_col) 



##### High Responder Make the MBP #####

# generate the plot
my_high_respon_plot <- matrix_barplot(data = dat_mbp, aes(x = neg_logp, y = Compound, color = median_log2FC), 
                          sig_lines_at = NA,
                          show_strip.x.text = TRUE, show_strip.y.text = TRUE, 
                          remove.x.labels = TRUE, remove.y.labels = FALSE, switch = FALSE,
                          cols = vars(time), rows = vars(Class), orientation = "portrait") +
  
  scale_color_gradient2(low = 'blue', mid = 'grey75', high = 'red', midpoint = 0, limits = c(-cap, cap))+
  
  # customize the theme
  theme(#legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black",  size = 4),
    strip.background = element_rect(fill = "black", size = 0.5))

# view the plot
plot(my_high_respon_plot)

# save the condensed plot for the main figure
ggsave(plot = my_high_respon_plot, paste(Sys.Date(), "Hourigansk_0005_MBP_High_Responders.svg", sep = "_"), width = 2.5, height = 6, units = "in")


# save the extended height plot for the supplement
ggsave(plot = my_high_respon_plot, paste(Sys.Date(), "Hourigansk_0005_MBP_High_Responders_EXT_Height.svg", sep = "_"), width = 2.5, height = 9, units = "in")





#### Low Responder MBP Analysis ####




##### Metab Low Responder Finalization #####
# pull results from above responder analysis
# contains metabolomics data, scfa data for total cohort.
dat_metab_final <- dat_final_non_responders %>% 
  rename(Compound = metabolite, median_log2FC = log2FC, Class = Pathway_adj) %>% 
  select(-Pathway) %>% 
  mutate(Order_col = Class)



##### Low Responder Bile Acid Data Coercion #####

dat_bile_t <- dat_bile %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bile_merged <- dat_bile_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) %>% 
  
  # responder analysis portion
  rename(Patient_id = Patient) %>% 
  left_join(dat_responders, join_by(Patient_id)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient_id) %>% 
  filter(Response == "Non-responder") %>% 
  select(-c(SubjectID, Response, Colonoscopy)) %>% 
  rename(Patient = Patient_id)


##### Low Responder Bile Acid Data Stats #####

dat_bile_fc_pre <- dat_bile_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bile_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bile_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bile <- dat_bile_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### Low Responder Bile Calcs for MBP #####
df <- dat_bile_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### Low Responder Bile Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.

# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
# contains metabolomics data, scfa data for total cohort.
dat_bile_final <- fc_results_bile %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% bile_high_responder_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
  
  left_join(dat_bla_conj, join_by(Compound)) %>% 
  
  mutate(Order_col = "zzz_bile",
         Class = paste("zz", Class, sep = "_")) %>% 
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)






##### Low Responder BLA Ratio Coercion #####

dat_bla_ratios_t <- dat_bla_ratios %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_bla_ratios_merged <- dat_bla_ratios_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  relocate(Patient, Time_point, .after = Sample_name) %>% 
  
  # responder analysis portion
  rename(Patient_id = Patient) %>% 
  left_join(dat_responders, join_by(Patient_id)) %>% 
  relocate(SubjectID, Response, Colonoscopy, .after = Patient_id) %>% 
  filter(Response == "Non-responder") %>% 
  select(-c(SubjectID, Response, Colonoscopy)) %>% 
  rename(Patient = Patient_id)


##### Low Responder BLA Ratio Stats #####

dat_bla_ratios_fc_pre <- dat_bla_ratios_merged %>% 
  rename(time = Time_point) %>% 
  pivot_longer(4:ncol(dat_bla_ratios_merged), names_to = "metabolite", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) %>% 
  filter(!str_detect(Sample_name, "QC"))


# this strategy pairs the fold changes then takes the median value. 
Day0_value <- dat_bla_ratios_fc_pre %>% 
  filter(time == "Day 0") %>% 
  group_by(metabolite, Patient) %>% 
  mutate(day0_value = Value) %>% 
  ungroup() %>% 
  select(Patient, metabolite, day0_value)

# calculate paired fold changes
fc_results_bla_ratios <- dat_bla_ratios_fc_pre %>% 
  
  # merging method to pair the time point data
  left_join(Day0_value, join_by(metabolite, Patient), relationship = "many-to-many") %>% 
  
  # calculate the fold changes and transform with log2
  group_by(Patient, metabolite, time) %>% 
  mutate(fold_change = Value/day0_value,
         log2FC = log2(fold_change)) %>%
  
  # calculate the median values
  group_by(time, metabolite) %>% 
  mutate(median_log2FC = median(log2FC)) %>% 
  distinct(metabolite, .keep_all = TRUE) %>% 
  ungroup() %>% 
  
  # select needed columns for downstream analysis
  select(metabolite, time, median_log2FC)



##### Low Responder BLA Ratio Calcs for MBP #####
df <- dat_bla_ratios_fc_pre %>% 
  select(-Sample_name)



# remove the time point that every other time point will be compared back to.
combinations <- df %>% 
  filter(time != "Day 0") %>% 
  distinct(metabolite, time) %>% 
  rename(timepoint = time)




##### Low Responder BLA Ratio Day2 Sig Only #####
# this section focuses specifically on significant differences between Day2 and Day0
# aimed to trim the data set further. Sig hits will be used in the MBP.


# calculate the paired p-values using the wilcox.test()
p_results <- combinations %>% 
  rowwise() %>% 
  mutate(p_value = wilcox_test(metabolite, timepoint),
         p_adj = p.adjust(p_value, method = "BH"),
         neg_logp = -log10(p_adj)) %>% 
  ungroup() %>% 
  rename(time = timepoint) %>% 
  select(metabolite, time, neg_logp)


# combine the results
dat_bla_ratios_final <- fc_results_bla_ratios %>% 
  left_join(p_results, join_by(metabolite, time)) %>% 
  filter(metabolite %in% bla_ratio_high_responder_sig_hits) %>%
  
  filter(time != "Day 0") %>% 
  
  rename(Compound = metabolite) %>% 
  
  left_join(dat_bla_ratios_conj, join_by(Compound)) %>% 
  
  mutate(Order_col = "zzz_bile_ratios",
         Class = paste("zzzz", Class, sep = "_")) %>%  
  
  select(Compound, time, Class, Order_col, median_log2FC, neg_logp)






##### Low Responder MBP prep ####

# merge the low responder data sets
dat_final <- dat_metab_final %>% 
  bind_rows(dat_bile_final, dat_bla_ratios_final)


# define a capping value
cap <- 2

# make final adjustments specific to the matrix barplot
dat_mbp_pre <- dat_final %>% 
  mutate(median_log2FC = case_when(median_log2FC > cap ~ cap,
                                   median_log2FC < -cap ~ -cap,
                                   TRUE~median_log2FC)) %>% 
  mutate(time = factor(time, levels = unique(time)))

# set manually to match the total cohort MBP
end <- 6

# this data frame is used to generate a significance legend in the matrix barplot
legend <- data.frame("time" = rep("Baseline", end),
                    "Compound" = paste0("zz", sprintf("%02d", 1:end)), 
                    "Order_col" = paste0("ZZZZZ", sprintf("%02d", 1:end)), 
                    "Class" = paste0("ZZZZZ", sprintf("%02d", 1:end)),
                    "median_log2FC" = rep(0, end) , 
                    "neg_logp" = c(1: end))


# combine the legend data frame with the real data
dat_mbp <- bind_rows(dat_mbp_pre, legend) %>% 
  drop_na(Class) %>% 
  mutate(time = factor(time, levels = c("Baseline", "Day 2", "Day 8", "Day 14", "Day 28"))) %>% 
  arrange(Order_col) 



##### Low Responder Make the MBP #####

# generate the plot
my_low_respon_plot <- matrix_barplot(data = dat_mbp, aes(x = neg_logp, y = Compound, color = median_log2FC), 
                                      sig_lines_at = NA,
                                      show_strip.x.text = TRUE, show_strip.y.text = TRUE, 
                                      remove.x.labels = TRUE, remove.y.labels = FALSE, switch = FALSE,
                                      cols = vars(time), rows = vars(Class), orientation = "portrait") +
  
  scale_color_gradient2(low = 'blue', mid = 'grey75', high = 'red', midpoint = 0, limits = c(-cap, cap))+
  
  # customize the theme
  theme(#legend.position = "none",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black",  size = 4),
    strip.background = element_rect(fill = "black", size = 0.5))

# view the plot
plot(my_low_respon_plot)

# save the condensed plot
ggsave(plot = my_low_respon_plot, paste(Sys.Date(), "Hourigansk_0005_MBP_Low_Responders.svg", sep = "_"), width = 2.5, height = 6, units = "in")


# save the extended height plot for the supplement
ggsave(plot = my_low_respon_plot, paste(Sys.Date(), "Hourigansk_0005_MBP_Low_Responders_EXT_Height.svg", sep = "_"), width = 2.5, height = 9, units = "in")





#### Display all MBPs in a Grid ####

# generate the extended grid
# relative heights are proportional
# will need to balance the relative heights with the plot height that is set in ggsave()
extended_grid <- cowplot::plot_grid(my_tot_cohort_plot, my_tot_cohort_plot, my_high_respon_plot, my_low_respon_plot, 
                                    ncol = 1, 
                                    align = "v", 
                                    axis = "lr", 
                                    rel_heights = c(1.1, 1.7, 1.7, 1.7))


# save the extended grid with legible metabolite labels
ggsave(plot = extended_grid, paste(Sys.Date(), "Hourigansk_0005_Grid_MBP_Long.svg", sep = "_"), width = 3, height = 34, units = "in")



