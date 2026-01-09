
#### SOP ####
# data pre-processing for pair wise correlation between omics layers


#### Required Libraries ####
library(tidyverse)


#### Pivot Longer ####

# pivot longer, without scaling
# remove ME0, which did not correspond to a specific cluster
eigengenes_longer <- eigengenes %>% 
  pivot_longer(2:ncol(eigengenes), names_to = "Feature", values_to = "Value") %>% 
  mutate(Layer = "eigengene") %>% 
  filter(Feature != "ME0")

# pivot longer and z-scale
# filter for high RSD metabolites
dat_metab_longer <- dat_metab %>% 
  pivot_longer(2:ncol(dat_metab), names_to = "Acquired_Sample_ID", values_to = "Value") %>% 
  filter(Metabolite %in% metab_rsd$Metabolite) %>% 
  group_by(Metabolite) %>% 
  mutate(Value = as.numeric(scale(Value, center = TRUE, scale = TRUE))) %>% 
  ungroup() %>% 
  mutate(Layer = "metabolomics") %>% 
  rename(Feature = Metabolite)

# pivot longer and z-scale
dat_bile_longer <- dat_bile %>% 
  pivot_longer(2:ncol(dat_bile), names_to = "Acquired_Sample_ID", values_to = "Value") %>% 
  group_by(Bile_acid) %>% 
  mutate(Value = as.numeric(scale(Value, center = TRUE, scale = TRUE))) %>% 
  ungroup() %>% 
  mutate(Layer = "bile acids") %>% 
  rename(Feature = Bile_acid) 



# pivot longer 
single_cell_longer <- single_cell %>% 
  pivot_longer(2:ncol(single_cell), names_to = "Separate_me", values_to = "Value") %>% 
  separate(Separate_me, into = c("Patient_id", "Time_point")) %>% 
  mutate(Layer = "SingleCell") %>% 
  filter(Time_point == "T2")

# determine the patient ids in the single cell data
single_cell_patients <- single_cell_longer %>% 
  distinct(Patient_id) %>% 
  pull()


# select needed columns for filtering
data_for_filtering_simplified <- data_for_filtering %>% 
  left_join(dat_prot_meta, by = join_by(EntrezGeneSymbol)) %>% 
  select(Target, EntrezGeneSymbol)

# pre-filter for proteins of interest as defined by Victor
dat_prot_filtered <- dat_prot %>% 
  select(SampleId, SampleType, any_of(data_for_filtering_simplified$Target)) %>% 
  filter(SampleType == 'Sample') %>% 
  filter(!str_detect(SampleId, "QC")) %>% 
  select(-SampleType)

# pivot the filtered data frame longer
dat_prot_longer <- dat_prot_filtered %>% 
  pivot_longer(2:ncol(dat_prot_filtered), names_to = "Target", values_to = "Value") %>% 
  mutate(Layer = "proteomics")


#### Merge -omics and Meta Data ####

dat_meta_trim <- dat_meta %>% 
  select(Submitted_Sample_ID, Submitted_Sample_Name, Patient_id, Time_point) %>% 
  distinct(Submitted_Sample_ID, .keep_all = TRUE)

eigengenes_merged <- eigengenes_longer %>% 
  left_join(dat_meta_trim, join_by(Submitted_Sample_Name)) %>% 
  relocate(any_of(meta_names), Layer, .after = Submitted_Sample_Name)

dat_metab_merged <- dat_metab_longer %>% 
  left_join(dat_meta, join_by(Acquired_Sample_ID), relationship = 'many-to-many') %>% 
  relocate(any_of(meta_names), .after = Acquired_Sample_ID) %>% 
  relocate(Layer, Feature, .before = Value) %>% 
  select(-Acquired_Sample_ID, -MS_method)

dat_bile_merged <- dat_bile_longer %>% 
  left_join(dat_meta, join_by(Acquired_Sample_ID), relationship = 'many-to-many') %>% 
  relocate(any_of(meta_names), .after = Acquired_Sample_ID) %>% 
  relocate(Layer, Feature, .before = Value)%>% 
  select(-Acquired_Sample_ID, -MS_method)


dat_prot_merged <- dat_prot_longer %>% 
  left_join(dat_samp_meta, join_by(SampleId)) %>%
  rename(Feature = Target) %>% 
  select(Patient_id, Time_point, Layer, Feature, Value)

dat_prot_diagnostic <- dat_prot_longer %>% 
  left_join(dat_samp_meta, join_by(SampleId)) %>% 
  filter(is.na(Patient_id) == TRUE) %>% 
  distinct(SampleId)

#write.csv(dat_prot_diagnostic, "CHI-25-007_Samples_Without_Meta_Data.csv", row.names = F)

length(unique(dat_prot_merged$SampleId))
length(unique(dat_samp_meta$SampleId))
single_cell_patients



#### Merge the -omics+Clin data sets ####

# merge the standardized omics data, relevant meta data included
dat_semi_final <- bind_rows(eigengenes_merged, dat_metab_merged, dat_bile_merged, dat_prot_merged) %>% 
  select(Patient_id, Time_point, Layer, Feature, Value) %>% 
  filter(Time_point == "Day 2") %>% 
  filter(Patient_id %in% unique(dat_prot_merged$Patient_id))


dat_final <- bind_rows(dat_semi_final, single_cell_longer)


# save the file. 
# this output will be called in script 2 for correlation matrix calculation
write.csv(dat_final, file = paste(Sys.Date(), "Cor_ready.csv", sep = "_"), row.names = F)




# prepare Immune Panel vs Metabolomics specific data frame
# this has to be done separately in order to determine specific patients that have both data layers

out_1 <- bind_rows(dat_metab_merged, dat_bile_merged) %>% 
  select(Patient_id, Time_point, Layer, Feature, Value) %>% 
  filter(Time_point == "Day 2") %>% 
  bind_rows(dat_semi_final, single_cell_longer) %>% 
  filter(Patient_id %in% single_cell_patients)


# save the file. 
# this output will be called in script 2 for correlation matrix calculation
write.csv(out_1, file = paste(Sys.Date(), "Immune_Metab_cor_ready.csv", sep = "_"), row.names = F)



# prepare Immune Panel vs Eigengene specific data frame
# this has to be done separately in order to determine specific patients that have both data layers

out_2 <- eigengenes_merged %>% 
  filter(Patient_id %in% single_cell_patients) %>% 
  filter(Time_point == "Day 2") %>% 
  bind_rows(single_cell_longer)


# save the file. 
# this output will be called in script 2 for correlation matrix calculation
write.csv(out_2, file = paste(Sys.Date(), "Immune_Eig_cor_ready.csv", sep = "_"), row.names = F)


