
#### SOP ####
# the purpose of this analysis was to z-scale all omics layers, subtract Day2-Day0
# then prepare export, which will be used for pair-wise correlations


#### Required Libraries ####
library(tidyverse)
library(tibble)
library(janitor)


##### TSN Metabolomics Data Coercion #####
dat_metab_t <- dat_metab_tsn %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  rownames_to_column(var = "Sample_name")


dat_metab_tsn_longer <- dat_metab_t %>% 
  pivot_longer(4:ncol(dat_metab_t), names_to = "Metabolite", values_to = "Value") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  relocate(SubjectID, .after = Patient) %>% 
  mutate_at(c("Value", "SubjectID"), as.numeric)  



  
paired_patients <- dat_metab_tsn_longer %>%
  filter(Group %in% c("Day 0", "Day 2")) %>%
  distinct(Patient, Group) %>%
  count(Patient) %>% 
  filter(n == 2) %>% 
  pull(Patient)


dat_metab_tsn_mod <- dat_metab_tsn_longer %>% 
  filter(Patient %in% paired_patients) %>% 
  filter(Group %in% c("Day 0", "Day 2")) %>% 
  group_by(Metabolite) %>% 
  mutate(Value_scaled = as.numeric(scale(Value, center = TRUE, scale = TRUE))) %>% 
  mutate(Value_scaled = min(Value_scaled)+Value_scaled) %>% 
  group_by(Patient, Metabolite) %>% 
  summarise(Value_diff = Value_scaled[Group == "Day 2"] - Value_scaled[Group == "Day 0"]) %>% 
  mutate(Layer = "Metabolomics") %>% 
  rename(Feature = Metabolite) %>% 
  select(Patient, Layer, Feature, Value_diff)



##### TSN Bile Acid Data Coercion #####
dat_bile_t <- dat_bile_tsn %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(row_number = 1) %>% 
  select(-Group) %>% 
  rename(Group = `Sample Type`, Sample_name = `Sample Name`)%>% 
  filter(!str_detect(Group, "QC"))

dat_bile_merged <- dat_bile_t %>% 
  left_join(dat_bile_meta, join_by(Sample_name)) %>% 
  select(-Time_point) %>% 
  relocate(Patient, .before = Sample_name)
  


dat_bile_tsn_longer <- dat_bile_merged %>% 
  pivot_longer(4:ncol(dat_bile_merged), names_to = "Bile_acid", values_to = "Value") %>% 
  mutate_at("Value", as.numeric) 



dat_bile_tsn_mod <- dat_bile_tsn_longer %>% 
  filter(Patient %in% paired_patients) %>% 
  filter(Group %in% c("Day 0", "Day 2")) %>% 
  group_by(Bile_acid) %>% 
  mutate(Value_scaled = as.numeric(scale(Value, center = TRUE, scale = TRUE))) %>% 
  mutate(Value_scaled = min(Value_scaled)+Value_scaled) %>% 
  group_by(Patient, Bile_acid) %>% 
  summarise(Value_diff = Value_scaled[Group == "Day 2"] - Value_scaled[Group == "Day 0"]) %>% 
  mutate(Layer = "Bile acids") %>% 
  rename(Feature = Bile_acid) %>% 
  select(Patient, Layer, Feature, Value_diff)


##### Eigengene Coercion #####

dat_meta_trim <- dat_meta %>% 
  select(Submitted_Sample_ID, Submitted_Sample_Name, Patient_id, Time_point) %>% 
  distinct(Submitted_Sample_ID, .keep_all = TRUE)


# pivot longer, without scaling
eigengenes_longer <- eigengenes %>% 
  pivot_longer(2:ncol(eigengenes), names_to = "Feature", values_to = "Value_scaled")


eigengenes_mod <- eigengenes_longer %>% 
  
  # merge trimmed meta data
  left_join(dat_meta_trim, join_by(Submitted_Sample_Name)) %>% 
  
  filter(Feature %in% c("ME01", "ME02", "ME03", "ME16", "ME20", "ME21", "ME24")) %>% 
  rename(Patient = Patient_id, Group = Time_point) %>% 
  mutate(Value_scaled = min(Value_scaled)+Value_scaled) %>% 
  group_by(Patient, Feature) %>% 
  summarise(Value_diff = Value_scaled[Group == "Day 2"] - Value_scaled[Group == "Day 0"]) %>% 
  mutate(Layer = "Eigengenes") %>% 
  select(Patient, Layer, Feature, Value_diff)


##### Data Export ####

out_isl <- bind_rows(dat_metab_tsn_mod, dat_bile_tsn_mod, eigengenes_mod)

write.csv(out_isl, paste(Sys.Date(), "Cor_Ready.csv", sep = "_"), row.names = F)
