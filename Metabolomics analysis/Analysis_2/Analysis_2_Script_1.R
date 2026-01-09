
#### SOP ####

# the purpose of this script is  to look at Pseudomonadota day 2 and correlate to Day 0 metabolites.
# this is likely the script you are searching for. As this analysis was solidified for the final manuscript.

# below are the coercion steps to prepare the data for spearman correlation. 
# then calculations are available for the spearman correlation.
# strategy for visualization of the correlation results using a U-plot are available as well.


#### Required Libraries ####
library(tidyverse)
library(ggpubr)


#### Meta Data Coercion ####

dat_clin_meta_trim <- dat_clin_meta %>% 
  select(Sample, SubjectID, Timepoint)

paired_patients <- dat_clin_meta_trim %>% 
  filter(all(c("Day 0", "Day 2") %in% Timepoint)) %>% 
  filter(Timepoint %in% c("Day 0", "Day 2"))



#### Metagenomics Coercion ####

dat_genomics_longer <- dat_genomics %>% 
  select(family, c(10:ncol(dat_genomics))) %>% 
  pivot_longer(2:ncol(pseudomonadota), names_to = "Sample", values_to = "Value_raw") %>% 
  left_join(dat_clin_meta_trim, join_by(Sample)) %>%
  mutate(Layer = "metagenomics") %>% 
  relocate(Timepoint, .after = Sample)

# Normalization and scaling
dat_genomics_scaled <- dat_genomics_longer %>% 
  
  # median normalize by sample
  group_by(Sample) %>% 
  mutate(Value_norm = (Value_raw +1) / median((Value_raw+1), na.rm = TRUE)) %>% 
  
  # log2 transform by family
  group_by(family) %>% 
  mutate(Value_scaled = as.numeric(log2(Value_norm)))


# Density plot
ggplot(dat_genomics_scaled, aes(x = Value_scaled, color = Sample)) +
  #scale_x_continuous(limits = c(-2,4))+
  facet_wrap(~Timepoint)+
  geom_density(alpha = 0.3, size = 1) +   # alpha for transparency
  theme_minimal() +
  theme(legend.position = "none")+
  labs(
    title = "Day0 Density Plot by Sample",
    x = "Z-scaled Value",
    y = "Density"
  )



#### Pseudomonadota Selection ####

pseudomonadota_day2 <- dat_genomics_scaled %>% 
  
  # filter for pseudomonadota measurements in patients with time paired samples only
  filter(family == "Pseudomonadota" & Sample %in% paired_patients$Sample & SubjectID != 163114) %>% 
  
  # this analysis is focused on Day 2
  filter(Timepoint == "Day 2") %>% 
  
  rename(Feature = family) %>%
  group_by(SubjectID) %>%
  distinct(SubjectID, .keep_all = TRUE) %>% 
  select(SubjectID, Feature, Value_scaled)



#### Bile Coercion ####

dat_bile_t <- dat_bile %>% 
  t() %>% 
  as.data.frame() %>% 
  janitor::row_to_names(row_number = 1) %>% 
  tibble::rownames_to_column(var = "Sample_name") %>% 
  rename(Timepoint = `Sample Type`)

dat_bile_longer <- dat_bile_t %>% 
  pivot_longer(5:ncol(dat_bile_t), names_to = "Feature", values_to = "Value_raw") %>% 
  filter(!str_detect(Group, "QC")) %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  relocate(SubjectID, .after = Patient) %>% 
  mutate_at(c("Value_raw"), as.numeric) %>% 
  left_join(dat_clin_meta_trim, join_by(SubjectID, Timepoint)) %>% 
  mutate(Layer = "metabolomics") %>% 
  relocate(Sample, Layer, .after = SubjectID)


dat_bile_scaled <- dat_bile_longer %>%
  group_by(Sample) %>% 
  mutate(Value_norm = Value_raw / median(Value_raw, na.rm = TRUE)) %>% 
  group_by(Feature) %>% 
  mutate(Value_scaled = as.numeric(scale(Value_raw, center = TRUE, scale = TRUE))) %>% 
  ungroup()


# Density plot
ggplot(dat_bile_scaled, aes(x = Value_scaled, color = Group)) +
  #scale_x_continuous(limits = c(-2,4))+
  facet_wrap(~Timepoint)+
  geom_density(alpha = 0.3, size = 1) +   # alpha for transparency
  theme_minimal() +
  theme(legend.position = "none")+
  labs(
    title = "Day0 Density Plot by Sample",
    x = "Z-scaled Value",
    y = "Density"
  )


dat_bile_day0 <- dat_bile_scaled %>% 
  filter(Timepoint == "Day 0") %>% 
  filter(SubjectID %in% paired_patients$SubjectID) %>% 
  filter(SubjectID != 163114) %>% 
  group_by(SubjectID, Feature) %>%
  distinct(SubjectID, .keep_all = TRUE) %>% 
  select(SubjectID, Feature, Value_scaled)




#### Metabolomics coercion ####

dat_metab_t <- dat_metab %>% 
  t() %>% 
  as.data.frame() %>% 
  janitor::row_to_names(row_number = 1) %>% 
  tibble::rownames_to_column(var = "Sample_name")


dat_metab_longer <- dat_metab_t %>% 
  pivot_longer(4:ncol(dat_metab_t), names_to = "Feature", values_to = "Value_raw") %>% 
  filter(Group != "QC") %>% 
  select(-Sample_name) %>% 
  mutate(SubjectID = gsub("pat", "", Patient)) %>% 
  relocate(SubjectID, .after = Patient) %>% 
  mutate_at(c("Value_raw"), as.numeric) %>% 
  rename(Timepoint = Group) %>% 
  left_join(dat_clin_meta_trim, join_by(SubjectID, Timepoint)) %>% 
  mutate(Layer = "metabolomics") %>% 
  relocate(Sample, Layer, .after = SubjectID)



dat_metab_scaled <- dat_metab_longer %>% 
  group_by(Feature) %>% 
  mutate(Value_scaled = as.numeric(scale(Value_raw, center = TRUE, scale = TRUE))) %>% 
  ungroup()



# Density plot
ggplot(dat_metab_scaled, aes(x = Value_scaled, color = Sample)) +
  #scale_x_continuous(limits = c(-2,4))+
  facet_wrap(~Timepoint)+
  geom_density(alpha = 0.3, size = 1) +   # alpha for transparency
  theme_minimal() +
  theme(legend.position = "none")+
  labs(
    title = "Density Plot by Sample",
    x = "Z-scaled Value",
    y = "Density"
  ) 



dat_metab_day0 <- dat_metab_scaled %>% 
  filter(Timepoint == "Day 0") %>% 
  filter(SubjectID %in% paired_patients$SubjectID) %>% 
  filter(SubjectID != 163114) %>% 
  group_by(SubjectID, Feature) %>%
  distinct(SubjectID, .keep_all = TRUE) %>% 
  select(SubjectID, Feature, Value_scaled)



#### Data merge and export ####

dat_combined_longer <- bind_rows(pseudomonadota_day2, dat_metab_day0, dat_bile_day0)


dat_ready_wide <- dat_combined_longer %>% 
  filter(SubjectID != 163114) %>% 
  pivot_wider(names_from = Feature, values_from = Value_scaled)


write.csv(dat_ready_wide, paste(Sys.Date(), "Hourigansk_0005_Cor_Ready_Wide.csv", sep = "_"), row.names = F)


dat_ready_long <- dat_combined_longer %>% 
  pivot_wider(names_from = Feature, values_from = Value_scaled) %>% 
  filter(SubjectID != 163114) %>% 
  pivot_longer(3:245, names_to = "Compound", values_to = "Value_scaled") 



write.csv(dat_ready_long, paste(Sys.Date(), "Hourigansk_0005_Cor_Ready_Long.csv", sep = "_"), row.names = F)




#### Statistical Analysis: Spearman Correlation ####

corr_p_values <- dat_ready_long %>% 
  
  group_by(Compound) %>% 
  summarize(
    rho = round(cor(Pseudomonadota, Value_scaled, method = "spearman"),3),
    p = cor.test(Pseudomonadota, Value_scaled, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  
  # adjust the raw p-value using Benjamini Hochberg method
  mutate(p_adj = p.adjust(p, method = "BH")) %>% 
  
  # transform the raw p-value
  mutate(neg_log_p = round(-log10(p), 3)) %>% 
  
  arrange(p) 


write.csv(corr_p_values, paste(Sys.Date(), "Hourigansk_0005_Psuedomondata_vs_metab.csv", sep = "_"), row.names = F)




#### U-Plot Preparation ####

my_threshold <- corr_p_values %>%
  filter(p_adj <= 0.1) %>%
  summarise(max_p = max(p, na.rm = TRUE)) %>%
  pull(max_p)

my_threshold <- 0.05

cor_df <- corr_p_values %>% 
  left_join(dat_pathways, join_by(Compound)) %>% 
  mutate(Significance = case_when(p <= my_threshold & rho > 0 ~ "positive significant",
                                  p <= my_threshold & rho < 0 ~ "negative significant",
                                  TRUE ~ "insignificant"),
         Color = case_when(p <= my_threshold & rho > 0 ~ Pathway_adj,
                           p <= my_threshold & rho < 0 ~ Pathway_adj,
                           TRUE ~ "insignificant"),
         Label = case_when(neg_log_p > 1 & rho > 0  ~ Compound,
                           neg_log_p > 1 & rho < 0  ~ Compound,
                           TRUE ~ "")) 



write.csv(cor_df, paste(Sys.Date(), "Hourigansk_0005_pseudomonadota_spearman_cor_results.csv", sep = "_"), row.names = F)

# count correlates

cor_df %>%
  filter(Significance == "positive significant") %>%
  nrow()

cor_df %>%
  filter(Significance == "negative significant") %>%
  nrow()


cor_df %>%
  filter(Significance == "insignificant") %>%
  nrow()



##### Color Palette Generation #####

path <- unique(dat_pathways$Pathway_adj)

# determine the number of distinct colors needed
print( paste("I need", length(unique(path)), "distinct colors"))

n <- length(unique(path))


# use Rcolor Brewer to generate base colors from Dark2 color palette
base_colors <-  RColorBrewer::brewer.pal(8, "Dark2")

# adjust to accomidate the number of colors required for this particular analysis
expanded_colors <- colorRampPalette(base_colors)(n)

# match with grouping variable
pathway_colors <- setNames(expanded_colors, path)

# specify that insignificant points will be grey
pathway_colors["insignificant"] <- "grey"




###### Correlation U-Plot ######


## Create scatterplot using the ggplot function ##
cor_plot <- ggplot(cor_df, aes(rho, neg_log_p, col = Color, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = pathway_colors)+
  
  #scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), limits = c(-0.5, 0.5))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,3))+ 
  
  xlab("Spearman Rho") + 
  ylab("-log10(p-value)") + 
  
  
  geom_abline(intercept = -log10(my_threshold), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8, color = "black"), 
        axis.title = element_text(size = 8, color = "black"), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.ticks = element_line(size = 0.25, color = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  box.padding = 0.7)


# view the plot
plot(cor_plot)

# save the plot
ggsave(paste(Sys.Date(), "Spearman_Corr_UPlot.svg", sep = "_"), bg = "transparent",
       plot = cor_plot, width = 2, height = 3.5, dpi = 300)






