
#### SOP ####
# perform pair-wise correlations between omics layers
# visualize correlation results using a sankey diagram


#### Package Installation ####

# run this line every time there is an update to the MetLipNet package
devtools::install_github("IDSS-NIAID/MetLipNet")


#### Required Libraries ####

# run this chunk to open libraries containing functions to execute the script
library(MetLipNet)
library(tidyverse)
library(corrr)
library(ggsankey)




##### Data pre processing #####
meta_data_names <- NULL

dat_ready <- dat_in %>% 
  select(-c(Layer, Time_point)) %>% 
  pivot_wider(id_cols = Feature, names_from = Patient_id, values_from = Value) 


##### Correlations #####
# This function calculates the pairwise Pearson correlations for metabolite intensities 
# and provides the correlation estimate and p-value for each pair.

# view the help page
?cal_met_cor

# The returned data frame will be ready for the igraph function.
dat_2_cor_result <- cal_met_cor(data = dat_ready, 
                                meta_cols = meta_data_names, 
                                identifier_col = "Feature", 
                                method = "spearman",
                                cor_threshold = 0.5,
                                p_threshold = 0.05) 





#### Immune to Metab Sankey Diagram Prep ####

# need to make a key that will be used to match 'from' and 'to' columns back to the layer
feature_layer <- dat_in %>% 
  select(Feature, Layer) %>% 
  mutate(Layer = case_when(Layer == "bile acids" ~ "metabolomics",
                           TRUE ~ Layer))

# tag each feature with its data layer
dat_2_cor_result_merged <- dat_2_cor_result %>% 
  left_join(feature_layer, by = c("from" = "Feature")) %>% 
  rename(from_Layer = Layer) %>% 
  left_join(feature_layer, by = c("to" = "Feature")) %>% 
  rename(to_Layer = Layer)

# remove within-layer correlations
corr_between_layers <- dat_2_cor_result_merged %>% 
  filter(from_Layer != to_Layer) %>% 
  
  # drop duplicates - each unique link should only appear once
  distinct(from, from_Layer, to, to_Layer, .keep_all = TRUE)





sankey_ready <- corr_between_layers %>% 
  filter(from_Layer == "SingleCell" & to_Layer == "metabolomics") %>%
  filter(p< 1e-2 & cor > 0.6) %>% 
  rename(singlecell = from, metab = to, corr_eig_metab = cor) 



sankey_long <- sankey_ready %>% 
  make_long(singlecell, metab) %>% 
  mutate(fill = case_when(x == "singlecell" ~ node, TRUE ~ "")) %>% 
  left_join(dat_pathway, join_by(node)) %>% 
  group_by(x) %>% 
  arrange(x==2, Pathway_adj, node) %>% 
  mutate(node = factor(node, levels = unique(node))) %>% 
  ungroup()



##### Generate the Sankey Diagram #####

# make the plot
my_sankey <- ggplot(sankey_long, aes(x = x, next_x = next_x, 
                                     node = node, next_node = next_node, 
                                     label = node, fill = fill))+
  geom_sankey(flow.alpha = 0.7, node.color = "gray30")+
  geom_sankey_label(size = 3, fill = "white", color = "black")+
  theme_sankey(base_size = 10)+
  theme(legend.position = "none")

# view the plot
plot(my_sankey)

# save the plot
ggsave(plot = my_sankey, paste(Sys.Date(), "Sankey.svg", sep = "_" ), width = 15, height = 18, units = "in")








#### Somalogic Proteomics Sankey Diagram Prep ####

# need to make a key that will be used to match 'from' and 'to' columns back to the layer
feature_layer <- dat_final %>% 
  select(Feature, Layer) %>% 
  mutate(Layer = case_when(Layer == "bile acids" ~ "metabolomics",
                           TRUE ~ Layer))

# tag each feature with its data layer
dat_2_cor_result_merged <- dat_2_cor_result %>% 
  left_join(feature_layer, by = c("from" = "Feature")) %>% 
  rename(from_Layer = Layer) %>% 
  left_join(feature_layer, by = c("to" = "Feature")) %>% 
  rename(to_Layer = Layer)

# remove within-layer correlations
corr_between_layers <- dat_2_cor_result_merged %>% 
  filter(from_Layer != to_Layer) %>% 
  
  # drop duplicates - each unique link should only appear once
  distinct(from, from_Layer, to, to_Layer, .keep_all = TRUE)




sankey_ready <- corr_between_layers %>% 
  filter(from_Layer == "proteomics" & to_Layer == "metabolomics") %>%
  filter(p< 1e-2 & cor > 0.6) %>% 
  rename(prot = from, metab = to, corr_eig_metab = cor) 



sankey_long <- sankey_ready %>% 
  make_long(prot, metab) %>% 
  mutate(fill = case_when(x == "prot" ~ node, TRUE ~ "")) %>% 
  left_join(dat_pathway, join_by(node)) %>% 
  group_by(x) %>% 
  arrange(x==2, Pathway_adj, node) %>% 
  mutate(node = factor(node, levels = unique(node))) %>% 
  ungroup()



##### Generate the Sankey Diagram #####
my_sankey <- ggplot(sankey_long, aes(x = x, next_x = next_x, 
                                     node = node, next_node = next_node, 
                                     label = node, fill = fill))+
  geom_sankey(flow.alpha = 0.7, node.color = "gray30")+
  geom_sankey_label(size = 3, fill = "white", color = "black")+
  #scale_fill_brewer(palette = "Dark2")+
  theme_sankey(base_size = 10)+
  theme(legend.position = "none")

plot(my_sankey)


# save the plot
ggsave(plot = my_sankey, paste(Sys.Date(), "Sankey.svg", sep = "_" ), width = 15, height = 15, units = "in")




#### Eigengene to Immune Specific Analysis ####

#### Read in Data ####

# Day 8 immune panel data, Day 2 Eigengene Module Data
dat_in <- read.csv("2026-01-05_Immune_Eig_cor_ready.csv")



#### Analysis: Correlations ####

##### Data prep #####
meta_data_names <- NULL

dat_ready <- dat_in %>% 
  select(-c(Layer, Time_point)) %>% 
  pivot_wider(id_cols = Feature, names_from = Patient_id, values_from = Value) 


##### Correlations #####
# This function calculates the pairwise Pearson correlations for metabolite intensities 
# and provides the correlation estimate and p-value for each pair.

# view the help page
?cal_met_cor

# The returned data frame will be ready for the igraph function.
dat_2_cor_result <- cal_met_cor(data = dat_ready, 
                                meta_cols = meta_data_names, 
                                identifier_col = "Feature", 
                                method = "spearman",
                                cor_threshold = 0.5,
                                p_threshold = 0.05) 



# need to make a key that will be used to match 'from' and 'to' columns back to the layer
feature_layer <- dat_in %>% 
  select(Feature, Layer)

# tag each feature with its data layer
dat_2_cor_result_merged <- dat_2_cor_result %>% 
  left_join(feature_layer, by = c("from" = "Feature")) %>% 
  rename(from_Layer = Layer) %>% 
  left_join(feature_layer, by = c("to" = "Feature")) %>% 
  rename(to_Layer = Layer)

# remove within-layer correlations
corr_between_layers <- dat_2_cor_result_merged %>% 
  filter(from_Layer != to_Layer) %>% 
  
  # drop duplicates - each unique link should only appear once
  distinct(from, from_Layer, to, to_Layer, .keep_all = TRUE)



sankey_ready <- corr_between_layers %>% 
  filter(from_Layer == "SingleCell" & to_Layer == "eigengene") %>%
  #filter(p< 0.05 & cor < -0.5) %>% 
  rename(SingleCell = from, Eig = to, corr_eig_metab = cor) 



sankey_long <- sankey_ready %>% 
  make_long(SingleCell, Eig) %>% 
  mutate(fill = case_when(x == "SingleCell" ~ node, TRUE ~ ""))



##### Generate the Sankey Diagram #####
my_sankey <- ggplot(sankey_long, aes(x = x, next_x = next_x, 
                                     node = node, next_node = next_node, 
                                     label = node, fill = fill))+
  geom_sankey(flow.alpha = 0.7, node.color = "gray30")+
  geom_sankey_label(size = 3, fill = "white", color = "black")+
  #scale_fill_brewer(palette = "Dark2")+
  theme_sankey(base_size = 10)+
  theme(legend.position = "none")

plot(my_sankey)


# save the plot
ggsave(plot = my_sankey, paste(Sys.Date(), "Sankey.svg", sep = "_" ), width = 8, height = 5, units = "in")





