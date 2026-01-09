
#### Required libraries ####
library(tidyverse)
library(ggsankey)


#### SOP ####
# this script generates a sankey diagram for correlation results.
# leverages z-scaled data for each omic layer. Then differences were calculated as Day2 - Day0.


##### Data Prep #####
# Eig -> Metab 
# need to make a key that will be used to match 'from' and 'to' columns back to the layer
feature_layer <- dat_2 %>% 
  select(Feature, Layer) %>% 
  mutate(Layer = case_when(Layer == "Bile acids" ~ "Metabolomics",
                           TRUE ~ Layer))

# tag each feature with its data layer
dat_3_cor_result_merged <- dat_3_cor_result %>% 
  left_join(feature_layer, by = c("from" = "Feature")) %>% 
  rename(from_Layer = Layer) %>% 
  left_join(feature_layer, by = c("to" = "Feature")) %>% 
  rename(to_Layer = Layer)

# remove within-layer correlations
corr_between_layers <- dat_3_cor_result_merged %>% 
  filter(from_Layer != to_Layer) %>% 
  
  # drop duplicates - each unique link should only appear once
  distinct(from, from_Layer, to, to_Layer, .keep_all = TRUE)


# prepare final filtering steps
sankey_ready <- corr_between_layers %>% 
  filter(from_Layer == "eigengene" & to_Layer == "metabolomics") %>%
  filter(p<= 5e-2 & cor >= 0.4) %>% 
  rename(eig = from, metab = to, corr_eig_metab = cor)



##### Formatting for ggsankey #####

# convert to sankey format
sankey_long <- sankey_ready %>% 
  make_long(eig, metab) %>% 
  mutate(fill = case_when(x == "eig" ~ node, TRUE ~ "")) %>% 
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
  scale_fill_brewer(palette = "Dark2")+
  theme_sankey(base_size = 10)+
  theme(legend.position = "none")

# view the plot
plot(my_sankey)

# save the plot
ggsave(plot = my_sankey, paste(Sys.Date(), "Sankey.svg", sep = "_" ), width = 4, height = 9, units = "in")

