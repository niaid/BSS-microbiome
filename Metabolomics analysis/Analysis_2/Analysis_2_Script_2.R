
#### SOP ###
# dimension reduction through PCA
# visualization of scores and loadings


#### Required Libraries ####
library(tidyverse)
library(RColorBrewer)
library(ggrepel)


#### Z-scale the TSN TCM data ####

dat_scaled <- dat %>% 
  pivot_longer(2:ncol(dat), names_to = "Sample", values_to = "Value_tsn") %>% 
  filter(! str_detect(Sample, "^QC")) %>% 
  group_by(Metabolite) %>% 
  mutate(Value_scaled = as.numeric(scale(Value_tsn, center = TRUE, scale = TRUE))) %>% 
  ungroup()

#### PCA Complete Panel of Metabolites ####

dat_pca <- dat_scaled %>%
  left_join(dat_meta, join_by(Sample)) %>% 
  left_join(dat_responders, join_by(Patient_id), relationship = "many-to-many") %>% 
  filter(Time_point %in% c("Day 0", "Day 2", "Day 8")) %>% 
  relocate(Submitted_Sample_Name,Submitted_Sample_ID, Patient_id,Time_point,MS_method, SubjectID, Response, Colonoscopy,  .after = Sample) %>% 
  relocate(Metabolite, .before = Value_tsn) %>% 
  pivot_wider(id_cols = 1:9, names_from = Metabolite, values_from = Value_scaled)

meta_names <- dat_pca %>% 
  select(c(1:9)) %>% 
  names()

dat_meta <- dat_pca %>% 
  select(all_of(meta_names))

dat_pca <- dat_pca %>% 
  select(-all_of(meta_names))

# perform PCA
pca_results <- prcomp(dat_pca, center = FALSE, scale. = FALSE) 

# extract PCA scores and merge meta data
pca_data <- as.data.frame(pca_results$x) %>% 
  bind_cols(dat_meta) 

###### Scores Plot #######
# calculate variance explained by PC1 and PC2
variance_explained <- summary(pca_results)$importance[2, 1:2] * 100
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

#Design PCA plot with appropriate parameters
pca_plot <- ggplot(pca_data, aes(x = PC1, y = -1*PC2, color = Time_point, fill = Time_point, shape = Response))+
  #Putting in x/y axes (lty = line type (2 = dashed); lwd = line width) #fill = column
  geom_hline(yintercept = 0, lty = 2, lwd = 0.25)+
  geom_vline(xintercept = 0, lty = 2, lwd = 0.25)+
  
  #Adjusting point size (alpha = transparency; shape = 21 refers to specific shape (21-25 contain black border with fill))
  geom_point(size = 4, alpha = 0.8)+
  #Labeling plot/coloring groups
  labs(#title = "PCA Scores Plot", 
       #subtitle = paste("Sample type: Human stool"), 
       x = paste("PC1 (", pc1_var, "% variance )"), 
       y = paste("PC2 (", pc2_var, "% variance )"))+
  
  scale_color_manual(values = c("#e6ab02","#e7298a","#7570b3", "#d95f02", "#666666", "#66A61E"))+
  scale_fill_manual(values = c("#e6ab02","#e7298a","#7570b3", "#d95f02", "#666666", "#66A61E"))+
  scale_shape_manual(values = c(21, 24, 22))+
  
  theme_minimal()+
  theme(#legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5), 
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA), 
    panel.grid = element_blank())  

# view the plot
plot(pca_plot)

# save the plot
ggsave(plot = pca_plot, filename = paste(Sys.Date(), "Hourigansk_0005_PCA_Scores_PC1_PC2_Responders.svg", sep = "_"), width = 5.5, height = 4, units = "in", dpi = 300)




###### Loadings Plot ######
loadings <- as.data.frame(pca_results$rotation[ ,1:4]) 
loadings$Variable <-  rownames(loadings)

loadings <- loadings %>% 
  left_join(dat_metabolite_pathway, join_by(Variable)) 


loadings_plot <- ggplot(loadings, aes(x = PC1, y = PC2, label = Variable))+
  
  geom_hline(yintercept = 0, lty = 2, lwd = 0.25)+
  geom_vline(xintercept = 0, lty = 2, lwd = 0.25)+
  
  geom_point(size = 2)+
  
  ggrepel::geom_text_repel(aes(label = Variable), size = 2.5, max.overlaps = 15)+
  
  labs(x = "Loadings 1", 
       y = "Loadings 2")+
  
  #scale_color_manual(values = c("orange","blue","red","green", "black"))+
  
  theme_minimal()+
  theme(#legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA), 
        panel.grid = element_blank())

plot(loadings_plot)
ggsave(plot = loadings_plot, filename = paste(Sys.Date(), "Hourigansk_0005_PCA_Loadings_Plot_Time_PC1_PC2.svg", sep = "_"), width = 5, height = 5, units = "in", dpi = 300)

