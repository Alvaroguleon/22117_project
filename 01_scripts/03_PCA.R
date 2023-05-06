library(tidyverse)

data <- read_csv("../00_datasets/processed/data.csv")
metadata <- read_csv("../00_datasets/processed/metadata.csv")

# perform PCA on the data ------
pca <- prcomp(data, scale. = TRUE)

# create a data frame of the first two principal components
pca_data <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Residue = metadata$Residue,
  Mutated_residue = metadata$Mutated_residue,
  Residue_type = metadata$residue_type,
  Mutant = metadata$MUT,
  Exposure = metadata$Exposure,
  Clin = metadata$Clinical_significance,
  Expo = metadata$Exposure,
  Chain = metadata$Chain
)

# WILD TYPE FULL
wt_full <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                color = Residue_type, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
  labs(color = "WT residue type") +
  theme_classic()+
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   

wt_full
ggsave("../02_figures/res_wt_full.png", wt_full, width = 24, height = 18, units = "cm") #+
# ggtitle("A)") + 



# WILD TYPE ZOOM
wt_zoom <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                color = Residue_type, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "WT residue type") +
  theme_classic()+
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))  +
  xlim(-NA, 5) + 
  #theme(legend.position = c(0.15, 0.8))+
  theme(legend.position = c(0.09, 0.85))+
  theme(legend.position = "none")


wt_zoom 
ggsave("../02_figures/res_wt_zoom.png", wt_zoom, width = 24, height = 18, units = "cm")
#+
# ggtitle("A)") 


wt_full
# EXPOSURE  FULL
expo_full <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = Expo, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
  labs(color = "WT residue type") +
  theme_classic()+
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   
ggsave("../02_figures/expo_full.png", expo_full, width = 24, height = 18, units = "cm")


# EXPOSURE ZOOM
expo_zoom <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = Expo, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "WT residue type") +
  theme_classic()+
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   +
  xlim(-NA, 5) + 
  #theme(legend.position = c(0.15, 0.8))+
  theme(legend.position = c(0.09, 0.85))+
  theme(legend.position = "none")


expo_zoom 
ggsave("../02_figures/expo_zoom.png", expo_zoom, width = 24, height = 18, units = "cm")


# MUTANT FULL
custom_colors <- c("#FF6C6A", "#C59B00", "#00C029", "#00C2C5", "#6698FF", "#FF4AE6", "grey")
mut_full <- pca_data %>% filter(Mutated_residue != "Undocumented") %>% 
  ggplot(aes(x = PC1, y = PC2, 
             color = Mutant, label = Mutated_residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "Mutant residue type") +
  scale_color_manual(values = c("#FF6C6A", "#C59B00", "#00C029", "#00C2C5", "#6698FF", "#FF4AE6", "grey"), 
                     labels = c("Acidic", "Aromatic", "Basic", 
                                "Hydrophobic", "Polar", "Special", "Undocumented")) +
  theme_classic()+
  # theme(legend.position = c(0.9, 0.23)) +
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   +
  ylim(-5.5, NA) 

ggsave("../02_figures/res_mut_full.png", mut_full, width = 24, height = 18, units = "cm") #+


wt_zoom
# MUTANT ZOOM
custom_colors <- c("#FF6C6A", "#C59B00", "#00C029", "#00C2C5", "#6698FF", "#FF4AE6", "grey")
mut_zoom <- pca_data %>% filter(Mutated_residue != "Undocumented") %>% 
  ggplot(aes(x = PC1, y = PC2, 
             color = Mutant, label = Mutated_residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "Mutant residue type") +
  scale_color_manual(values = c("#FF6C6A", "#C59B00", "#00C029", "#00C2C5", "#6698FF", "#FF4AE6", "grey"), 
                     labels = c("Acidic", "Aromatic", "Basic", 
                                "Hydrophobic", "Polar", "Special", "Undocumented")) +
  theme_classic()+
  # theme(legend.position = c(0.9, 0.23)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   +
  ylim(-5.5, NA)+
  xlim(-NA, 5) + 
  theme(legend.position = c(0.10, 0.85)) +
  theme(legend.position = "none")
mut_zoom
ggsave("res_mut_zoom.png", mut_zoom, width = 24, height = 18, units = "cm")

mut_zoom
clin_zoom
# Define custom colors for each group
custom_colors_marrones <- c("#0072B2", "#56B4E9", "#E69F00", "#D55E00", "#009E73", "#CC79A7")
custom_colors <- c("#0072B2", "#93B1FF", "#E69F00", "#D55E00", "#009E73", "#BFBFBF")

custom_colors <- c("#0072B2", "#01C620", "#B99F03", "#FF6068", "#009E73", "#BFBFBF")
custom_colors <- c("#0072B2", "#93B1FF", "#B99F03", "#FF6068", "#009E73", "#BFBFBF")

#Exportar 954-635
# Clinical FULL
custom_colors <- c("#0072B2", "#93B1FF", "#B99F03", "#FF6068", "#009E73", "#BFBFBF")
clin_full <- pca_data %>% filter(Mutated_residue != "Undocumented") %>% 
  ggplot(aes(x = PC1, y = PC2, 
             color = Clin, label = Mutated_residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6,show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "Clinical significance") +
  scale_color_manual(values = c("#0072B2", "#93B1FF", "#B99F03", "#FF6068", "#009E73", "#BFBFBF"), 
                     labels = c("Benign", "Likely benign", "Likely pathogenic", 
                                "Pathogenic", "Uncertain significance", "Undocumented")) +
  guides(label = FALSE) +
  theme_classic()+
  theme(legend.position = c(0.9, 0.2)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))   +
  ylim(-5.5, NA) 

ggsave("../02_figures/res_clin_full.png", clin_full, width = 24, height = 18, units = "cm")

# Clinical ZOOM
custom_colors <- c("#0072B2", "#93B1FF", "#B99F03", "#FF6068", "#009E73", "#BFBFBF")
clin_zoom <- pca_data %>% filter(Mutated_residue != "Undocumented") %>% 
  ggplot(aes(x = PC1, y = PC2, 
             color = Clin, label = Mutated_residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6,show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "Clinical significance") +
  scale_color_manual(values = c("#0072B2", "#93B1FF", "#B99F03", "#FF6068", "#009E73", "#BFBFBF"), 
                     labels = c("Benign", "Likely benign", "Likely pathogenic", 
                                "Pathogenic", "Uncertain significance", "Undocumented")) +
  guides(label = FALSE) +
  theme_classic()+
  #theme(legend.position = c(0.85, 0.2))+
  #theme(legend.text = element_text(size = 10), legend.title = element_text(size = 10))  +
  ylim(-5.5, NA)+
  xlim(-NA, 5) + 
  theme(legend.position = "none")
#theme(legend.position = c(0.13, 0.87)

clin_zoom
ggsave("../02_figures/res_clin_zoom.png", clin_zoom, width = 24, height = 18, units = "cm")






