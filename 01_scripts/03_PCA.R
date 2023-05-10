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

#write_csv(pca_data, "../00_datasets/processed/pca_data.csv")
# Create a data frame with the variance of each principal component
var_df <- data.frame(
  PC = paste0("PC", 1:ncol(pca$rotation)),
  Variance = pca$sdev^2 / sum(pca$sdev^2)
)

# Sort the data frame by variance in decreasing order
var_df <- var_df[order(-var_df$Variance),]

# Create a bar plot
ggplot(var_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Variance of Principal Components", x = "Principal Component", y = "Variance")

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
  theme(legend.position = c(0.88, 0.18)) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  ggtitle("B)") + theme(plot.title = element_text(size=20))

wt_full
ggsave("../02_figures/res_wt_full.png", wt_full, width = 24, height = 18, units = "cm") 



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
  theme(legend.position = c(0.09, 0.85))+
  theme(legend.position = "none")


wt_zoom 
ggsave("../02_figures/res_wt_zoom.png", wt_zoom, width = 24, height = 18, units = "cm")


wt_full
# EXPOSURE  FULL
expo_full <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = Expo, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
  labs(color = "Exposure class") +
  theme_classic()+
  theme(legend.position = c(0.88, 0.13)) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18))   +
  ggtitle("A)") + theme(plot.title = element_text(size=20))
ggsave("../02_figures/expo_full.png", expo_full, width = 24, height = 18, units = "cm")


# EXPOSURE ZOOM
expo_zoom <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                  color = Expo, label = Residue)) +
  geom_point(size = 3) + 
  geom_text(hjust=0.45, vjust=1.5, size = 6, show.legend = FALSE) + 
  labs(x = "PC1", y = "PC2") +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  labs(color = "Exposure class") +
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
  theme(legend.position = c(0.87, 0.18)) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18))   +
  ylim(-5.5, NA)  +
  ggtitle("C)") + theme(plot.title = element_text(size=20))

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

ggsave("../02_figures/res_mut_zoom.png", mut_zoom, width = 24, height = 18, units = "cm")



# Clinical FULL
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
  theme(legend.position = c(0.83, 0.18)) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18))   +
  ylim(-5.5, NA)  +
  ggtitle("D)") + theme(plot.title = element_text(size=20))
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




