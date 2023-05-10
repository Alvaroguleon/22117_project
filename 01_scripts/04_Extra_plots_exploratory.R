library(tidyv)
# Loading the difference of free energies -----
alpha <- read_csv("../00_datasets/raw/energies_alpha.csv")
beta <- read_csv("../00_datasets/raw/energies_beta.csv")
names(alpha)[2] <- 'Chain'
names(beta)[2] <- 'Chain'
beta$Chain[beta$Chain == FALSE] <- "F"


# Uniting the columns into one
alpha <- unite(alpha, residue, "WT residue type":"Residue #", sep = "")
# Remove second letter and add "A" to each element
alpha$residue <- paste0("A", substr(alpha$residue, 1, 1), substr(alpha$residue, 3, nchar(alpha$residue)))

beta <- unite(beta, residue, "WT residue type":"Residue #", sep = "")
# Remove second letter and add "B" to each element
beta$residue <- paste0("B", substr(beta$residue, 1, 1), substr(beta$residue, 3, nchar(beta$residue)))


# data <- data %>% separate(residue, into = c("WT", "MUT"), sep = "_")
alpha <- column_to_rownames(alpha, var = "residue")
beta <- column_to_rownames(beta, var = "residue")

# Merge datasets
data <- rbind(alpha, beta)

data <- data %>% rownames_to_column("residue")

#---- Hydrophobic

df <- data %>% filter(residue == "AF128") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BW15") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG of hydrophobic residues", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()




#---- arriba y abajo

df <- data %>% filter(residue == "BS89") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BW15") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG of arriba residues", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()



# ----- Left plot
#df <- data %>% filter(residue == "AT67") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AN97") %>% column_to_rownames("residue")
df3 <- data %>% filter(residue == "BN102") %>% column_to_rownames("residue")

df <- rbind(df2,df3) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG on the left", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()


#---- In the middle middle

df <- data %>% filter(residue == "BV54") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AA69") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG in the middle middle", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

#---- In the middle 2

df <- data %>% filter(residue == "AH87") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BH92") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG in the middle 2", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

#---- In the middle down

df <- data %>% filter(residue == "AF128") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BW15") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG in the middle down", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

# ----- Glycines
df <- data %>% filter(residue == "BG64") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AG59") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG on the right", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()



#---- In the axis of coordinates

df <- data %>% filter(residue == "AA19") %>% column_to_rownames("residue")


df <- rbind(df) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG in the axis", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()


#---- In the axis of coordinates symmetric

df <- data %>% filter(residue == "BT50") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AR92") %>% column_to_rownames("residue")
df3 <- data %>% filter(residue == "AP44") %>% column_to_rownames("residue")
df4 <- data %>% filter(residue == "BM55") %>% column_to_rownames("residue")


df <- rbind(df,df2,df3,df4) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "B)", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

# ---- Many
df <- data %>% filter(residue == "AN97") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AH87") %>% column_to_rownames("residue")
df3 <- data %>% filter(residue == "AP44") %>% column_to_rownames("residue")
df4 <- data %>% filter(residue == "BV54") %>% column_to_rownames("residue")
df5 <- data %>% filter(residue == "AG59") %>% column_to_rownames("residue")
df6 <- data %>% filter(residue == "BA140") %>% column_to_rownames("residue")


df <- rbind(df,df3,df5, df6) %>% rownames_to_column("Residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "Residue")

# Define the desired order of residues
residue_order <- c("AG59", "AP44", "BA140" ,"AN97")  # Replace these with the actual residue names

# Define corresponding colors and shapes
residue_colors <- c("AG59" = "#1F77B4", "AP44" = "black", "BA140" = "red", "AN97" = "#2CA02C") 
residue_shapes <- c("AG59" = 19, "AP44" = 17, "BA140" = 23, "AN97" = 15) 

plot <- ggplot(data = df_long, aes(x = variable, y = value, shape = Residue, color = Residue, fill = Residue)) +
  geom_point(size = 3) +
  geom_smooth(aes(group = Residue), method = "lm", se = FALSE) +
  labs(x = "Mutation", y = "ΔΔG (kcal/mol)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_classic() +
  theme(legend.position = c(0.1, 0.88)) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  scale_color_manual(values = residue_colors, breaks = residue_order) +
  scale_shape_manual(values = residue_shapes, breaks = residue_order) + 
  scale_fill_manual(values = residue_colors, breaks = residue_order) +  # add this line to fill the shapes with the desired color
  theme(axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

plot
ggsave("../02_figures/residues.png", plot, width = 24, height = 18, units = "cm") 


 




# -----
df <- data %>% filter(residue == "BY130") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BW15") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
plot <- ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(x = "Mutation", y = "ΔΔG (kcal/mol)", color = "Mutation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()
plot
ggsave("../02_figures/points.png", plot, width = 20, height = 14, units = "cm") 

# -----
df <- data %>% filter(residue == "BG74") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AS124") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
plot <- ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(x = "Mutation", y = "ΔΔG (kcal/mol)", color = "Mutation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()
plot
ggsave("../02_figures/points2.png", plot, width = 20, height = 14, units = "cm") 



# -----
df <- data %>% filter(residue == "BG74") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AF128") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
plot <- ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(x = "Mutation", y = "ΔΔG (kcal/mol)", color = "Mutation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()
plot
ggsave("../02_figures/points3.png", plot, width = 20, height = 14, units = "cm") 

