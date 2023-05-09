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
  labs(title = "DDG in the middle middle", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_classic()

# ---- Many
df <- data %>% filter(residue == "BS89") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "BW15") %>% column_to_rownames("residue")
df3 <- data %>% filter(residue == "BG64") %>% column_to_rownames("residue")
df4 <- data %>% filter(residue == "AN97") %>% column_to_rownames("residue")


df <- rbind(df,df2,df3,df4) %>% rownames_to_column("residue")

# Convert the data frame to a long format
df_long <- reshape2::melt(df, id.vars = "residue")

# Create a scatter plot
ggplot(data = df_long, aes(x = variable, y = value, color = residue)) +
  geom_point(size = 3) +
  labs(title = "DDG of different clusters", x = "Residue", y = "Value", color = "Residue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_classic()
