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
df <- data %>% filter(residue == "AT67") %>% column_to_rownames("residue")
df2 <- data %>% filter(residue == "AN97") %>% column_to_rownames("residue")

df <- rbind(df,df2) %>% rownames_to_column("residue")

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
