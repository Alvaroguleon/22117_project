# Processing of energies
library(tidyverse)


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



#### Metadata -----
metadata <- data.frame(
  id = rownames(data)
) 

metadata <- separate(metadata, id, into = c("chain", "residue"), sep = 1)
metadata <- separate(metadata, residue, into = c("residue", "position"), sep = 1)

# Define a dictionary of 3-letter to 1-letter codes
aa_dict <- c(A = "Hydrophobic", R = "Positive", N = "Polar", D = "Negative", C = "Sulfur",
             E = "Negative", Q = "Polar", G = "Hydrophobic", H = "Aromatic", I = "Hydrophobic",
             L = "Hydrophobic", K = "Positive", M = "Hydrophobic", "F" = "Aromatic", P = "Hydrophobic",
             S = "Polar", "T" = "Polar", "W" = "Aromatic", Y = "Aromatic", V = "Hydrophobic")

# Create a vector of residue types
residue_type <- ifelse (metadata$residue == "V", "Hydrophobic", 
                        ifelse (metadata$residue == "L", "Hydrophobic",
                                ifelse (metadata$residue == "S", "Polar",
                                        ifelse (metadata$residue == "P", "Special",
                                                ifelse (metadata$residue == "A", "Hydrophobic",
                                                        ifelse (metadata$residue == "I", "Hydrophobic",
                                                                ifelse (metadata$residue == "M", "Hydrophobic",
                                                                        ifelse (metadata$residue == "F", "Hydrophobic",
                                                                                ifelse (metadata$residue == "Y", "Aromatic",
                                                                                        ifelse (metadata$residue == "W", "Aromatic",
                                                                                                ifelse (metadata$residue == "G", "Special",
                                                                                                        ifelse (metadata$residue == "C", "Special",
                                                                                                                ifelse (metadata$residue == "N", "Polar",
                                                                                                                        ifelse (metadata$residue == "Q", "Polar",
                                                                                                                                ifelse (metadata$residue == "T", "Polar",
                                                                                                                                        ifelse (metadata$residue == "R", "Basic",
                                                                                                                                                ifelse (metadata$residue == "H", "Basic",
                                                                                                                                                        ifelse (metadata$residue == "K", "Basic",
                                                                                                                                                                ifelse (metadata$residue == "D", "Acidic",
                                                                                                                                                                        ifelse (metadata$residue == "E", "Acidic", NA))))))))))))))))))))

# Add the vector as a new column to the data frame
metadata <- cbind (metadata, residue_type)
metadata <- unite(metadata, Residue, "chain":"position", sep = "")



## Obtain surface accesibility
# Read text file into R
surface_E <- read.table("../00_datasets/raw/E_Sequence.netsurfp.txt", sep = "")%>% 
  select(V1, V2, V5) %>% 
  mutate(V2 = paste0(V2, seq_along(V2))) %>%
  mutate(V2 = str_sub(V2, start = 1, end = 1) %>%
           str_c("E", str_sub(V2, start = 2)))
surface_F <- read.table("../00_datasets/raw/F_Sequence.netsurfp.txt", sep = "") %>% 
  select(V1, V2, V5) %>% 
  mutate(V2 = paste0(V2, seq_along(V2))) %>%
  mutate(V2 = str_sub(V2, start = 1, end = 1) %>%
           str_c("F", str_sub(V2, start = 2)))

colnames(surface_E) <- c("Exposure", "Residue", "RSA") 


colnames(surface_F) <- c("Exposure", "Residue", "RSA")


surface <- rbind(surface_E, surface_F)

# Modify the Residue column based on conditions
surface$Residue <- ifelse(substr(surface$Residue, 2, 2) == "E",
                          paste0("A", substr(surface$Residue, 1, 1), substr(surface$Residue, 3, nchar(surface$Residue))),
                          ifelse(substr(surface$Residue, 2, 2) == "F",
                                 paste0("B", substr(surface$Residue, 1, 1), substr(surface$Residue, 3, nchar(surface$Residue))),
                                 surface$Residue))


# Obtain data from clinvar
metadata2 <- read_tsv("../00_datasets/processed/clinvar.txt") %>% filter(Clinical_significance != "Other") 

metadata2 <- metadata2 %>% select("Protein_change", "Clinical_significance", "MUT", "Mutated_residue") %>%
  mutate(Protein_change = str_sub(Protein_change, end = -2)) %>% 
  rename("Residue" = Protein_change)

# Define the order of priority for the clinical significance values
priority_order <- c("Pathogenic", "Likely pathogenic", "Benign", "Likely benign", "Uncertain significance")

# Remove duplicates based on Residue column and priority order of clinical significance
metadata2 <- metadata2 %>%
  arrange(Residue, match(Clinical_significance, priority_order)) %>%
  distinct(Residue, .keep_all = TRUE)

# Insert "F" in second position of cells with length < 3
metadata2$Residue <- ifelse(nchar(metadata2$Residue) < 3, 
                            paste0(substr(metadata2$Residue, 1, 1), "F", substr(metadata2$Residue, 2, nchar(metadata2$Residue))), 
                            metadata2$Residue)

metadata2$Residue[1] <- "AE143"


# Modify the Residue column based on conditions
metadata2$Residue <- ifelse(substr(metadata2$Residue, 2, 2) == "E",
                            paste0("A", substr(metadata2$Residue, 1, 1), substr(metadata2$Residue, 3, nchar(metadata2$Residue))),
                            ifelse(substr(metadata2$Residue, 2, 2) == "F",
                                   paste0("B", substr(metadata2$Residue, 1, 1), substr(metadata2$Residue, 3, nchar(metadata2$Residue))),
                                   metadata2$Residue))




# Remove first row because it is not in data
# metadata2 <- metadata2[!grepl("AE143", metadata2$Residue), ]

metadata <- left_join(metadata, metadata2, by=c("Residue")) %>%
  replace_na(list(Clinical_significance = "Undocumented"))%>%
  replace_na(list(MUT = "Undocumented")) %>% 
  replace_na(list(Mutated_residue = "Undocumented")) 

metadata <- left_join(metadata, surface, by = "Residue")

metadata <- metadata %>% 
  mutate(Exposure = str_replace_all(Exposure, c("B" = "Buried", "E" = "Exposed")))

metadata <- metadata[-149, ]



metadata <- metadata %>%
  mutate(Chain = ifelse(str_starts(Residue, "A"), "Alpha", "Beta")) %>% 
  mutate(
    Mutated_residue = if_else(
      Mutated_residue != "Undocumented",
      str_c(str_sub(Residue, 1, 1), Mutated_residue),
      Mutated_residue
    )
  )

write_csv(data, "../00_datasets/processed/data.csv")
write_csv(metadata, "../00_datasets/processed/metadata.csv")
