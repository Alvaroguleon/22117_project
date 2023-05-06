library(tidyverse)
library(data.table)

### Processing clinvar input ------
hba1 <- read_tsv("../00_datasets/raw/clinvar_result_hba1.txt") %>% select("Gene(s)", "Protein change", "Clinical significance (Last reviewed)")
hba2 <- read_tsv("../00_datasets/raw//clinvar_result_hba2.txt") %>% select("Gene(s)", "Protein change", "Clinical significance (Last reviewed)")
hbb <-  read_tsv("../00_datasets/raw/clinvar_result_hbb.txt") %>% select("Gene(s)", "Protein change", "Clinical significance (Last reviewed)")


hba <- rbind(hba1, hba2)
hba <- hba %>%
  mutate(Residue = str_extract(`Protein change`, "[A-Z]\\d+"))
hbb <- hbb %>%
  mutate(Residue = str_extract(`Protein change`, "[A-Z]\\d+"))


# Subtract 1 from the numbers in Protein_change to match with sequece
hba <- hba %>% rename(Protein_change = "Protein change") %>% 
  mutate(Protein_change = str_replace(Protein_change, "\\d+", function(x) as.numeric(x) - 1))%>% 
  mutate(Residue = str_replace(Residue, "\\d+", function(x) as.numeric(x) - 1))

hbb <- hbb %>% rename(Protein_change = "Protein change") %>% 
  mutate(Protein_change = str_replace(Protein_change, "\\d+", function(x) as.numeric(x) - 1))%>% 
  mutate(Residue = str_replace(Residue, "\\d+", function(x) as.numeric(x) - 1))



## Obtain surface accesibility
# Read text file into R
surface_E <- read.table("../00_datasets/raw/E_Sequence.netsurfp.txt", sep = "") %>% select(V1, V2, V5) %>% mutate(V2 = paste0(V2, seq_along(V2)))
surface_F <- read.table("../00_datasets/raw/F_Sequence.netsurfp.txt", sep = "") %>% select(V1, V2, V5) %>% mutate(V2 = paste0(V2, seq_along(V2)))

colnames(surface_E) <- c("Exposure", "Residue", "RSA")
colnames(surface_F) <- c("Exposure", "Residue", "RSA")


# Join the two datasets based on the "Residue" column
alpha_joined <- left_join(hba, surface_E, by = "Residue")
beta_joined <- left_join(hbb, surface_F, by = "Residue")

alpha_joined <- alpha_joined %>%
  mutate(Protein_change = str_sub(Protein_change, start = 1, end = 1) %>%
           str_c("E", str_sub(Protein_change, start = 2)))


beta_joined <- beta_joined %>%
  mutate(Protein_change = str_sub(Protein_change, start = 1, end = 1) %>%
           str_c("F", str_sub(Protein_change, start = 2)))

# Merging all into one
# positions <- rbind(alpha_joined, beta_joined) %>%  select(!Residue)
positions <- rbind(alpha_joined, beta_joined) %>%
  mutate(
    Residue = str_sub(Residue, 1), # Remove the first position in Residue
    Residue = str_replace(Residue, "^.", str_sub(Protein_change, -1)) # Replace the first position in Residue with the last position of Protein_change
  ) %>% 
  rename(Mutated_residue = Residue) 



positions <- positions %>%
  mutate(Gene = case_when(
    str_detect(`Gene(s)`, "HBB") ~ "HBB",
    str_detect(`Gene(s)`, "HBA1") & !str_detect(`Gene(s)`, "HBA2") ~ "HBA1",
    str_detect(`Gene(s)`, "HBA2") & !str_detect(`Gene(s)`, "HBA1") ~ "HBA2",
    str_detect(`Gene(s)`, "HBA1") & str_detect(`Gene(s)`, "HBA2") ~ "HBA1|HBA2",
    TRUE ~ ""
  )) %>% select(!"Gene(s)")


# Curating the data
positions <- positions %>%
  mutate(`Clinical significance (Last reviewed)` = str_replace(`Clinical significance (Last reviewed)`, "\\(.*\\)", ""))

colnames(positions) <- str_replace_all(colnames(positions), " ", "_")
positions <- positions %>% rename("Clinical_significance" = "Clinical_significance_(Last_reviewed)")

positions <- positions %>%
  separate_rows(Protein_change, sep = ",")

# positions <- positions %>%
#  separate_rows(Gene, sep = " | ")


# Create new columns for the first and last letters of the Protein_change column
positions <- positions %>%
  mutate(WT = str_sub(Protein_change, start = 1, end = 1), # extract first letter and assign to new column "WT"
         MUT = str_sub(Protein_change, start = -1)) # extract last letter and assign to new column "MUT"


# Define the amino acid dictionary with capitalized first letters
amino_acid_dict <- c(V = "Hydrophobic",
                  L = "Hydrophobic",
                  S = "Polar",
                  P = "Special",
                  A = "Hydrophobic",
                  I = "Hydrophobic",
                  M = "Hydrophobic",
                  F = "Hydrophobic",
                  Y = "Aromatic",
                  W = "Aromatic",
                  G = "Special",
                  C = "Special",
                  N = "Polar",
                  Q = "Polar",
                  T = "Polar",
                  R = "Basic",
                  H = "Basic",
                  K = "Basic",
                  D = "Acidic",
                  E = "Acidic")


# Apply the function to the column of 3-letter amino acid codes using the stringr package's str_replace_all function
positions <- positions %>%
  mutate(WT = recode(WT, !!!amino_acid_dict))

positions <- positions %>%
  mutate(MUT = recode(MUT, !!!amino_acid_dict))



# Renaming Clinical_significance cells
conditions <- list(
  str_detect(positions$Clinical_significance, "Conflicting interpretations of pathogenicity") ~ "Uncertain significance",
  positions$Clinical_significance == "not provided" ~ "Other",
  positions$Clinical_significance == "no interpretation for the single variant" ~ "Uncertain significance",
  positions$Clinical_significance == "other" ~ "Other",
  positions$Clinical_significance == "Pathogenic; other" ~ "Pathogenic",
  positions$Clinical_significance == "Benign/Likely benign" ~ "Likely benign",
  positions$Clinical_significance == "Pathogenic/Likely pathogenic" ~ "Likely pathogenic",
  positions$Clinical_significance == "Pathogenic/Likely pathogenic; other" ~ "Likely pathogenic"
)

# Apply the changes using case_when()
positions$Clinical_significance <- case_when(
  !!!conditions,
  TRUE ~ positions$Clinical_significance
)



positions <- select(positions, "Gene", "Protein_change", "WT", "MUT", "Mutated_residue", "Clinical_significance",
                    "Exposure", "RSA")


write_tsv(positions, "../00_datasets/processed/clinvar.txt")




