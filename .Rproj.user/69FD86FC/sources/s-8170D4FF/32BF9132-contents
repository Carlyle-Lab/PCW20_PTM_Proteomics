library(tidyverse)
library(openxlsx)

## Import non-modified protein data, v. basic QC, convert to long table

Non_mod_two_peps_data <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1a  non mod")
colnames(Non_mod_two_peps_data)

head(Non_mod_two_peps_data)

ID_conversion <- read_tsv("Input files/Human_reference_biomart.txt") %>%
  select(`UniProtKB/Swiss-Prot ID`, `Gene name`) %>%
  rename(Accession = `UniProtKB/Swiss-Prot ID`) %>%
  rename(Gene_symbol = `Gene name`)

contaminants <- read_csv("Input files/Contaminants_Accession.csv")

# Make a nice long dataframe using Normalized abundances for each protein, and convert Uniprot ID to Gene symbol.  Remove common contaminants. Convert zeroes to NAs
Non_mod_abundance <- Non_mod_two_peps_data %>%
  select(Accession, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-Accession, names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  select(Accession,Gene_symbol, Region, Norm_quant) %>%
  unique() %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Region = gsub("Perietal Cortex", "Parietal Cortex", Region)) %>%
  mutate(Region = gsub("Pia matter \\(anterior\\)", "Pia mater (anterior)", Region))



Protein_NAs <- Non_mod_abundance %>%
  group_by(Accession) %>%
  summarise(n_NA = sum(is.na(Norm_quant))) %>%
  filter(n_NA < 19)

All_missing <- Protein_NAs %>%
  filter(n_NA == 18) 

##Some proteins are entirely NAs - filter out

filt_non_mod_abundance <- Non_mod_abundance %>%
  filter(!(Accession %in% All_missing$Accession))
table(Protein_NAs$n_NA)

Protein_NAs <- filt_non_mod_abundance %>%
  group_by(Accession) %>%
  summarise(n_NA = sum(is.na(Norm_quant))) %>%
  filter(n_NA < 19)


NM_protein_NAs <- ggplot(Protein_NAs, aes(x = n_NA)) + geom_histogram(binwidth = 1) + xlab("Number of missing values per protein") + xlim(-1,18)
ggsave("Figures/NM_protein_NAs.pdf", width = 10, height = 10, units = "cm") 
## Full rows of NAs are present in the original spreadsheet - not an artefact of ID conversion
## Filter out proteins with empty rows

NM_all_protein_abundanc <- ggplot(filt_non_mod_abundance, aes(x = Region, y = log2(Norm_quant))) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("log2(Non-modified protein abundance)")


ggsave("Figures/NM_protein_abundance_boxplots.pdf", width = 10, height = 10, units = "cm")









## Import phospho and SIA peptide data, v. basic QC, convert to long table

phos_only_data <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1b Phospho")
SIA_only_data <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1c SIA")
both_data <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1d Phospho_SIA")


# Make a nice long dataframe using Normalized abundances for each protein, and convert Uniprot ID to Gene symbol.  Remove common contaminants. Convert zeroes to NAs
phos_only_abundance <- phos_only_data %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>%
  mutate(phos_number = str_extract(Modifications, pattern = ".{2}(?=Phospho)")) %>%
  mutate(phos_site = str_extract(Modifications, pattern = "Phospho \\[([^\\]]+)\\]")) %>%
  mutate(phos_site = str_remove(phos_site, "Phospho \\[")) %>%
  mutate(phos_site = str_remove(phos_site, "\\]")) %>%
  mutate(phos_site = str_replace_all(phos_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length, phos_number, phos_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Phos = "Yes") %>%
  mutate(SIA = "No")

SIA_only_abundance <- SIA_only_data %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>%
  mutate(SIA_number = str_extract(Modifications, pattern = ".{2}(?=SIA)")) %>%
  mutate(SIA_site = str_extract(Modifications, pattern = "SIA \\[([^\\]]+)\\]")) %>%
  mutate(SIA_site = str_remove(SIA_site, "SIA \\[")) %>%
  mutate(SIA_site = str_remove(SIA_site, "\\]")) %>%
  mutate(SIA_site = str_replace_all(SIA_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length,SIA_number, SIA_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Phos = "No") %>%
  mutate(SIA = "Yes")

Phos_and_SIA_abundance <- both_data %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>%
  mutate(phos_number = str_extract(Modifications, pattern = ".{2}(?=Phospho)")) %>%
  mutate(phos_site = str_extract(Modifications, pattern = "Phospho \\[([^\\]]+)\\]")) %>%
  mutate(phos_site = str_remove(phos_site, "Phospho \\[")) %>%
  mutate(phos_site = str_remove(phos_site, "\\]")) %>%
  mutate(phos_site = str_replace_all(phos_site, "; ", "_")) %>%
  mutate(SIA_number = str_extract(Modifications, pattern = ".{2}(?=SIA)")) %>%
  mutate(SIA_site = str_extract(Modifications, pattern = "SIA \\[([^\\]]+)\\]")) %>%
  mutate(SIA_site = str_remove(SIA_site, "SIA \\[")) %>%
  mutate(SIA_site = str_remove(SIA_site, "\\]")) %>%
  mutate(SIA_site = str_replace_all(SIA_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length,phos_number, SIA_number, phos_site, SIA_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Phos = "Yes") %>%
  mutate(SIA = "Yes")

Phos_SIA_all_sites <- phos_only_abundance %>%
  bind_rows(SIA_only_abundance) %>%
  bind_rows(Phos_and_SIA_abundance)

Protein_NAs <- Phos_SIA_all_sites %>%
  group_by(Gene_symbol) %>%
  summarise(n_NA = sum(is.na(Norm_quant))) %>%
  filter(n_NA < 19)

Phos_peptide_NAs <- ggplot(Protein_NAs, aes(x = n_NA)) + geom_histogram(binwidth = 1) + xlab("Number of missing values per peptide")
ggsave("Figures/Phospho_SIA_peptide_NAs.pdf", width = 10, height = 10, units = "cm")



proteins_to_keep <- Protein_NAs %>%
  filter(n_NA < 15)

filt_all_phos_SIA <- Phos_SIA_all_sites %>%
  filter(Gene_symbol %in% proteins_to_keep$Gene_symbol)

ggplot(filt_all_phos_SIA, aes(x = Region, y = log2(Norm_quant))) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("log2(Phos & SIA peptide abundance)")
ggsave("Figures/Phospho_SIA_peptide_Abundances.pdf", width = 10, height = 10, units = "cm")





## Import glycosylated peptide data, v. basic QC, convert to long table

Glyco_data_NXC <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1e NXC Glycos")
Glyco_data_NXST <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1f NXST Glycos")


# Make a nice long dataframe using Normalized abundances for each protein, and convert Uniprot ID to Gene symbol.  Remove common contaminants. Convert zeroes to NAs
Glyco_abundance_NXC <- Glyco_data_NXC %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>% 
  mutate(glyco_site = str_extract(Modifications, pattern = "\\[([^\\]]+)\\]")) %>%
  mutate(glyco_site = str_remove(glyco_site, "\\[")) %>%
  mutate(glyco_site = str_remove(glyco_site, "\\]")) %>%
  mutate(glyco_site = str_replace(glyco_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length, glyco_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Consensus_seq = "NXC")

Glyco_abundance_NXST <- Glyco_data_NXST %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>%
  mutate(glyco_site = str_extract(Modifications, pattern = "\\[([^\\]]+)\\]")) %>%
  mutate(glyco_site = str_remove(glyco_site, "\\[")) %>%
  mutate(glyco_site = str_remove(glyco_site, "\\]")) %>%
  mutate(glyco_site = str_replace(glyco_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length, glyco_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region)) %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Consensus_seq = "NXST")

Glyco_all_sites <- Glyco_abundance_NXC %>%
  bind_rows(Glyco_abundance_NXST)

Protein_NAs <- Glyco_all_sites %>%
  group_by(Gene_symbol) %>%
  summarise(n_NA = sum(is.na(Norm_quant))) %>%
  filter(n_NA < 19)

Glyco_Protein_NAs <- ggplot(Protein_NAs, aes(x = n_NA)) + geom_histogram(binwidth = 1) + 
  xlab("Number of missing values per peptide")
ggsave("Figures/Glyco_peptide_NAs.pdf", width = 10, height = 10, units = "cm")



proteins_to_keep <- Protein_NAs %>%
  filter(n_NA < 15)

filt_all_glycos <- Glyco_all_sites %>%
  filter(Gene_symbol %in% proteins_to_keep$Gene_symbol)



##Stability of acquisition sample by sample

glycosylated_peptides_abundance <- ggplot(filt_all_glycos, aes(x = Region, y = log2(Norm_quant))) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("log2(Glycosylated peptide abundance)")
ggsave("Figures/Glyco_peptide_abundances.pdf", width = 10, height = 10, units = "cm")




## Import acetylated peptide data, v. basic QC, convert to long table

acet_data <- read.xlsx("Supplementary tables/Table S1.xlsx", sheet = "Table S1g Acetylation")


# Make a nice long dataframe using Normalized abundances for each protein, and convert Uniprot ID to Gene symbol.  Remove common contaminants. Convert zeroes to NAs
acet_abundance <- acet_data %>%
  rename(Accession = "Master.Protein.Accessions") %>%
  filter(!(Accession %in% contaminants$Accession)) %>%
  separate(Accession, into = c("Accession", NA), sep = "-") %>%
  left_join(ID_conversion) %>%
  unique() %>%
  separate(Annotated.Sequence, into = c(NA, "Peptide", NA), sep = "\\.", remove = FALSE) %>%
  mutate(short_name = str_sub(Peptide, start = 1L, end = 4L)) %>%
  mutate(length = str_length(Peptide)) %>%
  mutate(acet_number = str_extract(Modifications, pattern = ".{2}(?=Acetyl)")) %>%
  mutate(acet_site = str_extract(Modifications, pattern = "Acetyl \\[([^\\]]+)\\]")) %>%
  mutate(acet_site = str_remove(acet_site, "Acetyl \\[")) %>%
  mutate(acet_site = str_remove(acet_site, "\\]")) %>%
  mutate(acet_site = str_replace_all(acet_site, "; ", "_")) %>%
  unite("Unique_ID", c(Gene_symbol,short_name, length, acet_number, acet_site), sep = "_", remove = FALSE) %>%
  select(Gene_symbol, Unique_ID, `Abundances.(Normalized):.F1:.127N,.Sample,.Cingulate.Cortex`: `Abundances.(Normalized):.F1:.132N,.Sample,.Thalamus.(posterior)`) %>%
  pivot_longer(-c(Gene_symbol, Unique_ID), names_to = "Region", values_to = "Norm_quant") %>%
  separate(Region, into = c(NA, NA, "Region"), sep = ",\\.") %>%
  mutate(Norm_quant = ifelse(Norm_quant == 0, NA, Norm_quant)) %>%
  mutate(Region = gsub("Perietal\\.Cortex", "Parietal.Cortex", Region)) %>%
  mutate(Region = gsub("Pia\\.matter\\.\\(anterior\\)", "Pia.mater.(anterior)", Region))  %>%
  mutate(Region = gsub("\\.", " ", Region)) %>%
  mutate(Acetylated = "Yes")


Protein_NAs <- acet_abundance %>%
  group_by(Gene_symbol) %>%
  summarise(n_NA = sum(is.na(Norm_quant))) %>%
  filter(n_NA < 19)

ggplot(Protein_NAs, aes(x = n_NA)) + geom_histogram(binwidth = 1) + xlab("Number of missing values per peptide")
ggsave("Figures/Acet_peptide_NAs.pdf", width = 10, height = 10, units = "cm")



proteins_to_keep <- Protein_NAs %>%
  filter(n_NA < 15)

filt_acet <- acet_abundance %>%
  filter(Gene_symbol %in% proteins_to_keep$Gene_symbol)



##Stability of acquisition sample by sample

ggplot(filt_acet, aes(x = Region, y = log2(Norm_quant))) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("log2(Acetylated peptide abundance")
ggsave("Figures/Acet_peptide_Abundances.pdf", width = 10, height = 10, units = "cm")









##Combine into one master long table

long_format <- filt_non_mod_abundance %>%
  rename(NM_protein_abundance = Norm_quant) 


phos_data <- filt_all_phos_SIA %>%
  rename(Phos_SIA_abundance = Norm_quant) %>%
  separate(Unique_ID, into = c(NA, "Peptide", "length", NA), sep = "_", remove = FALSE) %>%
  unite("Peptide_ID", Peptide:length, sep = "_") %>%
  rename(Phos_unique_ID = Unique_ID) %>%
  select(Gene_symbol, Peptide_ID, Phos_unique_ID, Phos, SIA, Region, Phos_SIA_abundance) %>%
  mutate(Phos_unique_ID = str_remove_all(Phos_unique_ID, pattern = "\\(([^\\]]+)\\)"))



acet_data <- filt_acet %>%
  rename(Acet_abundance = Norm_quant) %>%
  separate(Unique_ID, into = c(NA, "Peptide", "length", NA), sep = "_", remove = FALSE) %>%
  unite("Peptide_ID", Peptide:length, sep = "_") %>%
  rename(Acet_unique_ID = Unique_ID) %>%
  select(Gene_symbol, Peptide_ID, Acet_unique_ID, Acetylated, Region, Acet_abundance) %>%
  mutate(Acet_unique_ID = str_remove_all(Acet_unique_ID, pattern = "\\(([^\\]]+)\\)"))

glyco_data <- filt_all_glycos %>%
  rename(Glyco_abundance = Norm_quant) %>%
  separate(Unique_ID, into = c(NA, "Peptide", "length", NA), sep = "_", remove = FALSE) %>%
  unite("Peptide_ID", Peptide:length, sep = "_") %>%
  rename(Glyco_unique_ID = Unique_ID) %>%
  mutate(Glycosylated = "Yes") %>%
  select(Gene_symbol, Peptide_ID, Glyco_unique_ID, Glycosylated, Region, Glyco_abundance) %>%
  mutate(Glyco_unique_ID = str_remove_all(Glyco_unique_ID, pattern = "\\(([^\\]]+)\\)"))


## so that merge works correctly, first create a master list of all peptides and proteins in the study

nm_prots <- filt_non_mod_abundance %>%
  select(Gene_symbol) %>%
  unique()

phos_peptides <- phos_data %>%
  select(Gene_symbol, Peptide_ID) %>%
  unique()

acet_peptides <- acet_data %>%
  select(Gene_symbol, Peptide_ID) %>%
  unique()

glyco_peptides <- glyco_data %>%
  select(Gene_symbol, Peptide_ID) %>%
  unique()

all_peps <- phos_peptides %>%
  bind_rows(acet_peptides) %>%
  bind_rows(glyco_peptides) %>%
  unique()

##Unique call decreases rows from 21624 to 21589, so there are some peptides with more than one type of mod
##Create master list of proteins/peptides to merge

tmp <- which(!(nm_prots$Gene_symbol %in% all_peps$Gene_symbol))
prot_only <- nm_prots[tmp,] %>%
  mutate(Peptide_ID = NA)

master_list <- prot_only %>%
  bind_rows(all_peps)



master_table <- master_list %>%
  left_join(long_format) %>%
  left_join(phos_data, by = c("Gene_symbol", "Peptide_ID", "Region")) %>%
  mutate(Phos_SIA_prop = Phos_SIA_abundance/NM_protein_abundance) %>%
  left_join(acet_data) %>%
  mutate(Acet_prop = Acet_abundance/NM_protein_abundance) %>%
  left_join(glyco_data) %>%
  mutate(Glyco_prop = Glyco_abundance/NM_protein_abundance) %>%
  arrange(Gene_symbol) %>%
  mutate()





## Add RNA counts to the master table
## RNA matrix obtained as download from Brainspan https://www.brainspan.org/static/download.html

RNA_matrix <- read_csv("Input files/expression_matrix.csv", col_names = FALSE) %>%
  column_to_rownames(var = "X1")

col_data <- read_csv("Input files/columns_metadata.csv") %>%
  unite(Unique_ID, c(donor_id, structure_acronym), sep = "_")

matrix_colnames <- c(col_data$Unique_ID)
colnames(RNA_matrix) <- matrix_colnames

row_data <- read_csv("Input files/rows_metadata.csv") %>%
  unite(Unique_ID, c(gene_symbol, ensembl_gene_id), sep = "_")


matrix_rownames <- c(row_data$Unique_ID)
rownames(RNA_matrix) <- matrix_rownames



RNA_filt <- RNA_matrix %>%
  rownames_to_column(var = "Gene_ID") %>%
  separate(Gene_ID, into = c("Gene_symbol", "Ensembl_ID"), sep = "_") %>%
  filter(Gene_symbol %in% master_table$Gene_symbol)

freq_table <- data.frame(table(RNA_filt$Gene_symbol)) %>%
  filter(Freq == 2)

samples_to_keep <- c("19 pcw", "21 pcw")

RNA_filt_long <- RNA_filt %>%
  pivot_longer(-c(Gene_symbol, Ensembl_ID), names_to = "Unique_ID", values_to = "RPKM") %>%
  left_join(col_data) %>%
  filter(age %in% samples_to_keep)

matched_regions <- read_csv("Input files/Brainspan_matched_regions.csv")

RNA_filt_long <- RNA_filt_long %>%
  filter(structure_name %in% matched_regions$Brainspan_region)


to_correl <- RNA_filt_long %>%
  unite(Unique_gene_ID, c(Gene_symbol, Ensembl_ID), remove = FALSE, sep = "_") %>%
  select(Gene_symbol, Unique_gene_ID, RPKM, age) %>%
  pivot_wider(names_from = age, values_from = RPKM, values_fn = mean)

ggplot(to_correl, aes(x = log10(`19 pcw`), y = log10(`21 pcw`))) + geom_point()

cor(to_correl$`19 pcw`, to_correl$`21 pcw`, method = "pearson")
##0.91 pearson correlation - probably just pick one sample to use - pick 19 pcw as fewer zeroes?

RNA_19PCW <- RNA_filt_long %>%
  filter(age == "19 pcw") %>%
  select(Gene_symbol, RPKM, structure_name) %>%
  rename(Brainspan_region = structure_name) %>%
  left_join(matched_regions)

master_protein_RNA_table <- master_table %>%
  left_join(RNA_19PCW)


##Write master csv file


write_csv(master_protein_RNA_table, "Supplementary tables/S2 Master_table.csv")

##background set for GO

background_set <- master_protein_RNA_table %>%
  select(Gene_symbol) %>%
  unique()

write_csv(background_set, "Intermediate files/background_set.csv")


##Summary Upset plot

protein_present <- master_protein_RNA_table %>%
  select(Gene_symbol, NM_protein_abundance) %>%
  group_by(Gene_symbol) %>%
  summarise(NM_protein = mean(NM_protein_abundance)) 



for_Upset <- master_protein_RNA_table %>%
  select(Gene_symbol, Phos, SIA, Acetylated, Glycosylated) %>%
  left_join(protein_present)%>%
  unique() %>%
  mutate(Phos = ifelse(Phos == "No", NA, Phos)) %>%
  mutate(SIA  = ifelse(SIA == "No", NA, SIA)) %>%
  filter(!is.na(Gene_symbol))

list_gen <- function(col){
  for_Upset %>%
    filter(!is.na(col)) %>%
    pull(Gene_symbol)
}



protein_list <- apply(for_Upset[,-1], 2, list_gen)
library(UpSetR)

pdf(file="Figures/protein_level_upset.pdf", width = 8.5, height = 6)
upset(fromList(protein_list), order.by = "freq", nsets = 10, nintersects = 20)
dev.off() 

length(unique(protein_list$NM_protein))






