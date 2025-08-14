library(tidyverse)

master_table <- read_csv("Supplementary tables/S2 Master_table.csv", guess_max = 2000)

glimpse(master_table)

glycos_data <- master_table %>%
  filter(Glycosylated == "Yes")

length(unique(glycos_data$Gene_symbol))
length(unique(glycos_data$Glyco_unique_ID))

to_correl <- glycos_data %>%
  select(Region, Glyco_unique_ID, Glyco_abundance, NM_protein_abundance) %>%
  filter(Region == "Frontal Cortex")

ggplot(to_correl, aes(x = log2(NM_protein_abundance), y = log2(Glyco_abundance))) + geom_point() + 
  geom_smooth(method = "lm") + xlab("log2(Unmodified protein quant)") + ylab("log2(Glycosylated peptide abundance)")
ggsave("Figures/Glycos_NM_correl_FC.pdf", width = 10, height = 10, units = "cm")


all_glyco_IDs <- glycos_data %>%
  select(Gene_symbol) %>%
  unique()

write_csv(all_glyco_IDs, "Intermediate files/All_glyco_IDs.csv")

##Pie chart for lm
glyco_to_model <- glycos_data %>%
  select(Glyco_unique_ID, Region, Glyco_abundance, NM_protein_abundance)

glyco_lm <- glyco_to_model %>%
  split(.$Glyco_unique_ID) %>%
  map(~lm(log2(Glyco_abundance) ~ log2(NM_protein_abundance), data = .)) %>%
  tibble(Glyco_unique_ID = names(.),
         untidied = .) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-untidied) %>%
  filter(!(term == "(Intercept)")) %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  mutate(Significant = ifelse(padj < 0.05, "Yes", "No"))

lm_outcome <- data.frame(table(glyco_lm$Significant)) 

table(lm_outcome$Freq)

ggplot(lm_outcome, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  scale_fill_brewer(palette="RdYlBu") + theme(panel.background = element_blank())
ggsave("Figures/LM_predicts_glyco_pie.pdf", width = 10, height =10, units = "cm")

## dot plot for all glyco terms - most significant

MF_all_glyco <- read_tsv("Intermediate files/all_glyco_molecular_function.tsv")

MF_all_glyco_filt <- MF_all_glyco %>%
  filter(`false discovery rate` < 0.00001)

ggplot(MF_all_glyco_filt, aes(x = strength, y = reorder(`term description`, strength), size = `observed gene count`)) +
  geom_point() + theme_bw() + xlim(0,1.1) 
ggsave("Figures/All_glyco_GO.pdf", width = 25, height = 20, units = "cm") 


## pivot - small number of peptides have variable deamidation mods at the same site - mean them for now.  Most are very similar.
Glyco_for_PCA <- glycos_data %>%
  select(Glyco_unique_ID, Region, Glyco_abundance) %>%
  pivot_wider(names_from = Glyco_unique_ID, values_from = "Glyco_abundance", values_fn = mean) %>%
  column_to_rownames(var = "Region") %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) 


Glyco_PCA <- prcomp(Glyco_for_PCA, 
                    center = TRUE, 
                    scale. = TRUE)

plot(Glyco_PCA)

plot(Glyco_PCA, type = "l")

summary(Glyco_PCA)

codes <- read_csv("Input files/Region_codes.csv")
Glyco_PC_data <- data.frame(Glyco_PCA$x) %>%
  rownames_to_column("Region") %>%
  left_join(codes)



head(Glyco_PC_data)

Region_col = c("Cingulate Cortex" = "#08519c", "Frontal Cortex" = "#3182bd", "Hippocampus (anterior)" = "#feb24c", 
               "Hippocampus (posterior)" = "#fed976", "Hypothalamus" = "#999999", "Insular Cortex/Claustrum"  = "#6baed6",
               "Motor Cortex" = "#9ecae1", "Occipital Cortex" = "#c6dbef", "Parietal Cortex" = "#54278f",
               "Pia mater (medial)" = "#fd8d3c", "Pia mater (posterior)" = "#f03b20", "Pia mater (anterior)"= "#bd0026", 
               "Primary Visual Cortex" = "#9e9ac8" , "Somatosensory Cortex"  = "#756bb1", "Striatum" = "#252525", 
               "Thalamus (anterior)" = "#006837", "Thalamus (medial)" = "#d9f0a3", "Thalamus (posterior)" = "#78c679")


ggplot(Glyco_PC_data, aes(x = PC1, y = PC2, label = Short_name)) + 
  geom_point(aes(color = Region), size = 3, show.legend = FALSE) + 
  geom_text_repel() +
  scale_color_manual(values = Region_col) + theme_bw()

ggsave("Figures/GlycoPCA.pdf", width = 15, height = 15, units = "cm")


PCA_loadings <- Glyco_PCA$rotation 

PC1_loadings <- data.frame(PCA_loadings) %>%
  rownames_to_column("Unique_ID") %>%
  select(Unique_ID, PC1) %>%
  separate(Unique_ID, into = c("Gene_symbol", NA, NA, NA), remove = FALSE)



##GSEA doesn't show anything as so many proteins have multiple modified peptides with very different rotations on the first PCs - it looks like it just takes the first occurence in the list.
## So let's just take a look at the tails of the distribution and see what those peptides are doing across the samples

PC1_tails_top <- PC1_loadings %>%
  slice_min(order_by = PC1, n = nrow(PC1_loadings)/20) %>%
  mutate(Direction = "Bottom")

PC1_tails <- PC1_loadings %>%
  slice_max(order_by = PC1, n = nrow(PC1_loadings)/20) %>%
  mutate(Direction = "Top") %>%
  bind_rows(PC1_tails_top)

write_csv(PC1_tails, "Intermediate files/Glyco_PC1_tails.csv")

##Made figures for tails in string


SV2A_abundance <- glycos_data %>%
  filter(Gene_symbol == "SV2A") %>%
  select(Glyco_unique_ID, Region, Glyco_abundance) %>%
  pivot_wider(names_from = Glyco_unique_ID, values_from = Glyco_abundance) %>%
  column_to_rownames(var = "Region")

SV2A_protein <- glycos_data %>%
  filter(Gene_symbol == "SV2A") %>%
  select(Gene_symbol, Region, NM_protein_abundance) %>%
  unique() %>%
  pivot_wider(names_from = Gene_symbol, values_from = NM_protein_abundance, values_fn = mean) %>%
  column_to_rownames(var = "Region")

SV2A_names <- read.csv("Input files/SV2A_sites.csv") %>%
  select(Feature, Site)

SV2A_to_tile <- SV2A_abundance %>%
  bind_cols(SV2A_protein) %>%
  mutate(rank = rank(SV2A)) %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(-c(Region, rank), names_to = "Feature", values_to = "Abundance") %>%
  left_join(SV2A_names) %>%
  separate(Site, into = c("Glyco_site", "Non_glyco_deam"), sep = "_") %>%
  group_by(Region,rank, Glyco_site) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  left_join(codes)


ggplot(SV2A_to_tile, aes(x = reorder(Short_name, rank), y = Glyco_site, fill = log2(Abundance))) +
  geom_tile() + scale_fill_gradientn(colours = c("steelblue", "white", "darkred"),
                                   values = scales::rescale(c(-0.5, -0.25, 0, 0.25, 0.5))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("Figures/SV2A_sites.pdf", width = 18, height = 10, units = "cm") 

  
SV2A_to_correl<- SV2A_to_tile %>%
  ungroup() %>%
  select(-c(rank,Region)) %>%
  pivot_wider(names_from = Glyco_site, values_from = Abundance) %>%
  pivot_longer(-c(Short_name, Protein), names_to = "Site", values_to = "Abundance")

site_color <- c("Asn573" = "#7b3294",  "Asn498" = "#a6dba0")

ggplot(SV2A_to_correl, aes(x = log2(Protein), y = log2(Abundance), group = Site)) + geom_point(aes(color = Site)) + 
  geom_smooth(method = "lm", aes(color = "black")) + scale_color_manual(values = site_color) + theme_bw() + 
  xlab("log2(unmodified protein abundance)") + ylab("log2(modified peptide abundance)")
ggsave("Figures/SV2A_correlation.pdf", width = 8, height = 6, units = "cm") 




