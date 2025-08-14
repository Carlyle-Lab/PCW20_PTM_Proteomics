library(tidyverse)

master_table <- read_csv("Supplementary tables/S2 Master_table.csv", guess_max = 2000)
glimpse(master_table)

acet_data <- master_table %>%
  filter(Acetylated == "Yes")

all_acet_IDs <- acet_data %>%
  select(Gene_symbol) %>%
  unique()

write_csv(all_acet_IDs, "Intermediate files/All_acet_IDs.csv")

length(unique(acet_data$Gene_symbol))
length(unique(acet_data$Acet_unique_ID))

to_correl <- acet_data %>%
  select(Region, Acet_unique_ID, Acet_abundance, NM_protein_abundance) %>%
  filter(Region == "Frontal Cortex")

ggplot(to_correl, aes(x = log2(NM_protein_abundance), y = log2(Acet_abundance))) + geom_point() + 
  geom_smooth(method = "lm") + xlab("log2(Unmodified protein quant)") + ylab("log2(Acetylated peptide abundance)")
ggsave("Figures/Acet_NM_correl_FC.pdf", width = 10, height = 10, units = "cm")


##Acet predicted by total protein

acet_to_model <- acet_data %>%
  select(Acet_unique_ID, Region, Acet_abundance, NM_protein_abundance)

acet_lm <- acet_to_model %>%
  split(.$Acet_unique_ID) %>%
  map(~lm(log2(Acet_abundance) ~ log2(NM_protein_abundance), data = .)) %>%
  tibble(Acet_unique_ID = names(.),
         untidied = .) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-untidied) %>%
  filter(!(term == "(Intercept)")) %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  mutate(Significant = ifelse(padj < 0.05, "Yes", "No"))

lm_outcome <- data.frame(table(acet_lm$Significant)) 

ggplot(lm_outcome, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  scale_fill_brewer(palette="RdYlBu") + theme(panel.background = element_blank())
ggsave("Figures/LM_predicts_acet_pie.pdf", width = 10, height =10, units = "cm")


## dot plot for all acet terms - most significant

BP_all_acet <- read_tsv("Intermediate files/Acet_GO_BP.tsv")

BP_all_acet_filt <- BP_all_acet %>%
  filter(`false discovery rate` < 0.00001)

ggplot(BP_all_acet_filt, aes(x = strength, y = reorder(`term description`, strength), size = `observed gene count`)) +
  geom_point() + theme_bw() +xlim(0,0.9)
ggsave("Figures/All_acet_GO.pdf", width = 25, height = 20, units = "cm") 


## pivot - small number of peptides have variable deamidation mods at the same site - mean them for now.  Most are very similar.
Acet_for_PCA <- acet_data %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  pivot_wider(names_from = "Acet_unique_ID", values_from = "Acet_abundance", values_fn = mean) %>%
  column_to_rownames(var = "Region") %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) 



Acet_PCA <- prcomp(Acet_for_PCA, 
                    center = TRUE, 
                    scale. = TRUE)

plot(Acet_PCA)

plot(Acet_PCA, type = "l")

summary(Acet_PCA)

codes <- read_csv("Input files/Region_codes.csv")


Acet_PCA_data <- data.frame(Acet_PCA$x) %>%
  rownames_to_column("Region") %>%
  left_join(codes)
 
head(Acet_PCA_data)

Region_col = c("Cingulate Cortex" = "#08519c", "Frontal Cortex" = "#3182bd", "Hippocampus (anterior)" = "#feb24c", 
               "Hippocampus (posterior)" = "#fed976", "Hypothalamus" = "#999999", "Insular Cortex/Claustrum"  = "#6baed6",
               "Motor Cortex" = "#9ecae1", "Occipital Cortex" = "#c6dbef", "Parietal Cortex" = "#54278f",
               "Pia mater (medial)" = "#fd8d3c", "Pia mater (posterior)" = "#f03b20", "Pia mater (anterior)"= "#bd0026", 
               "Primary Visual Cortex" = "#9e9ac8" , "Somatosensory Cortex"  = "#756bb1", "Striatum" = "#252525", 
               "Thalamus (anterior)" = "#006837", "Thalamus (medial)" = "#d9f0a3", "Thalamus (posterior)" = "#78c679")
library(ggrepel)

ggplot(Acet_PCA_data, aes(x = PC1, y = PC2, label = Short_name)) + 
  geom_point(aes(color = Region), size = 3, show.legend = FALSE) + 
  geom_text_repel() +
  scale_color_manual(values = Region_col) + theme_bw()

ggsave("Figures/Acet_PCA.pdf", width = 15, height = 15, units = "cm")



PCA_loadings <- Acet_PCA$rotation 

PC1_loadings <- data.frame(PCA_loadings) %>%
  rownames_to_column("Unique_ID") %>%
  select(Unique_ID, PC1) %>%
  separate(Unique_ID, into = c("Gene_symbol", NA, NA, NA), remove = FALSE)


PC1_tails_top <- PC1_loadings %>%
  slice_min(order_by = PC1, n = nrow(PC1_loadings)/20) %>%
  mutate(Direction = "Bottom")

PC1_tails <- PC1_loadings %>%
  slice_max(order_by = PC1, n = nrow(PC1_loadings)/20) %>%
  mutate(Direction = "Top") %>%
  bind_rows(PC1_tails_top)

write_csv(PC1_tails, "Intermediate files/Acet_PC1_tails.csv")


##Made figures for tails in string

HIST <- grep("HIST", PC1_tails$Gene_symbol)
histones <- PC1_tails[HIST,]

table(histones$Gene_symbol)
##Histone H1B is the most commonly modified

H1B <-  histones %>%
  filter(Gene_symbol == "HIST1H1B")


Histone_protein_quants <- acet_data %>%
  filter(Gene_symbol == "HIST1H1B")

SETA <- Histone_protein_quants %>%
  filter(Acet_unique_ID == "HIST1H1B_SETA_20_1x_K16" | Acet_unique_ID == "HIST1H1B_SETA_20_1x_N-Term") %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  group_by(Region) %>%
  summarise(Acet_abundance = sum(Acet_abundance)) %>%
  mutate(Acet_unique_ID = "HIST1H1B_SETA_20_summed") %>%
  select(Acet_unique_ID, Region, Acet_abundance)

KPAA <- Histone_protein_quants %>%
  filter(Acet_unique_ID == "HIST1H1B_KPAA_9_1x_K8" | Acet_unique_ID == "HIST1H1B_PAAA_8_1x_K7") %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  group_by(Region) %>%
  summarise(Acet_abundance = sum(Acet_abundance)) %>%
  mutate(Acet_unique_ID = "HIST1H1B_KPAA_9_summed") %>%
  select(Acet_unique_ID, Region, Acet_abundance)

KATG <- Histone_protein_quants %>%
  filter(Acet_unique_ID == "HIST1H1B_KATG_19_1x_K13" | Acet_unique_ID == "HIST1H1B_ATGP_18_1x_K12") %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  group_by(Region) %>%
  summarise(Acet_abundance = sum(Acet_abundance)) %>%
  mutate(Acet_unique_ID = "HIST1H1B_KATG_19_summed") %>%
  select(Acet_unique_ID, Region, Acet_abundance)

KATK <- Histone_protein_quants %>%
  filter(Acet_unique_ID == "HIST1H1B_KATK_8_1x_K4" | Acet_unique_ID == "HIST1H1B_ATKS_7_1x_K3") %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  group_by(Region) %>%
  summarise(Acet_abundance = sum(Acet_abundance)) %>%
  mutate(Acet_unique_ID = "HIST1H1B_KATK_8_summed") %>%
  select(Acet_unique_ID, Region, Acet_abundance)

HIST1H1B_protein <- Histone_protein_quants %>%
  filter(Gene_symbol == "HIST1H1B") %>%
  select(Gene_symbol, Region, NM_protein_abundance) %>%
  unique() %>%
  mutate(Acet_unique_ID = "Non-modified protein") %>%
  rename(Acet_abundance = NM_protein_abundance) %>%
  select(Acet_unique_ID, Region, Acet_abundance) 

Histone_protein_quants_to_plot <- Histone_protein_quants %>%
  select(Acet_unique_ID, Region, Acet_abundance) %>%
  filter(Acet_unique_ID == "HIST1H1B_AAGA_9_1x_K8" | Acet_unique_ID == "HIST1H1B_NGLS_10_1x_K9" | Acet_unique_ID == "HIST1H1B_ALAA_15_1x_K11") %>%
  bind_rows(SETA) %>%
  bind_rows(KPAA) %>%
  bind_rows(KATG) %>%
  bind_rows(KATK) %>%
  bind_rows(HIST1H1B_protein) %>%
  mutate(Acet_unique_ID = gsub("HIST1H1B_", "", Acet_unique_ID)) %>%
  pivot_wider(names_from = Acet_unique_ID, values_from = Acet_abundance) %>%
  left_join(codes) %>%
  select(-Region) %>%
  column_to_rownames(var = "Short_name")




library(gplots)

pdf("Figures/Acet_HIST1H1B_heatmap.pdf", width=12, height=10)  

heatmap.2(as.matrix(log2(Histone_protein_quants_to_plot)), trace = "none", 
          scale="col", 
          col=colorRampPalette(c("steelblue","white","darkred"))(50))
dev.off()




Histone_protein_quants_correl_plot <- Histone_protein_quants_to_plot %>%
  rownames_to_column(var = "Region") %>%
  pivot_longer(-c(Region, `Non-modified protein`), names_to = "Acet_unique_ID", values_to = "Quant")
head(Histone_protein_quants_correl_plot)


ggplot(Histone_protein_quants_correl_plot, aes(x = log2(`Non-modified protein`), y = log2(Quant), group = Acet_unique_ID, color = Acet_unique_ID)) +
  geom_point() + geom_smooth(method = lm) + theme_bw() + scale_color_brewer(palette = "RdYlBu") + xlab("log2(Non-modified protein)") +
  ylab("log2(Modified peptide abundance)")
ggsave("Figures/histone_acet_correl.pdf", width = 15, height =15, units = "cm")


