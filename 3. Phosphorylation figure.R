library(tidyverse)

master_table <- read_csv("Supplementary tables/S2 Master_table.csv", guess_max = 2000)

phos_data <- master_table %>%
  filter(Phos == "Yes")

to_correl <- phos_data %>%
  select(Region, Phos_unique_ID, Phos_SIA_abundance, NM_protein_abundance) %>%
  filter(Region == "Frontal Cortex")

ggplot(to_correl, aes(x = log2(NM_protein_abundance), y = log2(Phos_SIA_abundance))) + geom_point() + 
  geom_smooth(method = "lm") + xlab("log2(Unmodified protein quant)") + ylab("log2(Phospho peptide abundance")
ggsave("Figures/Phos_NM_correl_FC.pdf", width = 10, height = 10, units = "cm")

phos_gene <- phos_data %>%
  select(Peptide_ID, Phos_unique_ID) %>%
  unique()

length(unique(phos_data$Phos_unique_ID))

number_summary <- data.frame(table(phos_gene$Peptide_ID)) %>%
  filter(Freq>1)



##Linear model for each phospho site - predicted abundance - then look for ones with biggest residuals to plot

##Phos_LM
phos_to_model <- phos_data %>%
  select(Phos_unique_ID, Region, Phos_SIA_abundance, NM_protein_abundance)

phos_lm <- phos_to_model %>%
  split(.$Phos_unique_ID) %>%
  map(~lm(log2(Phos_SIA_abundance) ~ log2(NM_protein_abundance), data = .)) %>%
  tibble(Phos_unique_ID = names(.),
         untidied = .) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-untidied) %>%
  filter(!(term == "(Intercept)")) %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  mutate(Significant = ifelse(padj < 0.05, "Yes", "No"))

lm_outcome <- data.frame(table(phos_lm$Significant)) 

ggplot(lm_outcome, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  scale_fill_brewer(palette="RdYlBu") + theme(panel.background = element_blank())
ggsave("Figures/LM_predicts_phos_pie.pdf", width = 10, height =10, units = "cm")

###SIA_LM
SIA_data <- master_table %>%
  filter(SIA == "Yes")

SIA_to_correl <- SIA_data %>%
  select(Region, Phos_unique_ID, Phos_SIA_abundance, NM_protein_abundance) %>%
  filter(Region == "Frontal Cortex")

ggplot(SIA_to_correl, aes(x = log2(NM_protein_abundance), y = log2(Phos_SIA_abundance))) + geom_point() + 
  geom_smooth(method = "lm") + xlab("log2(Unmodified protein quant)") + ylab("log2(SIA peptide abundance")
ggsave("Figures/SIA_NM_correl_FC.pdf", width = 10, height = 10, units = "cm")


SIA_to_model <- SIA_data %>%
  select(Phos_unique_ID, Region, Phos_SIA_abundance, NM_protein_abundance)

SIA_lm <- SIA_to_model %>%
  split(.$Phos_unique_ID) %>%
  map(~lm(log2(Phos_SIA_abundance) ~ NM_protein_abundance, data = .)) %>%
  tibble(Phos_unique_ID = names(.),
         untidied = .) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-untidied) %>%
  filter(!(term == "(Intercept)")) %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  mutate(Significant = ifelse(padj < 0.05, "Yes", "No"))

lm_outcome <- data.frame(table(SIA_lm$Significant)) 

ggplot(lm_outcome, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  scale_fill_brewer(palette="RdYlBu") + theme(panel.background = element_blank())
ggsave("Figures/LM_predicts_SIA_pie.pdf", width = 10, height =10, units = "cm")



  

##Plan for phospho page
## Outcome of lm - how many proteins have phos assoc with protein abundance, how many not
## Take the proteins where abundance doesn't predict phospho - PCA 
## If PCA looks interesting, then use STRING to look at tails
## If PCA looks boring, then just do STRING on all these proteins, select a category to heatmap and show individual sites

##PCA for non-assoc

phos_no_assoc <- phos_lm %>%
  filter(Significant == "No")

phos_to_PCA <- phos_to_model %>%
  filter(Phos_unique_ID %in% phos_no_assoc$Phos_unique_ID) %>%
  select(Phos_unique_ID, Region, Phos_SIA_abundance) %>%
  pivot_wider(names_from = Phos_unique_ID, values_from = Phos_SIA_abundance, values_fn = mean) %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) %>%
  column_to_rownames(var = "Region")


Phos_PCA <- prcomp(phos_to_PCA, 
                 center = TRUE, 
                 scale. = TRUE)

plot(Phos_PCA)

plot(Phos_PCA, type = "l")

summary(Phos_PCA)

codes <- read_csv("Input files/Region_codes.csv")

Phos_PC_data <- data.frame(Phos_PCA$x) %>%
  rownames_to_column("Region") %>%
  left_join(codes)
head(Phos_PC_data)

Region_col = c("Cingulate Cortex" = "#08519c", "Frontal Cortex" = "#3182bd", "Hippocampus (anterior)" = "#feb24c", 
               "Hippocampus (posterior)" = "#fed976", "Hypothalamus" = "#999999", "Insular Cortex/Claustrum"  = "#6baed6",
               "Motor Cortex" = "#9ecae1", "Occipital Cortex" = "#c6dbef", "Parietal Cortex" = "#54278f",
               "Pia mater (medial)" = "#fd8d3c", "Pia mater (posterior)" = "#f03b20", "Pia mater (anterior)"= "#bd0026", 
               "Primary Visual Cortex" = "#9e9ac8" , "Somatosensory Cortex"  = "#756bb1", "Striatum" = "#252525", 
               "Thalamus (anterior)" = "#006837", "Thalamus (medial)" = "#d9f0a3", "Thalamus (posterior)" = "#78c679")




ggplot(Phos_PC_data, aes(x = PC1, y = PC2, label = Short_name)) + geom_point(aes(color = Region, alpha = 0.6), size = 4) + 
  geom_text_repel(max.overlaps = 18) +
  scale_color_manual(values = Region_col) + theme_bw()
ggsave("Figures/Phos_PC1_PC2.pdf", width = 15, height =10, units = "cm")




long_PCA <- Phos_PC_data %>%
  pivot_longer(-c(Region, Short_name), names_to = "PC", values_to = "loading") 

ggplot(long_PCA, aes(x = Short_name, y = loading)) + geom_boxplot(aes(color = Region)) + 
  facet_wrap(~PC) + scale_color_manual(values = Region_col) 


## Looks very similar to total protein - PCs 3 & 4 look more interesting
## Do the tails then maybe select one or two proteins from the tail ends
 
PC_proteins <- data.frame(Phos_PCA$rotation) %>%
  rownames_to_column(var = "Phos_Unique_ID") %>%
  pivot_longer(-Phos_Unique_ID, names_to = "PC", values_to = "rotation")
head(PC_proteins)

PCs1_2_proteins_tails_max <- PC_proteins %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  separate(Phos_Unique_ID, into = c("Gene_symbol", NA, NA, NA), remove = FALSE) %>%
  group_by(PC) %>%
  slice_max(order_by = rotation, n = n_distinct(PC_proteins$Phos_Unique_ID)/20) %>%
  mutate(end = "max")

PCs1_2_proteins_tails_both <- PC_proteins %>%
  filter(PC == "PC1" | PC == "PC2") %>%
  separate(Phos_Unique_ID, into = c("Gene_symbol", NA, NA, NA), remove = FALSE) %>%
  group_by(PC) %>%
  slice_min(order_by = rotation, n = n_distinct(PC_proteins$Phos_Unique_ID)/20) %>%
  mutate(end = "min") %>%
  bind_rows(PCs1_2_proteins_tails_max)

write_csv(PCs1_2_proteins_tails_both, "Intermediate files/Phos_no_assoc_PC1_2_tails.csv") 

##  Plots for one direction for each PC - nothing exciting in other two halves
PC1_min_all <- read_tsv("Intermediate files/Phos_PC1_min_all.tsv") %>%
  mutate(direction = "min")
PC2_min_all <- read_tsv("Intermediate files/Phos_PC2_min_all.tsv") %>%
  mutate(direction = "min")
PC1_max_all <- read_tsv("Intermediate files/Phos_PC1_max_all.tsv") %>%
  mutate(direction = "max")
PC2_max_all <- read_tsv("Intermediate files/Phos_PC2_max_all.tsv") %>%
  mutate(direction = "max")


PC1_dot_color <- c("min" = "#e9a3c9", "max" = "#a1d76a")

PC1_max_all %>%
  bind_rows(PC1_min_all) %>%
  filter(`#category` == "UniProt Keywords") %>%
  ggplot(aes(x = `strength`, y = reorder(`term description`, `strength`), size = `observed gene count`, color = direction)) +
  geom_point()  + theme_bw() +scale_color_manual(values = PC1_dot_color)
ggsave("Figures/Phos_PC1_GO_Uniprot.pdf", width = 15, height = 10, units = "cm") 





PC2_dot_color <- c("max" = "#fc8d59", "min" = "#91bfdb")

PC2_min_all %>%
  bind_rows(PC2_max_all) %>%
  filter( `#category` == "UniProt Keywords" ) %>%
  ggplot(aes(x = `strength`, y = reorder(`term description`, `strength`), size = `observed gene count`, color = direction)) +
  geom_point() + theme_bw()  +scale_color_manual(values = PC2_dot_color)
ggsave("Figures/Phos_PC2_GO_Uniprot.pdf", width = 15, height = 10, units = "cm") 


##You're here - need to put these plots into the figure

to_plot <- c("KALRN_NEAT_7_1x_T4", "RALGAPA2_SSSP_10_1x_S3", "TSC2_SQSG_21_1x_S3")

phos_to_model %>%
  filter(Phos_unique_ID %in% to_plot) %>%
  filter(Region == "Thalamus (posterior)" | Region == "Hypothalamus" | Region == "Striatum") %>%
  ggplot(aes(x = Region, y = log2(Phos_SIA_abundance))) + geom_point()

## Map sites in study to known sites - interesting function?  
  

Rho_GTAPases_reg <- PC2_min_all %>%
  filter(`term description` == "Positive regulation of GTPase activity") %>%
  select(`matching proteins in your network (labels)`) %>%
  separate_rows(`matching proteins in your network (labels)`, sep = ",")


Rho_GTPases_reg_to_plot <- phos_to_model %>%
  separate(Phos_unique_ID, into = c("Gene_symbol", NA, NA, NA), sep = "_", remove = FALSE) %>%
  filter(Gene_symbol %in% Rho_GTAPases_reg$`matching proteins in your network (labels)`) %>%
  filter(Phos_unique_ID %in% PCs1_2_proteins_tails_both$Phos_Unique_ID) 

TSC2 <- phos_to_model %>%
  separate(Phos_unique_ID, into = c("Gene_symbol", NA, NA, NA), sep = "_", remove = FALSE) %>%
  filter(Gene_symbol == "TSC2")

write_csv(data.frame(unique(TSC2$Phos_unique_ID)), "Intermediate files/TSC2 peptides.csv")

TSC2_p_summary <- read_csv("Intermediate files/TSC2 peptides_phos_sites.csv") %>%
  rename(Phos_unique_ID = unique.TSC2.Phos_unique_ID.)

##some peptides are the same phospho site - sum these

TSC2_summed <- TSC2 %>%
  left_join(TSC2_p_summary) %>%
  select(-c(Phos_unique_ID,  Gene_symbol, `Phos predicted?`, `Known function`))%>%
  group_by(Phospho_site, Region, NM_protein_abundance) %>%
  summarise(summed_quant = sum(Phos_SIA_abundance)) %>%
  left_join(codes) %>%
  ungroup() %>%
  select(-Region) %>%
  pivot_wider(names_from = Phospho_site, values_from = summed_quant) %>%
  column_to_rownames(var = "Short_name")


library(gplots)

pdf("Figures/Phos_TSC2_heatmap.pdf", width=12, height=10)  

heatmap.2(as.matrix(log2(TSC2_summed)), trace = "none", 
          scale="col", 
          col=colorRampPalette(c("steelblue","white","darkred"))(100))

dev.off()

