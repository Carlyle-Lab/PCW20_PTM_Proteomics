library(tidyverse)
library(openxlsx)



filt_non_mod_abundance <- read_csv("Supplementary tables/S2 Master_table.csv", guess_max = 2000) %>%
  filter(!is.na(NM_protein_abundance)) %>%
  select(Gene_symbol, Accession, Region, NM_protein_abundance) %>%
  unique()

length(unique(filt_non_mod_abundance$Accession))

## Substitute tiny number of NAs with zeroes and pivot wider for PCA

NM_for_PCA <- filt_non_mod_abundance %>%
  select(-Gene_symbol) %>%
  pivot_wider(names_from = Accession, values_from = NM_protein_abundance, values_fn = sum) %>%
  mutate_all(~ifelse(is.na(.x), 0, .x)) %>%
  column_to_rownames("Region") 


NMQ_PCA <- prcomp(NM_for_PCA, 
                  center = TRUE, 
                  scale. = TRUE)

plot(NMQ_PCA)

plot(NMQ_PCA, type = "l")

summary(NMQ_PCA)

codes <- read_csv("Input files/Region_codes.csv")
PC_data <- data.frame(NMQ_PCA$x) %>%
  rownames_to_column("Region") %>%
  mutate(Region = factor(Region, levels = c("Frontal Cortex", "Cingulate Cortex", "Motor Cortex", "Insular Cortex/Claustrum",
                            "Parietal Cortex", "Somatosensory Cortex", "Primary Visual Cortex", "Occipital Cortex",
                            "`Pia mater (anterior`)", "Pia mater (medial)", "Pia mater (posterior)", "Hippocampus (anterior)",
                            "Hippocampus (posterior)", "Thalamus (anterior)", "Thalamus (medial)", "Thalamus (posterior)",
                            "Hypothalamus", "Striatum"))) %>%
  left_join(codes)
head(PC_data)



Region_col = c("Cingulate Cortex" = "#08519c", "Frontal Cortex" = "#3182bd", "Hippocampus (anterior)" = "#feb24c", 
               "Hippocampus (posterior)" = "#fed976", "Hypothalamus" = "#999999", "Insular Cortex/Claustrum"  = "#6baed6",
               "Motor Cortex" = "#9ecae1", "Occipital Cortex" = "#c6dbef", "Parietal Cortex" = "#54278f",
               "Pia mater (medial)" = "#fd8d3c", "Pia mater (posterior)" = "#f03b20", "Pia mater (anterior)"= "#bd0026", 
               "Primary Visual Cortex" = "#9e9ac8" , "Somatosensory Cortex"  = "#756bb1", "Striatum" = "#252525", 
               "Thalamus (anterior)" = "#006837", "Thalamus (medial)" = "#d9f0a3", "Thalamus (posterior)" = "#78c679")


library(ggrepel)
ggplot(PC_data, aes(x = PC1, y = PC2, label = Short_name)) + 
  geom_point(aes(color = Region), size = 3, show.legend = FALSE) + 
  geom_text_repel() +
  scale_color_manual(values = Region_col) + theme_bw()

ggsave("Figures/NM_PCA_all_regions.pdf", width = 15, height = 15, units = "cm")

##PC plot shows a few clusters and some outliers. Group into Pia (3 regions), Thalamus (3 regions), Cortex (8 regions)
## Hippocampus (2 regions), Striatum, Hypothalamus


PCA_loadings <- NMQ_PCA$rotation 

PC1_loadings <- data.frame(PCA_loadings) %>%
  rownames_to_column("Accession") %>%
  select(Accession, PC1) %>%
  arrange(desc(PC1))

write_csv(PC1_loadings, "Intermediate files/NM_PC1 loadings.csv")

STRING_GSEA_PC1 <- read_tsv("Intermediate files/NM PC1 GSEA BP.tsv")

##dotplot for top terms 

GSEA_PC1_most_sig <- STRING_GSEA_PC1 %>%
  slice_min(order_by = `false discovery rate`, n = 20)

dot_color <- c("top" = "#e9a3c9", "bottom" = "#a1d76a")

ggplot(GSEA_PC1_most_sig, aes(x = `enrichment score`, y = reorder(`term description`, `enrichment score`), size = `genes mapped`, color = direction)) +
  geom_point() + theme_bw() + xlim(0,8) + scale_color_manual(values = dot_color)
ggsave("Figures/NM_PC1_GSEA.pdf", width = 20, height = 10, units = "cm") 



##list tail genes for heatmap
list_tails <- PC1_loadings %>%
  mutate(rank = rank(PC1)) %>%
  filter(rank < 21 | rank > max(rank)-20 )

to_plot_tails <- list_tails %>%
  left_join(filt_non_mod_abundance)

  



PC1_tails_to_heatmap <- to_plot_tails %>%
  select(-c(PC1, rank, Accession)) %>%
  left_join(codes) %>%
  select(-Region) %>%
  pivot_wider(names_from = Short_name, values_from = NM_protein_abundance) %>%
  column_to_rownames("Gene_symbol")


library(gplots)

pdf("Figures/NM_PC1_heatmap.pdf", width=8, height=8)  
heatmap.2(t(as.matrix(log2(PC1_tails_to_heatmap))), trace = "none", 
          scale="col", 
          col=colorRampPalette(c("steelblue","white","darkred"))(50))

dev.off()





##How many transcription factors do we detect?

TF_true <- read_csv("Input files/TF_list.csv")

TFs_in_data <- filt_non_mod_abundance %>%
  filter(Gene_symbol %in% TF_true$Name)

unique(TFs_in_data$Gene_symbol)  
  
n_TFs <- length(unique(TF_true$Name))
n_in_data <- length(unique(TFs_in_data$Gene_symbol))

TF_stats <- tibble("Transcription factors" = c("All", "In brain data"), "Count" = c(n_TFs,n_in_data))

ggplot(TF_stats, aes(x="", y=Count, fill=`Transcription factors`)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + 
  scale_fill_brewer(palette="RdYlBu") + theme(panel.background = element_blank())
ggsave("Figures/TF_pie_chart.pdf", width = 10, height =10, units = "cm")



##  All region PCA - not shown in paper

TFs_to_PCA <- TFs_in_data %>%
  select(-Gene_symbol) %>%
  pivot_wider(names_from = Accession, values_from = NM_protein_abundance, values_fn = sum) %>%
  column_to_rownames("Region") %>%
  mutate_all(~ifelse(is.na(.x), 0, .x))


TF_PCA <- prcomp(TFs_to_PCA, 
                    center = TRUE, 
                    scale. = TRUE)

plot(TF_PCA)

plot(TF_PCA, type = "l")

summary(TF_PCA)

TF_PC_data <- data.frame(TF_PCA$x) %>%
  rownames_to_column("Region")
head(TF_PC_data)


ggplot(TF_PC_data, aes(x = PC1, y = PC2, label = Region)) + geom_point(aes(color = Region, alpha = 0.6), size = 5) + geom_text(hjust=0, vjust=0) +
  scale_color_manual(values = Region_col)

##Looks very similar with Pia and thalamus driving a lot of the changes - go down to just cortex


##Cortex only PCA

cortical_rows <- grep("Cortex", rownames(TFs_to_PCA))
cortical_regions_for_TF_PCA <- TFs_to_PCA[cortical_rows,]

cortical_TF_PCA <- prcomp(cortical_regions_for_TF_PCA, 
                 center = TRUE, 
                 scale. = TRUE)

plot(cortical_TF_PCA)


PCA_summary <- as.data.frame(t(summary(cortical_TF_PCA)$importance)) %>%
  rownames_to_column(var = "PC") %>%
  filter(PC != "PC8")
 
ggplot(PCA_summary, aes(x = PC, y = `Proportion of Variance`)) + geom_bar(stat = "identity") + theme_bw() +
  xlab("Principal Component")
ggsave("Figures/TF_PCA_variance.pdf", width = 15, height = 8, units = "cm") 


cortical_TF_PC_data <- data.frame(cortical_TF_PCA$x) %>%
  rownames_to_column("Region") %>%
  left_join(codes) %>%
  mutate(Region = factor(Region, levels = c("Frontal Cortex", "Cingulate Cortex", "Motor Cortex", "Insular Cortex/Claustrum",
                                            "Parietal Cortex", "Somatosensory Cortex", "Primary Visual Cortex", "Occipital Cortex")))
head(cortical_TF_PC_data)

##PC plot for transcription factors in cortex
ggplot(cortical_TF_PC_data, aes(x = PC1, y = PC2, label = Short_name)) + 
  geom_point(aes(color = Region), size = 3, show.legend = FALSE)  +
  scale_color_manual(values = Region_col) + geom_text_repel() + theme_bw()
ggsave("Figures/TF_NM_PCA_plot.pdf", width = 10, height = 10, units = "cm") 


cortical_TF_PCA_long <- cortical_TF_PC_data %>%
  pivot_longer(-c(Region, Short_name), names_to = "PC", values_to = "loadings") %>%
  mutate(Region = factor(Region, levels = c("Frontal Cortex", "Cingulate Cortex", "Motor Cortex", "Insular Cortex/Claustrum",
                                            "Parietal Cortex", "Somatosensory Cortex", "Primary Visual Cortex", "Occipital Cortex"))) %>%
  filter(PC != "PC8")

ggplot(cortical_TF_PC_data, aes(x = PC3, y = PC4, label = Short_name)) + 
  geom_point(aes(color = Region), size = 3, show.legend = FALSE)  +
  scale_color_manual(values = Region_col) + geom_text_repel() + theme_bw()
ggsave("Figures/TF_NM_PCA_plot_PC3_4.pdf", width = 10, height = 10, units = "cm") 





CTX_TFs_loadings <- data.frame(cortical_TF_PCA$rotation) %>%
  rownames_to_column("Accession") %>%
  pivot_longer(-Accession, names_to = "PC", values_to = "Rotation") %>%
  filter(PC != "PC8")

CTX_TFs_loading_tails_min <- data.frame(CTX_TFs_loadings) %>%
  group_by(PC) %>%
  slice_min(order_by = `Rotation`, n = 10)

tmp <- grep("Cortex", TFs_in_data$Region)
CTX_TF_quants <- TFs_in_data[tmp,]

CTX_TFs_loading_tails <- data.frame(CTX_TFs_loadings) %>%
  group_by(PC) %>%
  slice_max(order_by = `Rotation`, n = 10) %>%
  ungroup() %>%
  bind_rows(CTX_TFs_loading_tails_min) %>%
  filter(PC == "PC1" | PC == "PC2" | PC == "PC3" | PC == "PC4") %>%
  left_join(CTX_TF_quants) %>%
  unique()

CTX_TFs_PC1_to_heatmap <- CTX_TFs_loading_tails %>%
  filter(PC == "PC1") %>%
  select(-c(PC, Rotation)) %>%
  select(-Accession) %>%
  unique() %>%
  pivot_wider(names_from = Gene_symbol, values_from = NM_protein_abundance, values_fn = mean)%>%
  left_join(codes) %>%
  select(-Region) %>%
  column_to_rownames(var = "Short_name")


##PC1 TF heatmap

pdf("Figures/NM_TF_PC1_heatmap.pdf", width=10, height=10)  
heatmap.2(as.matrix(log2(CTX_TFs_PC1_to_heatmap)), trace = "none", 
          scale="col", 
          col=colorRampPalette(c("steelblue","white","darkred"))(50))
dev.off()

CTX_TFs_PC2_to_heatmap <- CTX_TFs_loading_tails %>%
  filter(PC == "PC2") %>%
  select(-c(PC, Rotation)) %>%
  select(-Accession) %>%
  unique() %>%
  pivot_wider(names_from = Gene_symbol, values_from = NM_protein_abundance)%>%
  left_join(codes) %>%
  select(-Region) %>%
  column_to_rownames(var = "Short_name")



##PC2 TF heatmap
pdf("Figures/NM_TF_PC2_heatmap.pdf", width=10, height=10)  

heatmap.2(as.matrix(log2(CTX_TFs_PC2_to_heatmap)), trace = "none", 
          scale="col", 
          col=colorRampPalette(c("steelblue","white","darkred"))(50))
dev.off()


  



##Line graph rather than heatmap?

CTX_TFs_PC3_4_to_line <- CTX_TFs_loading_tails %>%
  filter(PC == "PC3" | PC == "PC4") %>%
  mutate(Region = factor(Region, levels = c("Frontal Cortex", "Cingulate Cortex", "Motor Cortex", "Insular Cortex/Claustrum",
                                            "Parietal Cortex", "Somatosensory Cortex", "Primary Visual Cortex", "Occipital Cortex"))) %>%
  mutate(Direction = ifelse(Rotation > 0, "Top", "Bottom"))  %>%
  left_join(codes)


summary_line <- CTX_TFs_PC3_4_to_line %>%
  group_by(PC, Direction, Region) %>%
  summarise(mean_val = mean(NM_protein_abundance))  %>%
  left_join(codes)



ggplot(CTX_TFs_PC3_4_to_line, aes(x = Short_name, y = log2(NM_protein_abundance))) + 
  geom_line(aes(group = Gene_symbol, color = Gene_symbol)) +
  facet_grid(rows = vars(Direction), cols = vars(PC), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_smooth(data = summary_line, aes(x=Short_name, y=log2(mean_val), group = PC),color = "black",  method = "lm") +
  geom_text(data = subset(CTX_TFs_PC3_4_to_line, Short_name == "V1C"), aes(label = Gene_symbol, colour = Gene_symbol,  hjust = 1)) +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.3))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_smooth(method = "lm") + ylab("log2 Quant") + 
  xlab("Region")


ggsave("Figures/PC_3_4_line_plots.pdf", width = 20, height = 20, units = "cm")


