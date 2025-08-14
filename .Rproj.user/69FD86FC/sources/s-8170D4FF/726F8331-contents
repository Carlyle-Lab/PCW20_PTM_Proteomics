library(tidyverse)
library(openxlsx)

master_table <- read_csv("Supplementary tables/S2 Master_table.csv", guess_max = 2000)
glimpse(master_table)


RNA_protein_to_plot <- master_table %>%
  select(Gene_symbol, Region, NM_protein_abundance, RPKM) %>%
  unique()

ggplot(RNA_protein_to_plot, aes(x = log2(RPKM), y = log2(NM_protein_abundance))) + geom_point() + facet_wrap(~Region) 


zero_RNA <- RNA_protein_to_plot %>%
  filter(RPKM == 0)


zero_summary <- data.frame(table(zero_RNA$Gene_symbol, zero_RNA$Region)) %>%
  arrange(Var2)

##Create background list for enrichment throughout - all unique genes IDed in the study

protein_list<- data.frame(unique(master_table$Gene_symbol))

write_csv(protein_list, "Input files/background_list.csv")



##Original brainspan looked at fold changes between regions between the two metrics.  Let's try thalamus and prefrontal to start

MD_thal_Frontal_Protein <- RNA_protein_to_plot %>%
  select(-RPKM) %>%
  filter(Region == "Frontal Cortex" | Region == "Thalamus (medial)") %>%
  pivot_wider(names_from = Region, values_from =NM_protein_abundance, values_fn = mean) %>%
  mutate(log2FC_Protein = log2(`Frontal Cortex`) - log2(`Thalamus (medial)`)) %>%
  select(Gene_symbol, log2FC_Protein)



MD_thal_Frontal_RNA <- RNA_protein_to_plot %>%
  select(-NM_protein_abundance) %>%
  filter(Region == "Frontal Cortex" | Region == "Thalamus (medial)") %>%
  pivot_wider(names_from = Region, values_from =RPKM, values_fn = mean) %>%
  mutate(log2FC_RNA = log2(`Frontal Cortex`) - log2(`Thalamus (medial)`)) %>%
  select(Gene_symbol, log2FC_RNA) %>%
  left_join(MD_thal_Frontal_Protein) %>%
  mutate(RNA_diff = ifelse(log2FC_RNA < -1, "Down", ifelse(
    log2FC_RNA > 1, "Up", "No"
  ))) %>%
  mutate(Protein_diff = ifelse(log2FC_Protein < -1, "Down", ifelse(
    log2FC_Protein > 1, "Up", "No"
  ))) %>%
  mutate(RNA_protein_diff = log2FC_RNA - log2FC_Protein) %>%
  mutate(diff_colour = ifelse(RNA_protein_diff < -1, "Big", ifelse(
    RNA_protein_diff > 1, "Big", "Small"
  ))) %>%
  mutate(Final_colour = ifelse(RNA_diff == "No" & Protein_diff == "No", "grey", ifelse(
    RNA_diff == Protein_diff & diff_colour == "Big", "orange", ifelse(
      Protein_diff == "No" & (RNA_diff == "Up" | RNA_diff == "Down"), "blue", ifelse(
        RNA_diff == "No" & (Protein_diff == "Up" | Protein_diff == "Down"), "ltblue", ifelse(
          diff_colour == "Small" & (RNA_diff == "Up" & Protein_diff == "Up")|(RNA_diff == "Down" & Protein_diff == "Down"), "yellow", ifelse(
            diff_colour == "Big" & RNA_diff != Protein_diff, "red", NA
          )
        )
      )
    )
  )))


plot_colors <- c("grey" = "#FFFFFF", "yellow" = "#ffffbf", "blue" = "#2c7bb6", "ltblue" = "#abd9e9", "red" = "#d7191c", "orange" = "#fdae61")

ggplot(MD_thal_Frontal_RNA, aes(x = log2FC_Protein, y = log2FC_RNA)) + geom_point(aes(color=Final_colour), show.legend = FALSE) + 
  geom_abline(intercept = 1, slope = 1) +
  geom_abline(intercept = -1, slope = 1) + geom_vline(xintercept = -1) + geom_vline(xintercept = 1) +
  geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + xlim(-10,10) + ylim(-10,10) +
  scale_color_manual(values = plot_colors) +
  xlab("log2(fold change protein)") + ylab("log2(fold change RNA)") 
ggsave("Figures/Thal_Frontal_RNA_protein.pdf", width = 10, height = 10, units = "cm")

Color_freq <- data.frame(table(MD_thal_Frontal_RNA$Final_colour))

ggplot(Color_freq, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity", show.legend = FALSE) + 
  xlab("RNA protein outcome") + 
  scale_fill_manual(values = plot_colors) + xlab("RNA vs protein outcome") + ylab("Gene count")
ggsave("Figures/Thal_frontal_RNA_freq.pdf", width = 10, height = 10, units = "cm")




write_csv(MD_thal_Frontal_RNA, "Intermediate files/Thal_vs_frontal_for_GO.csv")

##GO analysis in STRING with defaults for enrichment and all proteins detected as background set


##Red category proteins - String output file for Uniprot keyword - only 4 sig terms

string_red_thal_RNA <- read_tsv("Intermediate files/Thal_FC_red.tsv") %>%
  mutate(colour = "red")

point_colour <- c("red" = "#d7191c")
ggplot(string_red_thal_RNA, aes(x = strength, y = reorder(`term description`, strength), size = `observed gene count`)) +
  geom_point(aes(color = colour)) + theme_bw() + xlim(0,0.5) + scale_color_manual(values = point_colour)
ggsave("Figures/Thal_FC_GO_red.pdf", width = 20, height = 10, units = "cm") 



##String output file for Molecular Function for blue, FC down

string_blue_thal_RNA <- read_tsv("Intermediate files/Thal_Blue_down_MF.tsv") %>%
  mutate(colour = "blue")
head(string_blue_thal_RNA)

point_colour <- c("blue" = "#2c7bb6")
ggplot(string_blue_thal_RNA, aes(x = strength, y = reorder(`term description`, strength), size = `observed gene count`)) +
  geom_point(aes(color = colour)) + xlim(0,1.2) + theme_bw()  + scale_color_manual(values = point_colour)
ggsave("Figures/Thal_FC_GO_blue.pdf", width = 20, height = 10, units = "cm") 




##Shiva's data 

RGC_pairs <- read_csv("Input files/thalamus_RGC_PFC_GW22.csv") %>%
  select(target, ligand.complex, receptor.complex,aggregate_rank)

Subplate_pairs <- read_csv("Input files/thalamus_subplate_PFC_GW22.csv") %>%
  select(target,ligand.complex, receptor.complex,aggregate_rank)


both_pairs <- RGC_pairs %>%
  full_join(Subplate_pairs) %>%
  filter(aggregate_rank < 0.05) %>%
  select(-c(aggregate_rank, target)) %>%
  pivot_longer(c(ligand.complex,receptor.complex), names_to = "Type", values_to = "Gene_symbol") %>%
  unique() %>%
  left_join(MD_thal_Frontal_RNA) 

library(ggrepel)
plot_colors <- c("ligand.complex" = "#af8dc3", "receptor.complex" = "#7fbf7b")
ggplot(both_pairs, aes(x = log2FC_Protein, y = log2FC_RNA)) + geom_point(aes(color=Type)) + geom_abline(intercept = 1, slope = 1) +
    geom_abline(intercept = -1, slope = 1) + geom_vline(xintercept = -1) + geom_vline(xintercept = 1) +
    geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + xlim(-5,5) + ylim(-5,5) +
    scale_color_manual(values = plot_colors) + 
    geom_text_repel(aes(label = Gene_symbol), max.overlaps = 5) 
ggsave("Figures/Receptor_ligand_all.pdf", width = 15, height = 10, units = "cm") 

## Looking at this plot and which ligands/receptors fall in the "right" places, focus in on a number

focus_ligand <- c( "L1CAM", "CNTN4",  "LGALS1")
focus_receptor <- c( "CD9","PTPRG",  "ITGB1" )
focus_pairs <- RGC_pairs %>%
  full_join(Subplate_pairs) %>%
  filter(aggregate_rank < 0.05) %>%
  filter(ligand.complex %in% focus_ligand) %>%
  filter(receptor.complex %in% focus_receptor) %>%
  mutate(pair_num = rank(receptor.complex))%>%
  select(-c(aggregate_rank, target)) %>%
  pivot_longer(c(ligand.complex,receptor.complex), names_to = "Type", values_to = "Gene_symbol") %>%
  unique() %>%
  left_join(MD_thal_Frontal_RNA) 

ggplot(focus_pairs, aes(x = log2FC_Protein, y = log2FC_RNA)) + geom_point(aes(color=as.factor(pair_num), shape = Type), size = 3) + geom_abline(intercept = 1, slope = 1) +
  geom_abline(intercept = -1, slope = 1) + geom_vline(xintercept = -1) + geom_vline(xintercept = 1) +
  geom_hline(yintercept = -1) + geom_hline(yintercept = 1) + xlim(-4,4) + ylim(-4,4) +
  scale_color_brewer(palette =  "RdYlBu") + 
  geom_text_repel(aes(label = Gene_symbol), max.overlaps = 20) 


ggsave("Figures/Receptor_ligand_pairs.pdf", width = 15, height = 10, units = "cm") 



