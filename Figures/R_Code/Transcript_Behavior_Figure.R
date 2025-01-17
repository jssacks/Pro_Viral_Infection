library(tidyverse)
library(viridis)
#library(vegan)
#library(ggdendro)

##
pal <- c("darkred", "gray", "steelblue")


###
transcript.file <- "Collaborator_Data/Transcripts/MED4_annot_DE+non-DE_082823_normalized.csv"
transcript.file.2 <- "Collaborator_Data/Transcripts/HV_LFC+pval_MED4norm.csv"

transcript.dat <- read_csv(transcript.file)
transcript.dat.2 <- read_csv(transcript.file.2)


t.sml <- transcript.dat.2 %>%
  select(GENE.NAME, KEGG, HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, PathwayID, Pathway_info, Details) %>%
  pivot_longer(cols = HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, names_to = "Samp_Name", values_to = "val") %>%
  mutate(param = case_when(str_detect(Samp_Name, "log2FoldChange") ~ "FC",
                           str_detect(Samp_Name, "padj") ~ "p_adj")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_log2FoldChange")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_padj")) %>%
  mutate(Samp_Name = as.numeric(str_remove(Samp_Name, "HV"))) %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) 


###
t.sml.sum <- t.sml %>%
  select(GENE.NAME, Samp_Name, val, param) %>%
  unique() %>%
  pivot_wider(id_cols = c(GENE.NAME, Samp_Name), names_from = param, values_from = val) %>%
  mutate(classification = case_when(p_adj < 0.05 & FC > 0 ~ "upregulated",
                                    p_adj < 0.05 & FC < 0 ~ "downregulated",
                                    TRUE ~ "unchanged")) %>%
  group_by(Samp_Name, classification) %>%
  reframe(count = n())


plot <- ggplot(t.sml.sum, aes(x = as.factor(Samp_Name), y = count, fill = classification)) +
  geom_col(color = "black", size = 0.1, width = 0.7) +
  scale_fill_manual(values = pal) +
  theme_test() +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  xlab("Timepoint (h)") +
  ylab("Number of genes")
plot


###export 
ggsave(plot, file = "Figures/Outputs/Transcript_behavior_plot.png",
       bg = "white", dpi = 600, units = "in", height = 5, width = 5)





