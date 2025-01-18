





library(tidyverse)


#define inputs:
transcript.file <- "Collaborator_Data/Transcripts/MED4_annot_DE+non-DE_082823_normalized.csv"
transcript.file.2 <- "Collaborator_Data/Transcripts/HV_LFC+pval_MED4norm.csv"


#load data:
transcript.dat <- read_csv(transcript.file)
transcript.dat.2 <- read_csv(transcript.file.2)

#organize data:
t.sml <- transcript.dat.2 %>%
  select(GENE.NAME, KEGG, HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, PathwayID, Pathway_info, Details) %>%
  pivot_longer(cols = HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, names_to = "Samp_Name", values_to = "val") %>%
  mutate(param = case_when(str_detect(Samp_Name, "log2FoldChange") ~ "FC",
                           str_detect(Samp_Name, "padj") ~ "p_adj")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_log2FoldChange")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_padj")) %>%
  mutate(Samp_Name = as.numeric(str_remove(Samp_Name, "HV"))) %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) %>%
  pivot_wider(id_cols = c(GENE.NAME, Samp_Name, KEGG, PathwayID, Pathway_info, Details), names_from = param, values_from = val) %>%
  unique() %>%
  rename("Time" = Samp_Name,
         "Gene_Name" = GENE.NAME,
         "KEGG_ID" = KEGG,
         "Pathway_ID" = PathwayID,
         "Pathway_details" = Details) %>%
  select(Time, everything())

#export supplemental table:
write_csv(t.sml, file = "Tables/Outputs/Transcripts_FC_Supplemental_Table.csv")















































