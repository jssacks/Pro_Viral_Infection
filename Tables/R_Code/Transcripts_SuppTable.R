





library(tidyverse)


#define inputs:
transcript.file <- "Collaborator_Data/Transcripts/MED4_annot_DE+non-DE_082823_normalized.csv"
transcript.file.2 <- "Collaborator_Data/Transcripts/HV_LFC+pval_MED4norm.csv"
transcript.count.file <-  "Collaborator_Data/Transcripts/Transcript_Counts.csv"
transcript.sample.names.file <- "Collaborator_Data/Transcripts/Transcript_Sample_Log.csv"
gene.ID.file <- "Collaborator_Data/Transcripts/ALl_MED4_ANNOT.csv"

#load data:
transcript.dat <- read_csv(transcript.file)
transcript.dat.2 <- read_csv(transcript.file.2)
transcript.count.dat <- read_csv(transcript.count.file)
transcript.sample.names <- read_csv(transcript.sample.names.file)
gene.ID.dat <- read_csv(gene.ID.file)


##get unique gene IDs:
gene.IDs <- gene.ID.dat %>%
  select(GENE.NAME, ID, start, end, X, strand) %>%
  rename("gene.unique.ID" = ID)


###Organize transcript counts and sample names
count.dat <- transcript.count.dat %>%
  pivot_longer(cols = S01:S48, names_to = "SampName", values_to = "Count") %>%
  left_join(transcript.sample.names) %>%
  filter(!is.na(Treatment)) %>%
  left_join(., gene.IDs) %>%
  select(GENE.NAME, gene.unique.ID, start, end, X, strand, SampName, Treatment, Timepoint, Replicate, Count) %>%
  filter(!is.na(gene.unique.ID)) %>%
  filter(!Treatment == "LGV")


###Transcript Fold Change data
transcript.fc.dat <- transcript.dat.2 %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) %>%
  left_join(., gene.IDs)  %>%
  filter(!is.na(gene.unique.ID))


#organize fold change data:
t.sml <- transcript.fc.dat %>%
  select(GENE.NAME, KEGG, HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, PathwayID, Pathway_info, Details, gene.unique.ID) %>%
  pivot_longer(cols = HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, names_to = "Samp_Name", values_to = "val") %>%
  mutate(param = case_when(str_detect(Samp_Name, "log2FoldChange") ~ "FC",
                           str_detect(Samp_Name, "padj") ~ "p_adj")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_log2FoldChange")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_padj")) %>%
  mutate(Samp_Name = as.numeric(str_remove(Samp_Name, "HV"))) %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) %>%
  pivot_wider(id_cols = c(GENE.NAME, Samp_Name, KEGG, PathwayID, Pathway_info, Details, gene.unique.ID), names_from = param, values_from = val) %>%
  unique() %>%
  rename("Time" = Samp_Name,
         "Gene_Name" = GENE.NAME,
         "Gene_Unique_ID" = gene.unique.ID,
         "KEGG_ID" = KEGG,
         "Pathway_ID" = PathwayID,
         "Pathway_details" = Details) %>%
  select(Time, Gene_Name, Gene_Unique_ID, everything())

#export supplemental table:
write_csv(t.sml, file = "Tables/Outputs/Transcripts_FC_Supplemental_Table.csv")















































