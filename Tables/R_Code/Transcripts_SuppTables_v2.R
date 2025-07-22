





library(tidyverse)

#inputs:
med4.transcript.file <- "Intermediates/MED4_DE_transcripts.csv"
gene.ID.file <- "Collaborator_Data/Transcripts/ALl_MED4_ANNOT.csv"
transcript.counts.file <- "Collaborator_Data/Transcripts/Transcript_Data/Counts_MED4.csv"


#load in data
med4.transcript.dat <- read_csv(med4.transcript.file)
gene.ID.dat <- read_csv(gene.ID.file)
transcript.count.dat <- read_csv(transcript.counts.file)



#Clean up each data frame

#fold change and p-values
t.sml <- med4.transcript.dat %>%
  #  select(GENE.NAME, KEGG,      HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, PathwayID, Pathway_info, Details) %>%
  pivot_longer(cols = HV0_LFC_nonDE_padj:HV36_LFC_nonDE_log2FoldChange, names_to = "Samp_Name", values_to = "val") %>%
  mutate(param = case_when(str_detect(Samp_Name, "log2FoldChange") ~ "FC",
                           str_detect(Samp_Name, "padj") ~ "p_adj")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_log2FoldChange")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_padj")) %>%
  mutate(Samp_Name = as.numeric(str_remove(Samp_Name, "HV"))) %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) 

#unique Gene IDs
gene.IDs <- gene.ID.dat %>%
  select(GENE.NAME, ID) %>%
  rename("Gene.ID" = ID) %>%
  mutate(Gene.ID = str_remove(Gene.ID, "cds-"))

#count data
count.dat <- transcript.count.dat %>%
  pivot_longer(cols = CTRL0_A:TRT36_C, names_to = "Sample", values_to = "Normalized.Counts") %>%
  separate(col = Sample, into = c("Sample", "Replicate"), sep = "_") %>%
  mutate("Time" = case_when(str_detect(Sample, "0") ~ 0,
                            str_detect(Sample, "12") ~ 12,
                            str_detect(Sample, "24") ~ 24,
                            str_detect(Sample, "36") ~ 36
                            )) %>%
  mutate(Treatment = case_when(str_detect(Sample, "CTRL") ~ "C",
                               TRUE ~ "HV")) %>%
  select(GENE.NAME, Time, Treatment, Replicate, Normalized.Counts)



#Organize Fold Change Data:
FC.supptable <- left_join(gene.IDs, t.sml) %>%
  select(Gene.ID, GENE.NAME, KEGG, PathwayID, Pathway_info, Samp_Name, val, param) %>%
  unique() %>%
  filter(!is.na(val)) %>%
  pivot_wider(id_cols = c(Gene.ID, GENE.NAME, KEGG, PathwayID, Pathway_info, Samp_Name), names_from = param, values_from = val) %>%
  rename(Gene.Name = GENE.NAME,
         KEGG.ID = KEGG,
         KEGG.Pathway.ID = PathwayID,
         KEGG.Pathway.Description = Pathway_info,
         Time = Samp_Name,
         Log2FC = FC,
         p.adj = p_adj)


#export FC supplemental table:
write_csv(FC.supptable, file = "Tables/Outputs/Transcripts_FC_Supplemental_Table_060625.csv")




#Organize Count Supplemental Table 
count.supptable <- left_join(gene.IDs, count.dat) %>%
  filter(Gene.ID %in% FC.supptable$Gene.ID) %>%
  rename(Gene.Name = GENE.NAME) %>%
  select(Gene.ID, Gene.Name, Time, Treatment, Replicate, Normalized.Counts)

#export count supplemental table:
write_csv(count.supptable, file = "Tables/Outputs/Transcripts_Count_Supplemental_Table_060625.csv")




#export supplemental table:
write_csv(t.sml, file = "Tables/Outputs/Transcripts_FC_Supplemental_Table.csv")










