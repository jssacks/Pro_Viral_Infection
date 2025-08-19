




library(tidyverse)





#Load in data:
untargeted.dat <- read_csv("Intermediates/ProMo_Particulate_QC2_data.csv")
sdt <- read_csv("Intermediates/MW_StudyDesignTable.csv")
targeted.dat <- read_csv("Intermediates/Final_Processed_Quantified_Data.csv")
standards <- read_csv("Meta_Data/Ingalls_Lab_Standards_03172023.csv")


#Organize Untargeted data for HILIC Pos, HILIC Neg, and RP:
qc2.dat.mod <- untargeted.dat %>%
  filter(!is.na(Time)) %>%
  unite(c("biorep", "treatment", "Time"), col = "Samp_Name", remove = FALSE) %>%
  unite(c("mz", "RT"), col = "mz_rt", remove = FALSE) %>%
  rename("rt_min" = RT,
         "Compound_Name" = Name) %>%
  mutate(Normalized_Area = Adjusted_Area/Vol.filt.mL) %>%
  select(fraction, mz_rt, mz, rt_min, Compound_Name, Samp_Name, Normalized_Area) %>%
  left_join(., sdt) %>%
  mutate(Compound_Name = str_remove(Compound_Name, "L-"))

#Make untargeted data tables: 
hpos.u.dat <- qc2.dat.mod %>%
  filter(fraction == "HILIC_Pos") %>%
  select(Samp_Name, mz_rt, Normalized_Area) %>%
  pivot_wider(id_cols = mz_rt, names_from = Samp_Name, values_from = Normalized_Area)

hneg.u.dat <- qc2.dat.mod %>%
  filter(fraction == "HILIC_Neg") %>%
  select(Samp_Name, mz_rt, Normalized_Area) %>%
  pivot_wider(id_cols = mz_rt, names_from = Samp_Name, values_from = Normalized_Area)

rp.u.dat <- qc2.dat.mod %>%
  filter(fraction == "RP") %>%
  select(Samp_Name, mz_rt, Normalized_Area) %>%
  pivot_wider(id_cols = mz_rt, names_from = Samp_Name, values_from = Normalized_Area)


#Organize Targeted Data for HILIC Pos, HILIC Neg, RP, and Vitamin data:
targeted.dat.mod <- targeted.dat %>%
  rename(Compound_Name = Name) %>%
  mutate(fraction = case_when(str_detect(MF, "Neg") ~ "HILIC_Neg",
                              str_detect(MF, "Pos") ~ "HILIC_Pos",
                              str_detect(MF, "RP") ~ "RP",
                              TRUE ~ "Vit"),
         ) %>%
  mutate(Compound_Name = str_remove(Compound_Name, "L-")) %>%
  mutate(Samp_Name = str_replace(SampID, "1C", "C"),
         Samp_Name = str_replace(Samp_Name, "3LV", "LV"),
         Samp_Name = str_replace(Samp_Name, "6HV", "HV"),
         Samp_Name = str_remove(Samp_Name, "T")
         ) %>%
  rename(Concentration_nM = nM.in.smp) %>%
  mutate(Concentration_nM = case_when(Concentration_nM < 0 ~ NA,
                                      TRUE ~ Concentration_nM)) %>%
  select(fraction, Samp_Name, Compound_Name, Concentration_nM) %>%
  left_join(., sdt) %>%
  filter(!is.na(Replicate))



#Make targeted data tables:
hpos.t.dat <- targeted.dat.mod %>%
  filter(fraction == "HILIC_Pos") %>%
  select(Samp_Name, Compound_Name, Concentration_nM) %>%
  pivot_wider(id_cols = Compound_Name, names_from = Samp_Name, values_from = Concentration_nM)

hneg.t.dat <- targeted.dat.mod %>%
  filter(fraction == "HILIC_Neg") %>%
  select(Samp_Name, Compound_Name, Concentration_nM) %>%
  pivot_wider(id_cols = Compound_Name, names_from = Samp_Name, values_from = Concentration_nM)

rp.t.dat <- targeted.dat.mod %>%
  filter(fraction == "RP") %>%
  select(Samp_Name, Compound_Name, Concentration_nM) %>%
  pivot_wider(id_cols = Compound_Name, names_from = Samp_Name, values_from = Concentration_nM)

vit.t.dat <- targeted.dat.mod %>%
  filter(fraction == "Vit") %>%
  filter(!Compound_Name %in% c("Pyridoxamine", "Vitamin B2", "Vitamin B5")) %>%
  select(Samp_Name, Compound_Name, Concentration_nM) %>%
  pivot_wider(id_cols = Compound_Name, names_from = Samp_Name, values_from = Concentration_nM) %>%
  mutate(Compound_Name = case_when(Compound_Name == "OH-B12" ~ "Hydroxocobalamin",
                                   Compound_Name == "OH-Pseudocob" ~ "Hyrodoxopseudocobalamin",
                                   Compound_Name == "CN-Pseudocob" ~ "Cyanopseudocobalamin"))



#Make metabolite metadata tables:
standards.mod <- standards %>%
  mutate(fraction = case_when(Column == "RP" ~ "RP",
                              Column == "HILIC" & z == -1 ~ "HILIC_Neg",
                              Column == "HILIC" & z == 1 ~ "HILIC_Pos")) %>%
  mutate(fraction = case_when(Compound_Name == "Hydroxocobalamin" ~ "Vit",
                              TRUE ~ fraction)) %>%
  mutate(Compound_Name = str_remove(Compound_Name, "L-"))

#Pseudocobalamin metadata
pb12.metadata <- tibble(
  Compound_Name = c("Hyrodoxopseudocobalamin", "Cyanopseudocobalamin"),
  fraction = c("Vit", "Vit"),
  mz = c(659.7750, 673.2805),
  RT_minute = c(3.50, 4.30)
)


#Organize all metadata:
all.metab.metadata <- targeted.dat.mod %>%
  filter(!Compound_Name %in% c("Pyridoxamine", "Vitamin B2", "Vitamin B5")) %>%
  mutate(Compound_Name = case_when(Compound_Name == "OH-B12" ~ "Hydroxocobalamin",
                                   Compound_Name == "OH-Pseudocob" ~ "Hyrodoxopseudocobalamin",
                                   Compound_Name == "CN-Pseudocob" ~ "Cyanopseudocobalamin",
                                   TRUE ~ Compound_Name)) %>%
  select(Compound_Name, fraction) %>%
  unique() %>%
  left_join(., standards.mod %>% full_join(., pb12.metadata)) %>%
  select(fraction, Compound_Name, mz, RT_minute, PubChem_Code, KEGG_Code)

#Make metadata table:
hpos.meta.dat <- all.metab.metadata %>%
  filter(fraction == "HILIC_Pos") %>%
  select(-fraction)

hneg.meta.dat <- all.metab.metadata %>%
  filter(fraction == "HILIC_Neg") %>%
  select(-fraction)

rp.meta.dat <- all.metab.metadata %>%
  filter(fraction == "RP") %>%
  select(-fraction)

vit.meta.dat <- all.metab.metadata %>%
  filter(fraction == "Vit") %>%
  select(-fraction)


##Write Everything to a folder

#untargeted:
write_tsv(hpos.u.dat, file = "Tables/Metabolomics_Workbench/HILICPos_Untargeted_Results.txt")
write_tsv(hneg.u.dat, file = "Tables/Metabolomics_Workbench/HILICNeg_Untargeted_Results.txt")
write_tsv(rp.u.dat, file = "Tables/Metabolomics_Workbench/RP_Untargeted_Results.txt")

#quantified:
write_tsv(hpos.t.dat, file = "Tables/Metabolomics_Workbench/HILICPos_Quantified_Results.txt")
write_tsv(hneg.t.dat, file = "Tables/Metabolomics_Workbench/HILICNeg_Quantified_Results.txt")
write_tsv(rp.t.dat, file = "Tables/Metabolomics_Workbench/RP_Quantified_Results.txt")
write_tsv(vit.t.dat, file = "Tables/Metabolomics_Workbench/Vitamin_Quantified_Results.txt")

#metadata:
write_tsv(hpos.meta.dat, file = "Tables/Metabolomics_Workbench/HILICPos_Metab_Metadata.txt")
write_tsv(hneg.meta.dat, file = "Tables/Metabolomics_Workbench/HILICNeg_Metab_Metadata.txt")
write_tsv(rp.meta.dat, file = "Tables/Metabolomics_Workbench/RP_Metab_Metadata.txt")
write_tsv(vit.meta.dat, file = "Tables/Metabolomics_Workbench/Vitamin_Metab_Metadata.txt")





