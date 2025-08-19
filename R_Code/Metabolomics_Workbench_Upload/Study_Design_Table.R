


library(tidyverse)
library(readr)

#define inputs

#HILIC
HILIC.Pos.MSDIAL.file <- "Data_Raw/MSDIAL/Pos/Area_1_20226141448.txt"

#RP
RP.MSDIAL.file <- "Data_Raw/MSDIAL/RP/Area_0_20226201356.txt"

#Vitamin
vit.file <- "Data_Raw/Skyline/ProMo_Vit_TQS_Results_Sept24.csv"


###Volumes filtered
vol.filt.file <- "Meta_Data/MEP_volume_filtered.csv"

####Source Functions
source("R_Code/Functions.R")


####Pull in MS-DIAL HILIC Positive Data:
MSDIAL.Pos.area <- MSDIAL_read(HILIC.Pos.MSDIAL.file, "HILICPos")  %>%
  select(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
         RT_matched, mz_matched, MS2_matched, SN_ave, `200929_Blk_DOM_blk_T0`:`200929_Smp_E_8G_T48`) %>%
  pivot_longer(-c(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
                  RT_matched, mz_matched, MS2_matched, SN_ave), names_to = "SampID", values_to = "Area")

####Pull in MS-DIAL RP Data:
##Pos
MSDIAL.RP.area <- MSDIAL_read(RP.MSDIAL.file, "HILICPos")  %>%
  select(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
         RT_matched, mz_matched, MS2_matched, SN_ave, starts_with("201006_")) %>%
  pivot_longer(-c(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
                  RT_matched, mz_matched, MS2_matched, SN_ave), names_to = "SampID", values_to = "Area")



#Pull in vitamin data:
vit.dat <- read_csv(vit.file)

#Pull in volume data:
vol.dat <- read_csv(vol.filt.file)










#Get Just sample names:
hilic.smp.names <- MSDIAL.Pos.area %>%
  select(SampID) %>%
  unique() %>%
  filter(!str_detect(SampID, "_E_"),
         !str_detect(SampID, "_4LGV_"),
         !str_detect(SampID, "DDA"),
         !str_detect(SampID, "_Filter_"),
         !str_detect(SampID, "_Media_"),
         !str_detect(SampID, "_8G_"),
         !str_detect(SampID, "_9V_"),
         !str_detect(SampID, "_Poo_"),
         !str_detect(SampID, "_Blk_")
         ) %>%
  mutate(Treatment = case_when(str_detect(SampID, "1C") ~ "C",
                               str_detect(SampID, "3LV") ~ "LV",
                               str_detect(SampID, "6HV") ~ "HV"),
         Replicate = case_when(str_detect(SampID, "_A_") ~ "A",
                               str_detect(SampID, "_B_") ~ "B",
                               str_detect(SampID, "_C_") ~ "C"),
         Timepoint = case_when(str_detect(SampID, "T0") ~ "0",
                               str_detect(SampID, "T6") ~ "6",
                               str_detect(SampID, "T12") ~ "12",
                               str_detect(SampID, "T18") ~ "18",
                               str_detect(SampID, "T24") ~ "24",
                               str_detect(SampID, "T30") ~ "30",
                               str_detect(SampID, "T36") ~ "36",
                               str_detect(SampID, "T42") ~ "42",
                               str_detect(SampID, "T48") ~ "48")) %>%
  unite(c("Replicate", "Treatment", "Timepoint"), col = "Samp_Name", remove = FALSE) %>%
  rename("HILIC_SampID" = SampID)

rp.smp.names <- MSDIAL.RP.area %>%
  select(SampID) %>%
  unique() %>%
  filter(!str_detect(SampID, "_E_"),
         !str_detect(SampID, "_4LGV_"),
         !str_detect(SampID, "DDA"),
         !str_detect(SampID, "_Filter_"),
         !str_detect(SampID, "_Media_"),
         !str_detect(SampID, "_8G_"),
         !str_detect(SampID, "_9V_"),
         !str_detect(SampID, "_Poo_"),
         !str_detect(SampID, "_blk_")
  ) %>%
  mutate(Treatment = case_when(str_detect(SampID, "1C") ~ "C",
                               str_detect(SampID, "3LV") ~ "LV",
                               str_detect(SampID, "6HV") ~ "HV"),
         Replicate = case_when(str_detect(SampID, "_A_") ~ "A",
                               str_detect(SampID, "_B_") ~ "B",
                               str_detect(SampID, "_C_") ~ "C"),
         Timepoint = case_when(str_detect(SampID, "T0") ~ "0",
                               str_detect(SampID, "T6") ~ "6",
                               str_detect(SampID, "T12") ~ "12",
                               str_detect(SampID, "T18") ~ "18",
                               str_detect(SampID, "T24") ~ "24",
                               str_detect(SampID, "T30") ~ "30",
                               str_detect(SampID, "T36") ~ "36",
                               str_detect(SampID, "T42") ~ "42",
                               str_detect(SampID, "T48") ~ "48")) %>%
  unite(c("Replicate", "Treatment", "Timepoint"), col = "Samp_Name", remove = FALSE)  %>%
  rename("RP_SampID" = SampID)

vit.smp.names <- vit.dat %>%
  rename("SampID" = Replicate) %>%
  select(SampID) %>%
  unique() %>%
  filter(!str_detect(SampID, "_E_"),
         !str_detect(SampID, "_4LGV_"),
         !str_detect(SampID, "DDA"),
         !str_detect(SampID, "_Filter_"),
         !str_detect(SampID, "_Media_"),
         !str_detect(SampID, "_8G_"),
         !str_detect(SampID, "_9V_"),
         !str_detect(SampID, "_Std_"),
         !str_detect(SampID, "_Poo_"),
         !str_detect(SampID, "_Blk_")
  ) %>%
  mutate(Treatment = case_when(str_detect(SampID, "1C") ~ "C",
                               str_detect(SampID, "3LV") ~ "LV",
                               str_detect(SampID, "6HV") ~ "HV"),
         Replicate = case_when(str_detect(SampID, "_A_") ~ "A",
                               str_detect(SampID, "_B_") ~ "B",
                               str_detect(SampID, "_C_") ~ "C"),
         Timepoint = case_when(str_detect(SampID, "T0") ~ "0",
                               str_detect(SampID, "T6") ~ "6",
                               str_detect(SampID, "T12") ~ "12",
                               str_detect(SampID, "T18") ~ "18",
                               str_detect(SampID, "T24") ~ "24",
                               str_detect(SampID, "T30") ~ "30",
                               str_detect(SampID, "T36") ~ "36",
                               str_detect(SampID, "T42") ~ "42",
                               str_detect(SampID, "T48") ~ "48")) %>%
  unite(c("Replicate", "Treatment", "Timepoint"), col = "Samp_Name", remove = FALSE)  %>%
  rename("Vitamin_SampID" = SampID)

##Volume data:
vol.smp.name <- vol.dat%>%
  rename("SampID" = Sample.Name) %>%
  select(SampID, Vol.filt.mL) %>%
  unique() %>%
  filter(!str_detect(SampID, "_E_"),
         !str_detect(SampID, "_4LGV_"),
         !str_detect(SampID, "DDA"),
         !str_detect(SampID, "_Filter_"),
         !str_detect(SampID, "_Media_"),
         !str_detect(SampID, "_8G_"),
         !str_detect(SampID, "_9V_"),
         !str_detect(SampID, "_Std_"),
         !str_detect(SampID, "_Poo_"),
         !str_detect(SampID, "_Blk_")
  ) %>%
  mutate(Treatment = case_when(str_detect(SampID, "1C") ~ "C",
                               str_detect(SampID, "3LV") ~ "LV",
                               str_detect(SampID, "6HV") ~ "HV"),
         Replicate = case_when(str_detect(SampID, "_A_") ~ "A",
                               str_detect(SampID, "_B_") ~ "B",
                               str_detect(SampID, "_C_") ~ "C"),
         Timepoint = case_when(str_detect(SampID, "T0") ~ "0",
                               str_detect(SampID, "T6") ~ "6",
                               str_detect(SampID, "T12") ~ "12",
                               str_detect(SampID, "T18") ~ "18",
                               str_detect(SampID, "T24") ~ "24",
                               str_detect(SampID, "T30") ~ "30",
                               str_detect(SampID, "T36") ~ "36",
                               str_detect(SampID, "T42") ~ "42",
                               str_detect(SampID, "T48") ~ "48")) %>%
  unite(c("Replicate", "Treatment", "Timepoint"), col = "Samp_Name", remove = FALSE) %>%
  select(-SampID)






##Combine all data together:
studytable <- left_join(hilic.smp.names, rp.smp.names) %>%
  left_join(., vit.smp.names) %>%
  left_join(., vol.smp.name) %>%
  rename("Volume_mL" = Vol.filt.mL,
         "Timepoint_h" = Timepoint) %>%
  mutate(Sample_Type = "Marine Microbe Culture") %>%
  select(Samp_Name, Sample_Type, Replicate, Treatment, Timepoint_h, Volume_mL, HILIC_SampID, RP_SampID, Vitamin_SampID)


#Save study design table:
write_csv(studytable, file = "Intermediates/MW_StudyDesignTable.csv")
write_tsv(studytable, file = "Tables/Metabolomics_Workbench/MW_StudyDesignTable.txt")
  
