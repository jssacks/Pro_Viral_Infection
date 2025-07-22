



library(tidyverse)



###define inputs

#hilic_data
hilic.BMIS.file <- "Intermediates/ProMo_Particulate_HILIC_BMISed_dat.csv"
hilic.QC1.file <- "Intermediates/ProMo_HILIC_QC1_output.csv"

#rp_data
rp.BMIS.file <- "Intermediates/ProMo_Particulate_RP_BMISed_dat.csv"
rp.QC1.file <- "Intermediates/ProMo_RP_QC1_output.csv"

#vit_data
#

#meta_data
vol.filt.file <- "Meta_Data/MEP_volume_filtered.csv"
abundance.file <- "Collaborator_Data/FCM/ProMo_Abundat_Combined.csv"



##read in area data
hilic.dat <- read_csv(hilic.BMIS.file)
rp.dat <- read_csv(rp.BMIS.file) 

#read in MF details data
hilic.qc1 <- read_csv(hilic.QC1.file) %>%
  select(ID, Name, RT, mz, MS2) %>%
  rename("MF" = ID) %>%
  unique()


rp.qc1 <- read_csv(rp.QC1.file) %>%
  select(ID, Name, RT, mz, MS2) %>%
  unique() %>%
  rename("MF" = ID) 

rp.blk.dat <- read_csv(hilic.QC1.file) %>%
  filter(str_detect(SampID, "DOM")) %>%
  rename("MF" = ID) %>%
  group_by(MF) %>%
  summarize(mean.blk.area = mean(Area))

##Add in MF information to hilic and rp data and combine into a single dataframe
hilic.info.dat <- left_join(hilic.dat, hilic.qc1)
rp.info.dat <- left_join(rp.dat, rp.qc1) %>%
  filter(!str_detect(SampID, "_blk"))

full.dat <- rbind(hilic.info.dat, rp.info.dat) %>%
  select(-runDate, -samp, -replicate)
full.dat <- full.dat %>%
  mutate(fraction = case_when(str_detect(.$MF, "HNeg") ~ "HILIC_Neg",
                              str_detect(.$MF, "HPos") ~ "HILIC_Pos",
                              str_detect(.$MF, "RP") ~ "RP",
                              str_detect(.$MF, "Vit") ~ "Vit"))




##Get just samples:
HILIC.samples <- full.dat %>%
  select(SampID, fraction) %>%
  unique() %>%
  filter(str_detect(fraction, "HILIC_Pos"))

RP.samples <- full.dat %>%
  select(SampID, fraction) %>%
  unique()

#Vit.samples <- XXX





#library(here)
library(readr)
library(tidyverse)

#Load Functions
source("R_Code/Functions.R")

## This script identifies targeted compounds present in the pooled samples, determines their batch- and sample-specific mass,
## accounting for mass defects, and their batch- and sample-specific retention time. This information is then exported as a 
## tab-separated .txt file that can be imported into MS DIAL to identify these targeted compounds in the untargeted dataset.


# Define files -----------------------------------------------------------
# Define your Skyline output and transition list filepaths from the first set of manual integrations.

## HILIC Positive
HPos.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_File1.csv"
HPos.tlist <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_TransitionList.csv"

## HILIC Negative
HNeg.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_File1.csv"
HNeg.tlist <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_TransitionList.csv"

## Reverse Phase
RP.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_RP_File1.csv"
RP.tlist <- "Data_Raw/Skyline/ProMo_Particulate_RP_TransitionList.csv"

## Min Area Threshold
min.area <- 40000

# Import defined files -----------------------------------------------------------

## HILIC Positive
HPos.samps <- skyline_read(HPos.file.1) %>%
  select(Rep) %>%
  unique() %>%
  filter(!str_detect(Rep, "4LGV"),
         !str_detect(Rep, "8G"),
         !str_detect(Rep, "9V"),
         !str_detect(Rep, "Std_"),
         !str_detect(Rep, "Blk_Filter"),
         !str_detect(Rep, "Blk_Media")
         ) %>%
  mutate(SampName = str_remove(Rep, "200929_"),
         SampName = str_remove(SampName, "Blk_"),
         SampName = str_remove(SampName, "Poo_"),
         SampName = str_remove(SampName, "Smp_")) %>%
  mutate(Treatment = case_when(str_detect(SampName, "1C") ~ "C",
                               str_detect(SampName, "3LV") ~ "LV",
                               str_detect(SampName, "6HV") ~ "HV",
                               TRUE ~ NA))
         
  
  

## HILIC Negative
HNeg.samps <- skyline_read(HNeg.file.1) %>%
  select(Rep) %>%
  unique() %>%
  filter(!str_detect(Rep, "4LGV"),
         !str_detect(Rep, "8G"),
         !str_detect(Rep, "9V"),
         !str_detect(Rep, "Std_"),
         !str_detect(Rep, "Blk_Filter"),
         !str_detect(Rep, "Blk_Media")
  )

## Reverse Phase
RP.samps <- skyline_read(RP.file.1) %>%
  select(Rep) %>%
  unique() %>%
  filter(!str_detect(Rep, "4LGV"),
         !str_detect(Rep, "8G"),
         !str_detect(Rep, "9V"),
         !str_detect(Rep, "Std_"),
         !str_detect(Rep, "Blk_Filter"),
         !str_detect(Rep, "Media_blk")
  )



















