

library(tidyverse)
library(readr)

#define inputs

#HILIC Positive
RP.MSDIAL.file <- "Data_Raw/MSDIAL/RP/Area_0_20226201356.txt"
Sky.RP.file1 <- "Data_Raw/Skyline/ProMo_Particulate_RP_File1.csv"
Sky.RP.file2 <- "Data_Raw/Skyline/ProMo_Particulate_RP_File2.csv"
Sky.RP.search.file <- "Intermediates/Sky_search_files/Particulate_RP_skysearch.txt"


####Adduct Lists
adduct.file <- "Meta_Data/FiehnLab_Adduct_List_modified.csv"


###Volumes filtered
vol.filt.file <- "Meta_Data/Sample_Lists/MEP_volume_filtered.csv"

####Source Functions
source("R_Code/Functions.R")

####Define QC parameters
int.diff.threshold <- 10 #Threshold of % difference between manual (Skyline) and automated (MSDIAL) integration deemed acceptable. Compounds above this threshold are manually integrated in all samples in skyline.
min.area <- 40000   #minimum peak area 
min.blk.ratio <- 3 #minimum ratio for a peak area to be above the blank
min.smp.number <-  12 #minimum numbers of samples an MF must be found in to pass QC
poo.smp.number <-  8   #number of pooled samples in dataset (both full and half poos)

##Define adduct detection parameters
adduct.error.ppm <- 5   ###PPM Error tolerance
RT_ad_tol <- 0.05   #Retention time tolerance
adduct.cor.threshold <- 0.95  #Adduct correlation threshold



####Pull in MS-DIAL HILIC Positive and HILIC Negative Data:
##Pos
MSDIAL.RP.area <- MSDIAL_read(RP.MSDIAL.file, "HILICPos")  %>%
  select(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
         RT_matched, mz_matched, MS2_matched, SN_ave, starts_with("201006_")) %>%
  pivot_longer(-c(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
                  RT_matched, mz_matched, MS2_matched, SN_ave), names_to = "SampID", values_to = "Area")

#####Pull in SKyline RP Data:

#RP
Sky.RP.area.1 <- sky_read(Sky.RP.file1) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.RP.area.2 <- sky_read(Sky.RP.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.RP <- rbind(Sky.RP.area.1, Sky.RP.area.2)
Sky.RP.poo <- Sky.RP %>%
  filter(str_detect(.$SampID, "Poo"))




####Calculate skyline area  correction to targeted MS-DIAL data to account for differences in baseline subtraction:
#RP.correction <- left_join(Sky.RP.poo, MSDIAL.RP.area) %>%
#  mutate(offset = Area - Sky.Area) %>%
 # group_by(ID, Name) %>%
#  mutate(mean.offset = mean(offset),
#         perc.offset = mean.offset/mean(Area)*100) %>%
#  select(Name, Area, Sky.Area, mean.offset, perc.offset) %>%
#  select(Name, mean.offset, perc.offset) %>%
#  unique()
#  ungroup() %>%
#  filter(!is.na(mean.offset))

  
RP.integration.check <- left_join(Sky.RP.poo, MSDIAL.RP.area) %>%
    mutate(offset = Area - Sky.Area) %>%
    group_by(ID, Name) %>%
    mutate(mean.offset = mean(offset),
           perc.offset = mean.offset/mean(Area)*100,
           RSD.offset = abs(sd(offset)/mean.offset*100)) %>%
    filter(abs(perc.offset) > int.diff.threshold) %>%
    select(Name, mean.offset, perc.offset, RSD.offset) %>%
    mutate(Fraction = "HILIC_Pos") %>%
    unique()

#Write to .csv
write_csv(RP.integration.check, file = "Intermediates/ProMo_Particulate_RP_Integration_Check.csv")


####Add in skyline compounds (Manually write in the names of compounds manually integrated in Skyline)
###Make compound lists
sky.RP.comps <- c("S-Adenosylhomocysteine", "S-Adenosylmethionine", "Agmatine",
                   "Pyridoxal", "Nicotinic acid", "Pantothenic acid", "Pyridoxine", "Xanthine",
                  "Vanillin")



###Get mass and RT information for Skyline comps
########load in files
sky.RP.details <- read_tsv(Sky.RP.search.file) %>%
  mutate(ID = paste("sky",row_number(), sep = ""))

###Add in mass + RT info and pull out manually integrated compounds from skyline files
Sky.RP.add <- left_join(Sky.RP, sky.RP.details)%>%
  filter(Name %in% sky.RP.comps) %>%
  rename("mz"= `Pos-m/z`,
         "Area" = Sky.Area) %>%
  select(-Adduct)


###Remove any skyline comps from MSDIAL data
MSDIAL.RP.area.2 <- MSDIAL.RP.area %>%
  filter(!Name %in% Sky.RP.add$Name)

####Add in skyline data and remove samples not relevant to this project (grazer treatments)
RP.area <- full_join(MSDIAL.RP.area.2, Sky.RP.add) %>%
  filter(!str_detect(SampID, "4LGV")) %>%
  filter(!str_detect(SampID, "LG_")) %>%
  filter(!str_detect(SampID, "8G"))



####Apply minimum area threshold to pooled samples
RP.min.area.qc <- RP.area %>%
  filter(str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "DDA")) %>%
  filter(Area > min.area) %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count == poo.smp.number)

RP.area.2 <- left_join(RP.min.area.qc, RP.area) %>%
  select(-count)



#####Apply qc requiring peak area in a sample to be 1) above minimum area threshold and 2) above blk threshold. MFs 
# will then be removed if not present in a minimum number of samples 

#Pos
RP.blk.dat <- RP.area.2 %>%
  filter(str_detect(.$SampID, "DOM")) %>%
  group_by(ID) %>%
  summarize(blk.av.area = mean(Area))

#Identify all MF-sample combinations that meet min.blk.ratio thresholds
RP.blk.qc.smp <- left_join(RP.area.2, RP.blk.dat) %>%
  filter(str_detect(.$SampID, "Smp")) %>%
  select(ID, SampID, Area, blk.av.area) %>%
  mutate(blk.ratio = Area/(blk.av.area + 1)) %>%
  filter(blk.ratio >= min.blk.ratio) %>%
  filter(Area >= min.area) %>%
  select(ID, SampID)

#Identify all MFs that are present above min.blk.ratio in minimum number of samples triplicates (treatmentXtimepoint combinations)
RP.blk.qc.MF <- RP.blk.qc.smp %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count >= min.smp.number) %>%
  select(-count)

RP.blk.qc <- left_join(RP.blk.qc.MF, RP.area.2)


#Add Pooled samples and blanks back into dataframe
RP.area.poo.blk <- left_join(RP.blk.qc.MF, RP.area.2, by = "ID") %>%
  filter(str_detect(.$SampID, "Poo") |
           str_detect(.$SampID, "Blk") |
           str_detect(.$SampID, "blk")) %>%
  filter(!str_detect(.$SampID, "DDA"))

RP.area.3 <- left_join(RP.blk.qc, RP.area.2)  %>%
  rbind(., RP.area.poo.blk) %>%
  mutate(polarity = 1,
         fraction = "RP") %>%
  mutate(ID = paste(fraction, ID, sep = "_"))


####Final MF
Area.dat <- RP.area.3

write_csv(Area.dat, file = "Intermediates/ProMo_RP_QC1_output.csv")






