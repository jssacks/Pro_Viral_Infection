

library(tidyverse)
library(readr)

#define inputs

#HILIC Positive
HILIC.Pos.MSDIAL.file <- "Data_Raw/MSDIAL/Pos/Area_1_20226141448.txt"
Sky.Pos.file1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_File1.csv"
Sky.Pos.file2 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_File2.csv"
Sky.Pos.search.file <- "Intermediates/Sky_search_files/Particulate_HPos_skysearch.txt"


#HILIC Negative
HILIC.Neg.MSDIAL.file <- "Data_Raw/MSDIAL/Neg/Area_1_2022616159.txt" #need to change
Sky.Neg.file1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_File1.csv"
Sky.Neg.file2 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_File2.csv"
Sky.Neg.search.file <- "Intermediates/Sky_search_files/Particulate_HNeg_skysearch.txt"


####Adduct Lists
adduct.file <- "Meta_Data/FiehnLab_Adduct_List_modified.csv"

####Source Functions
source("R_Code/Functions.R")

####Define QC parameters
int.diff.threshold <- 10 #Threshold of % difference between manual (Skyline) and automated (MSDIAL) integration deemed acceptable. Compounds above this threshold are manually integrated in all samples in skyline.
min.area <- 40000   #minimum peak area 
min.blk.ratio <- 3 #minimum ratio for a peak area to be above the blank
min.smp.number <-  12 #minimum numbers of samples an MF must be found in to pass QC
poo.smp.number <-  12   #number of pooled samples in dataset (both full and half poos)




####Pull in MS-DIAL HILIC Positive and HILIC Negative Data:
##Pos
MSDIAL.Pos.area <- MSDIAL_read(HILIC.Pos.MSDIAL.file, "HILICPos")  %>%
  select(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
         RT_matched, mz_matched, MS2_matched, SN_ave, `200929_Blk_DOM_blk_T0`:`200929_Smp_E_8G_T48`) %>%
  pivot_longer(-c(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
                  RT_matched, mz_matched, MS2_matched, SN_ave), names_to = "SampID", values_to = "Area")

##Neg
MSDIAL.Neg.area <-  MSDIAL_read(HILIC.Neg.MSDIAL.file, "HILICNeg")  %>%
  select(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
         RT_matched, mz_matched, MS2_matched, SN_ave, `200929_Blk_DOM_blk_T0`:`200929_Smp_E_8G_T48`) %>%
  pivot_longer(-c(ID, RT, mz, Name, Adduct, Note, Fill, MS2, Ref_RT, Ref_mz, 
                  RT_matched, mz_matched, MS2_matched, SN_ave), names_to = "SampID", values_to = "Area")



#####Pull in SKyline HILIC Positive and HILIC Negative Data:
#Pos
Sky.Pos.area.1 <- sky_read(Sky.Pos.file1) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.Pos.area.2 <- sky_read(Sky.Pos.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.Pos <- rbind(Sky.Pos.area.1, Sky.Pos.area.2)
Sky.Pos.poo <- Sky.Pos %>%
  filter(str_detect(.$SampID, "Poo")) %>%
  mutate(Poo = case_when(str_detect(.$SampID, "Half") ~ "half",
                         TRUE ~ "full"))

#Neg
Sky.Neg.area.1 <- sky_read(Sky.Neg.file1) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.Neg.area.2 <- sky_read(Sky.Neg.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound,
         "Sky.Area" = Area) %>%
  select(SampID, Name, Sky.Area)

Sky.Neg <- rbind(Sky.Neg.area.1, Sky.Neg.area.2)
Sky.Neg.poo <- Sky.Neg %>%
  filter(str_detect(.$SampID, "Poo")) %>% 
  mutate(Poo = case_when(str_detect(.$SampID, "Half") ~ "half",
                         TRUE ~ "full"))


####Calculate differences between Skyline and MSDIAL integrations in pooled samples to identify compounds to 
#   manually integrate (above a user defined average % difference threshold)
Pos.integration.check <- left_join(Sky.Pos.poo, MSDIAL.Pos.area) %>%
  mutate(offset = Area - Sky.Area) %>%
  group_by(ID, Name) %>%
  mutate(mean.offset = mean(offset),
         perc.offset = mean.offset/mean(Area)*100,
         RSD.offset = abs(sd(offset)/mean.offset*100)) %>%
  filter(abs(perc.offset) > int.diff.threshold) %>%
  select(Name, mean.offset, perc.offset, RSD.offset) %>%
  mutate(Fraction = "HILIC_Pos") %>%
  unique()

Neg.integration.check <- left_join(Sky.Neg.poo, MSDIAL.Neg.area) %>%
  mutate(offset = Area - Sky.Area) %>%
  group_by(ID, Name) %>%
  mutate(mean.offset = mean(offset),
         perc.offset = mean.offset/mean(Area)*100,
         RSD.offset = abs(sd(offset)/mean.offset*100)) %>%
  filter(abs(perc.offset) > int.diff.threshold) %>%
  select(Name, mean.offset, perc.offset, RSD.offset) %>%
  mutate(Fraction = "HILIC_Neg") %>%
  unique()

#combine HILIC_Pos and HILIC_Neg and export to csv
HILIC.integration.check <- rbind(Pos.integration.check, Neg.integration.check)
write_csv(HILIC.integration.check, file = "Intermediates/ProMo_Particulate_HILIC_Integration_Check.csv")


####Add in skyline compounds (Manually write in the names of compounds manually integrated in Skyline)

###Make compound lists
sky.pos.comps <- c("Creatine", "Glucosamine", "Hypotaurine", "L-Arginine", "L-Valine", "L-Serine", "L-Ornithine",
                   "L-Arginosuccinic acid", "L-Cystathionine", "L-Glutamic acid", "L-Isoleucine", "L-Homoserine",
                   "L-Methionine", "L-Methionine S-oxide", "Ophthalmic acid", "Proline betaine", "N6-Acetyl-L-lysine",
                   "Guanine", "Glycine betaine", "Glutatione disulfide", "Pipecolic acid", "beta-Alanine")

sky.neg.comps <- c("N-Acetyl-L-glutamic acid", "Adenosine monophosphate", "Citric acid", "D-Gluconic acid", "Malic acid",
                   "L-Cysteic acid", "Isethionic acid", "2-Ketoglutaric acid", "Salicylic acid", "Threonic acid",
                   "Succinic acid", "Taurine", "Thymine", "Uracil")

###Get mass and RT information for Skyline comps
########load in files
sky.pos.details <- read_tsv(Sky.Pos.search.file) %>%
  mutate(ID = paste("sky",row_number(), sep = ""))
sky.neg.details <- read_tsv(Sky.Neg.search.file) %>%
  mutate(ID = paste("sky",row_number(), sep = ""))

###Add in mass + RT info and pull out manually integrated compounds from skyline files
Sky.Pos.add <- left_join(Sky.Pos, sky.pos.details)%>%
  filter(Name %in% sky.pos.comps) %>%
  rename("mz"= `Pos-m/z`,
         "Area" = Sky.Area) %>%
  select(-Adduct) %>%
  filter(!is.na(ID))
  
Sky.Neg.add <- left_join(Sky.Neg, sky.neg.details) %>%
  filter(Name %in% sky.neg.comps) %>%
  rename("mz"= `Pos-m/z`,
         "Area" = Sky.Area) %>%
  select(-Adduct) %>%
  filter(!is.na(ID))

###Remove any skyline comps from MSDIAL data
MSDIAL.Pos.area.2 <- MSDIAL.Pos.area %>%
  filter(!Name %in% Sky.Pos.add$Name)

MSDIAL.Neg.area.2 <- MSDIAL.Neg.area %>%
  filter(!Name %in% Sky.Neg.add$Name)



####Add in skyline data and remove samples not relevant to this project (grazer treatments)
Pos.area <- full_join(MSDIAL.Pos.area.2, Sky.Pos.add) %>%
  filter(!str_detect(SampID, "4LGV")) %>%
  filter(!str_detect(SampID, "LG_")) %>%
  filter(!str_detect(SampID, "8G"))

Neg.area <- full_join(MSDIAL.Neg.area.2, Sky.Neg.add) %>%
  filter(!str_detect(SampID, "4LGV")) %>%
  filter(!str_detect(SampID, "LG_")) %>%
  filter(!str_detect(SampID, "8G"))




####Apply minimum area threshold to pooled samples, peak must be above this value in all pooled samples 
Pos.min.area.qc <- Pos.area %>%
  filter(str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "DDA")) %>%
  filter(Area > min.area) %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count == poo.smp.number)

Pos.area.2 <- left_join(Pos.min.area.qc, Pos.area) %>%
  select(-count)

Neg.min.area.qc <- Neg.area %>%
  filter(str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "DDA")) %>%
  filter(Area > min.area) %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count == poo.smp.number)

Neg.area.2 <- left_join(Neg.min.area.qc, Neg.area) %>%
  select(-count)

#####Apply qc requiring peak area in a sample to be 1) above minimum area threshold and 2) above blk threshold. MFs 
# will then be removed if not present in a minimum number of samples 

#Pos
Pos.blk.dat <- Pos.area.2 %>%
  filter(str_detect(.$SampID, "DOM")) %>%
  group_by(ID) %>%
  summarize(blk.av.area = mean(Area))

#Identify all MF-sample combinations that meet min.blk.ratio thresholds
Pos.blk.qc.smp <- left_join(Pos.area.2, Pos.blk.dat) %>%
  filter(str_detect(.$SampID, "Smp")) %>%
  select(ID, SampID, Area, blk.av.area) %>%
  mutate(blk.ratio = Area/(blk.av.area + 1)) %>%
  filter(blk.ratio >= min.blk.ratio) %>%
  filter(Area >= min.area) %>%
  select(ID, SampID)

#Identify all MFs that are present above min.blk.ratio in minimum number of samples 
Pos.blk.qc.MF <- Pos.blk.qc.smp %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count >= min.smp.number) %>%
  select(-count)

Pos.blk.qc <- left_join(Pos.blk.qc.MF, Pos.area.2)

#Add Pooled samples and blanks back into dataframe
Pos.area.poo.blk <- left_join(Pos.blk.qc.MF, Pos.area.2, by = "ID") %>%
  filter(str_detect(.$SampID, "Poo") |
           str_detect(.$SampID, "Blk")) %>%
  filter(!str_detect(.$SampID, "DDA"))

Pos.area.3 <- Pos.blk.qc %>%
  rbind(., Pos.area.poo.blk) %>%
  mutate(polarity = 1,
         fraction = "HPos") %>%
  mutate(ID = paste(fraction, ID, sep = "_"))


#Neg
Neg.blk.dat <- Neg.area.2 %>%
  filter(str_detect(.$SampID, "DOM")) %>%
  group_by(ID) %>%
  summarize(blk.av.area = mean(Area))

#Identify all MF-sample combinations that meet min.blk.ratio thresholds
Neg.blk.qc.smp <- left_join(Neg.area.2, Neg.blk.dat) %>%
  filter(str_detect(.$SampID, "Smp")) %>%
  select(ID, SampID, Area, blk.av.area) %>%
  mutate(blk.ratio = Area/(blk.av.area + 1)) %>%
  filter(blk.ratio >= min.blk.ratio) %>%
  filter(Area >= min.area) %>%
  select(ID, SampID)

#Identify all MFs that are present above min.blk.ratio in minimum number of samples 
Neg.blk.qc.MF <- Neg.blk.qc.smp %>%
  group_by(ID) %>%
  summarize(count = n()) %>%
  filter(count >= 12) %>%
  select(-count) 

Neg.blk.qc <- left_join(Neg.blk.qc.MF, Neg.area.2) 

#Add pooled samples and blanks back into dataframe 
Neg.area.poo.blk <- left_join(Neg.blk.qc.MF, Neg.area.2, by = "ID") %>%
  filter(str_detect(.$SampID, "Poo") |
        str_detect(.$SampID, "Blk")) %>%
  filter(!str_detect(.$SampID, "DDA"))

Neg.area.3 <- Neg.blk.qc %>%
  rbind(., Neg.area.poo.blk) %>%
  mutate(polarity = -1,
         fraction = "HNeg") %>%
  mutate(ID = paste(fraction, ID, sep = "_"))

#Join positive and negative data sets together:
Area.dat <- rbind(Pos.area.3, Neg.area.3)

write_csv(Area.dat, file = "Intermediates/ProMo_HILIC_QC1_output.csv")




