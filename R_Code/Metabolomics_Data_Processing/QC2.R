



#
#
library(tidyverse)
library(readr)

###define inputs

#hilic_data
hilic.BMIS.file <- "Intermediates/ProMo_Particulate_HILIC_BMISed_dat.csv"
hilic.QC1.file <- "Intermediates/ProMo_HILIC_QC1_output.csv"

#rp_data
rp.BMIS.file <- "Intermediates/ProMo_Particulate_RP_BMISed_dat.csv"
rp.QC1.file <- "Intermediates/ProMo_RP_QC1_output.csv"

#vit_data
#vit.QC1.file <- "Intermediates/ProMo_Particulate_Vit_QC1.csv"

#meta_data
vol.filt.file <- "Meta_Data/MEP_volume_filtered.csv"
abundance.file <- "Collaborator_Data/FCM/ProMo_Abundat_Combined.csv"

##Define cutoffs:
RSD.cutoff <- 0.30
min.reps <- 7
min.blk.ratio <- 3
min.area <- 40000

##read in area data
hilic.dat <- read_csv(hilic.BMIS.file)
rp.dat <- read_csv(rp.BMIS.file) 
#vit.dat <- read_csv(vit.QC1.file) %>%
#  select(-polarity, -fraction) %>%
#  rename("MF" = ID) %>% separate(SampID, 
#                                 c("runDate",
#                                   "type","samp","replicate"),"_", remove = FALSE) %>%
#  select(-runDate, -samp, -replicate) %>%
#  mutate(Adjusted_Area = case_when(str_detect(SampID, "Half") ~ Area*2,
#                                   TRUE ~ Area)) %>%
#  mutate(Name = MF)

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

full.sum <- full.dat %>%
  select(MF, fraction) %>%
  unique() %>%
  group_by(fraction) %>%
  summarize(count = n())

#####Run RSD test on Pooled samples to identify "high quality" mass features
rsd.test <- full.dat %>%
  select(MF, type, SampID, Adjusted_Area) %>%
  filter(type == "Poo") %>%
  group_by(MF) %>%
  mutate(Poo.RSD = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
  select(MF, type, Poo.RSD) %>%
  filter(Poo.RSD <= RSD.cutoff) %>%
  ungroup() %>%
  select(MF) %>%
  unique() 

##########
select.dat <- full.dat %>%
  filter(MF %in% rsd.test$MF) %>%
  mutate(SampID = str_replace(.$SampID, "201006_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "201013_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "200929_", "")) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  mutate(SampID = str_replace(.$SampID, "Smp_", "")) %>%
  select(-Orig_RSD, -type) 




##########
select.dat <- full.dat %>%
  filter(MF %in% rsd.test$MF) %>%
  mutate(SampID = str_replace(.$SampID, "201006_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "201013_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "200929_", "")) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  mutate(SampID = str_replace(.$SampID, "Smp_", "")) %>%
  select(-Orig_RSD, -type) 
  

####################

###
##Pull in abundance and volume filtered dat
abund.dat <- read_csv(abundance.file) %>%
  separate(Rep, c("biorep", "treatment_number")) %>%
  mutate(treatment = NA) %>%
  mutate(treatment = case_when(treatment_number == "1" & !biorep == "E" ~ "C",
                               treatment_number == "3" & !biorep == "E" ~ "LV",
                               treatment_number == "4" & !biorep == "E" ~ "LGV",
                               treatment_number == "6"& !biorep == "E" ~ "HV",
                               biorep == "E" ~ "LG")) %>%
  select(Time, biorep, treatment, everything()) %>%
  unite("Samp", c("treatment_number", "treatment"), sep = "", remove = FALSE) %>%
  mutate("T" = "T") %>%
  unite("time", c("T", "Time"), sep = "", remove = FALSE) %>%
  unite("SampID", c("biorep", "Samp", "time"), remove = FALSE) %>%
  select(-Samp)

vol.dat <- read_csv(vol.filt.file)  %>%
  rename("SampID" = Sample.Name) %>%
  mutate(SampID = str_replace(.$SampID, "201013_Smp_", ""))  

meta.dat <- left_join(vol.dat, abund.dat) %>%
  select(-T)

####combine meta data and peak list
dat <- full_join(meta.dat, select.dat)

########Run QC to require untargeted MFs (not annotated, name = "Unknown) be present in all three replicates of 
####     at least 7 treatment_timepoint combnations (ex. HV_T24)


#####Run min replicate number test to see if compounds are above minimum area and blk_poo ratio thresholds in a minimum number of 
# triplicates (treatment X timepoint combinations)

hilic.blk.dat <- read_csv(hilic.QC1.file) %>%
filter(str_detect(SampID, "DOM")) %>%
  rename("MF" = ID) %>%
  group_by(MF) %>%
  summarize(mean.blk.area = mean(Area))

rp.blk.dat <- read_csv(rp.QC1.file) %>%
  filter(str_detect(SampID, "DOM")) %>%
  rename("MF" = ID) %>%
  group_by(MF) %>%
  summarize(mean.blk.area = mean(Area))

#vit.blk.dat <- read_csv(vit.QC1.file) %>%
#  filter(str_detect(SampID, "DOM")) %>%
#  rename("MF" = ID) %>%
#  group_by(MF) %>%
#  summarize(mean.blk.area = mean(Area))

all.blk.dat <- rbind(hilic.blk.dat, rp.blk.dat)

####
blk.area.dat <- left_join(full.dat, all.blk.dat) %>%
  filter(str_detect(.$SampID, "Smp")) %>%
  select(MF, SampID, Area, mean.blk.area) %>%
  mutate(blk.ratio = Area/(mean.blk.area + 1)) %>%
  mutate(SampID = str_replace(.$SampID, "201006_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "201013_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "200929_", "")) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  mutate(SampID = str_replace(.$SampID, "Smp_", ""))

##Run QC to identify MFs meeting thresholds 
dat.rep.qc <- left_join(dat, blk.area.dat) %>%
  filter(Area >= min.area) %>%
  filter(blk.ratio >= min.blk.ratio) %>%
  group_by(MF, Name, treatment, Time) %>%
  summarize(count = n()) %>%
  filter(count == 3) %>%
  ungroup() %>%
  select(MF, Name) %>%
  group_by(MF, Name) %>%
  summarize(count = n()) %>%
  filter(count >= min.reps | !Name == "Unknown") %>%
  ungroup()



####Export final MF list
dat.export <- dat %>%
  filter(MF %in% dat.rep.qc$MF)

write_csv(dat.export, file = "Intermediates/ProMo_Particulate_QC2_data.csv")



####Export data for quantification

###write known MFxSample pairs that pass blk_qc to csv for use in quant script.
MF.samps.qc.for.quant <- left_join(dat.export, blk.area.dat) %>%
  filter(blk.ratio >= min.blk.ratio) %>%
  filter(!Name == "Unknown") %>%
  select(MF, Name, SampID)
  
write_csv(MF.samps.qc.for.quant, file = "Intermediates/QC2_MFSamps_forQuant.csv")

  
####Export final MF list
#dat.export <- dat %>%
 # filter(MF %in% dat.rep.qc$MF)

#test <- dat.export %>%
#  select(MF, fraction) %>%
#  unique() %>%
#  group_by(fraction) %>%
#  summarize(count = n())
  

#####Export data
#write_csv(dat.export, file = "Intermediates/ProMo_Particulate_QC2_data.csv")








#Supplementary Code


####Summarize total number of MFs
MF.sum <- dat.export %>%
  select(MF, fraction) %>%
  unique() %>%
  group_by(fraction) %>%
  summarize(count = n())

#Summarize total number of annotated MFs
MF.sum.ids <- dat.export %>%
  filter(!Name == "Unknown") %>%
  select(Name, fraction) %>%
  unique() %>% 
  group_by(fraction) %>%
  summarize(count = n())

###summarize total #s of MFs
dat.export %>%
  select(MF, fraction) %>%
  unique() %>%
  group_by(fraction) %>%
  summarize(count = n())
