
library(tidyverse)
library(broom)
####Source Functions
source("R_Code/Functions.R")


###Define inputs
vit.file <- "Data_Raw/Skyline/ProMo_Vit_TQS_Results_Sept24.csv"

#meta_data
vol.filt.file <- "Meta_Data/MEP_volume_filtered.csv"
abundance.file <- "Collaborator_Data/FCM/ProMo_Abundat_Combined.csv"
tqs.transitions.file <- "Meta_Data/TQS_Vit_Transitions.csv"

###Define Thresholds 
blk.ratio <- 3
min.area <- 40000
min.reps <-  7


###Quantification Parameters
#CN.B12.pB12.ratio <- ###????####
OH_B12.OH_pB12.ratio <- 1.675    #average of two calibration curves values (1.54 and 1.81) from KRH



###load vit monitored transitons
t.dat <- read_csv(tqs.transitions.file)


#Load in and QC datasets 
vit.dat <- read_csv(vit.file) %>%
  #filter(Peptide %in% c("OH-Pseudocob", "CN-Pseudocob", "OH-B12", "CN-B12",)) #%>%
  rename("Compound" = "Peptide",
         "SampID" = "Replicate",
         "Product_mz" = "Product Mz") %>%
  filter(!SampID == "230511_Poo_TruePooMortality_Full4") %>%
  filter(Product_mz %in% t.dat$Product_mz) %>%
  select(Compound, SampID, Area)


####Normalize B1 and B2 to internal standards

#Calculate IS normalization factor
vit.IS.dat <- vit.dat %>%
  filter(Compound %in% c("B1_IS", "B2_IS")) %>%
  group_by(Compound) %>%
  mutate(IS.norm = Area/mean(Area)) %>%
  mutate(match = case_when(Compound == "B1_IS" ~ "Vitamin B1",
                           Compound == "B2_IS" ~ "Vitamin B2")) %>%
  ungroup() %>%
  select(match, SampID, IS.norm) %>%
  rename("Compound" = match)
  
#apply IS normalization factor 
IS.norm <- left_join(vit.IS.dat, vit.dat) %>%
  mutate(norm.area = Area/IS.norm) %>%
  select(Compound, SampID, norm.area) %>%
  rename("Area" = norm.area)

#Add in IS normalized data back into dataset:
vit.dat.2 <- vit.dat %>%
  filter(!Compound %in% c("Vitamin B1", "Vitamin B2")) %>%
  rbind(., IS.norm)



###Run Quality Control 
Blk.dat <- vit.dat.2 %>%
  filter(str_detect(SampID, "DOM_blk")) %>%
  group_by(Compound) %>%
  summarize(mean.blk.area = mean(Area))

RSD.dat <- vit.dat %>%
  filter(str_detect(SampID, "Poo")) %>%
  group_by(Compound) %>%
  summarize(RSD.poo = sd(Area)/mean(Area))

##Add blk.qc flag, min.area.flag, RSD flag 
vit.qc <- left_join(vit.dat, Blk.dat) %>%
  mutate(min.area.flag = case_when(Area <= blk.ratio*mean.blk.area ~ T,
                                   TRUE ~ F)) %>% 
  mutate(blk.area.flag = case_when(Area <= min.area ~ T,
                                   TRUE ~ F)) %>%
  left_join(., RSD.dat) 

###########Perform min sample present comparison to only keep compounds that pass pass QC as being present in a 
# minimum number of triplicates (compound above min area and blk threshold in all three triplicates). Min number of triplicates is user defined. 
present.dat <- vit.qc %>%
  filter(min.area.flag == FALSE) %>%
  filter(blk.area.flag == FALSE) %>%
  filter(str_detect(SampID, "_1C") | str_detect(SampID, "_3LV") | str_detect(SampID, "6HV")) %>% 
  separate(SampID,  c("runDate", "type","samp","replicate", "Time"),"_", remove = FALSE) %>%
  group_by(Compound, replicate, Time) %>%
  summarize(count = n()) %>%
  filter(count == 3) %>%
  group_by(Compound) %>%
  summarize(numb.reps = n()) %>%
  filter(numb.reps > min.reps)

vit.qc.dat <- vit.qc %>%
  filter(Compound %in% present.dat$Compound) %>%
  select(Compound, SampID, Area) %>%
  mutate(SampID = str_replace(.$SampID, "230511_Smp_", "")) %>%
  filter(!str_detect(SampID, "230511_"))



# Pull in metadata and combine --------------------------------------------

######Pull in abundance and volume filtered dat
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

vol.abun.dat <- left_join(abund.dat, vol.dat) %>%
  select(-T)


#Join vit and meta data
vit.output <- left_join(vit.qc.dat, vol.abun.dat)

###Save to csv:
write_csv(vit.output, file = "Intermediates/Vit_TQS_QC2.dat")












# Quantification of Vitamins using Calibration Curve ----------------------

###Make calibration curves for CN-B12 and OH-B12
calib.dat <- vit.dat %>%
  filter(!str_detect(Compound, "Pseudo")) %>%
  filter(str_detect(SampID, "Std")) %>%
  separate(SampID,  c("runDate", "type","samp","run", "rep"),"_", remove = FALSE) %>%
  mutate(samp = str_remove(samp, "nM"),
         samp = as.numeric(str_replace(samp, "O", "0")))

###Visualize calibration curves:
ggplot(calib.dat, aes(x = samp, y = Area, color = rep)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(.~Compound, scales = "free")

###fit linear models to calibration curves and calculate slopes 
dat.cal.curve <- calib.dat %>%
  group_by(Compound)

calib.curves <- do(dat.cal.curve,
                             tidy(
                               lm(Area ~ samp, data = .)))

B12.cal.vals <- calib.curves %>%
  filter(term == "samp") %>%
  select(Compound, estimate) %>%
  rename("slope" = estimate) 

pseudo.cal.vals <- B12.cal.vals %>% 
  mutate(Compound = str_replace(Compound, "CN-B12", "CN-Pseudocob")) %>%
  mutate(Compound = str_replace(Compound, "OH-B12", "OH-Pseudocob")) %>%
  mutate(slope = slope/OH_B12.OH_pB12.ratio)

all.cal.vals <- rbind(B12.cal.vals, pseudo.cal.vals)


####Quantify compounds
quant.dat <- left_join(vit.qc.dat, all.cal.vals) %>%
  mutate(Conc_nM_vial = Area/slope) %>%
  mutate(SampID = str_replace(.$SampID, "230511_Smp_", "")) %>%
  filter(!str_detect(SampID, "230511_"))




###
viz.dat <- quant.dat %>%
  separate(SampID,  c("samp","replicate", "Time"),"_", remove = FALSE) %>%
  filter(replicate %in% c("1C", "3LV", "6HV"))

ggplot(viz.dat, aes(x = Time, y = Conc_nM_vial, fill = samp)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(Compound~replicate, scales = "free")




# Pull in metadata and combine --------------------------------------------

######Pull in abundance and volume filtered dat
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

vol.abun.dat <- left_join(abund.dat, vol.dat)

vit.output <- left_join(quant.dat, vol.abun.dat) %>%
  mutate(nM.in.smp = Conc_nM_vial*400e-6*(1/(Vol.filt.mL*10^-3))) %>%
  select(-T)

######
ggplot(vit.output, aes(x = Time, y = nM.in.smp, color = Compound)) +
  geom_point() +
  facet_wrap(treatment~biorep)


###Save to csv:
write_csv(vit.output, file = "Intermediates/Vit_TQS_QC2_Quant_Dat.csv")

