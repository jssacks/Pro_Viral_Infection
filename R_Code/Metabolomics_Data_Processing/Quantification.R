


#####
library(readr)
library(tidyverse)

#Inputs

#Stds and RFs and RFratios
stds.file <- "Meta_Data/Ingalls_Lab_Standards_03172023.csv"
rf.file <- "Intermediates/Stds_RFs_RFratios.csv"


#IS
hilic.is.file <- "Intermediates/ProMo_Particulate_HILIC_ISdat.csv"
rp.is.file <- "Intermediates/ProMo_Particulate_RP_ISdat.csv"


### MF dat
qc2.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"
hilic.qc1.file <- "Intermediates/ProMo_HILIC_QC1_output.csv"
rp.qc1.file <- "Intermediates/ProMo_RP_QC1_output.csv"


###QC check for MFs to be 3x blk values in Samp (quantifiable)
blk.qc.file <- "Intermediates/QC2_MFSamps_forQuant.csv"


####
stds.not.in.mix <- c("Choline sulfate", "Pipecolic acid", "Taurocyamine", "Glycine")



#Quantify values using standards 
rf.dat <- read_csv(rf.file)


##Calculate concentration in Vial using normal RF and RFratio approach, 
# also remove compounds that did not have a standard run in the stds mix for this run or that are very poor ionizers
vial.quant.dat <- read_csv(qc2.file) %>%
  filter(!Name == "Unknown") %>%
  left_join(., rf.dat) %>%
  select(MF, Name, mz, RT, fraction, z, SampID, Vol.filt.mL, Adjusted_Area, RF, RFmax, RFmin, RFratio) %>%
  mutate(umol.in.vial.ave = Adjusted_Area/RF/RFratio,
         umol.in.vial.max = Adjusted_Area/RFmin/RFratio,
         umol.in.vial.min = Adjusted_Area/RFmax/RFratio) %>%
  filter(!str_detect(.$SampID, "Blk")) %>%
  filter(!Name %in% stds.not.in.mix)

v.names <- vial.quant.dat %>%
  select(Name) %>%
  unique()

#Quantify compounds with matched IS
hilic.is.dat <- read_csv(hilic.is.file)
rp.is.dat <- read_csv(rp.is.file)

full.is.dat <- rbind(hilic.is.dat, rp.is.dat) %>%
  rename("IS" = MF,
         "IS.area" = Area) %>%
  unique()

full.is.num <- full.is.dat %>%
  select(IS) %>%
  unique()


#get internal standard concentration values 
is.std <- read_csv(stds.file) %>%
  mutate(Compound_Name = str_replace(Compound_Name, "Sucrose, 13C12", "Sucrose, 13C")) %>%
  mutate(Compound_Name = str_replace(Compound_Name, "Trehalose, 13C12", "Trehalose, 13C")) %>%
  mutate(Compound_Name = str_replace(Compound_Name, "Adenosine monophosphate, 15N5", "AMP, 15N5")) %>%
  mutate(Compound_Name = str_replace(Compound_Name, "Guanosine monophosphate, 15N5", "GMP, 15N5")) %>%
  mutate(Compound_Name = str_replace(Compound_Name, "Riboflavin-dioxopyrimidine, 13C4, 15N2", "Vitamin B2, 13C4, 15N2")) %>%
  filter(Compound_Name %in% full.is.dat$IS) %>%
  select(Compound_Name, Concentration_uM) %>%
  rename("IS" = Compound_Name,
         "IS.Conc.uM" = Concentration_uM) 


####Pull in QC1 data 
hilic.qc1.dat <- read_csv(hilic.qc1.file)
rp.qc1.dat <- read_csv(rp.qc1.file)

qc1.dat.full <- rbind(hilic.qc1.dat, rp.qc1.dat) %>%
  rename("MF" = ID) %>%
  filter(!Name == "Unknown") %>%
  select(MF, Name, SampID, Area) %>%
  full_join(., full.is.dat)  %>%
  mutate(match.name = str_extract(.$IS, Name)) %>%
  filter(!is.na(match.name)) %>%
  filter(!Name == "Glycine") %>%
  filter(MF %in% vial.quant.dat$MF) %>% 
  mutate(SampID = str_replace(.$SampID, "201006_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "201013_", "")) %>%
  mutate(SampID = str_replace(.$SampID, "200929_", "")) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  mutate(SampID = str_replace(.$SampID, "Smp_", "")) %>%
  filter(!str_detect(.$SampID, "Blk"))  %>%
  unique()

vial.is.quant.dat <- left_join(qc1.dat.full, is.std) %>%
  mutate(umol.in.vial.ave = Area*IS.Conc.uM/IS.area)

is.mfs <- vial.is.quant.dat %>%
  select(Name) %>%
  unique()

vq.mfs <- vial.quant.dat %>%
  select(Name) %>%
  unique()

same <- inner_join(is.mfs, vq.mfs)


#####Combine normal and IS quantification data
vol.filt <- read_csv(qc2.file) %>%
  select(SampID, Vol.filt.mL) %>%
  unique()

smp.quant.dat.all <- vial.quant.dat %>%
  select(MF, Name, SampID, umol.in.vial.ave) %>%
  filter(!MF %in% vial.is.quant.dat$MF) %>%
  rbind(., vial.is.quant.dat %>% 
          select(MF, Name, SampID, umol.in.vial.ave)) %>%
  filter(!Name == "Taurine, 2H4") %>%
  left_join(., vol.filt) %>%
  mutate(nM.in.smp = umol.in.vial.ave*10^-6*400/(Vol.filt.mL*10^-3)*1000)%>%
  unique()
  



#Calculate nM C and N per sample
std.formula <- read_csv(stds.file) %>%
  select(Compound_Name, Empirical_Formula) %>%
  rename("Name" = Compound_Name) %>%
  unique() %>%
  mutate(C = ifelse(is.na(str_extract(Empirical_Formula, "^C\\d\\d")),
                    str_extract(Empirical_Formula, "^C\\d"), 
                    str_extract(Empirical_Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Empirical_Formula, "N\\D"),
                    1, str_extract(Empirical_Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", "")))


######keep only compounds in each sample that passed blank QC 
blk.qc.dat <- read_csv(blk.qc.file)
smp.quant.dat.bqc <- left_join(blk.qc.dat, smp.quant.dat.all) %>%
  filter(!Name %in% stds.not.in.mix)


###Export quantified dat
samp.quant.dat <- left_join(smp.quant.dat.bqc, std.formula) %>%
  unique() %>%
  mutate(nM_C = nM.in.smp*C,
         nM_N = nM.in.smp*N) %>%
  select(MF, Name, SampID, Vol.filt.mL, umol.in.vial.ave, nM.in.smp, nM_C, nM_N)

write_csv(samp.quant.dat, file = "Intermediates/ProMo_Particulate_Quant_Output.csv")















#
