




library(readr)
library(tidyverse)



####Source Functions
source("R_Code/Functions.R")

#Define Inputs:


Sky.Pos.file2 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_File2.csv"
Sky.Neg.file2 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_File2.csv"
Sky.RP.file2 <- "Data_Raw/Skyline/ProMo_Particulate_RP_File2.csv"






Sky.Pos.area.2 <- sky_read(Sky.Pos.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound) %>%
  select(SampID, Name, Area) %>%
  mutate(fraction = "HPos",
         Column = "HILIC")

Sky.Neg.area.2 <- sky_read(Sky.Neg.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound) %>%
  select(SampID, Name, Area) %>%
  mutate(fraction = "HNeg",
         Column = "HILIC")

Sky.RP.area.2 <- sky_read(Sky.RP.file2) %>%
  filter(!str_detect(.$Rep, "DDA")) %>%
  rename("SampID" = Rep,
         "Name" = Compound) %>%
  select(SampID, Name, Area) %>%
  mutate(fraction = "RP",
         Column = "RP")




Sky.area <- rbind(Sky.Pos.area.2, Sky.Neg.area.2, Sky.RP.area.2) %>%
  filter(str_detect(.$SampID, "Std"))




#HILIC.data.file <- "Intermediates/ProMo_HILIC_QC1_output.csv"
Stds.info.file <- "Meta_Data/Ingalls_Lab_Standards_03172023.csv" 

####Load in standards
stds.dat <- Sky.area %>%
  mutate(z = case_when(fraction == "HNeg" ~ -1,
                       TRUE ~ 1)) %>%
  select(Name, SampID, z, Area, Column) %>%
  filter(!Name == "Unknown") %>%
  mutate(Mix = str_extract(SampID, "Mix\\d")) %>%
  mutate(Mix = case_when(Column == "RP" ~ "RP",
                         TRUE ~ Mix))


###Get Mix and Concentration info:
stds.info <- read_csv(Stds.info.file) %>%
  filter(Priority == TRUE) %>%
  select(Compound_Name, z, Column, HILIC_Mix, Concentration_uM) %>%
  rename("Name" = Compound_Name) %>%
  mutate(HILIC_Mix = case_when(Column == "RP" ~ "RP", TRUE ~ HILIC_Mix))

#HILIC.stds.info <- Stds.info %>%
 # filter(Column == "HILIC") %>%
 # select(Compound, z, HILICMix, Conc)

###Join stuff together + remove Matrix Samples
stds.dat.info <- left_join(stds.dat, stds.info) 

##Calculate RFs
RF.dat <- stds.dat.info %>%
  filter(Mix == HILIC_Mix) %>%
  select(-Mix, -HILIC_Mix) %>%
  filter(!str_detect(.$SampID, "Matrix")) %>%
  mutate(RF = as.numeric(Area)/Concentration_uM, NA) %>%
  group_by(Name, z) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF),
            RF = mean(RF, na.rm = TRUE))  %>%
  ungroup()

RFratio.dat <- stds.dat.info %>%
  filter(HILIC_Mix == Mix | is.na(Mix)) %>% 
  mutate(RunNumber = str_extract(SampID, "_\\d$")) %>%
  mutate(RunType = ifelse(str_detect(SampID, "StdsMix\\dInH2O")|
                            str_detect(SampID, "StdsInH2O"), "Std_in_h2O", 
                          ifelse(str_detect(SampID, "StdsMix\\dInMatrix") |
                                   str_detect(SampID, "StdsInMatrix"), "Std_in_matrix",
                                 "Matrix_in_h2O"))) %>%
  filter(HILIC_Mix == Mix | is.na(Mix)) %>% 
  select(-Mix, -HILIC_Mix, -SampID, -Concentration_uM) %>%
  spread(key = RunType, value = Area ) 

RF.ratios <- RFratio.dat %>%
  ungroup() %>%
  group_by(Name, RunNumber, z) %>%
  summarize("Std_in_matrix" = Std_in_matrix,
            "Matrix_in_h2o" = Matrix_in_h2O,
            "Std_in_h2o" = Std_in_h2O) %>%
  mutate(RFratio = (Std_in_matrix - Matrix_in_h2o)/ Std_in_h2o) %>%
  group_by(Name, z) %>%
  summarise(RFratio = mean(RFratio)) %>%
  ungroup()

###Join it all together
RF.RFratios <- left_join(RF.dat, RF.ratios)
write_csv(RF.RFratios, file = "Intermediates/Stds_RFs_RFratios.csv")
















