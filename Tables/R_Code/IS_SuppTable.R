





library(tidyverse)



#define inputs:
#Stds and RFs and RFratios
stds.file <- "Meta_Data/Ingalls_Lab_Standards_03172023.csv"
hilic.is.file <- "Intermediates/ProMo_Particulate_HILIC_ISdat.csv"
rp.is.file <- "Intermediates/ProMo_Particulate_RP_ISdat.csv"








#Make list of internal standards used
hilic.is.dat <- read_csv(hilic.is.file) %>%
  mutate(Fraction = "HILIC")
rp.is.dat <- read_csv(rp.is.file) %>%
  mutate(Fraction = "RP")

full.is.dat <- rbind(hilic.is.dat, rp.is.dat) %>%
  rename("IS" = MF,
         "IS.area" = Area) %>%
  unique()

full.is.num <- full.is.dat %>%
  select(IS, Fraction) %>%
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


#compbine and make table
is.table <- left_join(is.std, full.is.num)

#export table
write_csv(is.table, file = "Tables/Outputs/IS_Supplemental_Table.csv")









































