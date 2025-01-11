



library(tidyverse)



# Define inputs -----------------------------------------------------------

#QC2 File:
qc2.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"    

#Adduct Files:
adduct.remove.file <- "Intermediates/QC2_Adducts_to_Remove.csv"
adduct.id.file <- "Intermediates/QC2_Adducts_to_Annotate.csv"

#Annotation Files:
annotation.file <- "Intermediates/ProMo_Particulate_SIRIUS_MF_Annotations.csv"

#Quantification Files:
quant.file <- "Intermediates/ProMo_Particulate_Quant_Output.csv"

#Vitamin Files:
vit.file <- "Intermediates/Vit_TQS_QC2_Quant_Dat.csv"




# Create final untargeted dataset -----------------------------------------

#load in datasets:
qc2.dat <- read_csv(qc2.file)
adduct.remove.dat <- read_csv(adduct.remove.file)
adduct.id.dat <- read_csv(adduct.id.file)
annotation.dat <- read_csv(annotation.file)
vit.dat <- read_csv(vit.file)


#combine all datasets
untargeted.data <- qc2.dat %>%
  filter(!MF %in% adduct.remove.dat$MF) %>%
  left_join(., adduct.id.dat %>% rename("Predicted_Adduct_Ion" = Ion)) %>%
  left_join(., annotation.dat %>% rename("SIRIUS_molecularFormula" = molecularFormula,
                                         "SIRIUS_name" = name) %>%
                                  select(MF, SIRIUS_molecularFormula, SIRIUS_name)) %>%
  full_join(., vit.dat %>%
              rename("MF" = Compound) %>% 
              select(MF, SampID, Area, time, Time, biorep, treatment, treatment_number,
                                                 Pro, BacTot, BacSml, BacLrg, Grz, Phage, Vol.filt.mL) %>%
              mutate(Adjusted_Area = Area,
                     Name = MF,
                     fraction = "Vit")
              )



#summarize statistics of final dataset
untargeted.sum <- untargeted.data %>%
  mutate(known = case_when(Name == "Unknown" ~ "no",
                           TRUE ~ "yes")) %>%
  select(MF, fraction, known) %>%
  unique() %>%
  group_by(fraction, known) %>%
  reframe(count = n())


#export final dataset
write_csv(untargeted.data, file = "Intermediates/Final_Processed_Untargeted_Data.csv")


# Create final quantified targeted dataset --------------------------------
quant.dat <- read_csv(quant.file) %>%
  full_join(., vit.dat %>% select(Compound, SampID, Vol.filt.mL, nM.in.smp) %>%
              mutate(MF = Compound, Name = Compound) %>% select(-Compound)) 

#export final dataset:
write_csv(quant.dat, file = "Intermediates/Final_Processed_Quantified_Data.csv")















































































