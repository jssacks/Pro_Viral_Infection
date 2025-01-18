

###library
library(tidyverse)


#define inputs
pATP.file <- "Collaborator_Data/ATP/pATP_data.csv"
sample_meta_data_file <- "Intermediates/ProMo_Particulate_QC2_data.csv"

#define values
ATP.molar.mass <- 507.18
 


## Load in data, convert ATP to nM, and convert ATP to biomass
pATP.dat <- read_csv(pATP.file)

meta.dat <- read_csv(sample_meta_data_file) %>%
  select(SampID, Time, biorep, treatment, Vol.filt.mL) %>%
  rename("time" = Time) %>%
  unique()

##add in sample meta data and convert ATP in ug/L to ATP in nM/L
pATP.processed <- pATP.dat %>%
  left_join(., meta.dat) %>%
  mutate(nM.in.smp = PATP*1000/ATP.molar.mass,
         live.biomass.ug.L = PATP*250,
         live.biomass.fg.mL = PATP*250*(1E9)/(1E3)) 

#visualize data
ggplot(pATP.processed, aes(x = time, y = nM.in.smp)) + 
  geom_point() +
  geom_line() +
  facet_grid(treatment~biorep)


## Summarize ATP data into means and standard deviations_____________________
# Remove as they appear to be fliers:
# B_1C_T12_pseudorep2
# B_6HV_T48_pseudorep3

#remove fliers
pATP.keep <- pATP.processed %>%
  select(pseudo_rep, SampID, everything()) %>%
  unite("pseudo_SampID", pseudo_rep:SampID, remove = FALSE) %>%
  filter(!pseudo_SampID %in% c("2_B_1C_T12", "3_B_6HV_T48")) 


#summarize
pATP.sum <- pATP.keep %>%
  group_by(treatment, biorep, time) %>%
  summarize(mean.conc.nM = mean(nM.in.smp),
            sd.conc.nM = sd(nM.in.smp),
            mean.live.C.ug.L = mean(live.biomass.ug.L),
            sd.live.C.ug.L = sd(live.biomass.ug.L),
            mean.live.C.fg.mL = mean(live.biomass.fg.mL),
            sd.live.C.fg.mL = sd(live.biomass.fg.mL))

#write to csv
write_csv(pATP.sum, file = "Intermediates/pATP_nM_processed.csv")


#make supplemental table:
atp.sup <- pATP.sum %>%
  rename("Treatment" = treatment,
         "Time" = time,
         "Rep" = biorep,
         "Mean_Particulate_ATP_Conc_nM" = mean.conc.nM,
         "SD_Particulate_ATP_Conc_nM" = sd.conc.nM,
         "Mean_ATP_Estimated_Living_Biomass_fg_mL" = mean.live.C.fg.mL,
         "SD_ATP_Estimated_Living_Biomass_fg_mL" = sd.live.C.fg.mL) %>%
  select(Treatment, Time, Rep, Mean_Particulate_ATP_Conc_nM, SD_Particulate_ATP_Conc_nM,
         Mean_ATP_Estimated_Living_Biomass_fg_mL, SD_ATP_Estimated_Living_Biomass_fg_mL)

#save supplemental table
write_csv(atp.sup, file = "Tables/Outputs/ATP_Supplemental_Table.csv")







