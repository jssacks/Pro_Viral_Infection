


library(tidyverse)



#define inputs:
mpm.file <- "Collaborator_Data/MPM/MPM_results_final.csv"





##import data:
mpm.dat <- read_csv(mpm.file) 

##change data to factors:
mpm.dat.factor <- mpm.dat %>%
  mutate(treatment = as.factor(treatment),
         time = as.factor(time),
         replicate = as.factor(replicate)) %>%
  mutate(treatment = fct_relevel(treatment, c("Control", "Low virus", "High virus")))

##Make into supplemental table
mpm.supp.table <- mpm.dat.factor %>%
  rename("Treatment" = treatment,
         "Rep" = replicate,
         "Time" = time,
         "Mean_Division_Rate_per_day" = division,
         "SD_Division_Rate_per_day" = division_sd,
         "Mean_Carbon_Fixation_Rate_fg_C_per_day" = carbon_fixation,
         "SD_Carbon_Fixation_Rate_fg_C_per_day" = fixation_sd,
         "Mean_Carbon_Loss_Rate_fg_C_per_day" = carbon_loss,
         "SD_Carbon_Loss_Rate_fg_C_per_day" = loss_sd) %>%
  select(-experiment)

##export supplemental table:
write_csv(mpm.supp.table, file = "Tables/Outputs/MPM_Supplemental_Table.csv")










































