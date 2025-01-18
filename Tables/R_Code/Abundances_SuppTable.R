
library(tidyverse)



#Define inputs:
abu.file <- "Collaborator_Data/FCM/ProMo_Abundat_Combined.csv"
microb.c.file <- "Intermediates/microbial_C_content.csv"





#Abundances Supplemental Table
abu.dat <- read_csv(abu.file)  %>%
  separate(col = "Rep", into = c("Rep", "Treat.number")) %>%
  mutate(Treatment = case_when(Rep %in% c("A", "B", "C") & Treat.number == 1 ~ "C",
                               Rep %in% c("A", "B", "C") & Treat.number == 3 ~ "LV",
                               Rep %in% c("A", "B", "C") & Treat.number == 4 ~ "LGV", 
                               Rep %in% c("A", "B", "C") & Treat.number == 6 ~ "HV",
                               Rep %in% c("E") ~ "LG")) %>%
  filter(Treatment %in% c("C", "LV", "HV")) %>%
  mutate(Treatment = fct_relevel(Treatment, c("C", "LV", "HV"))) %>%
  select(Treatment, Time, Rep, Pro, BacSml, BacLrg, BacTot, Phage) %>%
  rename("Pro.count" = Pro,
         "BacSml.count" = BacSml,
         "BacLrg.count" = BacLrg,
         "BacTot.count" = BacTot)



c.dat <- read_csv(microb.c.file) %>%
  select(experiment, treatment, time, Pro.Qc, Pro.C, B1.Qc, B1.C, B2.Qc, B2.C, bacterial.C, total.C, Pro.C.perc) %>%
  rename("Rep" = experiment,
         "Treatment" = treatment,
         "Time" = time,
         "BacSml.Qc" = B1.Qc,
         "BacSml.C" = B1.C,
         "BacLrg.Qc" = B2.Qc,
         "BacLrg.C" = B2.C,
         "BacTot.C" = bacterial.C,
         "TotMicrobial.C" = total.C,
         "Pro.C.Fraction" = Pro.C.perc)


abu.c.dat <- left_join(abu.dat, c.dat) %>%
  select(Time, Treatment, Rep, Pro.count, Pro.Qc, Pro.C, BacSml.count, BacSml.Qc, BacSml.C, BacLrg.count, BacLrg.Qc, BacLrg.C, BacTot.count, BacTot.C, TotMicrobial.C, Pro.C.Fraction)



#Write Abundances Supplemental Table:
write_csv(abu.c.dat, file = "Tables/Outputs/FCM_Abundances_Supplemental_Table.csv")























































