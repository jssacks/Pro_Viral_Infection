

#
#
library(tidyverse)
library(rstatix)
library(viridis)
library(lme4)
#
#
#Define inputs
metab.file <- "Intermediates/Final_Processed_Untargeted_Data.csv"
fcm.carbon.file <- "Intermediates/microbial_C_content.csv"
atp.file <- "Intermediates/pATP_nM_processed.csv"
#




# Compile and tidy datasets -----------------------------------------------

#metabolite data::
metab.dat <- read_csv(metab.file) %>%
  select(SampID, Vol.filt.mL, Time, biorep, treatment, treatment_number, MF, Name, Adjusted_Area, Pro, BacTot, BacSml, BacLrg, Phage)  %>%
  rename("time" = Time)

###ATP data:

#pull in ATP data
atp.dat <- read_csv(atp.file) %>%
  select(-sd.conc.nM) 

#combine ATP data with metadata for each sample from metab.dat
atp.meta.dat <- metab.dat %>%
  select(SampID, Vol.filt.mL, time, biorep, treatment, treatment_number, Pro, BacTot, BacSml, BacLrg, Phage) %>%
  unique() %>%
  left_join(., atp.dat %>%
              select(time, biorep, treatment, mean.conc.nM)) %>%
  mutate(MF = "ATP",
         Name = "ATP") %>%
  mutate("Adjusted_Area" = mean.conc.nM*Vol.filt.mL) %>%
  select(-mean.conc.nM)


##Combine Metabolite and ATP Data and add in carbon data:
metab.atp.dat <- rbind(metab.dat, atp.meta.dat) %>%
  left_join(., read_csv(fcm.carbon.file) %>% rename("biorep" = experiment) %>%
              select(biorep, treatment, time, Pro.C, bacterial.C, total.C))






# Fold change analysis normalized to microbial biomass data  ---------------------
#####

#calculate fold change compared to control (normalized to volume filtered and biomass)
FC.Carbon.dat <- metab.atp.dat %>%
  mutate(vol.norm.area = (Adjusted_Area+1)/Vol.filt.mL,
         bio.norm.area = vol.norm.area/total.C) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  group_by(MF, time, biorep) %>%
  mutate(C.norm.area = bio.norm.area/(bio.norm.area[treatment == "C"]),
         log2.fc = log2(C.norm.area)) %>%
  ungroup() %>%
  filter(!SampID == "C_3LV_T48")


#prepare data for t-test
FC.Carbon.test <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_36", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, biorep, SampID, C.norm.area, log2.fc) %>%
  group_by(MF) %>%
  mutate(mean.log2.fc = mean(log2.fc)) %>%
  ungroup() 

#run t-test on data
FC.Carbon.results <- FC.Carbon.test %>%
  group_by(MF) %>%
  t_test(log2.fc ~1, mu = 0, p.adjust.method = "fdr") %>%
  select(MF, p)


###run Wilcoxon Sign test
FC.nonparametric <- FC.Carbon.test %>%
  group_by(MF) %>%
  wilcox_test(log2.fc~1, mu = 0, p.adjust.method = "fdr") %>%
  select(MF, p) %>%
  rename("wilcox_p" = p)


#combine FC data with t-test results
FC.Carbon.out <- left_join(FC.Carbon.test, FC.Carbon.results) %>%
  left_join(., FC.nonparametric) %>%
  mutate(signif = case_when(p <= 0.05 & wilcox_p <= 0.05 ~ TRUE,
                            TRUE ~ FALSE)) %>%
  unique() 

#####################MAKE SURE WILCOX TEST ANALYSIS IS INCORPORATED!!!!!!!!!!!!!!!




# Correlation analysis comparing bacterial carbon and metabolite abundance using linear models---------------------
#####
bact.corr.dat <- metab.atp.dat %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  mutate(vol.norm.area = (Adjusted_Area+1)/Vol.filt.mL) %>%
  group_by(MF) %>%
  mutate(bact.cor.pearson = cor(vol.norm.area, BacTot, method = "pearson"))

###data for linear model
lin.mod.dat <- bact.corr.dat %>%
  group_by(MF)

####Run linear model
lin.mod.out <- do(lin.mod.dat,
                  tidy(
                    lm(vol.norm.area ~ BacTot, data = .))) %>%
  filter(term == "BacTot") 

#adjust p-values using fdr correction
lin.mod.out$p.adj <- p.adjust(lin.mod.out$p.value, method = "fdr")

lin.mod.out <- lin.mod.out %>%
  select(MF, p.adj) %>%
  rename("bactlm.p.adj" = p.adj)

####combine correlation values with p-values
bact.corr.test.out <- bact.corr.dat %>%
  left_join(., lin.mod.out)






# Identify differentially abundant metabolites:

###Identification of confidently altered metabolites due to viral infection
fc.bcorr.dat <- left_join(FC.Carbon.out, bact.corr.test.out) %>%
  select(MF, Name, mean.log2.fc, p, wilcox_p, signif, bact.cor.pearson, bactlm.p.adj) %>%
  unique() 

#identify bact correlated MFs
bact.corr.mfs <- fc.bcorr.dat %>%
  filter(bact.cor.pearson > 0) %>%
  filter(bactlm.p.adj < 0.05)


#identify MFs meeting FC requirements (0.5 threshold)
fc.mfs <- fc.bcorr.dat %>%
  filter(mean.log2.fc >= 0.5 | mean.log2.fc <= -0.5) %>%
  filter(p <= 0.05)


#identify MFs signficantly and positively correlated with bacterial abundance. 

final.fc.mf.list <- fc.bcorr.dat %>%
  filter(MF %in% fc.mfs$MF) %>%
  filter(!MF %in% bact.corr.mfs$MF) %>%
  left_join(read_csv(metab.file) %>% select(MF, SIRIUS_molecularFormula, SIRIUS_name) %>% unique())


#save results:

#just virally altered MFs
write_csv(final.fc.mf.list, file =  "Intermediates/virally_altered_mfs.csv")

#All Mfs
write_csv(fc.bcorr.dat, file = "Intermediates/fc_analysis_all_mfs.csv")






