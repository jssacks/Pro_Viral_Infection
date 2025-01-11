


# This script has 3 sections:
# 1. Calculate Total Microbial Biomass and % Pro C for each sample estimated from flow cytometry
# 2. 





library(tidyverse)

##data files

#FCM
c.pro.file <- "Collaborator_Data/FCM/Prochloro_size_abundance_rates_C_content.csv"
c.bact.file <- "Collaborator_Data/FCM/Bacteria2pop_size_abundance_C_content.csv"

#Metabolomics
quant.metab.file <- "Intermediates/Final_Processed_Quantified_Data.csv"
metab.meta.file <- "Intermediates/Final_Processed_Untargeted_Data.csv"

#ATP
atp.file <- "Intermediates/pATP_nM_processed.csv"







# 1. Calculate Total Microbial Biomass and % Pro C for each sample -----------

##### organize Pro C-content data
pro.dat <- read_csv(c.pro.file) %>%
  mutate(Pro.C = abundance*Qc) %>%
  filter(treatment %in% c("Control", "High virus", "Low virus")) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV")) %>%
  mutate("pop" = "Pro") %>%
  select(-net, -loss, -pop) %>%
  rename("Pro.abu" = abundance,
         "Pro.Qc" = Qc)

####Organize bacteria C-content data 
bact.dat <- read_csv(c.bact.file) %>%
  filter(treatment %in% c("Control", "LV", "HV")) %>%
  mutate(Pop.C = abundance*Qc*1000)  %>%
  mutate(abundance = abundance*1000) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               TRUE ~ treatment)) %>%
  select(-series) 

bact1.dat <- bact.dat %>%
  filter(pop == "bacteria") %>%
  rename("B1.abu" = abundance,
         "B1.Qc" = Qc,
         "B1.C" = Pop.C) %>%
  select(-pop)

bact2.dat <- bact.dat %>%
  filter(pop == "bacteria2") %>%
  rename("B2.abu" = abundance,
         "B2.Qc" = Qc,
         "B2.C" = Pop.C) %>%
  select(-pop)

bact.dat.all <- left_join(bact1.dat, bact2.dat)

####Combine Pro and bacteria C-content data
all.C.dat <- left_join(pro.dat, bact.dat.all) %>%
  mutate(bacterial.C = B1.C+B2.C,
         total.C = bacterial.C + Pro.C) %>%
  mutate(Pro.C.perc = Pro.C/(total.C)) 


#export
write_csv(all.C.dat, file = "Intermediates/microbial_C_content.csv")



### 2. Calculate % of metabolite carbon in XYZ  _______________


#load in metabolite data:

#get meta data from untargeted data file
metab.meta.dat <- read_csv(metab.meta.file) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  select(SampID, Time, biorep, treatment) %>%
  rename("time" = Time) %>%
  unique()

#combine meta data with targeted data 
quant.metab.dat <- read_csv(quant.metab.file) %>%
  left_join(., metab.meta.dat) %>%
  filter(!Name == "Succinic acid") %>%
  filter(!Name == "Vit_OH-pB12") %>%
  filter(!Name == "Taurine, 2H4") %>%
  filter(!SampID == "C_3LV_T48") %>%
  filter(!is.na(nM_C))

#Calculate total metabolite C for each sample:
metab.tot.C <- quant.metab.dat %>%
  select(SampID, time, biorep, treatment, nM_C) %>%
  group_by(SampID, time, biorep, treatment) %>%
  reframe(tot_metab_C_nM = sum(nM_C))


#Combine metabolite C data and FCM C data 
metab.c.comparison <- metab.tot.C %>%
  left_join(., all.C.dat %>% rename("biorep" = experiment)) %>%
  select(SampID, time, biorep, treatment, tot_metab_C_nM, total.C, Pro.C.perc) %>%
  mutate(tot.C.nM = total.C*(1/12.011)*(1e-6)*1000,
         tot_metab_C_nM/tot.C.nM)

#export comparison data:
write_csv(metab.c.comparison, file = "Intermediates/Metabolite_Carbon_Comparison.csv")







### 3. Comparison of C estimated by ATP and FCM and total metabolite carbon _______________

#load in data
atp.dat <- read_csv(atp.file)

##combine atp + metab + c dat
c.dat <- left_join(atp.dat, metab.c.comparison) %>%
  mutate(microb.C.fg.mL = total.C) %>%
  mutate(metab.tot.C.fg.mL = tot_metab_C_nM*12.011/1000/(1e-6)) %>%
  mutate(metab.tot.C.fg.mL.10 = tot_metab_C_nM*12.011/1000/(1e-6)*10)

###calculate correlation coefficients between variables:
c.dat.cor <- c.dat %>%
  ungroup() %>%
  reframe(atp.fcm.cor = cor(mean.live.C.fg.mL, microb.C.fg.mL, use = "pairwise.complete.obs"),
            atp.fcm.r2 = atp.fcm.cor^2,
            atp.metab.cor = cor(mean.live.C.fg.mL, metab.tot.C.fg.mL.10, use = "pairwise.complete.obs"),
            atp.metab.r2 = atp.metab.cor^2,
            fcm.metab.cor = cor(microb.C.fg.mL, metab.tot.C.fg.mL.10, use = "pairwise.complete.obs"),
            fcm.metab.r2 = fcm.metab.cor^2)






























































