

#
#
library(tidyverse)
library(rstatix)
library(viridis)
#library(lme4)
#
#
#Define inputs
metab.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"
vit.file <- "Intermediates/Vit_TQS_QC2.dat"
pro.carbon.file <- "Meta_Data/Abundance_Dat/Prochloro_size_abundance_rates_C_content.csv"
bact.file <- "Meta_Data/Abundance_Dat/Bacteria2pop_size_abundance_C_content.csv"
atp.file <- "Intermediates/pATP_nM_processed.csv"
annotation.file <- "Intermediates/ProMo_Particulate_SIRIUS_MF_Annotations.csv"
adduct.file <- "Intermediates/QC2_Adducts_to_Remove.csv"
#




# Compile and tidy datasets -----------------------------------------------

#metabolite data::
metab.dat <- read_csv(metab.file) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  select(SampID, Vol.filt.mL, Time, biorep, treatment, treatment_number, MF, Name, Adjusted_Area, Pro, BacTot, BacSml, BacLrg, Phage) %>%
  rename("time" = Time) %>%
  filter(!str_detect(Name, "Vit_"))

vit.dat <- read_csv(vit.file) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  mutate("Adjusted_Area" = Area,
         "MF" = Compound,
         "Name" = Compound) %>%
  select(SampID, Vol.filt.mL, Time, biorep, treatment, treatment_number, MF, Name, Adjusted_Area, Pro, BacTot, BacSml, BacLrg, Phage) %>%
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
  # select(SampID, time, biorep, treatment, mean.conc.nM) %>%
  mutate(MF = "ATP",
         Name = "ATP") %>%
  # rename("vol.norm.area" = mean.conc.nM) %>%
  mutate("Adjusted_Area" = mean.conc.nM*Vol.filt.mL) %>%
  select(-mean.conc.nM)

metab.dat.mz <- read_csv(metab.file) %>%
  select(MF, Name, mz, RT) %>%
  unique()


#metabs to remove:
adduct.remove <- read_csv(adduct.file)



##Combine Metabolite, Vitamin, and ATP Data
metab.vit.dat <- rbind(metab.dat, vit.dat, atp.meta.dat) %>%
  filter(!MF %in% adduct.remove$MF)



###Carbon Data 
#####Pro Carbon
pro.carbon.dat <- read_csv(pro.carbon.file) %>%
  mutate(Pro.Carbon = abundance*Qc) %>%
  filter(!experiment == "E") %>%
  filter(treatment %in% c("Control", "High virus", "Low virus")) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV"
  )) %>%
  rename("biorep" = experiment,
         "Pro.abundance" = abundance,
         "Pro.Qc" = Qc)

###Bacteria Carbon

#Bacteria Population 1 (Small)
bact1.carbon.dat <- read_csv(bact.file) %>%
  filter(pop == "bacteria") %>%
  filter(!experiment == "E") %>%
  filter(treatment %in% c("Control", "HV", "LV")) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV",
                               TRUE ~ treatment
  )) %>%
  rename("biorep" = experiment,
         "Bact1.Qc" = Qc) %>%
  mutate("Bact1.abundance" = abundance*1000,
         "Bact1.Carbon" = Bact1.abundance*Bact1.Qc) %>%
  select(-series, -pop, - abundance)

#Bacteria Population 2 (Large)
bact2.carbon.dat <- read_csv(bact.file) %>%
  filter(pop == "bacteria2") %>%
  filter(!experiment == "E") %>%
  filter(treatment %in% c("Control", "HV", "LV")) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV",
                               TRUE ~ treatment
  )) %>%
  rename("biorep" = experiment,
         "Bact2.Qc" = Qc) %>%
  mutate("Bact2.abundance" = abundance*1000,
         "Bact2.Carbon" = Bact2.abundance*Bact2.Qc) %>%
  select(-series, -pop, -abundance)

#merge carbon datasets
pro.bact1.dat <- left_join(pro.carbon.dat, bact1.carbon.dat)

full.carbon.dat <- left_join(pro.bact1.dat, bact2.carbon.dat) %>%
  mutate(Pro.C.nM = Pro.Carbon*(1/12.011)*(1e-6)*1000,
         Bact.C.tot.nM = (Bact1.Carbon+Bact2.Carbon)*(1/12.011)*(1e-6)*1000,
         Tot.C.nM = Pro.C.nM + Bact.C.tot.nM)





###combine carbon and metabolite data
metab.meta.dat <- left_join(metab.vit.dat, full.carbon.dat) 




# Fold change analysis normalized to Microbial Carbon Biomass  ---------------------
#####

#calculate fold change
FC.Carbon.dat <- metab.meta.dat%>%
  mutate(vol.norm.area = (Adjusted_Area+1)/Vol.filt.mL,
         bio.norm.area = vol.norm.area/Pro.C.nM) %>%
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


###ATP FC dat
# ATP.Carbon.dat <- atp.meta.dat %>%
#   left_join(., full.carbon.dat) %>%
#   mutate(pro.norm.area = Adjusted_Area/Tot.C.nM) %>%
#   group_by(MF, time, biorep) %>%
#   mutate(C.norm.area = pro.norm.area/(bio.norm.area[treatment == "C"]),
#          log2.fc = log2(C.norm.area)) %>%
#   ungroup() %>%
#   filter(!SampID == "C_3LV_T48") %>%
#   filter(!treatment == "C") %>%
#   unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
#   filter(treatment_time %in% c("LV_36", "HV_12", "HV_24")) %>%
#   select(MF, Name, treatment, time, biorep, SampID, C.norm.area, log2.fc) %>%
#   group_by(MF) %>%
#   mutate(mean.log2.fc = mean(log2.fc)) %>%
#   ungroup() 


#run t-test on data
FC.Carbon.results <- FC.Carbon.test %>%
#  rbind(., ATP.Carbon.dat) %>%
  group_by(MF) %>%
  t_test(log2.fc ~1, mu = 0, p.adjust.method = "fdr") %>%
  select(MF, p)


###Try Wilcoxon Sign test
FC.nonparametric <- FC.Carbon.test %>%
#  rbind(., ATP.Carbon.dat) %>%
  #  filter(!MF %in% bact.corr.mfs$MF) %>%
  group_by(MF) %>%
  wilcox_test(log2.fc~1, mu = 0, p.adjust.method = "fdr") %>%
  select(MF, p) %>%
  rename("wilcox_p" = p)


#combine FC data with t-test results
FC.Carbon.out <- left_join(FC.Carbon.test, FC.Carbon.results) %>%
  left_join(., FC.nonparametric) %>%
  mutate(signif = case_when(p <= 0.05 ~ TRUE,
                            TRUE ~ FALSE)) %>%
  unique() 







# Correlation analysis comparing bacterial carbon and metabolite abundance using linear models---------------------
#####
bact.corr.dat <- metab.meta.dat %>%
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



###Identification of confidently altered metabolites due to viral infection
fc.bcorr.dat <- left_join(FC.Carbon.out, bact.corr.test.out) %>%
  select(MF, Name, mean.log2.fc, p, wilcox_p, signif, bact.cor.pearson, bactlm.p.adj) %>%
  unique() 

#identify MFs meeting FC requirements
# 0.25
fc.mfs <- fc.bcorr.dat %>%
  filter(mean.log2.fc >= log2(1.25) | mean.log2.fc <= log2(0.75)) %>%
  filter(p <= 0.05) %>%
  filter(wilcox_p <= 0.05)


#identify MFs signficantly and positively correlated with bacterial abundance. 
bact.corr.mfs <- fc.bcorr.dat %>%
  filter(bact.cor.pearson > 0) %>%
  filter(bactlm.p.adj < 0.05)

bcm.sum <- bact.corr.mfs %>%
  filter(MF %in% fc.mfs$MF) %>%
  group_by(Name) %>%
  summarize(count = n())

final.results <- fc.bcorr.dat %>%
  filter(MF %in% fc.mfs$MF) %>%
  filter(!MF %in% bact.corr.mfs$MF)

###write to csv
write_csv(final.results, file = "Intermediates/virally_altered_mfs.csv")

write_csv(fc.bcorr.dat, file = "Intermediates/fc_analysis_all_mfs.csv")




#########

comp.name <- final.results %>%
  filter(!Name == "Unknown")

ggplot(final.results, aes(x = mean.log2.fc, y = -log10(p))) +
  geom_point(data = FC.Carbon.out, aes(x = mean.log2.fc, y = -log10(p))) +
  geom_point(aes(color = mean.log2.fc)) +
  scale_color_viridis() +
  geom_vline(xintercept = log2(1.25)) +
  geom_vline(xintercept = log2(0.75)) +
  geom_hline(yintercept = -log10(0.05)) + 
  geom_label(data = comp.name, aes(x = mean.log2.fc, y = -log10(p), label = Name))


final.results.sum <- final.results %>%
  group_by()
  

