






library(tidyverse)
library(readr)

source("R_Code/Functions.R")

#define inputs


#QC2 MF Output      
MF.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"    

##Adduct list
adduct.file <- "Meta_Data/FiehnLab_Adduct_List_modified.csv"



##Define adduct detection parameters
adduct.error.ppm <- 5   ###PPM Error tolerance
RT_ad_tol <- 0.1   #Retention time tolerance
adduct.cor.threshold <- 0.95  #Adduct correlation threshold






#________________Adduct Section______________________
###run adduct/isotope/etc. detection function
adduct.list <- read_csv(adduct.file) %>%
  rename("ad.polarity" = polarity)

ID.dat <- read_csv(MF.file) %>%
  select(MF, mz, RT, fraction) %>%
  rename("ID" = MF) %>%
  unique() %>%
  mutate(polarity = case_when(str_detect(fraction, "HILIC_Neg") ~ -1,
                              TRUE ~ 1)) 

#%>%
#  select(ID, mz, RT, polarity) %>%
#  unique() 


adducts.dat <- tibble(
  adduct.ID = as.character(),
  RT = as.numeric(),
  mz = as.numeric(),
  Ion = as.character(),
  RT_low = as.numeric(),
  RT_high = as.numeric(),
  adduct.mass.low = as.numeric(),
  adduct.mass.high = as.numeric(),
  ID = as.character()
)

for (i in seq_along(ID.dat$ID)) {
  comp.ID <- ID.dat$ID[i]
  adduct.out <- find_adducts(comp.ID, ID.dat, adduct.list)
  adducts.dat <- rbind(adducts.dat, adduct.out)
}

#ad.dat.2 <- adducts.dat
ad.dat.test <- adducts.dat %>%
  rowwise() %>%
  filter(!adduct.ID == ID) %>%
  ungroup() %>%
  mutate(fraction = case_when(str_detect(.$ID, "Pos") ~ "Pos",
                              str_detect(.$ID, "Neg") ~ "Neg"))

ggplot(ad.dat.test, aes(y = Ion, fill = Ion)) +
  geom_bar() +
  facet_wrap(.~fraction)

ad.dat.4 <- ad.dat.test %>%
  filter(Ion %in% c("[M+H]+", "[M-H]-"))

ggplot(ad.dat.4, aes(x = Ion, fill = Ion)) +
  geom_bar() +
  facet_wrap(.~fraction)


#######
####Correlation among adducts and matched MFs:
Area.dat <- read_csv("Intermediates/ProMo_Particulate_QC2_data.csv")


corr.dat <- Area.dat %>%
  rename("ID" = MF) %>%
  select(ID, SampID, Area) # %>%
#  rename("Area" = Area.vol)

adduct.corr.dat <- corr.dat%>%
  rename("adduct.ID" = ID,
         "adduct_Area" = Area)

adduct.corr.test <- left_join(ad.dat.test, corr.dat)
adduct.corr.test.2 <- left_join(adduct.corr.test, adduct.corr.dat) %>%
  filter(!str_detect(.$SampID, "DDA")) %>%
  filter(!str_detect(.$SampID, "Std_")) %>%
  filter(!str_detect(.$SampID, "Poo")) %>%
  filter(!str_detect(.$SampID, "Blk")) %>%
  group_by(ID, adduct.ID) %>%
  mutate(adduct.corr = cor(Area, adduct_Area, use = "complete.obs"),
         na.num = sum(is.na(adduct_Area)),
         not.na.num = sum(!is.na(adduct_Area)),
         tot.num = sum(!is.na(adduct_Area))+na.num,
         perc.na = na.num/tot.num) %>%
  select(ID, adduct.ID, adduct.corr, not.na.num, perc.na, Ion) %>%
  unique() %>%
  filter(perc.na < 0.25) %>%
  filter(not.na.num > 20)

###Visualizations of adduct correlation by Ion and as a histogram
ggplot(adduct.corr.test.2, aes(y = Ion, x = adduct.corr)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(alpha = 0.5)

ggplot(adduct.corr.test.2, aes(adduct.corr)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.95, color = "red")



####select adducts with correlation values above a specific threshold. 
adduct.corr.test.3 <- adduct.corr.test.2 %>%
  filter(adduct.corr >= adduct.cor.threshold)

ggplot(adduct.corr.test.3, aes(y = Ion, fill = Ion)) +
  geom_bar() 




#Write Adducts to CSV

#rename data
ad.dat.final <- adduct.corr.test.3 %>%
  rename("MF" = adduct.ID,
         "adduct_of" = ID) %>%
  select(MF, adduct_of, Ion)


write_csv(ad.dat.final, file = "Intermediates/adduct_search_output.csv")




#######Decide which adducts/isotopes to remove:
ad.dat.2 <- read_csv("Intermediates/adduct_search_output.csv") 



###Decide between positive and negative mode for MFs detected in both sets, currently using signal to noise as the deciding factor
Neg.Hplus <- left_join(ad.dat.2, Area.dat) %>%
  ungroup() %>%
  filter(Ion %in% c("[M+H]+", "[M-H]-")) %>%
  filter(str_detect(.$adduct_of, "HPos")) %>%
  filter(str_detect(.$MF, "HNeg")) %>%
  select(MF, adduct_of, Name, FinalRSD) %>%
  unique() %>%
  rename("Neg.ID" = MF,
         "Pos.ID" = adduct_of,
         "Neg.Name" = Name,
         "Neg.RSD" = FinalRSD)


Pos.Hplus <- left_join(ad.dat.2, Area.dat) %>%
  ungroup() %>%
  filter(Ion %in% c("[M+H]+", "[M-H]-")) %>%
  filter(str_detect(.$adduct_of, "HNeg")) %>%
  filter(str_detect(.$MF, "HPos")) %>%
  select(MF, adduct_of, Name, FinalRSD) %>%
  unique() %>%
  rename("Pos.ID" = MF,
         "Neg.ID" = adduct_of,
         "Pos.Name" = Name,
         "Pos.RSD" = FinalRSD)



###Keep compounds that have IDs (Name does not = "Unknown")
### Then decide on compounds based on their RSD
Pos.v.neg <- full_join(Pos.Hplus, Neg.Hplus) 
  
PvN.Ions.to.keep <- Pos.v.neg %>%
  mutate(keep = case_when(!Pos.Name == "Unknown" ~ Pos.ID,
                          !Neg.Name == "Unknown" ~ Neg.ID)) %>%
  mutate(keep = case_when(!is.na(keep) ~ keep,
                          Pos.RSD < Neg.RSD ~ Pos.ID,
                          Neg.RSD < Pos.RSD ~ Neg.ID))%>%
  select(keep)


#### M+H Pos duplicates
Pos.v.Pos <- left_join(ad.dat.2, Area.dat) %>%
  ungroup() %>%
  filter(Ion %in% c("[M+H]+")) %>%
  filter(str_detect(.$adduct_of, "HPos")) %>%
  filter(str_detect(.$MF, "HPos")) %>%
  select(MF, adduct_of, Name, FinalRSD) %>%
  unique() %>%
  mutate(keep = case_when(!Name == "Unknown" ~ MF)) %>%
  select(keep)


### M-H Neg duplicates
Neg.v.Neg <- left_join(ad.dat.2, Area.dat) %>%
  ungroup() %>%
  filter(Ion %in% c("[M-H]-")) %>%
  filter(str_detect(.$adduct_of, "HNeg")) %>%
  filter(str_detect(.$MF, "HNeg")) %>%
  select(MF, adduct_of, Name, FinalRSD) %>%
  unique() %>%
  mutate(keep = case_when(!Name == "Unknown" ~ MF)) %>%
  select(keep)


###RP de-duplicating:
RP.pos.v.pos <- left_join(ad.dat.2, Area.dat) %>%
  ungroup() %>%
  filter(Ion %in% c("[M+H]+")) %>%
  filter(str_detect(.$adduct_of, "RP")) %>%
  filter(str_detect(.$MF, "RP")) %>%
  select(MF, adduct_of, Name, FinalRSD) %>%
  unique() %>%
  rename("RP_MF_1" = MF,
         "MF" = adduct_of,
         "MF1.RSD" = FinalRSD) %>%
  left_join(., Area.dat) %>%
  select(RP_MF_1, MF, Name, MF1.RSD, FinalRSD) %>%
  unique() %>%
  rename("RP_MF_2" = MF,
         "MF2.RSD" = FinalRSD) %>%
  mutate(keep = case_when(MF1.RSD <= MF2.RSD ~ RP_MF_1,
                          MF2.RSD <= MF1.RSD ~ RP_MF_2)) %>%
  select(keep)


#Finalize adduct "keep list"
final.ions.to.keep <- rbind(PvN.Ions.to.keep, Pos.v.Pos, Neg.v.Neg, RP.pos.v.pos) %>%
  unique() %>%
  filter(!is.na(keep)) %>%
  rename("MF" = keep)


####Finalize List of Adducts/Ions/Isotopes to remove:
ad.dat.remove <- ad.dat.2 %>%
  select(MF, Ion) %>%
  filter(!MF %in% final.ions.to.keep$MF) %>%
  filter(Ion %in% c("[13CM+H]+", "[M+H]+", "[M-H]-"))
write_csv(ad.dat.remove, file = "Intermediates/QC2_Adducts_to_Remove.csv")

ad.dat.annotate <- ad.dat.2 %>%
  select(MF, Ion) %>%
  filter(!MF %in% final.ions.to.keep$MF) %>%
  filter(!Ion %in% c("[13CM+H]+", "[M+H]+", "[M-H]-"))
write_csv(ad.dat.annotate, file = "Intermediates/QC2_Adducts_to_Annotate.csv")





