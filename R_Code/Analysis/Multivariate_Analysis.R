

############
#Multivariate analyses and figures

#Goals: 
# 1. NMDS of samples
# 4. ANOSIM + PerMANOVA (show treatments are significantly different)



###################
#Packages

library(tidyverse)
library(vegan)
library(viridis)


#Source Code
source("Source_Code/biostats.R")



####################
#Define inputs
metab.dat <- read_csv("Intermediates/Final_Processed_Untargeted_Data.csv")

MF.details <- metab.dat %>%
  select(MF, mz, RT) %>%
  unique()


####Sample Meta data
samp.meta.dat <- metab.dat %>%
  select(SampID, Time, biorep, treatment, Pro, BacTot, BacSml, BacLrg, Phage) %>%
  unique() %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  filter(!SampID == "C_3LV_T48")


###Generate Sample Subsets:
norm.dat <- metab.dat %>%
  select(MF, SampID, Adjusted_Area, Vol.filt.mL, Time, biorep, treatment, Pro, BacTot, BacSml, BacLrg, Phage, Name, fraction) %>%
  mutate(vol.norm.area = Adjusted_Area/Vol.filt.mL) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  filter(!MF == "Vit_OH-pB12") %>%
  filter(!SampID == "C_3LV_T48") %>%
  unique()# %>%
#  filter(Time %in% c(0, 12, 24, 36, 48))





####################
# NMDS Analysis of Samples 
nmds.dat <- norm.dat %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

nmds.dat.stand <- decostand(nmds.dat, method = "max", MARGIN = 2)
nmds.dist <- vegdist(x = nmds.dat.stand, method = "manhattan")

nmds.out <- metaMDS(nmds.dist, k=2, trymax = 50)

nmds.out.plot <- data.frame(nmds.out$points) %>%
  rownames_to_column(var = "SampID") %>%
  left_join(., samp.meta.dat) %>%
  filter(Time %in% c(0, 12, 24, 36, 48))


#install.packages()
hull <- nmds.out.plot %>%
  group_by(as.factor(Time), treatment) %>%
  slice(chull(MDS1, MDS2))


ggplot(nmds.out.plot, aes(x = MDS2, y = MDS1, shape = treatment)) + 
  geom_polygon(data = nmds.out.plot, aes(x = MDS2, y = MDS1, fill = factor(Time)), color = "gray", alpha = 0) +
 # geom_polygon(data = hull, aes(x = MDS2, y = MDS1), alpha = 0.2) +
  geom_point(size = 3.5, stroke = 1, alpha = 0.8, aes(fill = factor(Time)))  +
  theme_test() + 
  scale_shape_manual(values = c(22, 21, 24)) +
  scale_fill_viridis(option = "B", discrete = "TRUE", direction = 1) +
  scale_color_viridis(option = "B", discrete = "TRUE", direction = 1) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) 
  #+
  #scale_fill_gradient2(low = "red2", mid = "white",high = "black", midpoint = 24)
  #scale_color_viridis(option = "B", discrete = "TRUE") +
 # guides(fill = guide_legend(override.aes = list(fill = ))).      aes(fill = factor(Time)
  
nmds.scree(nmds.dat.stand, distance = 'manhattan', k = 10)

#scree plot suggests that 2-3 dimensions are effective at reducing stress, choosing 2

mc.out <- nmds.monte(nmds.dat.stand, distance = "manhattan", k=2)
mc.out






############################################
#ANOSIM and/or PerMANOVA 

samp.dat <- norm.dat %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")
samp.dat.stand <- decostand(samp.dat, method = "max", MARGIN = 2)
samp.dist <- vegdist(x = samp.dat.stand, method = "manhattan")



adonis2(nmds.dist ~ treatment*Time*biorep, data = samp.meta.dat, permutations = 999)


###Subset comparisons: 

#Impact of time and biorep within each treatment

detach("package:MASS", unload=TRUE)

#Control
C.dat <- norm.dat %>%
  filter(treatment == "C") %>%
  ungroup() %>%
  dplyr::select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

C.dat.stand <- decostand(C.dat, method = "max", MARGIN = 2)

C.dist <- vegdist(x = C.dat.stand, method = "manhattan")

C.meta.dat <- samp.meta.dat %>%
  filter(treatment == "C") 

adonis2(C.dist ~ Time*biorep, data = C.meta.dat, permutations = 999)


#Low Virus
LV.dat <- norm.dat %>%
  filter(treatment == "LV") %>%
  dplyr::select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

LV.dat.stand <- decostand(LV.dat, method = "max", MARGIN = 2)

LV.dist <- vegdist(x = LV.dat.stand, method = "manhattan")

LV.meta.dat <- samp.meta.dat %>%
  filter(treatment == "LV")

adonis2(LV.dist ~ Time*biorep, data = LV.meta.dat, permutations = 999)


#High Virus
HV.dat <- norm.dat %>%
  filter(treatment == "HV") %>%
  dplyr::select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

HV.dat.stand <- decostand(HV.dat, method = "max", MARGIN = 2)

HV.dist <- vegdist(x = HV.dat.stand, method = "manhattan")

HV.meta.dat <- samp.meta.dat %>%
  filter(treatment == "HV")

adonis2(HV.dist ~ Time*biorep, data = HV.meta.dat, permutations = 999)




#Comparison of treatments directly to one another:

##Control vs. Low Virus
CvLV.dat <- norm.dat %>%
  filter(treatment %in% c("C", "LV")) %>%
  dplyr::select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

CvLV.dat.stand <- decostand(CvLV.dat, method = "max", MARGIN = 2)

CvLV.dist <- vegdist(x = CvLV.dat.stand, method = "manhattan")

CvLV.meta.dat <- samp.meta.dat %>%
  filter(treatment %in% c("C", "LV"))

adonis2(CvLV.dist ~ treatment*Time*biorep, data = CvLV.meta.dat, permutations = 999)



##Control vs. High Virus
CvHV.dat <- norm.dat %>%
  filter(treatment %in% c("C", "HV")) %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

CvHV.dat.stand <- decostand(CvHV.dat, method = "max", MARGIN = 2)

CvHV.dist <- vegdist(x = CvHV.dat.stand, method = "manhattan")

CvHV.meta.dat <- samp.meta.dat %>%
  filter(treatment %in% c("C", "HV"))

adonis2(CvHV.dist ~ treatment*Time*biorep, data = CvHV.meta.dat, permutations = 999)



##Low Virus vs. High Virus 
LVvHV.dat <- norm.dat %>%
  filter(treatment %in% c("LV", "HV")) %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

LVvHV.dat.stand <- decostand(LVvHV.dat, method = "max", MARGIN = 2)

LVvHV.dist <- vegdist(x = LVvHV.dat.stand, method = "manhattan")

LVvHV.meta.dat <- samp.meta.dat %>%
  filter(treatment %in% c("LV", "HV"))

adonis2(LVvHV.dist ~ treatment*Time*biorep, data = LVvHV.meta.dat, permutations = 999)





####timepoints:
T24.dat <- norm.dat %>%
  filter(Time %in% c(48)) %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")

T24.dat.stand <- decostand(T24.dat, method = "max", MARGIN = 2)

T24.dist <- vegdist(x = T24.dat.stand, method = "manhattan")

T24.meta.dat <- samp.meta.dat %>%
  filter(Time %in% c(48))

adonis2(T24.dist ~ treatment, data = T24.meta.dat, permutations = 999)









