





library(tidyverse)
library(viridis)
library(vegan)
library(ggdendro)
library(ggpubr)




###______ Make Metabolite Figure ------------------------------------------------
metab.dat <- "Intermediates/FC_analysis_3_metab_output.csv"

pal.fc.plots <- c("aliceblue", "steelblue")


###
##Perform Fold change analysis using Pro Carbon data 
FC.Carbon.dat <- metab.meta.dat%>%
  mutate(vol.norm.area = (Adjusted_Area+1)/Vol.filt.mL,
         pro.norm.area = vol.norm.area/Pro.Carbon) %>%
  group_by(MF, time, biorep) %>%
  mutate(C.norm.area = pro.norm.area/(pro.norm.area[treatment == "C"]),
         log2.fc = log2(C.norm.area)) %>%
  ungroup() %>%
  filter(!SampID == "C_3LV_T48")


#Make L-Aspartic acid dataset
FC.AspA <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "L-Aspartic acid")

FC.AspA.means <- FC.AspA %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

AspA.plot <- ggplot(FC.AspA) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.AspA.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.AspA.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.AspA$Name)
AspA.plot



###########

#Make Adenine dataset
FC.Ade <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "Adenine")

FC.Ade.means <- FC.Ade %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

Ade.plot <- ggplot(FC.Ade) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.Ade.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.Ade.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.Ade$Name)
Ade.plot


#### Lysine
#Make Lysine dataset
FC.Lys <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "L-Lysine")

FC.Lys.means <- FC.Lys %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

Lys.plot <- ggplot(FC.Lys) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.Lys.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.Lys.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.Lys$Name)
Lys.plot

##############_____
#### UDP-N-acetylglucosamine
#Make dataset
FC.UDP <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "UDP-N-acetylglucosamine")

FC.UDP.means <- FC.UDP %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

UDP.plot <- ggplot(FC.UDP) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.UDP.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.UDP.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.UDP$Name)
UDP.plot



############
#Make dataset
FC.Des <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "Desthiobiotin")

FC.Des.means <- FC.Des %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

Des.plot <- ggplot(FC.Des) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.Des.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.Des.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.Des$Name)
Des.plot


####AMP figure
#Make dataset
FC.AMP <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("LV_0", "LV_12", "LV_24", "LV_36", "HV_0", "HV_12", "HV_24")) %>%
  select(MF, Name, treatment, time, treatment_time, biorep, SampID, C.norm.area, log2.fc) %>%
  mutate(High_Viral_Infection = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ TRUE,
                                          TRUE ~ FALSE)) %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  ungroup() %>%
  filter(Name == "Adenosine monophosphate")

FC.AMP.means <- FC.AMP %>%
  select(MF, Name, treatment, time, treatment_time, High_Viral_Infection, mean.log2.fc, sd.log2.fc) %>%
  unique()

AMP.plot <- ggplot(FC.AMP) +
  scale_fill_manual(values = pal.fc.plots) +
  geom_col(data = FC.AMP.means, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
  geom_errorbar(data = FC.AMP.means, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
  geom_vline(data=filter(FC.Carbon.fig, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
  theme_classic() +
  geom_hline(yintercept = 0) +
  xlab("Time (hr)") +
  ylab("Log2(FC) (Treatment/Control)") +
  theme(legend.position = "none") +
  labs(title = FC.AMP$Name)
AMP.plot



####Combine and assemble figures:
comb.fig <- ggarrange(AspA.plot, NA, Ade.plot,
                      NA, NA, NA,
                      Lys.plot, NA, UDP.plot,
                      NA, NA, NA,
                      Des.plot, NA, AMP.plot,
                      nrow = 5, ncol = 3, 
                      widths = c(0.47, 0.06, 0.47), heights = c(0.3, 0.05, 0.3, 0.05, 0.3),
                      labels = c("A", NA, "B",
                                 NA, NA, NA,
                                 "C", NA, "D",
                                 NA, NA, NA, 
                                 "E", NA, "F"))
comb.fig

ggsave(filename = "Figures/Supplemental_Metab_FC_Figure.png",
       dpi = 300, height = 7, width = 6, units = "in",
       bg = "white", scale = 1.25)









































































