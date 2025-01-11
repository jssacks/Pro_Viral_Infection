


library(tidyverse)
library(ggsci)
library(ggthemes)
library(vegan)
library(viridis)



###
quant.dat <- read_csv("Intermediates/ProMo_Particulate_Quant_Output.csv")
qc2.dat <- read_csv("Intermediates/ProMo_Particulate_QC2_data.csv")





# NMDS_Code ---------------------------------------------------------------
####################
#Define inputs

MF.details <- qc2.dat %>%
  select(MF, mz, RT) %>%
  unique()


####Sample Meta data
samp.meta.dat <- qc2.dat %>%
  select(SampID, Time, biorep, treatment, Pro, BacTot, BacSml, BacLrg, Phage) %>%
  unique() %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  filter(!SampID == "C_3LV_T48")


###Generate Sample Subsets:
norm.dat <- qc2.dat %>%
  select(MF, SampID, Adjusted_Area, Vol.filt.mL, Time, biorep, treatment, Pro, BacTot, BacSml, BacLrg, Phage, Name, fraction) %>%
  mutate(vol.norm.area = Adjusted_Area/Vol.filt.mL) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  filter(!MF == "Vit_OH-pB12") %>%
  filter(!SampID == "C_3LV_T48")# %>%
#  filter(Time %in% c(0, 12, 24, 36, 48))

nmds.dat <- norm.dat %>%
  select(MF, SampID, vol.norm.area) %>%
  pivot_wider(id_cols = SampID, names_from = MF, values_from = vol.norm.area) %>%
  column_to_rownames(var = "SampID")


####################
# NMDS Analysis of Samples 

norm.dat <- qc2.dat %>%
  select(MF, SampID, Adjusted_Area, Vol.filt.mL, Time, biorep, treatment, Pro, BacTot, BacSml, BacLrg, Phage, Name, fraction) %>%
  mutate(vol.norm.area = Adjusted_Area/Vol.filt.mL) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  filter(!MF == "Vit_OH-pB12") %>%
  filter(!SampID == "C_3LV_T48")# %>%
#  filter(Time %in% c(0, 12, 24, 36, 48))

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


nmds.fig <- ggplot(nmds.out.plot, aes(x = MDS2, y = MDS1, shape = treatment)) + 
  geom_polygon(data = nmds.out.plot, aes(x = MDS2, y = MDS1, fill = factor(Time)), color = "gray", alpha = 0) +
  # geom_polygon(data = hull, aes(x = MDS2, y = MDS1), alpha = 0.2) +
  geom_point(size = 3.5, stroke = 1, alpha = 0.8, aes(fill = factor(Time)))  +
  theme_test() + 
  scale_shape_manual(values = c(22, 21, 24)) +
  scale_fill_viridis(option = "B", discrete = "TRUE", direction = 1) +
  scale_color_viridis(option = "B", discrete = "TRUE", direction = 1) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  labs(fill = "Time") +
  theme(legend.position = "right")
nmds.fig





##################################
###Define inputs
microb.C.file <- "Intermediates/microbial_C_content.csv"
quant.metab.file <- "Intermediates/ProMo_Particulate_Quant_Output.csv"
metab.meta.file  <- "Intermediates/ProMo_Particulate_QC2_data.csv"


#metabolite data:
metab.meta.dat <- read_csv(metab.file) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  select(SampID, Time, biorep, treatment) %>%
  rename("time" = Time) %>%
  unique()

quant.metab.dat <- read_csv(quant.metab.file) %>%
  left_join(., metab.meta.dat) %>%
  filter(!Name == "Succinic acid") %>%
  filter(!Name == "Vit_OH-pB12") %>%
  filter(!Name == "Taurine, 2H4") %>%
  filter(!SampID == "C_3LV_T48")


##microbial C data:
m.C.dat <- read_csv(microb.C.file) %>%
  mutate(microb.C.nM = total.C*(1/12.011)*(1e-6)*1000) %>%
  select(experiment, treatment, time, microb.C.nM) %>%
  rename("biorep" = experiment)

###Combine metabolite and C data:
metab.c.dat <- left_join(quant.metab.dat, m.C.dat) %>%
  mutate(Metab_Perc_C = (nM_C/microb.C.nM)*100)


###Aggregate figure
agg.dat <- metab.c.dat %>%
  group_by(biorep, treatment, time) %>%
  summarize(metab.tot.C = sum(Metab_Perc_C)) %>%
  filter(!is.na(treatment)) %>%
  ungroup() %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV")))  

metab.c.perc.fig <- ggplot(agg.dat, aes(x = treatment, y = metab.tot.C)) +
  geom_boxplot(width = 0.5, alpha = 0.5) +
  geom_jitter(size = 2.5, width = 0.1, aes(shape = biorep, fill = treatment)) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("aliceblue", "steelblue", "darkblue")) +
  scale_y_continuous(limits = c(0,15)) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab("Percent of Microbial C in metabolites (%)")



################################# Quantified Metabolite Bar Plots

subset.dat <- quant.dat %>%
  filter(!str_detect(.$SampID, "LG_")) %>%
  # filter(!str_detect(.$SampID, "LGV")) %>%
  filter(!str_detect(.$SampID, "8G")) %>%
  filter(!str_detect(.$SampID, "9V")) %>%
  filter(!str_detect(.$SampID, "blk")) %>%
  separate(SampID, into = c("Rep", "Treatment", "Timepoint"), remove = FALSE) %>%
  filter(!Treatment == "4LGV") %>%
  mutate(Timepoint = str_replace(Timepoint, "T6", "T06")) %>%
  filter(!SampID == "C_3LV_T48")

#ggplot(subset.dat, aes(x = SampID, y = nM_C, fill = Name)) +
#  geom_col(color = "black")


###Generate figure
dat <- subset.dat %>%
  filter(!Name == "Succinic acid")

dat.2 <- dat %>%  
  rename(mean.conc = nM.in.smp,
         mean.Nmol.C = nM_C,
         mean.Nmol.N = nM_N) 
#dat.sum
###get 12 most abundant compounds
ma.dat <- dat %>%
  group_by(Name) %>%
  summarize(max.abu = max(nM_C)) %>%
  slice_max(max.abu, n = 19) %>%
  mutate(order = row_number()) %>%
  select(Name, order)

####Summarize the mean.conc, nmol.C, and nmol.N of the less abundant compounds
dat.other <- full_join(dat, ma.dat) %>%
  filter(is.na(order)) %>%
  group_by(SampID, Rep, Treatment, Timepoint) %>%
  summarise(mean.conc = sum(nM.in.smp, na.rm = TRUE),
            mean.Nmol.C = sum(nM_C, na.rm = TRUE),
            mean.Nmol.N = sum(nM_N, na.rm = TRUE)) %>%
  mutate(Name = "Other",
         order = 20)



###Put most abundant and other data together for figure  
ma.dat.fig <- full_join(dat.2, ma.dat) %>%
  filter(!is.na(order)) %>%
  select(Name, SampID, mean.conc, mean.Nmol.C, mean.Nmol.N, order, Rep, Treatment, Timepoint)

dat.fig <- rbind(ma.dat.fig, dat.other) %>%
  group_by(Name, order, Treatment, Timepoint) %>%
  reframe(nM.C = mean(mean.Nmol.C))


#ggplot(dat.fig, aes(x = Timepoint, y = mean.Nmol.C, fill = reorder(Name, order))) +
#  geom_col(color = "black", alpha = 0.9, size = 0.2) +
#  facet_grid(Rep~Treatment, scales = "free_x") + 
  # scale_fill_tableau() +
#  scale_fill_tableau(palette = "Tableau 20")+
  #scale_fill_manual(values = ggthemes_data)
  # scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
#  theme(axis.text.x = element_text(angle = 45))

barplot.fig <- ggplot(dat.fig, aes(x = Timepoint, y = nM.C, fill = reorder(Name, order))) +
  geom_col(color = "black", position = "fill", size = 0.2) +
  facet_grid(.~Treatment, scales = "free_x") + 
  scale_fill_tableau(palette = "Tableau 20")+
  theme_test() +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  ylab("Mean mol % C") +
  labs(fill = "Compound") +
  #scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
barplot.fig




###########Arrange into a nice plot:
metabolome.fig <- ggarrange(
  ggarrange(nmds.fig, NA, metab.c.perc.fig,
            nrow = 1, ncol = 3, widths = c(0.6, 0.05, 0.35),
            labels = c("A", NA, "B")), NA,
  barplot.fig, nrow = 3, ncol = 1, heights = c(0.4, 0.05, 0.55),
  labels = c(NA, NA, "C")
  )
metabolome.fig

ggsave(metabolome.fig, filename = "Figures/Figures/NMDS_QuantMetab_Fig.png",
       dpi = 800, scale = 1.3, bg = "white",
       units = "in", height = 8, width = 8)
































































































































