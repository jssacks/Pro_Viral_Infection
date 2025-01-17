



library(tidyverse)
library(ggsci)
library(ggthemes)
library(vegan)
library(viridis)
library(patchwork)




###load in datasets
metab.dat <- read_csv("Intermediates/Final_Processed_Untargeted_Data.csv")
quant.dat <- read_csv("Intermediates/Final_Processed_Quantified_Data.csv")

  
  
  


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
 # group_by(SampID, Rep, Treatment, Timepoint) %>% 
  #summarise(mean.conc = sum(nM.in.smp, na.rm = TRUE),
 #           mean.Nmol.C = sum(nM_C, na.rm = TRUE),
 #           mean.Nmol.N = sum(nM_N, na.rm = TRUE)) %>%
  mutate(Name = "Other",
         order = 20) %>%
  select(Name, SampID, nM.in.smp, nM_C, nM_N, order, Rep, Treatment, Timepoint)



###Put most abundant and other data together for figure  
ma.dat.fig <- full_join(subset.dat, ma.dat) %>%
  filter(!is.na(order)) %>%
  select(Name, SampID, nM.in.smp, nM_C, nM_N, order, Rep, Treatment, Timepoint)

dat.fig <- rbind(ma.dat.fig, dat.other) %>%
  mutate(Treatment = str_remove(Treatment, "1"), 
         Treatment = str_remove(Treatment, "3"),
         Treatment = str_remove(Treatment, "6")) %>%
  mutate(Treatment = as_factor(Treatment)) %>%
  mutate(Treatment = fct_relevel(Treatment, c("C", "LV", "HV")))
 # group_by(Name, order, Treatment, Timepoint) #%>%
 # reframe(nM.C = mean(mean.Nmol.C))

####absolute abundance fig
abs.fig <- ggplot(dat.fig, aes(x = Timepoint, y = nM_C, fill = reorder(Name, order))) +
  geom_col(color = "black", position = "stack", size = 0.2, width = 0.6) +
  facet_grid(Rep~Treatment, scales = "free_x") + 
  scale_fill_tableau(palette = "Tableau 20")+
  theme_test() +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  ylab("nM C") +
  labs(fill = "Compound") +
  #scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
abs.fig

####relative abundance fig
rel.fig <- ggplot(dat.fig, aes(x = Timepoint, y = nM_C, fill = reorder(Name, order))) +
  geom_col(color = "black", position = "fill", size = 0.2, width = 0.6) +
  facet_grid(Rep~Treatment, scales = "free_x") + 
  scale_fill_tableau(palette = "Tableau 20")+
  theme_test() +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  ylab("Mole Fraction C (%)") +
  labs(fill = "Compound") +
  #scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
rel.fig

####Make combined plot:
comb.plot <- abs.fig/rel.fig +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
comb.plot

##export plot:
ggsave(comb.plot, file = "Figures/Outputs/Supplemental_all_stacked_barcharts.png",
       bg = "white", dpi = 600, units = "in", height = 8, width = 8, scale = 1.3)



















































