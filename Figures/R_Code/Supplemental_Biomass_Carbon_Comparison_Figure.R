



library(tidyverse)
library(patchwork)


#define inputs:
atp.file <- "Intermediates/pATP_nM_processed.csv"
microb.c.file <- "Intermediates/Metabolite_Carbon_Comparison.csv"

#palette
pal <- c("#465f61", "#dd9914", "#73a596")






#load in data
metab.c.comparison <- read_csv(microb.c.file)
atp.dat <- read_csv(atp.file)

##combine atp + metab + c dat
c.dat <- left_join(atp.dat, metab.c.comparison) %>%
  mutate(microb.C.fg.mL = total.C) %>%
  mutate(metab.tot.C.fg.mL = tot_metab_C_nM*12.011/1000/(1e-6)) %>%
  mutate(metab.tot.C.fg.mL.10 = tot_metab_C_nM*12.011/1000/(1e-6)*10)


# Make plot of correlations between prameters -----------------------------

####
atp.fcm.plot <- ggplot(c.dat, aes(y = mean.live.C.fg.mL, x = microb.C.fg.mL)) +
  geom_smooth(method = "lm", color = "#dd9914") +
  geom_point(fill = "#465f61", shape = 21, stroke = 0.1, size = 3, alpha = 0.7) +
  geom_abline(slope = 1) +
  xlim(0, 2.25e9) +
  ylim(0, 2.25e9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab(expression(FCM~biomass~(fg~C~mL^-1))) +
  ylab(expression(ATP~biomass~(fg~C~mL^-1))) 
atp.fcm.plot
####
atp.metab.plot <- ggplot(c.dat, aes(y = mean.live.C.fg.mL, x = metab.tot.C.fg.mL)) +
  geom_smooth(method = "lm", color = "#dd9914") +
  geom_point(fill = "#465f61", shape = 21, stroke = 0.1, size = 3, alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab(expression(Metab~C~(fg~mL^-1))) +
  ylab(expression(ATP~biomass~(fg~C~mL^-1))) 
atp.metab.plot
###
fcm.metab.plot <- ggplot(c.dat, aes(y = microb.C.fg.mL, x = metab.tot.C.fg.mL)) +
  geom_smooth(method = "lm", color = "#dd9914") +
  geom_point(fill = "#465f61", shape = 21, stroke = 0.1, size = 3, alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab(expression(Metab~C~(fg~mL^-1))) +
  ylab(expression(FCM~biomass~(fg~C~mL^-1)))
fcm.metab.plot






# Make biomass plot by treatment and time ---------------------------------

###organize Means
c.dat.sum.means <- c.dat %>%
  group_by(treatment, time) %>%
  reframe(mean.fcm.c = mean(microb.C.fg.mL, na.rm = TRUE),
           # sd.fcm.c = sd(microb.C.fg.mL, na.rm = TRUE),
            mean.atp.c = mean(mean.live.C.fg.mL, na.rm = TRUE),
          #  sd.atp.c = sd(mean.live.C.fg.mL, na.rm = TRUE),
            mean.metab.c.10 = mean(metab.tot.C.fg.mL.10, na.rm = TRUE),
           # sd.metab.c.10 = sd(metab.tot.C.fg.mL.10, na.rm = TRUE)) %>%
  ) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) %>%
  pivot_longer(cols = mean.fcm.c:mean.metab.c.10) %>%
  rename("Mean_fg_C" = value) %>%
  mutate(parameter = case_when(str_detect(name, "atp") ~ "ATP",
                               str_detect(name, "fcm") ~ "FCM",
                               str_detect(name, "metab") ~ "metabX10")) %>%
  select(treatment, time, Mean_fg_C, parameter) 

#organize standard deviations
c.dat.sum.sds <- c.dat %>%
  group_by(treatment, time) %>%
  reframe(sd.fcm.c = sd(microb.C.fg.mL, na.rm = TRUE),
          sd.atp.c = sd(mean.live.C.fg.mL, na.rm = TRUE),
          sd.metab.c.10 = sd(metab.tot.C.fg.mL.10, na.rm = TRUE)
  ) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) %>%
  pivot_longer(cols = sd.fcm.c:sd.metab.c.10) %>%
  rename("SD_fg_C" = value) %>%
  mutate(parameter = case_when(str_detect(name, "atp") ~ "ATP",
                               str_detect(name, "fcm") ~ "FCM",
                               str_detect(name, "metab") ~ "metabX10")) %>%
  select(treatment, time, SD_fg_C, parameter) 

###merge datasets
c.dat.sum <- left_join(c.dat.sum.means, c.dat.sum.sds)

biomass.plot <- ggplot(c.dat.sum, aes(x = time, y = Mean_fg_C)) +
  geom_errorbar(aes(ymin = Mean_fg_C-SD_fg_C, ymax = Mean_fg_C+SD_fg_C), width = 0.7, size = 0.3) +
  geom_line(aes(color = parameter)) +
  geom_point(aes(fill = parameter), shape = 21, size = 3) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  facet_wrap(.~treatment) +
  theme_test() +
  xlab("Time (h)") + 
  ylab(expression(Biomass~(fg~C~mL^-1)))
biomass.plot




###combine plots___________
comb.plot <- (atp.fcm.plot | atp.metab.plot | fcm.metab.plot)/biomass.plot + plot_annotation(tag_levels = 'A')
comb.plot

##Export Plot:
ggsave(comb.plot, file = "Figures/Outputs/Biomass_FCM_ATP_Metab_Comparison_Supplemental_Fig.png", dpi = 600,
       bg = "white", units = "in", height = 6, width = 8, scale = 1.3)




