




library(tidyverse)
library(ggpubr)






#Inputs
atp.file <- "Intermediates/pATP_nM_processed.csv"
microb.per.c.file <- "Intermediates/Metab_Perc_Microb_C.csv"



#load in data
atp.dat <- read_csv(atp.file)
metab.c.dat <- read_csv(microb.per.c.file)


##combine atp + metab + c dat

c.dat <- left_join(atp.dat, metab.c.dat) %>%
  mutate(microb.C.fg.mL = microb.C.nM*12.011/1000/(1e-6)) %>%
  mutate(metab.tot.C.fg.mL = metab.tot.nM.C*12.011/1000/(1e-6)) %>%
  mutate(metab.tot.C.fg.mL.10 = metab.tot.nM.C*12.011/1000/(1e-6)*10)




####
atp.fcm.plot <- ggplot(c.dat, aes(y = mean.live.C.fg.mL, x = microb.C.fg.mL)) +
  geom_smooth(method = "lm") +
  geom_point() +
  geom_abline(slope = 1) +
  xlim(0, 2.25e9) +
  ylim(0, 2.25e9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("FCM Estimated C (fg/mL)") +
  ylab("ATP Estaimted C (fg/mL") 
atp.fcm.plot
####
atp.metab.plot <- ggplot(c.dat, aes(y = mean.live.C.fg.mL, x = metab.tot.C.fg.mL)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("Metab C (fg/mL)") +
  ylab("ATP Estaimted C (fg/mL")
atp.metab.plot
###
fcm.metab.plot <- ggplot(c.dat, aes(y = microb.C.fg.mL, x = metab.tot.C.fg.mL)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("Metab C (fg/mL)") +
  ylab("FCM Estaimted C (fg/mL")
fcm.metab.plot

####
ggplot(c.dat, aes(x = time)) +
  geom_point(aes(y = microb.C.fg.mL), color = "red") +
  geom_line(aes(y = microb.C.fg.mL), color = "red") +
  geom_point(aes(y = mean.live.C.fg.mL), color = "blue") +
  geom_line(aes(y = mean.live.C.fg.mL), color = "blue") +
  geom_point(aes(y = metab.tot.C.fg.mL.10), color = "black") +
  geom_line(aes(y = metab.tot.C.fg.mL.10), color = "black") +
  facet_grid(treatment~biorep)


###calculate correlation coefficients between variables:
c.dat.cor <- c.dat %>%
  ungroup() %>%
  summarize(atp.fcm.cor = cor(mean.live.C.fg.mL, microb.C.fg.mL, use = "pairwise.complete.obs"),
            atp.fcm.r2 = atp.fcm.cor^2,
            atp.metab.cor = cor(mean.live.C.fg.mL, metab.tot.C.fg.mL.10, use = "pairwise.complete.obs"),
            atp.metab.r2 = atp.metab.cor^2,
            fcm.metab.cor = cor(microb.C.fg.mL, metab.tot.C.fg.mL.10, use = "pairwise.complete.obs"),
            fcm.metab.r2 = fcm.metab.cor^2)



###Combine replicates:
c.dat.sum <- c.dat %>%
  group_by(treatment, time) %>%
  summarize(mean.fcm.c = mean(microb.C.fg.mL, na.rm = TRUE),
            sd.fcm.c = sd(microb.C.fg.mL, na.rm = TRUE),
            mean.atp.c = mean(mean.live.C.fg.mL, na.rm = TRUE),
            sd.atp.c = sd(mean.live.C.fg.mL, na.rm = TRUE),
            mean.metab.c.10 = mean(metab.tot.C.fg.mL.10, na.rm = TRUE),
            sd.metab.c.10 = sd(metab.tot.C.fg.mL.10, na.rm = TRUE),
            mean.metab.c.nM = mean(metab.tot.nM.C, na.rm = TRUE),
            sd.metab.c.nM = sd(metab.tot.nM.C, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV")))  

treat.plot <- ggplot(c.dat.sum, aes(x = time)) +
  geom_point(aes(y = mean.fcm.c), color = "red", shape = 17, size = 2.5) +
  geom_line(aes(y = mean.fcm.c), color = "red") +
  geom_errorbar(aes(ymin = mean.fcm.c-sd.fcm.c, ymax = mean.fcm.c+sd.fcm.c), color = "red", width = 0.2) +
  geom_point(aes(y = mean.atp.c), color = "blue", shape = 15, size = 2.5) +
  geom_line(aes(y = mean.atp.c), color = "blue") +
  geom_errorbar(aes(ymin = mean.atp.c-sd.atp.c, ymax = mean.atp.c+sd.atp.c), color = "blue", width = 0.2) +
  geom_point(aes(y = mean.metab.c.10), color = "black", shape = 19, size = 2.5) +
  geom_line(aes(y = mean.metab.c.10), color = "black") +
  geom_errorbar(aes(ymin = mean.metab.c.10-sd.metab.c.10, ymax = mean.metab.c.10+sd.metab.c.10), color = "black", width = 0.2) +
  facet_wrap(.~treatment) +
  ylab("Carbon (fg/mL)") +
  xlab("Time (hr)") +
  theme_bw()
treat.plot



x <- ggarrange(ggarrange(atp.fcm.plot, NA, atp.metab.plot, NA, fcm.metab.plot,
                         ncol = 5, nrow = 1, labels = c("A", NA, "B", NA, "C"), widths = c(0.31, 0.035, 0.31, 0.035, 0.31)),
               NA, treat.plot, ncol = 1, nrow = 3, labels = c(NA, NA, "D"), heights = c(0.47, 0.06, 0.47))
x

ggsave(filename = "Figures/ATP_FCM_Metab_Figure.png",
       height = 4, width = 6, dpi = 300,
       unit = "in", scale = 1.5, bg = "white")













##Just metabolites
ggplot(c.dat.sum, aes(x = time, color = treatment)) +
  geom_line(aes(y = mean.metab.c.nM), size = 1) +
  geom_errorbar(aes(ymin = mean.metab.c.nM-sd.metab.c.nM, ymax = mean.metab.c.nM+sd.metab.c.nM), width = 0.2) +
  geom_point(aes(y = mean.metab.c.nM, shape = treatment, fill = treatment), size = 3.5, color = "black") +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("aliceblue", "steelblue3", "darkblue")) +
  scale_color_manual(values = c("gray20", "steelblue3", "darkblue")) +
  theme_bw() +
  ylab("Mean Metabolite C (nM)")

ggsave(filename = "Figures/Metabolite_C_over_time.png", 
       height = 3, width = 4, dpi = 300, scale = 1.2)













