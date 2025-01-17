


# Matrix_Population_Model_Cfixation_Results -------------------------------


library(tidyverse)
library(ggpubr)





###define inputs:
mpm.file <- "Collaborator_Data/MPM/MPM_results_final.csv"
mpm.anova.groups.file <- "Intermediates/MPM_ANOVA_Groups.csv"
pal <- c("#cfd4d8", "#73a596", "#365759")

##load in dataframe:
mpm.dat <- read_csv(mpm.file) 

mpm.fig.means <- mpm.dat %>%
  ungroup() %>%
  group_by(time, treatment) %>%
  mutate(mean.div = mean(division),
         sd.div = sd(division),
         mean.c.fix = mean(carbon_fixation),
         sd.c.fix = sd(carbon_fixation),
         mean.c.loss = mean(carbon_loss),
         sd.c.loss = sd(carbon_loss)) %>%
  select(treatment, time, mean.div, mean.c.fix, mean.c.loss, sd.div, sd.c.fix, sd.c.loss) %>%
  unique() %>%
  ungroup() %>%
  left_join(., read_csv(mpm.anova.groups.file) %>%
              rename("time" = Time) %>%
              select(-alpha) %>%
              pivot_wider(id_cols = c(treatment, time), names_from = parameter, values_from = groups) %>%
              rename("c.fix.groups" = C_fixation,
                     "c.loss.groups" = c_loss,
                     "div.groups" = division)) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV")) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) %>%
  unique() %>%
  mutate(time = case_when(time == 24 ~ "0-24",
                          time == 48 ~ "24-48"))

mpm.fig.all <- mpm.dat %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV")) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) %>%
  mutate(time = case_when(time == 24 ~ "0-24",
                          time == 48 ~ "24-48")) 



#division
div.fig <- ggplot(data = mpm.fig.means, aes(x = as.factor(time), y = mean.div)) +
  geom_col(aes(y = mean.div, fill = treatment), width = 0.5, 
           alpha = 1, position=position_dodge(width = 0.7), color = "black", size = 0.25) +
  geom_point(data = mpm.fig.all, aes(x = as.factor(time), y = division, group = treatment), 
             shape = 21, stroke = 0.5, size = 2, fill = "white", position = position_dodge(width = 0.7)) + 
  geom_errorbar(aes(group = treatment, ymax = mean.div + sd.div, ymin = mean.div - sd.div), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black", alpha = 0.7) + 
  geom_text(aes(group = treatment, x = as.factor(time), y = mean.div + sd.div+0.05, label = div.groups),
            position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = pal) +
  theme_test() +
  ylab(expression(Division~rate~(day^-1))) +
#  ylab("Division rate (1/day)") +
  xlab("Time window (h)") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 0.65)) +
  labs(fill = "Treatment")
div.fig

#C-fixation
cfix.fig <- ggplot(data = mpm.fig.means, aes(x = as.factor(time), y = mean.c.fix)) +
  geom_col(aes(y = mean.c.fix, fill = treatment), width = 0.5, alpha = 1, 
           position=position_dodge(width = 0.7), color = "black", size = 0.25) +
  geom_point(data = mpm.fig.all, aes(x = as.factor(time), y = carbon_fixation, group = treatment), 
             shape = 21, stroke = 0.5, size = 2, fill = "white", position = position_dodge(width = 0.7)) + 
  geom_errorbar(aes(group = treatment, min = mean.c.fix + sd.c.fix, ymax = mean.c.fix - sd.c.fix), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black", alpha = 0.7) + 
  geom_text(aes(group = treatment, x = as.factor(time), y = mean.c.fix + sd.c.fix+2.5, label = c.fix.groups),
            position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = pal) +
  theme_test() +
  ylab(expression(Carbon~fixation~rate~(fg~C~cell^-1~day^-1))) +
 # ylab("Carbon Fixation rate (fg C/cell/day)") +
  xlab("Time window (h)") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 42))+
  labs(fill = "Treatment")
cfix.fig

#C-loss
closs.fig <- ggplot(data = mpm.fig.means, aes(x = as.factor(time), y = mean.c.loss)) +
  geom_col(aes(y = mean.c.loss, fill = treatment), width = 0.5, alpha = 1, 
           position=position_dodge(width = 0.7), color = "black", size = 0.25) +
  geom_point(data = mpm.fig.all, aes(x = as.factor(time), y = carbon_loss, group = treatment), 
             shape = 21, stroke = 0.5, size = 2, fill = "white", position = position_dodge(width = 0.7)) + 
  geom_errorbar(aes(group = treatment, min = mean.c.loss + sd.c.loss, ymax = mean.c.loss - sd.c.loss), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black", alpha = 0.7) + 
  geom_text(aes(group = treatment, x = as.factor(time), y = mean.c.loss + sd.c.loss+2, label = c.loss.groups),
            position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = pal) +
  theme_test() +
  ylab(expression(Carbon~loss~rate~(fg~C~cell^-1~day^-1))) +
 # ylab("Carbon Loss rate (fg C/cell/day)") +
  xlab("Time window (h)") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 25))+
  labs(fill = "Treatment")
closs.fig


###########
mpm.group.fig <- ggarrange(
  NA, NA, NA, NA, NA, NA, NA,
  NA,  cfix.fig, NA, closs.fig, NA, div.fig, NA,
  NA, NA, NA, NA, NA, NA, NA,
  nrow = 3, ncol = 7, 
  widths = c(0.02, 0.3, 0.03, 0.3, 0.03, 0.3, 0.02),
  heights = c(0.05, 0.9, 0.05),
  common.legend = TRUE, legend = "bottom",
  labels = c(NA, NA, NA, NA, NA, NA, NA,
             NA, "A", NA, "B", NA, "C", NA,
             NA, NA, NA, NA, NA, NA, NA))
mpm.group.fig


#export figure
ggsave(mpm.group.fig, filename = "Figures/Outputs/MPM_Results.png", dpi = 800, scale = 1.4,
       units = "in", height = 3, width = 6, bg = "white"
       )






































































