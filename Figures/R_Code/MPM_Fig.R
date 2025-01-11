


# Matrix_Population_Model_Cfixation_Results -------------------------------


###define inputs:
mpm.file <- "Data_Raw/Collaborator_Contributed_Data/MPM_results_final.csv"


##load in dataframe:
mpm.dat <- read_csv(mpm.file) 

mpm.means <- mpm.dat %>%
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
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV")) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) 


#division
div.fig <- ggplot(data = mpm.means, aes(x = as.factor(time), y = mean.div)) +
  geom_col(aes(y = mean.div, fill = treatment), width = 0.5, alpha = 1, position=position_dodge(width = 0.7), color = "black") +
  geom_errorbar(aes(group = treatment, min = mean.div + sd.div, ymax = mean.div - sd.div), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black") + 
  scale_fill_manual(values = c("aliceblue",  "steelblue", "darkblue")) +
  theme_test() +
  ylab("Division rate (1/day)") +
  xlab("Time window") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 0.65))
div.fig

#C-fixation
cfix.fig <- ggplot(data = mpm.means, aes(x = as.factor(time), y = mean.c.fix)) +
  geom_col(aes(y = mean.c.fix, fill = treatment), width = 0.5, alpha = 1, position=position_dodge(width = 0.7), color = "black") +
  geom_errorbar(aes(group = treatment, min = mean.c.fix + sd.c.fix, ymax = mean.c.fix - sd.c.fix), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black") + 
  scale_fill_manual(values = c("aliceblue",  "steelblue", "darkblue")) +
  theme_test() +
  ylab("Carbon Fixation rate (fg C/cell/day)") +
  xlab("Time window") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 42))
cfix.fig

#C-loss
closs.fig <- ggplot(data = mpm.means, aes(x = as.factor(time), y = mean.c.loss)) +
  geom_col(aes(y = mean.c.loss, fill = treatment), width = 0.5, alpha = 1, position=position_dodge(width = 0.7), color = "black") +
  geom_errorbar(aes(group = treatment, min = mean.c.loss + sd.c.loss, ymax = mean.c.loss - sd.c.loss), 
                width = 0.3, position=position_dodge(width = 0.7), color = "black") + 
  scale_fill_manual(values = c("aliceblue",  "steelblue", "darkblue")) +
  theme_test() +
  ylab("Carbon Loss rate (fg C/cell/day)") +
  xlab("Time window") +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 25))
closs.fig


###########
mpm.group.fig <- ggarrange(
  NA, div.fig, NA, cfix.fig, NA, closs.fig, NA,
  nrow = 1, ncol = 7, widths = c(0.02, 0.3, 0.03, 0.3, 0.03, 0.3, 0.02),
  common.legend = TRUE, legend = "right",
  labels = c(NA, "A", NA, "B", NA, "C", NA))
mpm.group.fig

ggsave(mpm.group.fig, filename = "Figures/Figures/MPM_Results.png", dpi = 800, scale = 1.2,
       units = "in", height = 3, width = 8, bg = "white"
       )






































































