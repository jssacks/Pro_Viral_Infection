##


install.packages("multcomp")

#load packages:
library(tidyverse)
library(multcomp)
library(agricolae)


#define inputs:
mpm.file <- "Collaborator_Data/MPM/MPM_results_final.csv"





##import data:
mpm.dat <- read_csv(mpm.file) 

##change data to factors:
mpm.dat.factor <- mpm.dat %>%
  mutate(treatment = as.factor(treatment),
         time = as.factor(time),
         replicate = as.factor(replicate)) %>%
  mutate(treatment = fct_relevel(treatment, c("Control", "Low virus", "High virus")))

#Visualize Datasets:

#carbon fixation
ggplot(mpm.dat.factor, aes(x = as.factor(time), y = carbon_fixation, color = treatment)) +
# # geom_col(position = position_dodge(width = 0.7), width = 0.4, color = "black", size = 0.1)+ #  position_dodge(preserve = "single")) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
 # geom_point(shape = 21, position = position_dodge(width = 0.72)) +
#  scale_fill_manual(values = pal) +
  ylim(0, 40) 

#cell division
ggplot(mpm.dat, aes(x = as.factor(time), y = division, color = treatment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  ylim(0,NA)

#carbon loss
ggplot(mpm.dat, aes(x = as.factor(time), y = carbon_loss, color = treatment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  ylim(0, 30)


#stats on carbon fixation

##2 way anova:
cfix.anova <- aov(carbon_fixation ~ treatment*time+replicate, mpm.dat.factor)
summary(cfix.anova)

##Post-Hoc Test:
cfix.hsd <- TukeyHSD(cfix.anova, which = "treatment:time")
cfix.hsd

cfix.hsd.groups <- HSD.test(cfix.anova,
               trt = c("treatment", "time"), 
               alpha = 0.01,
               console = TRUE)

cfix.groups <- as.data.frame(cfix.hsd.groups$groups) %>%
  rownames_to_column(var = "treatment_timepoint") %>%
  separate(col = treatment_timepoint, into = c("treatment", "Time"), sep = ":") %>%
  mutate(parameter = "C_fixation") %>%
  dplyr::select(treatment, Time, groups, parameter)





################
#stats on cell division

##2 way anova:
cdiv.anova <- aov(division ~ treatment*time+replicate, mpm.dat.factor)
summary(cdiv.anova)

##Post-Hoc Test:
cdiv.hsd <- TukeyHSD(cdiv.anova, which = "treatment:time")
cdiv.hsd

cdiv.hsd.groups <- HSD.test(cdiv.anova,
                            trt = c("treatment", "time"), 
                            alpha = 0.01,
                            console = TRUE)
cdiv.hsd.groups
cdiv.groups <- as.data.frame(cdiv.hsd.groups$groups) %>%
  rownames_to_column(var = "treatment_timepoint") %>%
  separate(col = treatment_timepoint, into = c("treatment", "Time"), sep = ":") %>%
  mutate(parameter = "division") %>%
  dplyr::select(treatment, Time, groups, parameter)
#plot(cdiv.hsd.groups)



################
#stats on carbon loss

##2 way anova:
closs.anova <- aov(carbon_loss ~ treatment*time+replicate, mpm.dat.factor)
summary(closs.anova)

##Post-Hoc Test:
closs.hsd <- TukeyHSD(closs.anova, which = "treatment:time")
closs.hsd

closs.hsd.groups <- HSD.test(closs.anova,
                            trt = c("treatment", "time"), 
                            alpha = 0.01,
                            console = TRUE)
closs.hsd.groups
closs.groups <- as.data.frame(closs.hsd.groups$groups) %>%
  rownames_to_column(var = "treatment_timepoint") %>%
  separate(col = treatment_timepoint, into = c("treatment", "Time"), sep = ":") %>%
  mutate(parameter = "c_loss") %>%
  dplyr::select(treatment, Time, groups, parameter)




###Combine and export all group information
all.anova.groups <- rbind(closs.groups, cdiv.groups, cfix.groups) %>%
  mutate(alpha = 0.01)

write_csv(all.anova.groups, file = "Intermediates/MPM_ANOVA_Groups.csv")






































































































