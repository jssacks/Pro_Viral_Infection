##


install.packages("multcomp")

#load packages:
library(tidyverse)
library(multcomp)


#define inputs:
mpm.file <- "Data_Raw/Collaborator_Contributed_Data/MPM_results_final_raw.csv"



##import data:
mpm.dat <- read_csv(mpm.file) 

##change data to factors:
mpm.dat.factor <- mpm.dat %>%
  mutate(treatment = as.factor(treatment),
         time = as.factor(time),
         replicate = as.factor(replicate))

#Visualize Datasets:

#carbon fixation
ggplot(mpm.dat, aes(x = as.factor(time), y = carbon_fixation, color = treatment)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  ylim(0, 40)

#cell division
ggplot(mpm.dat, aes(x = as.factor(time), y = cell_division, color = treatment)) +
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
cfix.hsd.groups


################
#stats on cell division

##2 way anova:
cdiv.anova <- aov(cell_division ~ treatment*time+replicate, mpm.dat.factor)
summary(cdiv.anova)

##Post-Hoc Test:
cdiv.hsd <- TukeyHSD(cdiv.anova, which = "treatment:time")
cdiv.hsd

cdiv.hsd.groups <- HSD.test(cdiv.anova,
                            trt = c("treatment", "time"), 
                            alpha = 0.01,
                            console = TRUE)
cdiv.hsd.groups
plot(cdiv.hsd.groups)



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
plot(closs.hsd.groups)










#cfix.anova.2 <- aov(carbon_fixation ~ treatment*time + replicate, mpm.dat.factor)

#summary(cfix.anova.2)

#TukeyHSD(cfix.anova, which = "treatment")

res_Tukey <- glht(
  aov(carbon_fixation ~ treatment*time+replicate, data = mpm.dat.factor),
  linfct = mcp(treatment = "Tukey")
)
summary(res_Tukey)


xx <- TukeyHSD(cfix.anova, which = "treatment:time")
xx









#install.packages("agricolae")

library(agricolae)

x2 <- HSD.test(cfix.anova,
               trt = c("treatment", "time"), 
               alpha = 0.01,
               console = TRUE)
x2


plot(x2, las = 1)

##############

library(nlme)
install.packages("emmeans")

mod <- lme(carbon_fixation ~ treatment*time, random = ~1 | replicate, data = mpm.dat.factor)


































































































