

#
#
library(tidyverse)
library(patchwork)




#Define inputs
metab.file <- "Intermediates/Final_Processed_Untargeted_Data.csv"
fcm.carbon.file <- "Intermediates/microbial_C_content.csv"
atp.file <- "Intermediates/pATP_nM_processed.csv"
v.mf.file = "Intermediates/virally_altered_mfs.csv"

pal.fc.plots <- c("#73a596", "aliceblue")
#




# Compile and tidy datasets -----------------------------------------------

#metabolite data::
metab.dat <- read_csv(metab.file) %>%
  select(SampID, Vol.filt.mL, Time, biorep, treatment, treatment_number, MF, Name, Adjusted_Area, Pro, BacTot, BacSml, BacLrg, Phage)  %>%
  rename("time" = Time)

#make list of known virally altered MFs:
vfc.known.dat <- read_csv(v.mf.file) %>%
  filter(!Name == "Unknown")
  


###ATP data:

#pull in ATP data
atp.dat <- read_csv(atp.file) %>%
  select(-sd.conc.nM) 

#combine ATP data with metadata for each sample from metab.dat
atp.meta.dat <- metab.dat %>%
  select(SampID, Vol.filt.mL, time, biorep, treatment, treatment_number, Pro, BacTot, BacSml, BacLrg, Phage) %>%
  unique() %>%
  left_join(., atp.dat %>%
              select(time, biorep, treatment, mean.conc.nM)) %>%
  mutate(MF = "ATP",
         Name = "ATP") %>%
  mutate("Adjusted_Area" = mean.conc.nM*Vol.filt.mL) %>%
  select(-mean.conc.nM)


##Combine Metabolite and ATP Data and add in carbon data:
metab.atp.dat <- rbind(metab.dat, atp.meta.dat) %>%
  left_join(., read_csv(fcm.carbon.file) %>% rename("biorep" = experiment) %>%
              select(biorep, treatment, time, Pro.C, bacterial.C, total.C))






# Fold change analysis normalized to microbial biomass data  ---------------------
#####

#calculate fold change compared to control (normalized to volume filtered and biomass)
FC.Carbon.dat <- metab.atp.dat %>%
  mutate(vol.norm.area = (Adjusted_Area+1)/Vol.filt.mL,
         bio.norm.area = vol.norm.area/total.C) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  group_by(MF, time, biorep) %>%
  mutate(C.norm.area = bio.norm.area/(bio.norm.area[treatment == "C"]),
         log2.fc = log2(C.norm.area)) %>%
  ungroup() %>%
  filter(!SampID == "C_3LV_T48")


#prepare data for t-test
FC.Carbon.test <- FC.Carbon.dat %>%
  filter(!treatment == "C") %>%
  unite(c("treatment", "time"), col = "treatment_time", remove = FALSE) %>%
  filter(treatment_time %in% c("HV_0", "LV_0", "LV_12", "LV_24", "LV_36", "HV_12", "HV_24")) %>%
  mutate(condition = case_when(treatment_time %in% c("LV_36", "HV_12", "HV_24") ~ "High Viral Infection",
                               TRUE ~ "No")) %>%
  select(MF, Name, treatment, time, biorep, SampID, C.norm.area, log2.fc, condition) %>%
  filter(Name %in% vfc.known.dat$Name)%>%
  unique() %>%
  group_by(MF, treatment, time) %>%
  mutate(mean.log2.fc = mean(log2.fc),
         sd.log2.fc = sd(log2.fc)) %>%
  select(-biorep, - SampID, -C.norm.area, -log2.fc) %>%
  unique()

#Make a dataset for each compound:
suc.dat <- FC.Carbon.test %>% filter(Name == "Sucrose")
UDP.dat <- FC.Carbon.test %>% filter(Name == "UDP-N-acetylglucosamine")
AMP.dat <- FC.Carbon.test %>% filter(Name == "Adenosine monophosphate")
GBB.dat <- FC.Carbon.test %>% filter(Name == "(3-Carboxypropyl)trimethylammonium")
lys.dat <- FC.Carbon.test %>% filter(Name == "L-Lysine")
bet.dat <- FC.Carbon.test %>% filter(Name == "Betonicine")
GG.dat <- FC.Carbon.test %>% filter(Name == "2-O-alpha-D-Glucosylglycerol")
GSSG.dat <- FC.Carbon.test %>% filter(Name == "Glutathione disulfide")
asp.dat <- FC.Carbon.test %>% filter(Name == "L-Aspartic acid")
Ad.dat <- FC.Carbon.test %>% filter(Name == "Adenine")
OA.dat <- FC.Carbon.test %>% filter(Name == "Ophthalmic acid")
Des.dat <- FC.Carbon.test %>% filter(Name == "Desthiobiotin")
SAM.dat <- FC.Carbon.test %>% filter(Name == "S-Adenosylmethionine")
pb12.dat <- FC.Carbon.test %>% filter(Name == "OH-Pseudocob")

# #Make test plot for sucrose:
# fig.dat <- suc.dat
# suc.plot <- ggplot(fig.dat) +
#   scale_fill_manual(values = pal.fc.plots) +
#   geom_col(aes(x = as.factor(time), y = mean.log2.fc, fill = condition),  alpha = 0.6, width = 0.60, color = "gray40") +
#  # geom_col(data = suc.dat, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
#   geom_errorbar(data = fig.dat, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
#   facet_grid(.~treatment, scales = "free", space = "free") +
#   geom_vline(data=filter(fig.dat, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
#   geom_vline(data=filter(fig.dat, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
#   theme_classic() +
#   geom_hline(yintercept = 0) +
#   xlab("Time (hr)") +
#   ylab("Log2(FC) (Treatment/Control)") +
#   theme(legend.position = "none") 
# suc.plot

#make function:

fc_fig_plotter <- function(input.dat) {
  fig.dat <- input.dat
  fig <- ggplot(fig.dat) +
    scale_fill_manual(values = pal.fc.plots) +
    geom_col(aes(x = as.factor(time), y = mean.log2.fc, fill = condition),  alpha = 0.6, width = 0.60, color = "gray40") +
    # geom_col(data = suc.dat, aes(x = as.factor(time), y = mean.log2.fc, fill = High_Viral_Infection), alpha = 0.6, width = 0.60, color = "gray40") +
    geom_errorbar(data = fig.dat, aes(x = as.factor(time), ymin = mean.log2.fc-sd.log2.fc, ymax = mean.log2.fc+sd.log2.fc), width = 0.1, alpha = 0.6) +
    facet_grid(.~treatment, scales = "free", space = "free") +
    geom_vline(data=filter(fig.dat, treatment == "HV"), aes(xintercept = 1.5), color = "gray", linetype = "dashed") +
    geom_vline(data=filter(fig.dat, treatment == "LV"), aes(xintercept = 3.5), color = "gray", linetype = "dashed") +
    theme_classic() +
    geom_hline(yintercept = 0) +
    xlab("Time (hr)") +
    ylab("Log2(FC) (Treatment/Control)") +
    theme(legend.position = "none") +
    labs(title = fig.dat$Name)
  fig
}




#make all plots:
suc.fig <- fc_fig_plotter(suc.dat)
UDP.fig <- fc_fig_plotter(UDP.dat)
AMP.fig <- fc_fig_plotter(AMP.dat)
GBB.fig <- fc_fig_plotter(GBB.dat)
lys.fig <- fc_fig_plotter(lys.dat)
bet.fig <- fc_fig_plotter(bet.dat)
GG.fig <- fc_fig_plotter(GG.dat)
GSSG.fig <- fc_fig_plotter(GSSG.dat)
asp.fig <- fc_fig_plotter(asp.dat)
Ad.fig <- fc_fig_plotter(Ad.dat)
OA.fig <- fc_fig_plotter(OA.dat)
Des.fig <- fc_fig_plotter(Des.dat)
SAM.fig <- fc_fig_plotter(SAM.dat)
pb12.fig <- fc_fig_plotter(pb12.dat)

#combine all plots
comb.fig <- suc.fig + UDP.fig + AMP.fig + GBB.fig + lys.fig + bet.fig + GG.fig +
  GSSG.fig + asp.fig + Ad.fig + OA.fig + Des.fig + SAM.fig + pb12.fig +
  plot_layout(ncol = 3)
comb.fig

#save figure:
ggsave(comb.fig, file = "Figures/Outputs/Supplemental_Fold_Change_Figure.png", dpi = 600, 
       height = 10, width = 8, units = "in",
       bg = "white", scale = 1.5)






















