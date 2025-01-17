





library(tidyverse)
source("Source_Code/biostats.R")
library(dendextend)
library(ggdendro)

#load data
metab.file <- "Intermediates/Final_Processed_Untargeted_Data.csv"


#organize data
metab.dat <- read_csv(metab.file) %>%
  filter(!MF == "Vit_OH-pB12") %>%
  filter(!SampID == "C_3LV_T48") %>%
  mutate(vol.norm.area = Adjusted_Area/Vol.filt.mL) %>%
  select(MF, Name, Time, biorep, treatment, SampID, vol.norm.area) %>%
  filter(treatment %in% c("C", "LV", "HV")) %>%
  rename("time" = Time) %>%
  group_by(MF) %>%
  mutate(max.norm.area = vol.norm.area/max(vol.norm.area)) 

metab.dat.check <- metab.dat %>%
  select(MF) %>%
  unique()


####Heirarchical clustering:
metab.clust.dat <- metab.dat %>%
  select(MF, SampID, max.norm.area) %>%
  unique() %>%
  pivot_wider(id_cols = MF, names_from = SampID, values_from = max.norm.area) %>%
  column_to_rownames(var = "MF")

library(vegan)

metab.dist <- vegdist(metab.clust.dat, method = "euclidean")

clust.out <- hclust(metab.dist, method = "average")
dend <- as.dendrogram(clust.out)
dend.dat <- dendro_data(dend)
dend.order <- dend.dat$labels %>%
  rename("order" = x, 
         "MF" = label) %>%
  select(MF, order)

metab.fig <- left_join(metab.dat, dend.order) %>%
  filter(time %in% c(0, 12, 24, 36, 48)) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV")))



###make heatmap
heatmap.fig <- ggplot(metab.fig, aes(x = reorder(SampID, time), y = reorder(MF, order), fill = max.norm.area)) +
  geom_tile() +
  geom_vline(xintercept = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 3.5, alpha = 0.5) +
  geom_vline(xintercept = 6.5, alpha = 0.5) +
  geom_vline(xintercept = 9.5, alpha = 0.5) +
  geom_vline(xintercept = 12.5, alpha = 0.5) +
  geom_vline(xintercept = 15.5, alpha = 0.5) +
  geom_vline(xintercept = 18.5, alpha = 0.5) +
  geom_vline(xintercept = 21.5, alpha = 0.5) +
  geom_vline(xintercept = 24.5, alpha = 0.5) +
  geom_vline(xintercept = 27.5, alpha = 0.5) +
  facet_grid(.~treatment, scales = "free", space = "free") +
  xlab("Sample") +
  ylab("Mass Feature") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "white", mid = "#73a596", high = "#365759", midpoint = 0.5) +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 0.5, xmax = 3.5), fill = "white", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 3.5, xmax = 6.5), fill = "black", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 6.5, xmax = 9.5), fill = "white", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 9.5, xmax = 12.5), fill = "black", color = "black") +
  geom_rect(data = filter(metab.fig, treatment %in% c("C", "HV")), mapping = aes(ymin = -7, ymax = 0, xmin = 12.5, xmax = 15.5), fill = "white", color = "black")+
  geom_rect(data = filter(metab.fig, treatment %in% c("LV")), mapping = aes(ymin = -7, ymax = 0, xmin = 12.5, xmax = 14.5), fill = "white", color = "black")

heatmap.fig

#save heatmap
ggsave(heatmap.fig, filename = "Figures/Outputs/heatmap.png", dpi = 800, scale = 1.2,
       units = "in", height = 8, width = 7, bg = "white")


