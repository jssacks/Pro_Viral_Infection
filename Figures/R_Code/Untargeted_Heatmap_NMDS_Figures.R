





library(tidyverse)
source("Source_Code/biostats.R")
library(dendextend)
library(ggdendro)
library(vegan)
library(viridis)
library(ggpubr)

metab.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"


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
  rename("Normalized Area" = max.norm.area)  %>%
  mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) 


###Figure out signfiance of clusters through bootstrapping:

#hclus.scree(clust.out)
#install.packages("pvclust")
#library(pvclust)

#metab.clust.dat.trans <- metab.dat %>%
#  select(MF, SampID, max.norm.area) %>%
#  unique() %>%
#  pivot_wider(id_cols = SampID, names_from = MF, values_from = max.norm.area) %>%
#  column_to_rownames(var = "SampID")

#clus.stab <- pvclust(metab.clust.dat.trans, 
#                     method.hclust="average",
#                     method.dist="euclidean",
#                     nboot=50)
#plot(clus.stab) 
#pvrect(clus.stab, alpha = 0.90)

#library(NbClust)

#x <- NbClust(metab.clust.dat, distance = "euclidean", min.nc = 2, max.nc = 15, method = "average")
#x
#clust.numb <- as.data.frame(cutree(clust.out, k = 5)) %>%
#  rownames_to_column(var = "MF") %>%
#  rename("cluster" = "cutree(clust.out, k = 7)") 

#dend.order.2 <- left_join(dend.order, clust.numb)
#metab.clust.dat.2 <- cbind(metab.clust.dat, clust.numb) %>%
#  select()




###heatmap
hm.fig <- ggplot(metab.fig, aes(x = reorder(SampID, time), 
                                y = reorder(MF, order), 
                                fill = `Normalized Area`)) +
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
#  geom_hline(yintercept = 186.5, alpha = 0.5) +
#  geom_hline(yintercept = 567.5, alpha = 0.5)+ 
#  geom_hline(yintercept = 570.5, alpha = 0.5)+
#  geom_hline(yintercept = 611.5, alpha = 0.5)+
  facet_grid(.~treatment, scales = "free", space = "free") +
  xlab("Sample") +
  ylab("Mass Feature") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "aliceblue", mid = "steelblue3", high = "darkblue", midpoint = 0.5) +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 0.5, xmax = 3.5), fill = "black", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 3.5, xmax = 6.5), fill = "white", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 6.5, xmax = 9.5), fill = "black", color = "black") +
  geom_rect(mapping = aes(ymin = -7, ymax = 0, xmin = 9.5, xmax = 12.5), fill = "white", color = "black") +
  geom_rect(data = filter(metab.fig, treatment %in% c("C", "HV")), mapping = aes(ymin = -7, ymax = 0, xmin = 12.5, xmax = 15.5), fill = "black", color = "black")+
  geom_rect(data = filter(metab.fig, treatment %in% c("LV")), mapping = aes(ymin = -7, ymax = 0, xmin = 12.5, xmax = 14.5), fill = "black", color = "black")# +


hm.fig





# NMDS_Code ---------------------------------------------------------------
####################
#Define inputs
qc2.dat <- read_csv("Intermediates/ProMo_Particulate_QC2_data.csv")

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
         shape = guide_legend(override.aes = list(fill = "black")))
nmds.fig



# Combine_figures ---------------------------------------------------------

combine.plot <- ggarrange(nmds.fig, hm.fig,
                          nrow = 2, ncol = 1, heights = c(0.35, 0.65),
                          labels = c("A", "B"))

combine.plot









ggsave(filename = "Figures/Main_Text_Figure_2.png", dpi = 300, height = 10, width = 6, scale = 1.2)
