library(tidyverse)
library(viridis)
library(vegan)
library(ggdendro)


###
transcript.file <- "Collaborator_Data/Transcripts/MED4_annot_DE+non-DE_082823_normalized.csv"
transcript.file.2 <- "Collaborator_Data/Transcripts/HV_LFC+pval_MED4norm.csv"

transcript.dat <- read_csv(transcript.file)
transcript.dat.2 <- read_csv(transcript.file.2)


t.sml <- transcript.dat.2 %>%
  select(GENE.NAME, KEGG, HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, PathwayID, Pathway_info, Details) %>%
  pivot_longer(cols = HV0_LFC_nonDE_log2FoldChange:HV36_LFC_nonDE_padj, names_to = "Samp_Name", values_to = "val") %>%
  mutate(param = case_when(str_detect(Samp_Name, "log2FoldChange") ~ "FC",
                           str_detect(Samp_Name, "padj") ~ "p_adj")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_log2FoldChange")) %>%
  mutate(Samp_Name = str_remove(Samp_Name, "_LFC_nonDE_padj")) %>%
  mutate(Samp_Name = as.numeric(str_remove(Samp_Name, "HV"))) %>%
  mutate(GENE.NAME = str_remove(GENE.NAME, "MED4_")) 

##Assign significance levels to p-values 
t.sml.p <- t.sml %>%  
  filter(param == "p_adj") %>%
  mutate(signif = case_when(val <= 0.05 & val > 0.01 ~ "*",
                            val <= 0.01 & val > 0.005 ~ "**",
                            val <= 0.005 ~ "***",
                            TRUE ~ "")) %>%
  select(GENE.NAME, Samp_Name, signif)

t.sml.dat <- t.sml %>%
  filter(param == "FC") %>%
  left_join(., t.sml.p) %>%
  unique() %>%
  rename("FC" = val) %>%
  select(-param)



###### Sucrose and GG Transcripts:
gg.t.dat <- t.sml.dat %>%
  filter(KEGG == "K05978")

suc.t.dat <- t.sml.dat %>%
  filter(str_detect(Details, "sucrose")) %>% 
  rbind(gg.t.dat)

suc.t.plots <- ggplot(suc.t.dat, aes(x = as.factor(Samp_Name), y = GENE.NAME, fill = FC)) +
  geom_tile(color = "black", height = 0.4, width = 1) +
  theme_test() +
  geom_text(aes(label = signif), color = "white", size = 5.5, vjust = 0.8) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2)) +
  facet_wrap(.~GENE.NAME, scales = "free") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 14, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.spacing.x = unit(1, 'cm')) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  xlab("Time") +
  ylab("Gene Name") 

suc.t.plots


ggsave(suc.t.plots, filename = "Figures/Outputs/Sucrose_Transcripts_Separate.png",
       height = 5, width = 10, scale = 1.2, dpi = 1200)



###### Aspartic acid Transcripts:
asp.t.dat <- t.sml.dat %>%
  filter(str_detect(Details, "aspartate"))

asp.t.plots <- ggplot(asp.t.dat, aes(x = as.factor(Samp_Name), y = GENE.NAME, fill = FC)) +
  geom_tile(color = "black", height = 0.4, width = 1) +
  theme_test() +
  geom_text(aes(label = signif), color = "white", size = 5.5, vjust = 0.8) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2)) +
  facet_wrap(.~GENE.NAME, scales = "free") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 14, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.spacing.x = unit(1, 'cm')) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  xlab("Time") +
  ylab("Gene Name") 

asp.t.plots

ggsave(asp.t.plots, filename = "Figures/Outputs/Aspartate_Transcripts_Separate.png",
       height = 7, width = 10, scale = 1.2, dpi = 1200)

####### Antioxidant Transcripts:

gssg.t.dat <- t.sml.dat %>%
  filter(str_detect(Details, "Glutathione")) 

gssg.t.plots <- ggplot(gssg.t.dat, aes(x = as.factor(Samp_Name), y = GENE.NAME, fill = FC)) +
  geom_tile(color = "black", height = 0.4, width = 1) +
  theme_test() +
  geom_text(aes(label = signif), color = "white", size = 5.5, vjust = 0.8) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2)) +
  facet_wrap(.~GENE.NAME, scales = "free") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 14, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.spacing.x = unit(1, 'cm')) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  xlab("Time") +
  ylab("Gene Name") 

gssg.t.plots


ggsave(gssg.t.plots, filename = "Figures/Outputs/GSSG_Transcripts_Separate.png",
       height = 5, width = 10, scale = 1.2, dpi = 1200)



####### SAM Transcripts:
sam.t.dat <- t.sml.dat %>%
  filter(str_detect(Details, "methionine")) 

sam.t.plots <- ggplot(sam.t.dat, aes(x = as.factor(Samp_Name), y = GENE.NAME, fill = FC)) +
  geom_tile(color = "black", height = 0.4, width = 1) +
  theme_test() +
  geom_text(aes(label = signif), color = "white", size = 5.5, vjust = 0.8) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2)) +
  facet_wrap(.~GENE.NAME, scales = "free") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 14, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.spacing.x = unit(1, 'cm')) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  xlab("Time") +
  ylab("Gene Name") 

sam.t.plots

ggsave(sam.t.plots, filename = "Figures/Outputs/SAM_Transcripts_Separate.png",
       height = 5, width = 10, scale = 1.2, dpi = 1200)





####### pB12 Synthesis Transcripts: 
cob.t.dat <- t.sml.dat %>%
  filter(str_detect(Details, "Porphyrin")) %>%
  filter(GENE.NAME %in% c("hemA", "hemL", "hemB", "hemC", "hemD", "cobA",
                          "cobI", "cobG", "cobJ", "cobM", "cobF", "cobK", 
                          "cobL", "cobH", "cobB", "cobN", "CobS", "cobT", 
                          "cobO", "cobQ", "cobC", "cobD", "cobP", "cobU", 
                          "cobV", "cobO-1", "cysG"))

cob.t.clust.dat <- cob.t.dat %>%
  select(GENE.NAME, Samp_Name, FC) %>%
  pivot_wider(id_cols = GENE.NAME, names_from = Samp_Name, values_from = FC) %>%
  column_to_rownames(var = "GENE.NAME")

cob.dist <- vegdist(cob.t.clust.dat, method = "euclidean")

clust.out <- hclust(cob.dist, method = "average")
dend <- as.dendrogram(clust.out)
dend.dat <- dendro_data(dend)
dend.order <- dend.dat$labels %>%
  rename("order" = x, 
         "GENE.NAME" = label) %>%
  select(GENE.NAME, order)


cob.t.fig <- left_join(cob.t.dat, dend.order)

pb12.t.fig <- ggplot(cob.t.fig, aes(x = as.factor(Samp_Name), y = reorder(GENE.NAME, -order), fill = FC)) + 
  geom_tile(color = "black") +
  theme_minimal() +
  geom_text(aes(label = signif), color = "white", size = 5.5, vjust = 0.8) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2)) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  xlab("Time") +
  ylab("Gene Name")
pb12.t.fig

#export pb12 transcript plot:
ggsave(pb12.t.fig, file = "Figures/Outputs/cob_transcripts_heatmap.png",
       dpi = 800, units = "in", height = 7, width = 4, bg = "white")





