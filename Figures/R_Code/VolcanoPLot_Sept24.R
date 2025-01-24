



library(tidyverse)
library(ggrepel)
library(scales)
library(ggpubr)
library(patchwork)


###define inputs
v.mf.file = "Intermediates/virally_altered_mfs.csv"
all.mf.file = "Intermediates/fc_analysis_all_mfs.csv"

all.mf.sum <- read_csv(v.mf.file) %>%
  filter(!Name == "Unknown")
#  mutate()
# ###Make figures for Cell size, Cellular carbon content, and infection rates
# ###C-content files
# c.pro.file <- "Meta_Data/Abundance_Dat/Prochloro_size_abundance_rates_C_content.csv"
# c.bact.file <- "Meta_Data/Abundance_Dat/Bacteria2pop_size_abundance_C_content.csv"
# 
# ####Pro carbon data
# pro.dat <- read_csv(c.pro.file) %>%
#   mutate(Pop.C = abundance*Qc) %>%
#   filter(treatment %in% c("Control", "High virus", "Low virus")) %>%
#   mutate(treatment = case_when(treatment == "Control" ~ "C",
#                                treatment == "High virus" ~ "HV",
#                                treatment == "Low virus" ~ "LV")) %>%
#   mutate("pop" = "Pro") %>%
#   select(-net, -loss) %>%
#   filter(!experiment == "E") %>%
#   mutate(treatment = fct_relevel(treatment, c("C", "LV", "HV"))) 
# 
# hvi.dat <- pro.dat %>%
#   mutate(keep = case_when(treatment == "HV" & time %in% c(12, 24) ~ "Keep",
#                           treatment == "LV" & time %in% c(36) ~ "Keep",
#                           TRUE ~ "FALSE")) %>%
#   filter(keep == "Keep")
#   
# 
# #Pro Carbon fig:
# Pro.C.fig <- ggplot(pro.dat, aes(x = time, y = Pop.C, fill = treatment, shape = experiment)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, alpha = 0.5) +
#   geom_point(data = hvi.dat, size = 2, aes(x = time, y = Pop.C, shape = experiment), stroke = 1, fill = NA, color = "red") +
#   scale_shape_manual(values = c(21, 22, 23)) +
#   scale_fill_manual(values = c("aliceblue", "steelblue", "darkblue")) +
#   theme_classic() +
#   guides(fill = guide_legend(override.aes = list(shape=21)),
#          shape = guide_legend(override.aes = list(color = "black", stroke = 0.5))) +
#   ylab("Prochlorococcus carbon (fgC/mL)") +
#   xlab("Time (hr)") +
#   theme(legend.position = "bottom")
# Pro.C.fig
# 
# 
# 
# 

  


########Make bar chart of differentially abundant metabolites:
barplot.dat <- read_csv(all.mf.file) %>%
  mutate(behavior = case_when(mean.log2.fc > 0.5 & signif == TRUE ~ "Increased",
                              mean.log2.fc < -0.5 & signif == TRUE ~ "Decreased",
                              TRUE ~ "Not Changed")) %>%
  mutate(behavior = as.factor(behavior)) %>%
  mutate(behavior = fct_relevel(behavior, c("Decreased", "Not Changed", "Increased"))) %>%
  mutate(sig.bact.cor = case_when(bact.cor.pearson > 0 & bactlm.p.adj < 0.05 ~ "Yes",
                                  TRUE ~ "No")) %>%
  mutate(sig.bact.cor = as.factor(sig.bact.cor)) %>%
  mutate(sig.bact.cor= fct_relevel(sig.bact.cor, c("Yes", "No"))) 

barplot.dat.sum <- barplot.dat %>%
  filter(!sig.bact.cor == "Yes") %>%
  group_by(behavior) %>%
  reframe(count = n())


count.plot <- ggplot(barplot.dat, aes(x = behavior, fill = sig.bact.cor)) + 
  geom_bar(width = 0.4, color = "black", size = 0.3) +
  scale_fill_manual(values = c("#cfd4d8", "#73a596")) +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0,370)) +
  labs(fill = "Correlated with\n Bacteia") +
  theme_test() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.title.x = element_blank())
count.plot




#######__________________ Make Volcano Plot_______________________
###Make datasets for volcano plot
sig.mfs <- read_csv(v.mf.file) %>%
  mutate(mean.log2.fc = case_when(mean.log2.fc > 5 ~ 5,
                                  mean.log2.fc < -5 ~ -5,
                                  TRUE ~ mean.log2.fc))


not.sig.mfs <- read_csv(all.mf.file) %>%
  filter(!MF %in% sig.mfs$MF) %>%
  mutate(mean.log2.fc = case_when(mean.log2.fc > 5 ~ 5,
                                  mean.log2.fc < -5 ~ -5,
                                  TRUE ~ mean.log2.fc)) %>%
  mutate(p = case_when(p < 3.1E-8 ~ 3.1E-8,
                       TRUE ~ p))


label.dat.neg <- sig.mfs %>%
  filter(!Name == "Unknown") %>%
  mutate(Name = case_when(Name == "Adenosine monophosphate" ~ "AMP",
                          Name == "S-Adenosylmethionine" ~ "SAM",
                          Name == "OH-Pseudocob" ~ "pB12",
                          Name == "Ophthalmic acid" ~ "OA",
                          Name == "Glutathione disulfide" ~ "GSSG",
                          Name == "UDP-N-acetylglucosamine" ~ "UDP-GlcNAc",
                          Name == "L-Aspartic acid" ~ "Aspartic acid",
                          TRUE ~ Name)) %>%
  filter(mean.log2.fc < 0)

label.dat.pos <- sig.mfs %>%
  filter(!Name == "Unknown") %>%
  mutate(Name = case_when(Name == "2-O-alpha-D-Glucosylglycerol" ~ "GG",
                          Name == "L-Lysine" ~ "Lysine",
                          Name == "(3-Carboxypropyl)trimethylammonium" ~ "gamma-BB",
                          TRUE ~ Name)) %>%
  filter(mean.log2.fc > 0) 

log2(1.25)
log2(0.75)
0.25

###plot volcano plot
volcano.plot <- ggplot() +
  geom_vline(xintercept = 0.5, size = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.5, size = 0.25, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), size = 0.25, linetype = "dashed") +
  geom_point(data = not.sig.mfs, aes(x = mean.log2.fc, y = -log10(p)), shape = 21, size = 3.5, fill = "darkgray", color = "black", stroke = 0.1, alpha = 0.65) +
  geom_point(data = sig.mfs, aes(x = mean.log2.fc, y = -log10(p), fill = mean.log2.fc), size = 3.5, color = "black", stroke = 0.1, shape = 21) +
  geom_label_repel(data = label.dat.pos, aes(x = mean.log2.fc, y = -log10(p), label = Name),
                   label.size = 0.1,
                   size = 3,
                   min.segment.length = 0.01, 
                   point.padding = 0, 
                   box.padding = 1,
                   segment.size = 0.25,
                   force = 5, 
                 #  fontface = "bold",
                   xlim = c(2, NA),
                   max.overlaps = 40) +
  geom_label_repel(data = label.dat.neg, aes(x = mean.log2.fc, y = -log10(p), label = Name),
                   label.size = 0.1,
                   size = 3,
                   min.segment.length = 0.01, 
                   point.padding = 0, 
                   box.padding = 1,
                   force = 5, 
                 #  fontface = "bold",
                   xlim = c(NA, -1.5),
                   segment.size = 0.25,
                   segment.angle = 30,
                   max.overlaps = 40) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2.5, 2.5), oob = squish) +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  xlab("Mean log2FC (HVI/C)") +
  ylab("-log(p)") +
  labs(fill = "Log2FC") +
  theme_test()

volcano.plot 

#ggsave(volcano.plot, filename = "Figures/Figures/Metab_Volcano_Plot.png", height = 5, width = 7, units = "in",
#       scale = 1, dpi = 600)

###############


####Arrange Plot:
comb.plot <- volcano.plot + inset_element(count.plot, left = 0.01, bottom = 0.71, right = 0.31, top = 0.98) + 
  plot_layout(guides = "collect")
comb.plot

#save plot
ggsave(comb.plot, filename = "Figures/Outputs/Volcano_Inset_Plot.png", 
       dpi = 800, bg = "white", scale = 1.3,
       units = "in", height = 6, width = 7)
#





#Arrange everything
group.plot <-  ggarrange(
  NA,
  ggarrange(NA, Pro.C.fig, NA, count.plot, NA,
            nrow = 1, ncol = 5, widths = c(0.04, 0.44, 0.03, 0.46, 0.02),
            labels = c("A", NA, NA, "B", NA)), 
  NA,
  volcano.plot, ncol = 1, nrow = 4, heights = c(0.03, 0.32, 0.03, 0.62),
  labels = c(NA, NA, NA, "C")
)
group.plot

ggsave(group.plot, filename = "Figures/Figures/Volcano_Full_Group_Plot.png", 
       dpi = 800, bg = "white", scale = 1.4,
       units = "in", height = 8, width = 7)






