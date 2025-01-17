


library(tidyverse)


###define inputs
v.mf.file = "Intermediates/virally_altered_mfs.csv"
all.mf.file = "Intermediates/fc_analysis_all_mfs.csv"


###Make datasets for volcano plot
sig.mfs <- read_csv(v.mf.file) %>%
  mutate(mean.log2.fc = case_when(mean.log2.fc > 5 ~ 5,
                                  mean.log2.fc < -5 ~ -5,
                                  TRUE ~ mean.log2.fc)) %>%
  select(-SIRIUS_molecularFormula, -SIRIUS_name)


not.sig.mfs <- read_csv(all.mf.file) %>%
  filter(!MF %in% sig.mfs$MF) %>%
  mutate(mean.log2.fc = case_when(mean.log2.fc > 5 ~ 5,
                                  mean.log2.fc < -5 ~ -5,
                                  TRUE ~ mean.log2.fc)) 


  
  
####Fix names and select only metabs I want:
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
                          TRUE ~ Name)) %>%
  filter(mean.log2.fc > 0) 

label.dat.other <- not.sig.mfs %>%
  filter(!Name == "Unknown") %>%
  filter(Name %in% c("L-Asparagine", "L-Methionine", "ATP",
                     "S-Adenosylhomocysteine", "Methylthioadenosine")) %>%
  mutate(Name = case_when(Name == "L-Asparagine" ~ "Asparagine",
                          Name == "L-Methionine" ~ "Methionine",
                          Name == "Methylthioadenosine" ~ "MTA",
                          Name == "S-Adenosylhomocysteine" ~ "SAH",
                          TRUE ~ Name)) 





####


###Combine 
known.dat <- rbind(label.dat.neg, label.dat.pos, other.metabs.dat) %>%
  mutate(sig.level = case_when(p <= 0.05 & p > 0.01 & abs(mean.log2.fc) > 0.5 ~ "*",
                            p <= 0.01 & p > 0.005 & abs(mean.log2.fc) > 0.5 ~ "**",
                            p <= 0.005 & abs(mean.log2.fc) > 0.5 ~ "***",
                            TRUE ~ "")) 


############
circle.metab.plots <- ggplot(known.dat, aes(x = 1, y = Name)) +
  geom_point(aes(fill = mean.log2.fc), shape = 22, size = 12, stroke = 0.2) +
  geom_text(aes(label = sig.level), color = "black", size = 5.5, vjust = 0.8) +
  theme_test() +
  facet_wrap(.~Name, scales = "free") +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  theme(legend.key.height= unit(1.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) +
  scale_fill_steps2(low = "steelblue3", mid = "white", high = "red2", n.breaks = 12, limits = c(-2, 2))

circle.metab.plots

ggsave(circle.metab.plots, filename = "Figures/Outputs/metab_metabolic_maps_circles.png",
       height = 5, width = 10, scale = 1.2, dpi = 1200)



