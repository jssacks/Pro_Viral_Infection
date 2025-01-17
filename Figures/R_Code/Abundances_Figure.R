library(tidyverse)
library(readr)
library(patchwork)
library(ggpubr)
#
#
#Define inputs
abu.file <- "Collaborator_Data/FCM/ProMo_Abundat_Combined.csv"
pal <- c("#cfd4d8", "#73a596", "#365759")





#_______Organize Data__________________________

#Pull in abundance Data
abu.dat <- read_csv(abu.file) %>%
  separate(col = "Rep", into = c("Rep", "Treat.number")) %>%
  mutate(Treatment = case_when(Rep %in% c("A", "B", "C") & Treat.number == 1 ~ "C",
                               Rep %in% c("A", "B", "C") & Treat.number == 3 ~ "LV",
                               Rep %in% c("A", "B", "C") & Treat.number == 4 ~ "LGV", 
                               Rep %in% c("A", "B", "C") & Treat.number == 6 ~ "HV",
                               Rep %in% c("E") ~ "LG")) %>%
  filter(Treatment %in% c("C", "LV", "HV")) %>%
  mutate(Treatment = fct_relevel(Treatment, c("C", "LV", "HV"))) %>%
  group_by(Treatment, Time) %>%
  reframe(Pro.avg = mean(Pro, na.rm = TRUE),
          Pro.sd = sd(Pro, na.rm = TRUE),
          Bact.tot.avg = mean(BacTot, na.rm = TRUE),
          Bact.tot.sd = sd(BacTot, na.rm = TRUE),
          Phage.avg = mean(Phage, na.rm = TRUE),
          Phage.sd = sd(Phage, na.rm = TRUE))


hvi.dat <- abu.dat %>%
  mutate(keep = case_when(Treatment == "HV" & Time %in% c(12, 24) ~ "Keep",
                          Treatment == "LV" & Time %in% c(36) ~ "Keep",
                          TRUE ~ "FALSE")) %>%
  filter(keep == "Keep") %>%
  mutate(Condition = "High viral infection treatment-timepoint")
  
  



#_______Make component figures__________________________

##Pro fig
pro.fig <- ggplot(abu.dat, aes(x = Time, y = Pro.avg, fill = Treatment, shape = Treatment)) +
  geom_errorbar(aes(ymax = Pro.avg+Pro.sd, ymin = Pro.avg-Pro.sd), size = 0.25, width = 0.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = pal) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  geom_point(data = hvi.dat, stroke = 0.75, fill = "#dd9914", color = "black", size = 2.5) +
  theme_classic() +
#  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab(expression(Prochlorococcus~(cells~mL^-1))) +
  xlab("Time (h)") + 
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))  +
  theme(legend.position = "none") 
  
pro.fig





##Bact fig
bact.fig <- ggplot(abu.dat%>%filter(!is.na(Bact.tot.avg)), aes(x = Time, y = Bact.tot.avg, fill = Treatment, shape = Treatment)) +
  geom_errorbar(aes(ymax = Bact.tot.avg+Bact.tot.sd, ymin = Bact.tot.avg-Bact.tot.sd), size = 0.25, width = 0.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = pal) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  
  geom_point(data = hvi.dat, stroke = 0.75, fill = "#dd9914", color = "black", size = 2.5) +
  theme_classic() +
  #  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab(expression(Heterotrophic~bacteria~(cells~mL^-1))) +
  xlab("Time (h)") + 
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  theme(legend.position = "none") 
bact.fig






##Bact fig
phage.fig <- ggplot(abu.dat%>%filter(!is.na(Phage.avg)), aes(x = Time, y = Phage.avg, fill = Treatment, shape = Treatment)) +
  geom_errorbar(aes(ymax = Phage.avg+Phage.sd, ymin = Phage.avg-Phage.sd), size = 0.25, width = 0.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = pal) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  theme_classic() +
  geom_point(data = hvi.dat, stroke = 0.75, fill = "#dd9914", color = "black", size = 2.5) +
  #  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab(expression(Phage~(copies~mL^-1))) +
  xlab("Time (h)") + 
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) +
  theme(legend.position = "none") 
phage.fig



##Make legends:

#normal legend
pro.fig.legend <- ggplot(abu.dat, aes(x = Time, y = Pro.avg, fill = Treatment, shape = Treatment)) +
  geom_errorbar(aes(ymax = Pro.avg+Pro.sd, ymin = Pro.avg-Pro.sd), size = 0.25, width = 0.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = pal) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  theme_classic() +
  theme(legend.position = "bottom") 
pro.fig.legend

legend.normal <- get_legend(pro.fig.legend)


#hvi legend
hvi.fig.legend <- ggplot(hvi.dat, aes(x = Time, y = Pro.avg, fill = Condition)) +
  geom_point(stroke = 0.75, color = "black", shape = 23, size = 2.5) +
  scale_fill_manual(values = c("#dd9914")) + 
  theme_classic() +
  theme(legend.position = "bottom") 

hvi.fig.legend

legend.hvi <- get_legend(hvi.fig.legend)







###put plots together:
combined.fig <- (pro.fig | plot_spacer() | bact.fig | plot_spacer() | phage.fig) +
  plot_layout(guides = 'collect', widths = c(1, 0.1, 1, 0.1, 1)) +
  theme(legend.position = "none") 
combined.fig

#export plots:
ggsave(combined.fig, file = "Figures/Outputs/Abundances_Fig.png", dpi = 800, scale = 1.3,
          units = "in", height = 2.5, width = 7, bg = "white")

#export normal legend:
ggsave(legend.normal, file = "Figures/Outputs/Abundances_Legend_1.png", dpi = 800, scale = 1,
       units = "in", height = 3, width = 4, bg = "white")

#export HVI legend:
ggsave(legend.hvi, file = "Figures/Outputs/Abundances_Legend_2.png", dpi = 800, scale = 1,
       units = "in", height = 3, width = 4, bg = "white")
