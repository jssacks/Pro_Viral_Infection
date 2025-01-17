




library(tidyverse)

####load sample info
sample.dat <- read_csv("Meta_Data/Sampling_Scheme.csv")


#define palette
pal.full <- c("#73a596", "#cfd4d8", "#dd9914", "#365759", "#4cb8da", "#465f61")



###define order
x <- data.frame(
  measurement = c("Prochlorococcus", "Bacteria", "Phage", "Metabolites", "Transcripts", "ATP"),
  order = c(1, 2, 3, 4, 5, 6, 7)
)


###make data frame longer
sample.long <- sample.dat %>%
  pivot_longer(cols = Prochlorococcus:Cellular_Buoyancy, names_to = "measurement", values_to = "time") %>%
  filter(!measurement == "Cellular_Buoyancy") %>%
  left_join(., x) %>%
  filter(!is.na(time))


metab.control.samps <- tibble(
  measurement = c("Metabolites", "Metabolites", "Metabolites", "Metabolites"),
  time = c(6, 18, 30, 42),
  order = c(4, 4, 4, 4)
)
  


#install.packages("devtools")
#library(devtools)
#devtools::install_github("johannesbjork/LaCroixColoR")

#library(LaCroixColoR)

library(ggthemes)

###define day/night scheme
samp.scheme <- ggplot(sample.long, aes(x = time, y = reorder(measurement, -order))) +
  geom_rect(xmin = 14, xmax = 24, ymin = 0, ymax = 10, fill = "black", alpha = 0.005) +
  geom_rect(xmin = 38, xmax = 48, ymin = 0, ymax = 10, fill = "black", alpha = 0.005) +
  geom_point(
    aes(fill = reorder(measurement, -order)), 
    size = 3,
    shape = 21, color = "black") +
  geom_point(data = metab.control.samps, aes(fill = reorder(measurement, -order)),
             size = 3, shape = 23, color = "black") +
  #  scale_fill_continuous(low = "aliceblue", mid = "steelblue", high = "darkblue") +
  scale_fill_manual(values = pal.full, ) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("Measurement") +
  theme(legend.position = "blank")
samp.scheme

ggsave(samp.scheme, filename = "Figures/Outputs/Samp_Scheme.png",
       dpi = 800, height = 3, width = 5, units = "in", scale = 1.2)

