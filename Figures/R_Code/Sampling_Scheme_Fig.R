




library(tidyverse)

####load sample info
sample.dat <- read_csv("Meta_Data/Sampling_Scheme.csv")

###define order
x <- data.frame(
  measurement = c("Prochlorococcus", "Bacteria", "Phage", "Metabolites", "Transcripts", "ATP", "Cellular_Buoyancy"),
  order = c(1, 2, 3, 4, 5, 6, 7)
)

###make data frame longer
sample.long <- sample.dat %>%
  pivot_longer(cols = Prochlorococcus:Cellular_Buoyancy, names_to = "measurement", values_to = "time") %>%
  left_join(., x) %>%
 # mutate(measurment = case_when(measurement ))
  mutate(measurement = str_replace(measurement, "Cellular_Buoyancy", "Buoyancy"))

#install.packages("devtools")
#library(devtools)
#devtools::install_github("johannesbjork/LaCroixColoR")

#library(LaCroixColoR)

library(ggthemes)

###define day/night scheme
samp.scheme <- ggplot(sample.long, aes(x = time, y = reorder(measurement, -order))) +
  geom_rect(xmin = 14, xmax = 24, ymin = 0, ymax = 10, fill = "lightgray") +
  geom_rect(xmin = 38, xmax = 48, ymin = 0, ymax = 10, fill = "lightgray") +
  geom_point(
    aes(fill = measurement), 
    size = 3,
    shape = 21, color = "black") +
#  scale_fill_continuous(low = "aliceblue", mid = "steelblue", high = "darkblue") +
  scale_fill_tableau(palette = "Hue Circle") +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("Measurement") +
  theme(legend.position = "blank")

ggsave(samp.scheme, filename = "Figures/Figures/Samp_Scheme_10032024.png",
       dpi = 800, height = 3, width = 5, units = "in", scale = 1.1)

