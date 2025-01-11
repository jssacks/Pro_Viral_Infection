


###Make figures for Cell size, Cellular carbon content, and infection rates
###C-content files
c.pro.file <- "Meta_Data/Abundance_Dat/Prochloro_size_abundance_rates_C_content.csv"
c.bact.file <- "Meta_Data/Abundance_Dat/Bacteria2pop_size_abundance_C_content.csv"


##Viral infection rate data:
inf.file <- "Meta_Data/Abundance_Dat/Modeled_Infection_Data_Combined.csv"





####Pro data
pro.dat <- read_csv(c.pro.file) %>%
  mutate(Pop.C = abundance*Qc) %>%
  filter(treatment %in% c("Control", "High virus", "Low virus")) %>%
  mutate(treatment = case_when(treatment == "Control" ~ "C",
                               treatment == "High virus" ~ "HV",
                               treatment == "Low virus" ~ "LV")) %>%
  mutate("pop" = "Pro") %>%
  select(-net, -loss) %>%
  filter(!experiment == "E")



pro.size.dat <- pro.dat %>%
  filter(experiment == "A")

#Pro Size fig:
ggplot(pro.size.dat, aes(x = time, y = Qc, fill = treatment)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5, shape = 21) +
  scale_fill_manual(values = c("aliceblue", "darkblue", "steelblue")) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  xlab("Time (hr)")
  
  
  geom_line(size = 2) +
  scale_fill_manual(values = c("aliceblue", "darkblue", "steelblue")) 




#Pro Carbon fig:
ggplot(pro.dat, aes(x = time, y = Pop.C, fill = treatment, shape = experiment)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("aliceblue", "darkblue", "steelblue")) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab("Prochlorococcus carbon (fgC/mL)") +
  xlab("Time (hr)")



###Infection dat
inf.dat <- read_csv(inf.file)

ggplot(inf.dat, aes(x = time, y = per_Inf, fill = treatment, shape = biorep)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("darkblue", "steelblue")) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  ylab("Modeled Percent Infection (%)")
