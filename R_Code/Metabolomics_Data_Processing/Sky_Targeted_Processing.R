





#library(here)
library(readr)
library(tidyverse)

#Load Functions
source("R_Code/Functions.R")

## This script identifies targeted compounds present in the pooled samples, determines their batch- and sample-specific mass,
## accounting for mass defects, and their batch- and sample-specific retention time. This information is then exported as a 
## tab-separated .txt file that can be imported into MS DIAL to identify these targeted compounds in the untargeted dataset.


# Define files -----------------------------------------------------------
# Define your Skyline output and transition list filepaths from the first set of manual integrations.

## HILIC Positive
HPos.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_File1.csv"
HPos.tlist <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Pos_TransitionList.csv"

## HILIC Negative
HNeg.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_File1.csv"
HNeg.tlist <- "Data_Raw/Skyline/ProMo_Particulate_HILIC_Neg_TransitionList.csv"

## Reverse Phase
RP.file.1 <- "Data_Raw/Skyline/ProMo_Particulate_RP_File1.csv"
RP.tlist <- "Data_Raw/Skyline/ProMo_Particulate_RP_TransitionList.csv"

## Min Area Threshold
min.area <- 40000

# Import defined files -----------------------------------------------------------

## HILIC Positive
HPos.data <- skyline_read(HPos.file.1)
HPos.tlist <- tlist_read(HPos.tlist) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

## HILIC Negative
HNeg.data <- skyline_read(HNeg.file.1)
HNeg.tlist <- tlist_read(HNeg.tlist) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

## Reverse Phase
RP.data <- skyline_read(RP.file.1)
RP.tlist <- tlist_read(RP.tlist) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

# Calculate exact masses -----------------------------------------------------------
## Combine pooled data and transition list information. Calculate batch-specific exact mass 
## by incorporating batch-specific mass defect.

HPos.exact <- calculate_exact_mass(HPos.data, HPos.tlist)
HNeg.exact <- calculate_exact_mass(HNeg.data, HNeg.tlist)
RP.exact <- calculate_exact_mass(RP.data, RP.tlist)

# Remove flagged compounds ------------------------------------------------
## Remove compounds that do not pass QC, or were not detected in all pooled samples.
## Change the count filter to equal the total number of pooled samples

HPos.qc.results <- QC_filter(HPos.data, pooled.count = 12, filter.compound = "Picolinic_Acid")
HNeg.qc.results <- QC_filter(HNeg.data, 12)
RP.qc.results <- QC_filter(RP.data, 8)

# Filter data by qc results for pooled sample compounds -------------------
HPos.data.export <- data_export(HPos.qc.results, HPos.exact, "[M+H]+")
HNeg.data.export <- data_export(HNeg.qc.results, HNeg.exact, "[M-H]-")
RP.data.export <- data_export(RP.qc.results, RP.exact, "[M+H]+")


# Export data -------------------------------------------------------------
write_tsv(HPos.data.export, file = "Intermediates/Sky_search_files/Particulate_HPos_skysearch.txt")

write_tsv(HNeg.data.export, file = "Intermediates/Sky_search_files/Particulate_HNeg_skysearch.txt")

write_tsv(RP.data.export, file = "Intermediates/Sky_search_files/Particulate_RP_skysearch.txt")
