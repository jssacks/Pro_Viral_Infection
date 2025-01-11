

library(tidyverse)

####Source Functions
source("R_Code/Functions.R")


###inputs

#QC2 MF Output      
MF.file <- "Intermediates/ProMo_Particulate_QC2_data.csv"    

#SIRIUS HILIC Pos Outputs
hpos.canopus.comp.file <- "Data_Raw/SIRIUS/HILIC_Pos/canopus_compound_summary.tsv"
hpos.compound.id.file <- "Data_Raw/SIRIUS/HILIC_Pos/compound_identifications.tsv"
hpos.formula.id.file <- "Data_Raw/SIRIUS/HILIC_Pos/formula_identifications.tsv"

#SIRIUS HILIC Neg Outputs
hneg.canopus.comp.file <- "Data_Raw/SIRIUS/HILIC_Neg/canopus_compound_summary.tsv"
hneg.compound.id.file <- "Data_Raw/SIRIUS/HILIC_Neg/compound_identifications.tsv"
hneg.formula.id.file <- "Data_Raw/SIRIUS/HILIC_Neg/formula_identifications.tsv"


#SIRIUS RP Outputs
RP.canopus.comp.file <- "Data_Raw/SIRIUS/RP/canopus_compound_summary.tsv"
RP.compound.id.file <- "Data_Raw/SIRIUS/RP/compound_identifications.tsv"
RP.formula.id.file <- "Data_Raw/SIRIUS/RP/formula_identifications.tsv"


####Define mass feature matching parameters
RT.tol.val <- 20
search.error.ppm.val <- 5





###MFs 
MF.dat <- read_csv(MF.file) %>%
  filter(MS2 == "TRUE") %>%
  select(MF, RT, mz, Name, MS2, fraction) %>%
  unique() %>%
  mutate(RT.seconds = RT*60)


HPos.dat <- MF.dat %>%
  filter(fraction == "HILIC_Pos")

HNeg.dat <- MF.dat %>%
  filter(fraction == "HILIC_Neg")

RP.dat <- MF.dat %>%
  filter(fraction == "RP")





##Read in and Combine HPos Sirius Output

#Read in files
hpos.c.comp.dat <- read_tsv(hpos.canopus.comp.file)
hpos.comp.id.dat <- read_tsv(hpos.compound.id.file)
hpos.form.id.dat <- read_tsv(hpos.formula.id.file)

##combine into one file containing Predicted Formula, Classification, and Identity 
hpos.s.dat.1 <- full_join(hpos.c.comp.dat, hpos.comp.id.dat) %>%
  select(-rank, -molecularFormula, -ionMass, -retentionTimeInSeconds, -SiriusScore, -adduct)
hpos.s.dat.2 <- hpos.form.id.dat %>%
  select(id, adduct, ionMass, retentionTimeInSeconds, molecularFormula, precursorFormula, SiriusScore, TreeScore, IsotopeScore)

hpos.sirius.dat <- left_join(hpos.s.dat.2, hpos.s.dat.1)

#________________________________
##Read in and Combine HNeg Sirius Output

#Read in files
hneg.c.comp.dat <- read_tsv(hneg.canopus.comp.file)
hneg.comp.id.dat <- read_tsv(hneg.compound.id.file)
hneg.form.id.dat <- read_tsv(hneg.formula.id.file)

##combine into one file containing Predicted Formula, Classification, and Identity 
hneg.s.dat.1 <- full_join(hneg.c.comp.dat, hneg.comp.id.dat) %>%
  select(-rank, -molecularFormula, -ionMass, -retentionTimeInSeconds, -SiriusScore, -adduct)
hneg.s.dat.2 <- hneg.form.id.dat %>%
  select(id, adduct, ionMass, retentionTimeInSeconds, molecularFormula, precursorFormula, SiriusScore, TreeScore, IsotopeScore)

hneg.sirius.dat <- left_join(hneg.s.dat.2, hneg.s.dat.1)


#________________
##Read in and Combine RP Sirius Output

#Read in files
RP.c.comp.dat <- read_tsv(RP.canopus.comp.file)
RP.comp.id.dat <- read_tsv(RP.compound.id.file)
RP.form.id.dat <- read_tsv(RP.formula.id.file)

##combine into one file containing Predicted Formula, Classification, and Identity 
RP.s.dat.1 <- full_join(RP.c.comp.dat, RP.comp.id.dat) %>%
  select(-rank, -molecularFormula, -ionMass, -retentionTimeInSeconds, -SiriusScore, -adduct)
RP.s.dat.2 <- RP.form.id.dat %>%
  select(id, adduct, ionMass, retentionTimeInSeconds, molecularFormula, precursorFormula, SiriusScore, TreeScore, IsotopeScore)

RP.sirius.dat <- left_join(RP.s.dat.2, RP.s.dat.1)

#####



#####Use find_sirius_match() function to identify matches between SIRIUS annotations and quality controlled MFs,
# in some cases, there may be more than one annotations per MF 

hpos.match.out <- find_sirius_match(HPos.dat, hpos.sirius.dat, RT.tol.val, search.error.ppm.val) %>%
  mutate(fraction = "HPos")
hneg.match.out <- find_sirius_match(HNeg.dat, hneg.sirius.dat, RT.tol.val, search.error.ppm.val)%>%
  mutate(fraction = "HNeg")
rp.match.out <- find_sirius_match(RP.dat, RP.sirius.dat, RT.tol.val, search.error.ppm.val)%>%
  mutate(fraction = "RP")


#Combine All 3 fractions into a single dataframe and clean up for export
all.match.dat <- rbind(hpos.match.out, hneg.match.out, rp.match.out)

write_csv(all.match.dat, file = "Intermediates/ProMo_Particulate_SIRIUS_MF_Annotations.csv")



##summarize number of MFs annotated with formulas and/or IDs and/or classifications
all.match.sum <- all.match.dat %>%
  select(MF, fraction) %>%
  unique() %>%
  group_by(fraction) %>%
  summarize(count = n())
  



















