###
#
#
#
#
#Functions for ProMo Analysis
#
#
library(tidyverse)
library(readr)
#
#
#

#Functions for Sky_Targeted_Processing.R Script

# Define functions --------------------------------------------------------

calculate_exact_mass <- function(skyline.df, tlist) {
  # Calculate batch-specific exact mass. 
  #
  # Args
  #   skyline.df: Skyline dataframe
  #   tlist: Skyline transition list
  #
  exact.mass.df <- left_join(skyline.df, tlist) %>%
    filter(str_detect(.$Rep, "Poo")) %>%
    mutate(Exact_Mass = Mass + (Mass_Error_PPM/10^6)*Mass) %>%
    group_by(Compound) %>%
    summarize(Mean.Exact.Mass = mean(Exact_Mass), 
              Mean.RT = mean(RT)) %>%
    rename("Exact_Mass" = Mean.Exact.Mass,
           "RT" = Mean.RT)
}

data_export <- function(qc.results, df.exact, adduct) {
  # Calculate batch-specific exact mass. 
  #
  # Args
  #   qc.results: Skyline dataframe filtered for compounds.
  #   df.exact: Skyline dataframe with exact masses.
  #   adduct: User-entered string of correct adduct; must be in "[M+H]+" format.
  #
  export <- left_join(qc.results, df.exact) %>%
    select(Compound, RT, Exact_Mass) %>%
    mutate("Adduct" = adduct) %>%
    rename("Name" = Compound,
           "Pos-m/z" = Exact_Mass) %>%
    select(Name, `Pos-m/z`, RT, Adduct)
}

skyline_read <- function(file) {
  # Read and import files from Skyline.
  output <- read_csv(file) %>%
    mutate(Rep =`Replicate Name`,
           Compound = `Precursor Ion Name`,
           RT = as.numeric(`Retention Time`),
           Area = as.numeric(`Area`),
           Background = as.numeric(`Background`),
           Height = as.numeric(`Height`),
           Mass_Error_PPM = as.numeric(`Mass Error PPM`)
    ) %>%
    select(Rep, Compound, RT, Area, Background, Height, Mass_Error_PPM)
  print(output)
}

QC_filter <- function(skyline.df, pooled.count, filter.compound) {
  # Apply QC filter to imported Skyline file, filter user-defined 
  # compounds that don't pass QC, and change the count filter to equal 
  # the total number of pooled samples.
  qc.results <- skyline.df %>%
    filter(str_detect(.$Rep, "Poo")) %>%
    filter(!is.na(Area)) %>%
    filter(!Area < min.area) %>%
    group_by(Compound) %>%
    summarise(count = n()) %>%
    filter(count == pooled.count) %>%
    select(Compound)
  
  if(missing(filter.compound)) {
    
    return(qc.results)
  } else {
    qc.results.cmpd <- qc.results %>%
      filter(!Compound == filter.compound)
    return(qc.results.cmpd)
  }
}

tlist_read <- function(file) {
  # Read and import Skyline transition lists.
  output <- read_csv(file, col_names = FALSE) %>%
    mutate(Mass = as.numeric(X2),
           Compound = X4) %>%
    select(Compound, Mass)
}




#Functions for data processing:

#####Define Skyline read data read in function
sky_read <- function(file) {
  output <- read_csv(file) %>%
    mutate(Rep =`Replicate Name`,
           Compound = `Precursor Ion Name`,
           RT = as.numeric(`Retention Time`),
           Area = as.numeric(`Area`),
           Background = as.numeric(`Background`),
           Height = as.numeric(`Height`),
           Mass_Error_PPM = as.numeric(`Mass Error PPM`)
    ) %>%
    select(Rep, Compound, RT, Area, Background, Height, Mass_Error_PPM)
  print(output)
}

####Define Skyline transition list read in function
tl_read <- function(file) {
  output <- read_csv(file, col_names = FALSE) %>%
    mutate(Mass = as.numeric(X2),
           Compound = X4) %>%
    select(Compound, Mass)
}

###Define function to read in MS-DIAL Data
#mode should be the analytical fraction ("HILICPos", "HILICNeg", or "RP")

MSDIAL_read <- function(file1, Mode) {
  read_delim(file1,
             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = Mode) %>% 
  mutate("ID" = as.character(`Alignment ID`)) %>%
  rename("RT" = `Average Rt(min)`,
         'mz' = `Average Mz`,
         "Name" = `Metabolite name`, 
         "Adduct" = `Adduct type`,
         "Note" = `Post curation result`,
         "Fill" = `Fill %`,
         "MS2" = `MS/MS assigned`,
         "Ref_RT" = `Reference RT`,
         "Ref_mz" = `Reference m/z`,
         "RT_matched" = `RT matched`,
         "mz_matched" = `m/z matched`,
         "MS2_matched" = `MS/MS matched`,
         "SN_ave" = `S/N average`)
}




####Adduct, isotope, and pos/neg match finder function
###adduct finding function
find_adducts <- function(ID_num, dataset, adduct_list) {
  
  dataset.2 <- dataset %>%
    mutate(merge.key = "x")  
  
  MF.limits <- dataset %>%
    filter(ID == ID_num) %>%
    mutate(RT_low = RT - RT_ad_tol,
           RT_high = RT + RT_ad_tol) 
  ad.limits <- cbind(MF.limits, adduct_list) %>%
    mutate(adduct.mass = ((mz) - (1.007276*polarity))/abs(Charge) + Mass.change,
           adduct.mass.high = adduct.mass+(adduct.error.ppm/10^6)*adduct.mass,
           adduct.mass.low = adduct.mass-(adduct.error.ppm/10^6)*adduct.mass) %>%
    mutate(merge.key = "x") %>%
    select(merge.key, Ion, RT_low, RT_high, adduct.mass.low, adduct.mass.high, ad.polarity) %>%
    rename("polarity" = ad.polarity)
  
  ad.detect <- full_join(dataset.2, ad.limits) %>%
    select(!merge.key) %>%
    rowwise() %>%
    filter(RT >= RT_low & RT <= RT_high) %>%
    filter(mz >= adduct.mass.low & mz <= adduct.mass.high) %>%
    rename("adduct.ID" = ID) %>%
    mutate(ID = ID_num) 
  
  print(ad.detect)
}

###Function
find_sirius_match <- function(MF.file, SIRIUS.file, RT_tol, search.error.ppm) {
  
  MF.limits <- MF.file %>%
    mutate(RT_low = RT.seconds - RT_tol,
           RT_high = RT.seconds + RT_tol) %>%
    mutate(search.mass = mz,
           search.mass.high = search.mass+(search.error.ppm/10^6)*search.mass,
           search.mass.low = search.mass-(search.error.ppm/10^6)*search.mass) %>%
    mutate(match.key = "x")
  
  sirius.MFs <- SIRIUS.file %>%
    select(id, ionMass, retentionTimeInSeconds) %>%
    mutate(match.key = "x")
  
  MF.match <- full_join(MF.limits, sirius.MFs) %>%
    rowwise() %>%
    filter(retentionTimeInSeconds >= RT_low & retentionTimeInSeconds <= RT_high) %>%
    filter(ionMass >= search.mass.low & ionMass <= search.mass.high) 
  
  MF.annotate.full <- MF.match %>%
    select(MF, Name, mz, RT, RT.seconds, id) %>%
    left_join(., SIRIUS.file)
  
  print(MF.annotate.full)
  
}


#################Functions for interfacing with KEGG 

#function to get all reaction IDs for a set of compound IDs
get_reactionIDs <- function(CompID) {
dat.react <- keggGet(CompID)
df.react <- data.frame(x = print(dat.react[[1]]$REACTION)) %>%
  separate(x, into = c("1", "2", "3", "4", "5", "6", "7", "8"), sep = " ") %>%
  mutate(Comp_ID = CompID) %>%
  pivot_longer(cols = "1":"8", names_to = "numb", values_to = "Reaction_ID") %>%
  select(-numb) %>%
  filter(!is.na(Reaction_ID))
df.react
}



#function to get gene IDs (KOs) for a set of reaction IDs 
get_gene_IDs <- function(Reaction_ID) {
  kegg.out <- keggGet(Reaction_ID)
  ko.df <- data.frame(gene_name = kegg.out[[1]]$ORTHOLOGY) %>%
    rownames_to_column(var = "KEGG_gene_ID") %>%
    mutate(Reaction_ID = Reaction_ID) 
  ko.df
}










