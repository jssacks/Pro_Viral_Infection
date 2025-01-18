







#define inputs:
mf.file <-  "Intermediates/Final_Processed_Untargeted_Data.csv"
annotation.file <- "Intermediates/ProMo_Particulate_SIRIUS_MF_Annotations.csv"
quant.file <- "Intermediates/Final_Processed_Quantified_Data.csv"
metab.c.file <- "Intermediates/Metabolite_Carbon_Comparison.csv"
fc.file <- "Intermediates/fc_analysis_all_mfs.csv"


### Make All MF supplemental Table________________________

#load data:
mf.dat <- read_csv(mf.file) %>%
  mutate(vol.norm.area = Adjusted_Area/Vol.filt.mL) %>%
  select(SampID, Time, biorep, treatment, MF, Name, FinalBMIS, fraction, RT, mz, vol.norm.area, Predicted_Adduct_Ion) %>%
  unique() %>%
  rename("Rep" = biorep,
         "Treatment" = treatment,
         "Fraction" = fraction,
         "Volume_Normalized_Area" = vol.norm.area)

#Export Table
write_csv(mf.dat, file = "Tables/Outputs/Supplemental_All_MF_Table.csv")



### Make MF annotation supplemental table
anot.dat <- read_csv(annotation.file) %>%
  select(MF, Name, mz, RT, fraction, adduct, molecularFormula, `ClassyFire#all classifications`, `ClassyFire#most specific class Probability`, name,  ConfidenceScore) %>%
  rename("Fraction" = fraction,
         "SIRUS_adduct" = adduct,
         "SIRIUS_molecular_formula" = molecularFormula,
         "Classyfire_all_classifications" = `ClassyFire#all classifications`,
         "Classyfire_most_specific_classification_probability" = `ClassyFire#most specific class Probability`,
         "SIRIUS_predicted_identity" = name,
         "SIRIUS_identity_confidence_score" = ConfidenceScore) %>%
  unique()

#Export annotations:
write_csv(anot.dat, file = "Tables/Outputs/Supplemental_MF_Annotation_Table.csv")








### Make Quant Metab supplemental Table________________________
quant.dat <- read_csv(quant.file) %>%
  left_join(., read_csv(metab.c.file)) %>%
  select(SampID, time, biorep, treatment, MF, Name, nM.in.smp, nM_C, nM_N, tot.C.nM) %>%
  mutate(metab_C_perc = nM_C/tot.C.nM*100) %>%
  rename("Time" = time,
         "Rep" = biorep,
         "Treatment" = treatment,
         "Metab_nM_in_sample" = nM.in.smp,
         "Metab_nM_C_in_sample" = nM_C,
         "Metab_nM_N_in_sample" = nM_N,
         "Total_FCM_biomass_C_in_sample_nM" = tot.C.nM,
         "Percent_metab_C" = metab_C_perc)

#Export quant table:
write_csv(quant.dat, file = "Tables/Outputs/Supplemental_Metab_Quant_Table.csv")







####_Make Fold Change Analysis Table ___________________
fc.dat <- read_csv(fc.file) %>%
  mutate(bact_cor_determination = case_when(bact.cor.pearson > 0 & bactlm.p.adj < 0.05 ~ "Bact_correlated",
                                                  TRUE ~ NA),
         final_fc_determination = case_when(mean.log2.fc > 0.5 & is.na(bact_cor_determination) & signif == TRUE ~ "Increased",
                                            mean.log2.fc < 0.5 & is.na(bact_cor_determination) & signif == TRUE ~ "Decreased",
                                            TRUE ~ "Not Changed")) %>%
  rename("mean_log2_fc" = mean.log2.fc,
         "t_test_p" = p,
         "bact_cor_pearson" = bact.cor.pearson,
         "bact_cor_lm_p" = bactlm.p.adj) %>%
  select(MF, Name, mean_log2_fc, t_test_p, wilcox_p, bact_cor_pearson, bact_cor_lm_p, bact_cor_determination, final_fc_determination)


#export fold change table:
write_csv(fc.dat, file = "Tables/Outputs/FC_Analysis_Supplemental_Table.csv")

















































































































