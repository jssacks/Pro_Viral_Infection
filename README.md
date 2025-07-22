# Pro_Viral_Infection

# Pro_Viral_Infection

Remodeling of Prochlorococcus metabolism during viral infection

A repository for the data and code for Sacks et al. (in prep). Written in R. 

___ 
# Repo structure:

# Data:
## Folder 1: Data_Raw: 
* This folder contains the integrated metabolomics data output files from Skyline and MSDIAL, and the predicted annotation outputs from SIRIUS 5 in three separate folders. Each folder is separated into the three metabolomics analytical fractions: HILIC positive (Pos), HILIC Negative (Neg), and Reverse Phase (RP) outputs. 

### Skyline
* This folder contains two skyline output .csv files that collectively contain all samples and compounds integrated in skyline for each analytical fraction (HILIC_Pos, HILIC_Neg, RP). The transition list (m/z and retention times for each compound) used to search for compounds in Skyline is also provided for each analytical fraction. Finally, the skyline output for the targeted vitamin analysis performed on the TQS mass spectrometer is also in this folder. 

### MSDIAL
* This folder contains the MSDIAL output for each analytical fraction (HILIC_Pos, HILIC_Neg, RP) for peak areas (Area files) and signal to noise (SN files).



### SIRIUS
* This folder contains all SIRIUS 5 output for unknown compounds for each analytical fraction



## Folder 2: Collaborator_Data
* This folder contains non-metabolomics data included in this paper that was produced outside of the Ingalls Lab. Specifically, it includes ATP data from the Karl Lab at the University of Hawaii, flow cytometry data from the Lindell Lab at the Technion and the Armbrust Lab at the University of Washington, the output of the matrix population model used to predict carbon fixation, growth, and loss rates for Prochlorococcocus from the Armbrust Lab at the University of Washington, and Transcript data produced by the Dyhrman Lab at Columbia University and the Chisholm Lab at MIT.

### ATP
* This folder contains particulate ATP data as a csv file.

### FCM
* This folder contains abundances, estimated mean sizes, carbon quotas for Prochlorococcus and two sizes of bacteria (large and small) estimated using flow cytometry. This folder also contains phage abundances estimated using qPCR.

### MPM
* This folder contains the output of the matrix population model for both all model runs (MPM_results_final_raw) and summarized for each time window and treatment (MPM_results_final.csv). 

### Transcripts
* This folder contains two folders: Transcript_Data and Archive. The Transcript_Data folder contains the normalized transcripts counts for MED4 genes, the associated gene annotations and functions from KEGG, and the sample metadata. The Archive folder contains older versions of transcript analysis outputs and raw data.


## Folder 3: Meta_Data
* This folder contains information related to metabolomics data processing including sample lists, the list of adducts searched for, the volumes filtered for each sample, the Ingalls lab standards list and compound information from when the samples were run, the overall sampling scheme, and the specific transitions monitored for the TQS vitamins.


___ 
# Data Processing, Analysis, and Visualization:

# Folder 1: Source_Code
* This folder contains the biostats.R code and documentation used for the multivariate analysis of the metabolomics data.

# Folder 2: R_Code
* This folder contains three folder and one script (detailed below). The folder “Metabolomics_Data_Processing” contains scripts for processing, normalizing, and quantifying the metabolomics data to produce final peak lists and quantified metabolite data tables. The folder “Analysis” contains scripts for all statistical analysis and comparisons on metabolite, transcript, flow cytometry, ATP, and matrix population model data. The folder “Metabolomics_Workbench_Upload” contains scripts to produce required tables and documents.

## Script 2.0: Functions.R
* This script contains functions used for reading and organizing skyline, MSDIAL, and SIRIUS output, finding adducts, matching SIRIUS annotations to MFs defined by MSDIAL, performing quality control, and matching transcripts to KEGG IDs and pathway information.

## Metabolomics_Data_Processing
### Script 2.1.1
#### Sky_Targeted_Processing.R
* This script identifies targeted compounds present in the pooled samples, determines their batch- and sample-specific mass, accounting for mass defects, and their batch- and sample-specific retention time. This information is then exported as a tab-separated .txt file that can be imported into MS DIAL to identify these targeted compounds in the untargeted dataset.

### Script 2.1.2
#### HILIC_QC1.R
* This script imports MSDIAL and Skyline output of HILIC Positive and HILIC Negative data, compares the quality of MSDIAL and Skyline integrations and identifies compounds that need to be integrated manually in skyline due to poor annotations, replaces those compounds with better manual skyline integrations, performs QC on MFs based on area of peaks above blanks in a minimum number of samples, and exports preliminary QCed data to Intermediates folder.

### Script 2.1.3
#### HILIC_BMIS.R
* Sources organized and combined HILIC Positive and HILIC Negative metabolite data exported from HILIC QC1 from Intermediates folder, runs best-matched internal standard normalization (BMIS) from Boysen et al. (2018), and exports BMISed data to Intermediates folder. 


### Script 2.1.4
#### RP_QC1.R
* This script imports MSDIAL and Skyline output of RP data, compares the quality of MSDIAL and Skyline integrations and identifies compounds that need to be integrated manually in skyline due to poor annotations, replaces those compounds with better manual skyline integrations, and then performs QC on MFs based on area of peaks above blanks in a minimum number of samples.


### Script 2.1.5
#### RP_BMIS.R
* Sources organized and combined RP metabolite data exported from RP QC1 from Intermediates folder, runs best-matched internal standard normalization (BMIS) from Boysen et al. (2018), and exports BMISed data to Intermediates folder.

### Script 2.1.6
#### QC2.R
* This script sources QC1 and BMIS output for HILIC and RP data from the intermediates folder, microbial abundance data from collaborator data folder, and volume filtered data from the meta data folder. This script combines RP and HILIC data, performs quality control based on presence in a minimum number of triplicates, relative standard deviation in a pooled sample, and minimum peak area, combines data with abundance and volume filtered data, and exports quality controlled data tables for all MFs and quantifiable MFs to Intermediates.

### Script 2.1.7
#### RF_RFratio_calcs.R
* Sources organized raw data from Intermediates folder and Ingalls lab standards data from Meta_Data folder, calculates response factors (RFs) and response factor ratios (RFratios) as in Boysen et al. (2018), and exports RFs and RFratios as a csv to Intermediates folder. 


### Script 2.1.8
#### Quantification.R
* Sources BMISed data, organized raw data, RFs and RFratios from Intermediates folder and Ingalls Lab standards information from the meta data folder, runs quantification based on RFs and internal standards as described in Boysen et al. 2018, and calculates concentrations in samples in nmol/L and nmol C/L and nmol N/L.


### Script 2.1.9
#### Vitamin_QC_Quant.R
* Sources skyline vitamin data output from raw data folder, sample volume and monitored transition data from meta data folder, QCs vitamin data based on minimum area and having peak area above the blank in a minimum number of replicates, quantifies vitamins as in Heal et al. 2015 and 2017, and exports vitamin data to Intermediates.

### Script 2.1.10
#### Annotation.R
* Sources final QCed peak list (“Intermediates/ProMo_Particulate_QC2_data.csv”) as well as all SIRIUS output for HILIC Pos, HILIC Neg, and RP fractions for formula predictions, compound identification predictions, and canopus compound predictions. This script then matches the tentative identification information from SIRIUS with the MSDIAL and Skyline identified mass features by exact mass and retention time with user specified retention time error and mass error (in ppm) tolerances. Matched annotations are exported to the Intermediates folder.

### Script 2.1.11
#### Adduct_Detection.R
* This script sources the final QCed peak list (“Intermediates/ProMo_Particulate_QC2_data.csv”) and a modified version of the  Fiehn Lab adduct list, looks for and annotates potential adducts, isotopes, and duplicates between HILIC positive and negative modes by exact mass and retention time (within user defined tolerance windows) and that are highly correlated (correlation greater than or equal to 0.95) within the peak list. Adducts are identified and reported and the correct duplicated mass features to keep are decided on based on which feature has a lower relative standard deviation in repeat injections of the pooled sample. The script exports a list of mass features to annotate as ions and to remove from downstream analyses.  

### Script 2.1.12
#### Final_Processed_Data_Organization.R
* This script sources the final QCed peak list, list of adducts to remove and annotate, matched SIRIUS annotations, quantified concentrations of targeted compounds, and quantified vitamin data. These files are combined into a final untargeted data table alongside sample metadata and abundance data and a final quantified, targeted data table and exported for analyses.




## Analysis

### Script 2.2.1
#### PATP_processing.R
* This script sources raw ATP data provided by Karl lab, tidys the ATP data into an easy to work with format with metadata consistent with other metabolomics data in this repository and calculates carbon biomass using an established conversion ratio from Karl and Holm-Hansen, 1978.


### Script 2.2.2
#### Carbon_Content_Analysis_and_Comparisons.R
* This script sources flow cytometry estimates of Prochlorococcus and bacterial abundances and carbon quotas, quantified metabolite data, and atp concentrations, converts abundances and carbon quotas to carbon biomass, calculates the percentage of microbial biomass that is in Prochlorococcus, calculates the percentage of microbial carbon in metabolites, and then compares the carbon estimates from flow cytometry, atp estimated carbon biomass, and metabolites. 


### Script 2.2.3
#### Multivariate_Analysis.R
* This script sources untargeted metabolomics data and functions from biostats.R, performs an NMDS analysis, a monteccarlo analysis to determine the significance of the projection, and a  PerMANOVA on all three independent variables ("treatment*timepoint*replicate") as well as each possible subset of variables (ex. treatment*replicate for each timepoint separately). 

### Script 2.2.4
#### FC_Analysis.R
* This script sources untargeted metabolite, atp, and flow cytometry data, calculates biomass normalized fold changes in peak areas for each metabolite in the high viral infection samples relative to the replicate matched controls, performs a t-test and Wilcoxon sign test to test for significant differences, tests for linear correlations between metabolites and bacterial biomass, and identifies a final set of metabolites that are altered under viral infection and not correlated with bacterial growth. The final metabolite data table of virally altered compounds and all fold changes are exported.



### Script 2.2.5
#### Transcript_DESeq_Analysis.R
* This script sources transcript counts, sample metadata, and Prochlorococcus strain MED4 gene annotations, organizes data and performs a fold change analysis using the DESeq2 package, annotates genes with gene names, and exports DESeq2 results for interpretation.

### Script 2.2.6
#### MPM_ANOVAs.R
* This script sources matrix population model results, visualizes carbon fixation, cell division, and carbon loss rates for each treatment and time window, performs  2-way ANOVAs to look for differences by treatment and time for each model parameter with biological replicate included as a random variable. A Tukey honest significance difference test is then performed and significantly different treatment and timepoint groups are identified. The discrete groups are exported for incorporation into figures downstream. 









___
#Folder 3: Figures:
##R_Code:

### Script 3.1
#### Sampling_Scheme_Fig.R
* This script produces figure 1b, an overview of all samples collected in this experiment.

### Script 3.2
#### Abundances.R
* This script produces figures 1c, 1d, and 1e which show the abundances of Prochlorococcus, bacteria, and cyanophage P-SSP7 over the course of the experiment.

### Script 3.3
#### Main_Text_Fig2_NMDS_QuantMetabs.R
* This script produces main text figure 2 which includes an NMDS plot of each sample based on its metabolome, a boxplot showing the percentage of total biomass in metabolites for each treatment, and a stacked bar chart showing the relative abundance of quantified, targeted metabolites at each timepoint in each treatment.

### Script 3.3
#### MainTextFig3_VolcanoPLot.R 
* This script produces main text figure 3 which includes a volcano plot showing the fold change under viral infection for all mass features and an inset plot showing the number of mass features that are correlated with bacterial abundance.

### Script 3.4
#### MPM_Fig.R 
* This script produces main text figure 4 which includes bar plots showing the matrix population model predicted growth rates, carbon fixation rate, and carbon loss rate for each treatment and experiment time window.

### Script 3.5
#### Metabolic_Maps_Metabs.R 
* This script produces the fold change heat maps for each metabolite for use in constructing the metabolic maps for main text figures 5–7. 

### Script 3.6
#### Transcripts_for_Metabolic_Maps.R 
* This script produces the fold change heat maps for each relevant transcript for use in constructing the metabolic maps for main text figures 5–7 as well as the transcript heatmap for supplemental figure 6.

### Script 3.7
#### Transcript_Behavior_Figure.R 
* This script produces supplemental figure 1 which summarizes the expression of MED4 transcripts relative to the control (increased, decreased, unchanged) in a filled bar chart.

### Script 3.8
#### Supplemental_Biomass_Carbon_Comparison_Figure.R 
* This script produces supplemental figure 2 which compares amount of biomass estimated by flow cytometry and ATP and compares the amount of metabolite carbon to ATP and flow cytometry estimated carbon in scatter plots and plots over time for the experiment for each treatment. 

### Script 3.9
#### SuppFig3_All_MF_Heatmap.R  
* This script produces supplemental figure 3, a heatmap showing all mass features, hierarchically clustered, for each sample arranged by treatment, timepoint, and replicate.

### Script 3.10
#### SuppFig4_Stacked_Barcharts.R  
* This script produces supplemental figure 4 which depicts stacked and filled stacked bar charts of quantified metabolites for each treatment and replicate over the course of the experiment. 

### Script 3.11
#### SuppFig5_Metab_FC_Figure.R  
* This script produces supplemental figure 5 which shows the fold change relative to the control of all significantly different, known metabolites across treatments and timepoints.

## AD_Figures:
* This folder contains the affinity designer documents and outputs for making main text figures 1, 5, 6, and 7 and supplemental figure 6.


# Folder 4: Tables:

## R_Code:
### Script 4.1
#### Abundances_SuppTable.R
* This script compiles and produces supplemental table 1 containing the abundances, carbon quotas, and biomass estimates from flow cytometry for Prochlorococcus and the size classes of bacteria.

### Script 4.2
#### MF_and_metab_Quant_SuppTables.R
* This script compiles and produces supplemental tables 3, 5, 6, and 7 containing information on all mass features, quantified metabolite concentrations, metabolites fold changes, and untargeted MF SIRIUS annotations. 

### Script 4.3
#### Transcripts_SuppTable_v2.R
* This script produces supplemental tables 2 and 11 containing the normalized transcript counts, fold changes analysis and DESeq2 fold changes and p-values. 

### Script 4.4
#### IS_SuppTable.R
* This script organizes and produces supplemental table 9 containing information on the internal standards used for metabolomics analysis. 

### Script 4.5
#### MPM_SuppTable.R
* This script compiles and produces supplemental table 8 contianing the output of the matrix population model for each timepoint and treatment.









___
## Output Folders:

## Intermediates:
* This folder contains all intermediates files used in workflow. 

## Tables/Outputs:
* This folder contains all final table outputs produced by this workflow. 

## Figures:
* This folder contains all final figure outputs produced by this workflow.
___
## Citation:
## TBD:









