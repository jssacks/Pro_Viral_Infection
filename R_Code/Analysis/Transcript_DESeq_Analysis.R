





####################################################################################################

# Differential Gene Expression Analysis of MED4 (Prochlorococcus) - Viral Infection vs. Control at the Same Time point
# May 21, 2024
# Provided Tables: (1) MED4 mapped counts (Counts_MED4.csv), (2) Metadata (Metadata.csv), (3) MED4 KEGG Annotations (MED4_KEGG.csv)

# Note on GENE.NAME: When the first sequence of MED4 was submitted to NCBI the organism prefix was "PMM." As genes were further identified or discovered in other Prochlorococcus clades they were replaced with new or orthologous gene names. However, there are still some genes unnamed. For this study, the "PMM" prefix remains for these unnamed genes.

####################################################################################################

# 1 - Set working environment and upload files

####################################################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
### Packages Used
library(DESeq2)
library(tidyverse)
library(dplyr)
library(data.table)


# Load tables
counts <- read.csv("Collaborator_Data/Transcripts/Transcript_Data/Counts_MED4.csv")
metadata <- read.csv("Collaborator_Data/Transcripts/Transcript_Data/Metadata.csv")
MED4_annot <- read.csv("Collaborator_Data/Transcripts/Transcript_Data/MED4_KEGG.csv")

####################################################################################################

# 2 - Prepare Comparisons for DESeq2
# The goal here is to isolate associated data from the counts and metadata for each comparison to input into DESeq2
# We start by defining the comparisons to input into DESeq2 (Control vs. Viral) for all 4 time points (0, 12, 24, and 36 hours)
# The final output are two lists (metadata and counts) for all four time point comparisons

####################################################################################################

# Time points to iterate over
comparisons <- list(c("Control 0", "Viral 0"), c("Control 12", "Viral 12"), c("Control 24", "Viral 24"), c("Control 36", "Viral 36"))

# Lists
counts_list <- list()
#subset_dataframes <- list()
metadata_list <- list()

for (comp in comparisons) {
  # Create df names based on comparison (Control Time vs. Viral Time)
  df_name <- paste(comp, collapse = "_")
  
  # Filter metadata based on comparison criteria
  comparison_meta <- metadata %>%
    filter(Condition %in% comp) 
  
  # Get unique sample IDs (SID)
  treatment.SID <- unique(comparison_meta$SID)
  
  # Select columns dynamically based on SID
  comparison_counts <- counts %>%
    # !! unquotes
    dplyr::select(GENE.NAME, !!treatment.SID)
  
  # Store the data frame in the list
  counts_list[[df_name]] <-  comparison_counts
  
  # Subset metadata rows based on SID, then store in metadata_list
  metadata_subset <- metadata %>%
    filter(SID %in% treatment.SID)
  
  metadata_list[[df_name]] <- metadata_subset
}

####################################################################################################

# 3 - (A) Perform DESeq Analysis Using Lists from (#2) and extract logFoldChange and padj data for each comparison, (B) Combine into single data frame, (C) Annotate with KEGG

# Important functions used: (A) DESeqDataSetFromMatrix() - prepares data to be input into next function, (B) DESeq() - estimates size factors, dispersion, and fits data to negative binomial GLM and applies Wald statistics, (C) results() - returns results table of previous function including log2 fold change, standard error, p-values, adjusted p-values, and more. For more information on each function review the respective R Documentation.

# Note on results() function setup: contrast = c("condition", "treated","untreated")

####################################################################################################

# Part A - Perform DE Analysis and Extract Data of Interest
# Get the names of the lists from Step 2
df_names <- names(counts_list)
meta_names <- names(metadata_list)

# Iterate over the names and match them
results_list <- list()
extracted_results_list<-list()

for (name in intersect(df_names, meta_names)) {
  # Access the data frame and metadata based on the matched name
  df <- counts_list[[name]]
  meta <- metadata_list[[name]]
  meta$Time.point <- as.factor(meta$Time.point) 
  print(name)
  
  # Create a nested list to store multiple sets of results
  # Results {SETUP: contrast = c("condition", "treated","untreated")}
  name_results <- list()
  
  # Comparison (Viral v. Control)
  cat("Performing DE for", name, "\n") # If you want it to read out what comparison it is performing
  dds <- DESeqDataSetFromMatrix(countData=df, colData=meta, design=~Treatment, tidy = TRUE)
  dds <- DESeq(dds, minReplicatesForReplace = Inf)
  res <-results(dds, contrast = c("Treatment" ,"Treatment", "Untreated"), cooksCutoff = FALSE, independentFiltering = FALSE)
  
  # Extract padj and log2FoldChange columns
  extracted_res <- as.data.frame(res[, c("padj", "log2FoldChange")])
  #extracted_res$comparison <- name
  
  # Store the results in a nested list with the sample name
  name_results$DE <- list(sample_name = name, results = res, extracted = extracted_res)
  cat("Completed DE for", name, "\n")
  
  # Store the results in the list
  results_list[[name]] <- name_results
  extracted_results_list[[name]] <- extracted_res
}


# Part B - Combine all lists into 1 dataframe for downstream analysis
combined_df_list <- lapply(names(extracted_results_list), function(name) {
  df <- extracted_results_list[[name]]
  # Add suffix to column names
  colnames(df) <- paste(colnames(df), name, sep = "_")
  return(df)
})

# Combine all data frames into one
combined_df <- do.call(cbind, combined_df_list) %>% rownames_to_column(., var = "GENE.NAME")

# Part C - Annotate genes with KEGG

DE_annotated <- merge(combined_df, MED4_annot, by = "GENE.NAME")

# BOOM! Ready for the fun part!

#organize and tidy for Figures + Tables:
transcript_dat <- DE_annotated %>%
  rename("HV0_LFC_nonDE_padj" = `padj_Control 0_Viral 0` ,
        "HV12_LFC_nonDE_padj" =  `padj_Control 12_Viral 12`,
        "HV24_LFC_nonDE_padj" = `padj_Control 24_Viral 24`,
        "HV36_LFC_nonDE_padj" = `padj_Control 36_Viral 36`) %>%
  rename("HV0_LFC_nonDE_log2FoldChange" = `log2FoldChange_Control 0_Viral 0`,
         "HV12_LFC_nonDE_log2FoldChange" = `log2FoldChange_Control 12_Viral 12`,
         "HV24_LFC_nonDE_log2FoldChange" = `log2FoldChange_Control 24_Viral 24`,
         "HV36_LFC_nonDE_log2FoldChange" = `log2FoldChange_Control 36_Viral 36`)

#export transcript data:
write_csv(transcript_dat, "Intermediates/MED4_DE_transcripts.csv")
































































