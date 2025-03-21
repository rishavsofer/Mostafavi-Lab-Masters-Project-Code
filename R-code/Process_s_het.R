#!/usr/bin/env Rscript

# Script to:
# 1. Merge all_gene_features.tsv with prior_mean from s_het_info.xlsx
# 2. Keep only the first 4 columns
# 3. Binarize prior_mean (top 10% as 1, others as 0)
# 4. Filter to keep only rows with prior_mean=1
# 5. Save as s_het_filtered_features.tsv

library(readxl)

library(readr)
library(dplyr)

# File paths
base_dir <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/"
tsv_file <- paste0(base_dir, "all_gene_features.tsv")
xlsx_file <- paste0(base_dir, "s_het_info.xlsx")
output_file <- paste0(base_dir, "s_het_filtered_features.tsv")

# Read the TSV file
message("Reading TSV file: ", tsv_file)
tsv_data <- read_tsv(tsv_file, show_col_types = FALSE)
message("TSV data dimensions: ", nrow(tsv_data), " rows x ", ncol(tsv_data), " columns")

# Read the XLSX file - specifically the second tab "Supplementary Table 1"
message("Reading XLSX file (Supplementary Table 1): ", xlsx_file)
xlsx_data <- read_excel(xlsx_file, sheet = "Supplementary Table 1")
message("XLSX data dimensions: ", nrow(xlsx_data), " rows x ", ncol(xlsx_data), " columns")

# Extract only the necessary columns from the XLSX file
xlsx_subset <- xlsx_data %>% 
  select(ensg, prior_mean)  # No need to rename as the column names match

# Check for column name matches before joining
message("TSV file column names: ", paste(colnames(tsv_data)[1:5], collapse=", "))
message("XLSX subset column names: ", paste(colnames(xlsx_subset), collapse=", "))

# Merge the datasets based on ensg column
message("Merging datasets...")
merged_data <- tsv_data %>%
  left_join(xlsx_subset, by = "ensg")  # Using the same column name in both files

# Check for NAs after joining
na_count <- sum(is.na(merged_data$prior_mean))
message("Number of rows with NA prior_mean after merging: ", na_count)

# Reorder columns to put prior_mean as the 4th column
message("Reordering columns...")
col_names <- colnames(merged_data)
prior_mean_col <- which(col_names == "prior_mean")
col_order <- c(1:3, prior_mean_col, setdiff(4:(ncol(merged_data)-1), prior_mean_col))
merged_data <- merged_data[, col_order]

# Keep only the first 4 columns
message("Keeping only first 4 columns...")
merged_data <- merged_data[, 1:4]

# Binarize the prior_mean column (top 10% as 1, others as 0)
message("Binarizing prior_mean column...")
threshold <- quantile(merged_data$prior_mean, 0.9, na.rm = TRUE)
message("Threshold for top 10%: ", threshold)
merged_data$prior_mean <- ifelse(merged_data$prior_mean >= threshold, 1, 0)

# Filter to keep only rows with prior_mean = 1
message("Filtering to keep only rows with prior_mean = 1...")
filtered_data <- merged_data %>% 
  filter(prior_mean == 1)

message("Final data dimensions: ", nrow(filtered_data), " rows x ", ncol(filtered_data), " columns")

# Save the result
message("Saving result to: ", output_file)
write_tsv(filtered_data, output_file)

message("Processing complete!")

