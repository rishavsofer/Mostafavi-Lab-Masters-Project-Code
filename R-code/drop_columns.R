#!/usr/bin/env Rscript

# Script to drop a specific column from the data

# Load readr
library(readr)

# File path
file_path <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/select_gene_features_binarized.tsv"
cat("Reading:", file_path, "\n")

# Read the data
data <- read_tsv(file_path, show_col_types = FALSE)
cat("Loaded", ncol(data), "columns and", nrow(data), "rows\n")

# Column to drop
column_to_drop <- "human_brain3_adult_projected_pcaloadings_clusters_pre_def.125"

# Check if the column exists
if (column_to_drop %in% names(data)) {
  # Drop the specific column
  data <- data[, !names(data) %in% column_to_drop]
  cat("Dropped column:", column_to_drop, "\n")
  
  # Save result
  out_file <- gsub("\\.tsv$", "_modified.tsv", file_path)
  write_tsv(data, out_file)
  cat("Saved result to:", out_file, "\n")
} else {
  cat("Column '", column_to_drop, "' not found in the data\n")
}

cat("Done\n")
