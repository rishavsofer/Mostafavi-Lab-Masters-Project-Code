#!/usr/bin/env Rscript

# Script to:
# 1. Binarize numerical features (top 10% as 1, bottom 90% as 0)
# 2. Select specific columns (1-5, 100, 10000)
# 3. Save to a new TSV file
if (!requireNamespace("purrr", quietly = TRUE)) {
  install.packages("purrr", repos = "https://cran.r-project.org")
}
library(purrr)

# Load required libraries
library(readr)
library(dplyr)
library(purrr)

# File paths
input_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/all_gene_features.tsv"
output_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gene_features_binarized_selected.tsv"

# Log start time and information
start_time <- Sys.time()
cat("Starting processing at:", format(start_time), "\n")
cat("Reading input file:", input_file, "\n")

# Read the data
data <- read_tsv(input_file, show_col_types = FALSE)
cat("Loaded", ncol(data), "columns and", nrow(data), "rows\n")

# Function to check if a column is numeric
is_numeric_column <- function(col) {
  all(is.na(col) | (is.numeric(col) | !is.na(as.numeric(col))))
}

# Function to check if a column is binary (contains only 0s and 1s, possibly with NAs)
is_binary_column <- function(col) {
  all(is.na(col) | col %in% c(0, 1))
}

# Function to binarize a numeric column (top 10% = 1, others = 0)
binarize_column <- function(col) {
  if (all(is.na(col))) {
    return(col)  # Return as is if all values are NA
  }
  
  # Calculate the 90th percentile (excluding NAs)
  threshold <- quantile(col, 0.9, na.rm = TRUE)
  
  # Binarize: values above threshold (top 10%) = 1, others = 0
  ifelse(col > threshold, 1, 0)
}

# Categorize columns
cat("Analyzing column types...\n")
col_categories <- list(
  character = character(0),
  binary = character(0),
  numeric = character(0),
  other = character(0)
)

for (col_name in names(data)) {
  col <- data[[col_name]]
  
  if (is.character(col)) {
    col_categories$character <- c(col_categories$character, col_name)
  } else if (is_binary_column(col)) {
    col_categories$binary <- c(col_categories$binary, col_name)
  } else if (is_numeric_column(col)) {
    col_categories$numeric <- c(col_categories$numeric, col_name)
  } else {
    col_categories$other <- c(col_categories$other, col_name)
  }
}

# Log column categories
cat("Found", length(col_categories$character), "character columns\n")
cat("Found", length(col_categories$binary), "binary columns\n")
cat("Found", length(col_categories$numeric), "numeric columns\n")
cat("Found", length(col_categories$other), "other columns\n")

# Process data: keep character and binary columns as is, binarize numeric columns
cat("Binarizing numeric columns...\n")
processed_data <- data

# Binarize all numeric columns
for (col_name in col_categories$numeric) {
  processed_data[[col_name]] <- binarize_column(data[[col_name]])
  # Rename the column to indicate binarization
  names(processed_data)[names(processed_data) == col_name] <- paste0(col_name, "_binarized")
}

# Verify binarization
verify_binary <- function(col) {
  all(is.na(col) | col %in% c(0, 1))
}

verification_results <- map_lgl(processed_data, verify_binary)
non_binary_cols <- names(processed_data)[!verification_results]

if (length(non_binary_cols) > 0) {
  cat("Warning: The following columns are not properly binarized:", 
      paste(non_binary_cols, collapse=", "), "\n")
} else {
  cat("Verification successful: All numeric columns have been properly binarized\n")
}

# Select specific columns: 1-5, 100, and 10000
cat("Selecting columns...\n")

# Get total number of columns
total_cols <- ncol(processed_data)

# Validate column indices
selected_indices <- c(1:5, 100, 10000)
valid_indices <- selected_indices[selected_indices <= total_cols]

if (length(valid_indices) < length(selected_indices)) {
  cat("Warning: Some requested columns do not exist. Maximum column index is", total_cols, "\n")
  cat("Will select columns:", paste(valid_indices, collapse=", "), "\n")
}

# Select columns
selected_data <- processed_data[, valid_indices]
cat("Selected", ncol(selected_data), "columns\n")

# Save the result
cat("Saving to output file:", output_file, "\n")
write_tsv(selected_data, output_file)

# Calculate and log elapsed time
end_time <- Sys.time()
elapsed <- end_time - start_time
cat("Processing completed in", format(elapsed), "\n")
cat("Done\n")
