#!/usr/bin/env Rscript

# Check if required packages are installed, install if needed
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr", repos = "https://cran.r-project.org")
}

# Load necessary libraries
library(readr)

# Set file path
file_path <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/select_gene_features.tsv"

# Read the TSV file
cat("Reading file:", file_path, "\n")
data <- read_tsv(file_path, show_col_types = FALSE)
cat("File loaded. Number of rows:", nrow(data), "Number of columns:", ncol(data), "\n")

# Function to check if a column is binary (contains only 0 and 1)
is_binary <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  all(x %in% c(0, 1, NA))
}

# Function to check if a column is numerical
is_numerical <- function(x) {
  is.numeric(x)
}

# Identify which columns are numerical but not binary
numerical_not_binary <- sapply(data, function(x) is_numerical(x) && !is_binary(x))
numerical_not_binary_cols <- names(data)[numerical_not_binary]

# Print the columns that are numerical but not binary
cat("\nColumns that are numerical but not binary:\n")
cat(paste(numerical_not_binary_cols, collapse = ", "), "\n")
cat("Total count:", length(numerical_not_binary_cols), "\n")

# Process all numerical non-binary columns
if (length(numerical_not_binary_cols) > 0) {
  # Create a copy of the original data
  binarized_data <- data
  
  # Process each column
  for (col in numerical_not_binary_cols) {
    column_data <- data[[col]]
    
    # Calculate the threshold for the top 10%
    threshold <- quantile(column_data, 0.9, na.rm = TRUE)
    
    # Create a new binarized column
    new_column_name <- paste0(col, "_binarized")
    binarized_data[[new_column_name]] <- ifelse(column_data > threshold, 1, 0)
    
    # Print summary for this column
    cat("\nBinarized column:", col, "\n")
    cat("  Threshold value (90th percentile):", threshold, "\n")
    cat("  Number of values set to 1:", sum(binarized_data[[new_column_name]] == 1, na.rm = TRUE), "\n")
    cat("  Number of values set to 0:", sum(binarized_data[[new_column_name]] == 0, na.rm = TRUE), "\n")
  }
  
  # Write the updated data back to a new file
  output_file <- gsub("\\.tsv$", "_binarized.tsv", file_path)
  write_tsv(binarized_data, output_file)
  cat("\nUpdated data written to:", output_file, "\n")
} else {
  cat("\nNo numerical non-binary columns found in the dataset.\n")
}

cat("\nScript completed.\n")
