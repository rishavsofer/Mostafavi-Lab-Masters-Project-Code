#!/usr/bin/env Rscript

# Script to:
# 1. Filter TSV file to drop rows with 0 values in numerical feature columns (4-7)
# 2. Merge with GENCODE GTF file based on 'ensg' to gene_id matching
# 3. Add annotation source, genomic start and end from GTF

# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# File paths
tsv_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gene_features_binarized_selected.tsv"
gtf_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gencode.v47lift37.basic.annotation.gtf"
output_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gene_features_filtered_annotated.tsv"

# Function to log with timestamp
log_message <- function(message) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), message, "\n")
}

# Start processing
log_message("Starting processing")

# Read the TSV file
log_message(paste("Reading TSV file:", tsv_file))
tsv_data <- read_tsv(tsv_file, show_col_types = FALSE)
log_message(paste("Loaded", nrow(tsv_data), "rows and", ncol(tsv_data), "columns from TSV file"))

# Print column names to confirm
log_message(paste("TSV file columns:", paste(names(tsv_data), collapse=", ")))

# Step 1: Filter rows based on numerical feature columns (4-7)
log_message("Filtering rows based on numerical feature columns (4-7)...")

# Check if columns 4-7 exist
if(ncol(tsv_data) >= 7) {
  # Get the numerical feature columns (4-7)
  feature_cols <- names(tsv_data)[4:7]
  
  # Create a filter condition: keep rows where ANY of columns 4-7 has a value of 1
  filter_condition <- rowSums(tsv_data[, feature_cols] == 1, na.rm = TRUE) > 0
  
  # Apply the filter
  filtered_data <- tsv_data[filter_condition, ]
  log_message(paste("After filtering: kept", nrow(filtered_data), "rows out of", nrow(tsv_data), 
                    "(", round(nrow(filtered_data)/nrow(tsv_data)*100, 2), "%)"))
} else {
  log_message("Warning: TSV file does not have enough columns (4-7). Continuing without filtering.")
  filtered_data <- tsv_data
}

# Check if 'ensg' column exists for merging
if(!"ensg" %in% names(filtered_data)) {
  log_message("Error: 'ensg' column not found in the TSV file. Cannot proceed with merging.")
  quit(status = 1)
}

# Step 2: Process GTF file to extract gene information
log_message(paste("Reading GTF file:", gtf_file))

# For GTF files, we'll read line by line and extract gene entries
gtf_lines <- readLines(gtf_file)
log_message(paste("Read", length(gtf_lines), "lines from GTF file"))

# Initialize data frame for GTF gene data
gtf_gene_data <- data.frame(
  gene_id = character(),
  annotation_source = character(),
  genomic_start = integer(),
  genomic_end = integer(),
  stringsAsFactors = FALSE
)

# Counter for progress reporting
line_count <- 0
gene_count <- 0

# Process GTF lines
log_message("Processing GTF file to extract gene information...")
for (line in gtf_lines) {
  line_count <- line_count + 1
  
  # Skip comment lines
  if (substr(line, 1, 1) == "#") {
    next
  }
  
  # Split the line by tabs
  columns <- strsplit(line, "\t")[[1]]
  
  # We only want gene entries
  if (length(columns) >= 9 && columns[3] == "gene") {
    gene_count <- gene_count + 1
    
    # Get the required columns
    annotation_source <- columns[2]
    genomic_start <- as.integer(columns[4])
    genomic_end <- as.integer(columns[5])
    
    # Parse the 9th column which contains key-value pairs
    attributes <- columns[9]
    
    # Extract gene_id from attributes
    gene_id_match <- str_match(attributes, 'gene_id "([^"]+)"')
    if (!is.na(gene_id_match[1])) {
      gene_id <- gene_id_match[2]
      
      # Add to our data frame
      gtf_gene_data <- rbind(gtf_gene_data, 
                            data.frame(
                              gene_id = gene_id,
                              annotation_source = annotation_source,
                              genomic_start = genomic_start,
                              genomic_end = genomic_end,
                              stringsAsFactors = FALSE
                            ))
    }
    
    # Report progress every 5,000 genes
    if (gene_count %% 5000 == 0) {
      log_message(paste("Processed", gene_count, "genes..."))
    }
  }
  
  # Report progress every 100,000 lines
  if (line_count %% 100000 == 0) {
    log_message(paste("Processed", line_count, "lines..."))
  }
}

log_message(paste("Extracted information for", nrow(gtf_gene_data), "genes from GTF"))

# Clean gene_id to match 'ensg' format (extract ENSG ID portion)
log_message("Cleaning gene IDs to match ensg format...")
gtf_gene_data$gene_id_clean <- str_extract(gtf_gene_data$gene_id, "ENSG\\d+")

# Step 3: Merge the datasets
log_message("Merging filtered TSV data with GTF gene information...")
merged_data <- left_join(filtered_data, gtf_gene_data, by = c("ensg" = "gene_id_clean"))

# Report on merge results
matched_count <- sum(!is.na(merged_data$annotation_source))
log_message(paste("Matched", matched_count, "out of", nrow(filtered_data), "rows (",
                 round(matched_count/nrow(filtered_data)*100, 1), "%)"))

# Show a preview of columns
if (nrow(merged_data) > 0) {
  log_message("Column names in merged data:")
  log_message(paste(names(merged_data), collapse=", "))
}

# Step 4: Write the merged data to output file
log_message(paste("Writing merged data to:", output_file))
write_tsv(merged_data, output_file)

log_message("Processing completed successfully")
