
# Script to merge select_gene_features_modified.tsv with gencode GTF file
# Matching 'ensg' column to annotation source and including genomic start/end locations
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cran.r-project.org")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr", repos = "https://cran.r-project.org")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr", repos = "https://cran.r-project.org")
}


# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# File paths
tsv_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/select_gene_features_binarized_modified.tsv"
gtf_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gencode.v47lift37.basic.annotation.gtf"
output_file <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gene_features_gencode_merged.tsv"

cat("Reading TSV file:", tsv_file, "\n")
# Read the TSV file
tsv_data <- read_tsv(tsv_file, show_col_types = FALSE)
cat("Loaded", nrow(tsv_data), "rows and", ncol(tsv_data), "columns from TSV file\n")

# Print column names to confirm
cat("TSV file columns:", paste(names(tsv_data), collapse=", "), "\n")

# Check if 'ensg' column exists
if (!"ensg" %in% names(tsv_data)) {
  stop("Error: 'ensg' column not found in the TSV file")
}

cat("Reading GTF file:", gtf_file, "\n")
# For GTF files, we need a custom approach as they have a specific format
# Read GTF file - this will take some time for large files
gtf_lines <- readLines(gtf_file)
cat("Read", length(gtf_lines), "lines from GTF file\n")

# Process only gene entries to extract gene IDs, annotation sources, and genomic locations
cat("Processing GTF file to extract gene information...\n")
gtf_gene_data <- data.frame(gene_id = character(), 
                           annotation_source = character(),
                           genomic_start = integer(),
                           genomic_end = integer(),
                           stringsAsFactors = FALSE)

# Counter for progress reporting
line_count <- 0
gene_count <- 0

# Process GTF lines
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
                             data.frame(gene_id = gene_id, 
                                       annotation_source = annotation_source,
                                       genomic_start = genomic_start,
                                       genomic_end = genomic_end,
                                       stringsAsFactors = FALSE))
    }
    
    # Report progress every 5,000 genes
    if (gene_count %% 5000 == 0) {
      cat("Processed", gene_count, "genes...\n")
    }
  }
  
  # Report progress every 100,000 lines
  if (line_count %% 100000 == 0) {
    cat("Processed", line_count, "lines...\n")
  }
}

cat("Extracted information for", nrow(gtf_gene_data), "genes from GTF\n")

# Clean gene_id to match 'ensg' format
gtf_gene_data$gene_id_clean <- str_extract(gtf_gene_data$gene_id, "ENSG\\d+")

# Merge the data
cat("Merging datasets...\n")
merged_data <- left_join(tsv_data, gtf_gene_data, by = c("ensg" = "gene_id_clean"))

# Report on merge results
matched_count <- sum(!is.na(merged_data$annotation_source))
cat("Matched", matched_count, "out of", nrow(tsv_data), "rows (", 
    round(matched_count/nrow(tsv_data)*100, 1), "%)\n")

# Show a preview of the first few rows of the final data
if (nrow(merged_data) > 0) {
  cat("\nPreview of first few rows and column names in final data:\n")
  cat("Column names:", paste(names(merged_data), collapse=", "), "\n\n")
}

# Write the merged data to a new file
cat("Writing merged data to:", output_file, "\n")
write_tsv(merged_data, output_file)

cat("Done\n")
