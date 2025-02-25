library(data.table)

# Start timing
start_time <- Sys.time()

# File paths
input_path <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/all_gene_features.tsv"
output_path <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/select_gene_features.tsv"

# Read the full TSV file
cat("Reading the file...\n")
gene_data <- fread(input_path, sep = "\t")

# Select specific columns
cat("Extracting columns...\n")
selected_columns <- c("ensg", "hgnc", "chrom", 
                      "human_brain3_adult_projected_pcaloadings_clusters_pre_def.125", 
                      "GO:0061154")

# Check if all columns exist in the dataset
missing_cols <- selected_columns[!selected_columns %in% names(gene_data)]
if(length(missing_cols) > 0) {
  warning("The following columns were not found: ", paste(missing_cols, collapse = ", "))
  # Continue with columns that do exist
  selected_columns <- selected_columns[selected_columns %in% names(gene_data)]
}

# Extract the selected columns
selected_data <- gene_data[, ..selected_columns]

# Write the selected data to a new TSV file
cat("Writing to output file...\n")
fwrite(selected_data, output_path, sep = "\t")

# End timing
end_time <- Sys.time()
run_time <- end_time - start_time
cat(paste("Runtime:", run_time, "seconds\n"))

# Print information about the output
cat(paste("Selected", ncol(selected_data), "columns and", nrow(selected_data), "rows\n"))
cat(paste("Output saved to:", output_path, "\n"))

