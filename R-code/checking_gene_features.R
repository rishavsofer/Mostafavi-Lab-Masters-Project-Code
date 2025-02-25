# Load the data.table package
library(data.table)

# Start timing
start_time <- Sys.time()

# Path to the file
file_path <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/all_gene_features.tsv"

# Read only the first 10 rows using fread
gene_data <- fread(file_path, nrows = 10, sep = "\t")

# Get the column names
col_names <- names(gene_data)

#Print gene data
print(gene_data)

# Print column names
print(col_names)

# End timing
end_time <- Sys.time()
run_time <- end_time - start_time
print(paste("Runtime:", run_time, "seconds"))

# Output folder
output_folder <- "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/R-outputs"

# Save the column names to a file
column_output_path <- file.path(output_folder, "all_gene_features_columns.txt")
write.table(col_names, file = column_output_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save the sample data
sample_output_path <- file.path(output_folder, "all_gene_features_sample.tsv")
fwrite(gene_data, sample_output_path, sep = "\t")
