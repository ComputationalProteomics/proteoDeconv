# Load the pure cell type data
pure_data_file <- system.file(
  "extdata",
  "cd8_mono.rds",
  package = "proteoDeconv"
)
pure_samples <- readRDS(pure_data_file)

pure_samples_matrix <- as.matrix(pure_samples[, -1])
rownames(pure_samples_matrix) <- pure_samples[[1]]
pure_samples <- pure_samples_matrix

# Load the mixed sample data
mix_data_file <- system.file(
  "extdata",
  "mixed_samples_matrix.rds",
  package = "proteoDeconv"
)
mixed_samples <- readRDS(mix_data_file)

mixed_samples_matrix <- as.matrix(mixed_samples[, -1])
rownames(mixed_samples_matrix) <- mixed_samples[[1]]
mixed_samples <- mixed_samples_matrix

# Load the precomputed signature matrix
signature_file <- system.file(
  "extdata",
  "cd8_mono_signature.rds",
  package = "proteoDeconv"
)
signature_matrix_df <- readRDS(signature_file)

signature_matrix <- as.matrix(signature_matrix_df[, -1])
rownames(signature_matrix) <- signature_matrix_df[[1]]

rename_columns <- function(colnames_vec) {
  # Extract the meaningful part (Mix, CD8T, or Mono) and the replicate number
  new_names <- sub(".*230516_(CD8T|Mono|Mix)-", "\\1_", colnames_vec) # Extract cell type
  new_names <- sub("\\.mzML$", "", new_names) # Remove file extension

  # Replace short names with full names
  new_names <- gsub("CD8T", "CD8+ T cells", new_names)
  new_names <- gsub("Mono", "Monocytes", new_names)

  # Keep "Genes" column unchanged
  new_names[colnames_vec == "Genes"] <- "Genes"

  return(new_names)
}

# Apply to mixed_samples and pure_samples
colnames(mixed_samples) <- rename_columns(colnames(mixed_samples))
colnames(pure_samples) <- rename_columns(colnames(pure_samples))

# Define the output directory inside the package
output_dir <- file.path("inst", "extdata")

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the matrices
saveRDS(
  mixed_samples,
  file = file.path(output_dir, "mixed_samples_matrix.rds")
)
saveRDS(
  pure_samples,
  file = file.path(output_dir, "pure_samples_matrix.rds")
)
saveRDS(
  signature_matrix,
  file = file.path(output_dir, "cd8t_mono_signature_matrix.rds")
)
