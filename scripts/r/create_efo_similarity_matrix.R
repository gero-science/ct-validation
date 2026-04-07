#!/usr/bin/env Rscript
# Create full EFO semantic similarity matrix using Lin + Resnik methods.
#
# Usage:
#   Rscript create_efo_similarity_matrix.R
#   Rscript create_efo_similarity_matrix.R --efo-path /path/to/efo.obo --output /path/to/output.parquet
#
# Dependencies: ontologySimilarity, ontologyIndex, arrow
# Install via: conda env create -f environment.yaml

library(optparse)

# Parse arguments
option_list <- list(
  make_option("--efo-path",
    type = "character",
    default = "/data/klim/mappings/efo_v3.84.0.obo",
    help = "Path to EFO ontology OBO file"
  ),
  make_option("--output",
    type = "character",
    default = "/data/klim/mappings/efo_full_similarity_matrix.parquet",
    help = "Output parquet path for full matrix"
  ),
  make_option("--lookup-output",
    type = "character",
    default = "/data/klim/mappings/efo_similarity_lookup_0.5.parquet",
    help = "Output parquet path for lookup table"
  ),
  make_option("--lookup-threshold",
    type = "double", default = 0.5,
    help = "Minimum similarity threshold for lookup table"
  ),
  make_option("--exclude-prefixes",
    type = "character", default = "PR,OBA",
    help = "Comma-separated ontology prefixes to exclude"
  ),
  make_option("--temp-dir",
    type = "character", default = "/data/klim/mappings",
    help = "Directory for temporary files"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Parse exclude prefixes
exclude_prefixes <- if (nchar(opt$`exclude-prefixes`) > 0) {
  strsplit(opt$`exclude-prefixes`, ",")[[1]]
} else {
  character(0)
}

# Setup packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ontologySimilarity", quietly = TRUE)) {
  BiocManager::install("ontologySimilarity", update = FALSE, ask = FALSE)
}

# Load required packages
suppressPackageStartupMessages({
  library(ontologySimilarity)
  library(ontologyIndex)
  library(arrow)
})

temp_lin_path <- file.path(opt$`temp-dir`, "temp_lin_similarity.rds")
temp_resnik_path <- file.path(opt$`temp-dir`, "temp_resnik_normalized.rds")

# Convert efo:EFO_1234567 -> EFO:1234567, keep others as-is
convert_from_ontology_format <- function(term_ids) {
  result <- character(length(term_ids))
  for (i in seq_along(term_ids)) {
    term <- term_ids[i]
    if (is.na(term)) {
      result[i] <- NA
    } else if (startsWith(term, "efo:EFO_")) {
      result[i] <- paste0("EFO:", substr(term, 9, nchar(term)))
    } else {
      result[i] <- term
    }
  }
  result
}

# Load ontology
cat("Loading EFO ontology from", opt$`efo-path`, "...\n")
efo <- get_ontology(opt$`efo-path`)
cat("Loaded", length(efo$id), "total terms\n")

# Calculate information content
cat("Calculating information content...\n")
ic <- descendants_IC(efo)

# Filter to valid terms
all_terms_raw <- efo$id
has_name <- all_terms_raw %in% names(efo$name)
is_not_obsolete <- !(all_terms_raw %in% efo$obsolete)
has_colon <- grepl(":", all_terms_raw, fixed = TRUE)

valid_prefix <- sapply(all_terms_raw, function(term) {
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return(FALSE)
  }
  prefix <- parts[1]
  if (nchar(prefix) > 20 || nchar(prefix) < 2) {
    return(FALSE)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", prefix)) {
    return(FALSE)
  }
  is_all_caps <- grepl("^[A-Z0-9_]+$", prefix)
  has_multiple_caps <- length(gregexpr("[A-Z]", prefix)[[1]]) >= 2
  has_underscore_or_number <- grepl("[_0-9]", prefix)
  is_efo <- prefix == "efo"
  is_orphanet <- prefix == "Orphanet"
  return(is_all_caps || has_multiple_caps || has_underscore_or_number || is_efo || is_orphanet)
})

all_terms <- all_terms_raw[has_name & is_not_obsolete & has_colon & valid_prefix]
cat("Filtered to", length(all_terms), "valid terms\n")

# Exclude specified prefixes
if (length(exclude_prefixes) > 0) {
  cat("Excluding prefixes:", paste(exclude_prefixes, collapse = ", "), "\n")
  prefixes <- sapply(strsplit(all_terms, ":"), function(x) x[1])
  all_terms <- all_terms[!(prefixes %in% exclude_prefixes)]
  cat("After exclusion:", length(all_terms), "terms\n")
}

# Calculate Lin similarity
cat("\nCalculating Lin similarity matrix (", length(all_terms), "x", length(all_terms), ")...\n")
sim_matrix_lin <- get_term_sim_mat(
  ontology = efo, information_content = ic, method = "lin",
  row_terms = all_terms, col_terms = all_terms
)
diag(sim_matrix_lin) <- 1.0
cat("Saving Lin matrix...\n")
saveRDS(sim_matrix_lin, temp_lin_path, compress = "gzip")
rm(sim_matrix_lin)
gc()

# Calculate Resnik similarity (normalized)
cat("Calculating Resnik similarity matrix...\n")
sim_matrix_resnik <- get_term_sim_mat(
  ontology = efo, information_content = ic, method = "resnik",
  row_terms = all_terms, col_terms = all_terms
)
max_ic <- max(ic, na.rm = TRUE)
sim_matrix_resnik <- sim_matrix_resnik / max_ic
diag(sim_matrix_resnik) <- 1.0
cat("Saving normalized Resnik matrix...\n")
saveRDS(sim_matrix_resnik, temp_resnik_path, compress = "gzip")
rm(sim_matrix_resnik, ic, efo)
gc()

# Average matrices
cat("\nAveraging Lin and Resnik matrices...\n")
sim_matrix_lin <- readRDS(temp_lin_path)
sim_matrix_resnik <- readRDS(temp_resnik_path)
sim_matrix_avg <- (sim_matrix_resnik + sim_matrix_lin) / 2
diag(sim_matrix_avg) <- 1.0
rm(sim_matrix_lin, sim_matrix_resnik)
gc()

# Convert labels
rownames(sim_matrix_avg) <- convert_from_ontology_format(rownames(sim_matrix_avg))
colnames(sim_matrix_avg) <- convert_from_ontology_format(colnames(sim_matrix_avg))

# Summary
cat("\nMatrix:", nrow(sim_matrix_avg), "x", ncol(sim_matrix_avg), "\n")
cat("Similarity range: [", sprintf("%.3f", min(sim_matrix_avg, na.rm = TRUE)),
  ", ", sprintf("%.3f", max(sim_matrix_avg, na.rm = TRUE)), "]\n",
  sep = ""
)
cat("Pairs >= 0.5:", sum(sim_matrix_avg >= 0.5, na.rm = TRUE), "\n")

# Save full matrix as parquet
cat("\nSaving full matrix to", opt$output, "...\n")
df <- as.data.frame(sim_matrix_avg)
write_parquet(df, opt$output)

# Create lookup table (pairs >= threshold, both directions)
cat("\nCreating lookup table (threshold >=", opt$`lookup-threshold`, ")...\n")
threshold <- opt$`lookup-threshold`

# Extract upper triangle pairs >= threshold (including diagonal for self-matches)
n <- nrow(sim_matrix_avg)
indices <- which(sim_matrix_avg >= threshold & row(sim_matrix_avg) <= col(sim_matrix_avg), arr.ind = TRUE)
cat("Found", format(nrow(indices), big.mark = ","), "pairs in upper triangle + diagonal\n")

lookup <- data.frame(
  efo_id_1 = rownames(sim_matrix_avg)[indices[, 1]],
  efo_id_2 = colnames(sim_matrix_avg)[indices[, 2]],
  similarity = sim_matrix_avg[indices]
)

# Add reverse direction for bidirectional lookup (skip diagonal - already symmetric)
is_diagonal <- lookup$efo_id_1 == lookup$efo_id_2
lookup_reverse <- data.frame(
  efo_id_1 = lookup$efo_id_2[!is_diagonal],
  efo_id_2 = lookup$efo_id_1[!is_diagonal],
  similarity = lookup$similarity[!is_diagonal]
)
lookup_full <- rbind(lookup, lookup_reverse)

cat("Lookup table:", format(nrow(lookup_full), big.mark = ","), "rows (both directions)\n")
cat("Unique terms:", length(unique(c(lookup_full$efo_id_1, lookup_full$efo_id_2))), "\n")

# Save lookup
write_parquet(lookup_full, opt$`lookup-output`)
cat("Saved lookup to", opt$`lookup-output`, "\n")

# Cleanup temp files
file.remove(temp_lin_path, temp_resnik_path)
cat("Done.\n")
