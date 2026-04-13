#!/usr/bin/env Rscript
# Create EFO semantic similarity matrix and lookup table using Lin + Resnik methods.
#
# Computes in row-chunks so the full dense matrix is never held in RAM.
# Peak memory ≈ 2 × chunk_size × n_terms × 8 bytes (Lin + Resnik chunk pair).
# Full matrix parquet is written incrementally (one row-group per chunk).
#
# Run from the repo root (ct-validation/):
#   Rscript scripts/r/create_efo_similarity_matrix.R
#   Rscript scripts/r/create_efo_similarity_matrix.R --efo-path data/mappings/efo_v3.84.0.obo
#
# Dependencies: ontologySimilarity, ontologyIndex, arrow
# Install via: conda env create -f environment.yaml

library(optparse)

option_list <- list(
  make_option("--efo-path",
    type = "character",
    default = "data/mappings/efo_v3.84.0.obo",
    help = "Path to EFO ontology OBO file"
  ),
  make_option("--output",
    type = "character",
    default = "data/mappings/efo_full_similarity_matrix.parquet",
    help = "Output parquet path for full matrix (written incrementally)"
  ),
  make_option("--lookup-output",
    type = "character",
    default = "data/mappings/efo_similarity_lookup_0.5.parquet",
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
  make_option("--chunk-size",
    type = "integer", default = 1000L,
    help = "Row-chunk size (tune to fit available RAM)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

exclude_prefixes <- if (nchar(opt$`exclude-prefixes`) > 0) {
  strsplit(opt$`exclude-prefixes`, ",")[[1]]
} else {
  character(0)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ontologySimilarity", quietly = TRUE)) {
  BiocManager::install("ontologySimilarity", update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(ontologySimilarity)
  library(ontologyIndex)
  library(arrow)
})

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

cat("Loading EFO ontology from", opt$`efo-path`, "...\n")
efo <- get_ontology(opt$`efo-path`)
cat("Loaded", length(efo$id), "total terms\n")

cat("Calculating information content...\n")
ic <- descendants_IC(efo)

all_terms_raw   <- efo$id
has_name        <- all_terms_raw %in% names(efo$name)
is_not_obsolete <- !(all_terms_raw %in% efo$obsolete)
has_colon       <- grepl(":", all_terms_raw, fixed = TRUE)

valid_prefix <- sapply(all_terms_raw, function(term) {
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  prefix <- parts[1]
  if (nchar(prefix) > 20 || nchar(prefix) < 2) return(FALSE)
  if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", prefix)) return(FALSE)
  is_all_caps <- grepl("^[A-Z0-9_]+$", prefix)
  has_multi_cap <- length(gregexpr("[A-Z]", prefix)[[1]]) >= 2
  has_us_num <- grepl("[_0-9]", prefix)
  is_efo <- prefix == "efo"
  is_orphanet <- prefix == "Orphanet"
  is_all_caps || has_multi_cap || has_us_num || is_efo || is_orphanet
})

all_terms <- all_terms_raw[has_name & is_not_obsolete & has_colon & valid_prefix]
cat("Filtered to", length(all_terms), "valid terms\n")

if (length(exclude_prefixes) > 0) {
  cat("Excluding prefixes:", paste(exclude_prefixes, collapse = ", "), "\n")
  prefixes  <- sapply(strsplit(all_terms, ":"), function(x) x[1])
  all_terms <- all_terms[!(prefixes %in% exclude_prefixes)]
  cat("After exclusion:", length(all_terms), "terms\n")
}

col_terms_conv <- convert_from_ontology_format(all_terms)

n          <- length(all_terms)
chunk_size <- opt$`chunk-size`
n_chunks   <- ceiling(n / chunk_size)
threshold  <- opt$`lookup-threshold`
max_ic     <- max(ic, na.rm = TRUE)

cat(sprintf(
  "\nChunked similarity: %d terms, chunk=%d, %d chunks, threshold>=%.2f\n",
  n, chunk_size, n_chunks, threshold
))

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`lookup-output`), recursive = TRUE, showWarnings = FALSE)

# Schema: efo_id (string) + one float64 column per term (matches original output dtype).
full_col_names <- c("efo_id", col_terms_conv)
full_schema <- arrow::schema(
  c(
    list(efo_id = arrow::field("efo_id", arrow::utf8())),
    setNames(
      lapply(col_terms_conv, function(nm) arrow::field(nm, arrow::float64())),
      col_terms_conv
    )
  )
)
full_props  <- arrow::ParquetWriterProperties$create(column_names = full_col_names)
full_sink   <- arrow::FileOutputStream$create(opt$output)
full_writer <- arrow::ParquetFileWriter$create(full_schema, full_sink, properties = full_props)

lookup_list <- vector("list", n_chunks)
t_start     <- proc.time()["elapsed"]

for (ci in seq_len(n_chunks)) {
  row_idx   <- seq((ci - 1) * chunk_size + 1, min(ci * chunk_size, n))
  row_terms <- all_terms[row_idx]
  row_conv  <- col_terms_conv[row_idx]

  lin    <- get_term_sim_mat(efo, ic, "lin",    row_terms = row_terms, col_terms = all_terms)
  resnik <- get_term_sim_mat(efo, ic, "resnik", row_terms = row_terms, col_terms = all_terms)

  resnik <- resnik / max_ic
  avg    <- (lin + resnik) / 2
  rm(lin, resnik)

  for (ri in seq_along(row_idx)) {
    avg[ri, row_idx[ri]] <- 1.0
  }

  colnames(avg) <- col_terms_conv

  # Write rows to full matrix parquet (one row-group per chunk)
  chunk_tbl <- arrow::as_arrow_table(
    cbind(data.frame(efo_id = row_conv, stringsAsFactors = FALSE),
          as.data.frame(avg)),
    schema = full_schema
  )
  full_writer$WriteTable(chunk_tbl, chunk_size = nrow(chunk_tbl))
  rm(chunk_tbl)

  # Accumulate lookup pairs (upper triangle + diagonal >= threshold)
  idx <- which(avg >= threshold & row(avg) + (row_idx[1] - 1) <= col(avg), arr.ind = TRUE)
  if (nrow(idx) > 0) {
    lookup_list[[ci]] <- data.frame(
      efo_id_1   = row_conv[idx[, 1]],
      efo_id_2   = col_terms_conv[idx[, 2]],
      similarity = avg[idx],
      stringsAsFactors = FALSE
    )
  }

  rm(avg)
  gc(verbose = FALSE)

  elapsed  <- proc.time()["elapsed"] - t_start
  per_chunk <- elapsed / ci
  eta       <- per_chunk * (n_chunks - ci)
  cat(sprintf(
    "  chunk %d/%d  (%.0fs elapsed, ~%.0fs remaining)\n",
    ci, n_chunks, elapsed, eta
  ))
}

full_writer$Close()
full_sink$close()
cat("Saved full matrix to", opt$output, "\n")

cat("\nMerging lookup chunks...\n")
lookup <- do.call(rbind, lookup_list)
rm(lookup_list); gc(verbose = FALSE)

is_diag        <- lookup$efo_id_1 == lookup$efo_id_2
lookup_reverse <- data.frame(
  efo_id_1   = lookup$efo_id_2[!is_diag],
  efo_id_2   = lookup$efo_id_1[!is_diag],
  similarity = lookup$similarity[!is_diag],
  stringsAsFactors = FALSE
)
lookup_full <- rbind(lookup, lookup_reverse)

cat("Lookup table:", format(nrow(lookup_full), big.mark = ","), "rows (both directions)\n")
cat("Unique terms:", length(unique(c(lookup_full$efo_id_1, lookup_full$efo_id_2))), "\n")

write_parquet(lookup_full, opt$`lookup-output`)
cat("Saved lookup to", opt$`lookup-output`, "\n")

total_elapsed <- proc.time()["elapsed"] - t_start
cat(sprintf("Done. Total time: %.0f min\n", total_elapsed / 60))
