# 01_meta_matrix.R — Build merged meta-analysis input matrix

suppressPackageStartupMessages({
  library(tidyverse)
})

if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("Cannot find config.R")
}

out_dir <- paste0(META_DIR, "meta_matrix/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dataset_list <- list()
merged_rows  <- list()

for (id in COHORTS) {
  in_path <- paste0(DESEQ_DIR, id, "_deseq2_results.csv")
  if (!file.exists(in_path)) {
    message("  SKIPPING ", id, " — DESeq2 results file not found")
    next
  }

  df <- read.csv(in_path)
  if (nrow(df) == 0) {
    message("  SKIPPING ", id, " — empty results file")
    next
  }

  # Map SYMBOL to Symbol if necessary, though the original script used Symbol
  # Based on my read_file, it's uppercase SYMBOL
  required <- c("SYMBOL", "log2FoldChange", "lfcSE", "pvalue", "padj")
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    message("  SKIPPING ", id, " — missing columns: ", paste(missing_cols, collapse = ", "))
    next
  }

  df_clean <- df |>
    filter(!is.na(SYMBOL), SYMBOL != "",
           !is.na(log2FoldChange), !is.na(lfcSE),
           !is.na(pvalue), lfcSE > 0, pvalue > 0, pvalue <= 1) |>
    mutate(
      Gene_Symbol    = SYMBOL,
      cohort         = id
    ) |>
    distinct(SYMBOL, .keep_all = TRUE) |>
    select(Gene_Symbol, log2FoldChange, lfcSE, pvalue, padj, cohort,
           any_of(c("baseMean", "stat", "GENENAME", "GENEBIOTYPE")))

  message("  ", id, ": ", nrow(df_clean), " genes")

  dataset_list[[id]] <- df_clean
  merged_rows[[id]]  <- df_clean
}

if (length(dataset_list) == 0) {
  stop("No datasets found in ", DESEQ_DIR)
}

merged_df <- bind_rows(merged_rows)
message("\nTotal rows in merged matrix: ", nrow(merged_df))
message("Unique symbols: ", n_distinct(merged_df$Gene_Symbol))
message("Cohorts represented: ", paste(names(dataset_list), collapse = ", "))

write.csv(merged_df, paste0(out_dir, "meta_matrix_merged.csv"), row.names = FALSE)
saveRDS(dataset_list, paste0(out_dir, "meta_matrix_list.rds"))
