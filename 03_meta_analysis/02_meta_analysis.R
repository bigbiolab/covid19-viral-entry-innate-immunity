# 02_meta_analysis.R — Random Effects Meta-analysis

suppressPackageStartupMessages({
  library(metafor)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("Cannot find config.R")
}

output_dir <- META_DIR
out_fig    <- paste0(FIGURES_DIR, "meta_analysis/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig,    showWarnings = FALSE, recursive = TRUE)

# ── Load meta-matrix 
in_path <- paste0(META_DIR, "meta_matrix/meta_matrix_list.rds")
if (!file.exists(in_path)) {
  stop("Meta-matrix not found — run 01_meta_matrix.R first")
}

dataset_list <- readRDS(in_path)
n_studies <- length(dataset_list)

# ── Helper: REM Analysis ──────────────────────────────────────────────────────
run_rem <- function(data_list, prefix) {
  gene_data <- bind_rows(lapply(names(data_list), function(id) {
    data_list[[id]] |>
      dplyr::select(Symbol = Gene_Symbol, log2FoldChange, lfcSE, pvalue, padj) |>
      dplyr::filter(!is.na(Symbol), Symbol != "", !is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0) |>
      dplyr::mutate(cohort = id)
  }))
  
  all_symbols <- unique(gene_data$Symbol)
  message("  [", prefix, "] Unique genes: ", length(all_symbols))
  
  # Use parallel processing if many genes, but for now standard lapply
  reml_rows <- lapply(all_symbols, function(sym) {
    sub <- gene_data[gene_data$Symbol == sym, ]
    k   <- nrow(sub)
    if (k >= 2) {
      fit <- tryCatch(
        metafor::rma(yi = sub$log2FoldChange, sei = sub$lfcSE, method = "DL", verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        return(data.frame(
          Symbol         = sym,
          n_studies      = k,
          pooled_log2FC  = as.numeric(fit$beta),
          pooled_lfcSE   = fit$se,
          rem_pvalue     = fit$pval,
          I2             = fit$I2,
          stringsAsFactors = FALSE
        ))
      }
    }
    # If only 1 study or fit fails, return single study data or NA
    data.frame(
      Symbol         = sym,
      n_studies      = k,
      pooled_log2FC  = sub$log2FoldChange[1],
      pooled_lfcSE   = sub$lfcSE[1],
      rem_pvalue     = sub$pvalue[1],
      I2             = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- bind_rows(reml_rows) |>
    dplyr::mutate(
      rem_padj = p.adjust(rem_pvalue, method = "BH"),
      multi_study = n_studies >= 2
    ) |>
    dplyr::arrange(rem_pvalue)
  
  return(results_df)
}

# ── 1. Global Meta-Analysis ───────────────────────────────────────────────────
message("Running Global Meta-Analysis...")
global_results <- run_rem(dataset_list, "Global")
global_results$is_robust <- global_results$n_studies >= ROBUST_N_STUDIES

# Identification of Robust DEGs
global_results <- global_results |>
  mutate(
    Confidence_Level = case_when(
      is_robust & rem_padj < PADJ_CUTOFF & abs(pooled_log2FC) >= LFC_CUTOFF ~ "High",
      is_robust & rem_padj < PADJ_CUTOFF ~ "Medium",
      rem_padj < PADJ_CUTOFF ~ "Low",
      TRUE ~ "NS"
    )
  )

# Export
saveRDS(global_results, paste0(output_dir, "meta_analysis_results.rds"))
write.csv(global_results, paste0(output_dir, "meta_analysis_results.csv"), row.names = FALSE)

# Robust Table
robust_table <- global_results |> 
  filter(Confidence_Level %in% c("High", "Medium")) |>
  arrange(rem_pvalue)

write.csv(robust_table, paste0(output_dir, "robust_degs.csv"), row.names = FALSE)

# Summary of conclusions
n_high <- sum(global_results$Confidence_Level == "High")
message("Conclusion: Found ", n_high, " high-confidence robust DEGs.")

# Volcano Plot
p_volcano <- ggplot(global_results, aes(x = pooled_log2FC, y = -log10(rem_pvalue), color = Confidence_Level)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c(High = "red", Medium = "orange", Low = "blue", NS = "grey")) +
  theme_minimal() +
  labs(title = "Global Meta-Analysis Volcano Plot",
       x = "Pooled Log2 Fold Change",
       y = "-log10(Meta P-value)")

ggsave(paste0(out_fig, "meta_volcano.png"), p_volcano, width = 8, height = 6)

# Print top 10
print(head(robust_table, 10))
