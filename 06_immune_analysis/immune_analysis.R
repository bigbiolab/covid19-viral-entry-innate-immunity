# 06_immune_analysis.R — Immune Deconvolution Meta-Analysis

suppressPackageStartupMessages({
  library(GSVA)
  library(pheatmap)
  library(ggplot2)
  library(tidyverse)
  library(metafor)
})

if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("Cannot find config.R")
}

out_dir <- paste0(OUT_DIR, "immune_analysis/")
out_fig <- paste0(FIGURES_DIR, "immune_analysis/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

message("=== COVID-19 Immune Meta-Analysis ===")

# Immune gene sets
signatures <- list(
  T_cell_CD8  = c("CD8A","CD8B","GZMB","PRF1","IFNG","GZMA","CD3E"),
  T_cell_CD4  = c("CD4","IL2RA","FOXP3","IL7R","TIGIT","CTLA4"),
  NK_cell     = c("NCAM1","KLRB1","KLRD1","NKG7","GNLY","FCGR3A"),
  Macrophage  = c("CD68","CD163","MRC1","FCGR3A","MSR1","CD14"),
  B_cell      = c("CD19","MS4A1","CD79A","CD79B","IGHM","IGKC"),
  Monocyte    = c("CD14","LYZ","S100A12","FCN1","VCAN","CSF1R"),
  Neutrophil  = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9","CEACAM8")
)

load_counts_robust <- function(path) {
  mat <- as.matrix(read.csv(path, check.names = FALSE, row.names = 1))
  mat
}

# Run ssGSEA per cohort
immune_results <- list()

for (id in COHORTS) {
  message("  Processing ", id, "...")
  counts_file <- paste0(OUT_DIR, "counts_data/tpm_counts/", id, "_tpm_counts.csv")
  meta_file   <- paste0(OUT_DIR, "metadata/", id, "_metadata.csv")
  
  if (!file.exists(counts_file) || !file.exists(meta_file)) {
    message("    Missing files for ", id)
    next
  }
  
  expr_mapped <- load_counts_robust(counts_file)
  if (max(expr_mapped) > 50) expr_mapped <- log2(expr_mapped + 1)
  
  # Run ssGSEA
  tryCatch({
    param <- GSVA::ssgseaParam(expr_mapped, signatures, minSize = 1L)
    ss_res <- GSVA::gsva(param, verbose = FALSE)
    
    # Differential analysis of immune scores
    meta_df <- read.csv(meta_file)
    common <- intersect(colnames(ss_res), meta_df$sample)
    ss_res <- ss_res[, common, drop = FALSE]
    meta_df <- meta_df[match(common, meta_df$sample), ]
    
    # Identify control and treatment
    conds <- unique(meta_df$condition)
    ctrl_names <- c("control", "Mock", "vector", "WT", "healthy", "wild-type")
    ctrl <- intersect(conds, ctrl_names)
    if (length(ctrl) == 0) ctrl <- conds[1] else ctrl <- ctrl[1]
    treat <- setdiff(conds, ctrl)
    if (length(treat) == 0) {
        message("    Only one condition found in ", id)
        next
    }
    treat <- treat[1]
    
    message("    Contrast: ", treat, " vs ", ctrl)
    
    cond_factor <- factor(meta_df$condition, levels = c(ctrl, treat))
    
    cohort_stats <- list()
    for (cell in rownames(ss_res)) {
      fit <- lm(ss_res[cell, ] ~ cond_factor)
      s <- summary(fit)
      if (nrow(s$coefficients) < 2) next
      cohort_stats[[cell]] <- data.frame(
        Cell_Type = cell,
        diff      = s$coefficients[2, 1],
        se        = s$coefficients[2, 2],
        pvalue    = s$coefficients[2, 4],
        cohort    = id
      )
    }
    immune_results[[id]] <- bind_rows(cohort_stats)
  }, error = function(e) message("    Error in ", id, ": ", e$message))
}

# ── Meta-analysis of Immune Scores ─────────────────────────────────────────────
all_immune_stats <- bind_rows(immune_results)
cells <- unique(all_immune_stats$Cell_Type)

immune_meta <- list()
for (cell in cells) {
  sub <- all_immune_stats |> filter(Cell_Type == cell)
  if (nrow(sub) < 2) {
    # If only 1 study, just use that
    if (nrow(sub) == 1) {
        immune_meta[[cell]] <- data.frame(
          Cell_Type = cell,
          pooled_diff = sub$diff,
          pooled_se   = sub$se,
          pvalue      = sub$pvalue,
          n_studies   = 1
        )
    }
    next
  }
  
  fit <- tryCatch(rma(yi = diff, sei = se, data = sub, method = "DL"), error = function(e) NULL)
  if (!is.null(fit)) {
    immune_meta[[cell]] <- data.frame(
      Cell_Type = cell,
      pooled_diff = as.numeric(fit$beta),
      pooled_se   = fit$se,
      pvalue      = fit$pval,
      n_studies   = nrow(sub)
    )
  }
}

immune_meta_df <- bind_rows(immune_meta) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  arrange(pvalue)

write.csv(immune_meta_df, paste0(out_dir, "immune_meta_analysis.csv"), row.names = FALSE)

# Plot
p_immune <- ggplot(immune_meta_df, aes(x = reorder(Cell_Type, pooled_diff), y = pooled_diff, fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "COVID-19 Immune Infiltration Meta-Analysis",
       subtitle = "Positive = Increased in Infection/Treatment",
       x = "Cell Type", y = "Pooled Score Difference") +
  theme_bw()

ggsave(paste0(out_fig, "immune_meta_barplot.png"), p_immune, width = 8, height = 6)

message("Immune Meta-Analysis Complete.")
print(immune_meta_df)
