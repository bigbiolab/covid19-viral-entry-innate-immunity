# 05_enrichment_analysis.R — functional enrichment (ORA + GSEA)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyverse)
})

if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("Cannot find config.R")
}

out_dir <- paste0(OUT_DIR, "enrichment/")
out_fig <- paste0(FIGURES_DIR, "enrichment/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

message("=== COVID-19 Functional Enrichment ===")

meta_path <- paste0(META_DIR, "meta_analysis_results.rds")
if (!file.exists(meta_path)) {
  stop("Meta-analysis results not found — run 03_meta_analysis/02_meta_analysis.R first")
}
meta <- readRDS(meta_path)

# ── 1. Over-Representation Analysis (ORA) ─────────────────────────────────────
# Up/Down on robust genes
up_genes <- meta |> filter(Confidence_Level %in% c("High", "Medium"), pooled_log2FC >= LFC_CUTOFF) |> pull(Symbol)
dn_genes <- meta |> filter(Confidence_Level %in% c("High", "Medium"), pooled_log2FC <= -LFC_CUTOFF) |> pull(Symbol)

# Universe: all genes in meta-analysis
universe <- meta |> pull(Symbol) |> unique()

# Conversion helper
convert_to_entrez <- function(symbols) {
  res <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(res$ENTREZID)
}

universe_entrez <- convert_to_entrez(universe)
up_entrez <- convert_to_entrez(up_genes)
dn_entrez <- convert_to_entrez(dn_genes)

# ORA GO BP
message("Running ORA for up-regulated genes...")
ego_up <- enrichGO(gene = up_entrez, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE, universe = universe_entrez)
message("Running ORA for down-regulated genes...")
ego_dn <- enrichGO(gene = dn_entrez, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE, universe = universe_entrez)

write.csv(as.data.frame(ego_up), paste0(out_dir, "ORA_GO_BP_up.csv"), row.names = FALSE)
write.csv(as.data.frame(ego_dn), paste0(out_dir, "ORA_GO_BP_down.csv"), row.names = FALSE)

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
    dotplot(ego_up, showCategory = 20, title = "Up-regulated GO BP (Robust)")
    ggsave(paste0(out_fig, "ORA_GO_BP_up_dotplot.png"), width = 8, height = 7)
}

if (!is.null(ego_dn) && nrow(as.data.frame(ego_dn)) > 0) {
    dotplot(ego_dn, showCategory = 20, title = "Down-regulated GO BP (Robust)")
    ggsave(paste0(out_fig, "ORA_GO_BP_down_dotplot.png"), width = 8, height = 7)
}

# ── 2. Gene Set Enrichment Analysis (GSEA) ────────────────────────────────────
message("Running GSEA...")
meta$rank_score <- sign(meta$pooled_log2FC) * -log10(pmax(meta$rem_pvalue, 1e-300))
# Optionally weight by n_studies
meta$rank_score <- meta$rank_score * (1 + log2(meta$n_studies))

ranked_list <- meta |> 
  filter(!is.na(Symbol), !is.na(rank_score)) |>
  arrange(desc(rank_score)) |>
  distinct(Symbol, .keep_all = TRUE)
  
gene_list <- setNames(ranked_list$rank_score, ranked_list$Symbol)

gse_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", verbose = FALSE)

write.csv(as.data.frame(gse_go), paste0(out_dir, "GSEA_GO_BP.csv"), row.names = FALSE)

if (!is.null(gse_go) && nrow(as.data.frame(gse_go)) > 0) {
    dotplot(gse_go, showCategory = 20, title = "GSEA GO BP")
    ggsave(paste0(out_fig, "GSEA_GO_BP_dotplot.png"), width = 9, height = 7)
}

message("Enrichment analysis complete.")
