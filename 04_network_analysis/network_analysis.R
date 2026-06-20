# 04_network_analysis.R — PPI Network Analysis

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(tidyverse)
})

if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("Cannot find config.R")
}

out_csv <- paste0(OUT_DIR, "network/")
out_fig <- paste0(FIGURES_DIR, "network/")
dir.create(out_csv, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ── Load meta-analysis results ────────────────────────────────────────────────
meta_path <- paste0(META_DIR, "meta_analysis_results.rds")
if (!file.exists(meta_path)) {
  stop("Meta-analysis results not found — run 03_meta_analysis/02_meta_analysis.R first")
}
meta <- readRDS(meta_path)

# Use High and Medium confidence genes for the network
sig <- meta |>
  filter(Confidence_Level %in% c("High", "Medium")) |>
  select(Symbol, pooled_log2FC, Confidence_Level)

message("Input genes for network (High/Med Confidence): ", nrow(sig))

if (nrow(sig) < 10) {
  message("  Expanding to all robust DEGs (rem_padj < 0.1)...")
  sig <- meta |>
    filter(is_robust, rem_padj < 0.1) |>
    select(Symbol, pooled_log2FC, Confidence_Level)
}

if (nrow(sig) < 2) {
    stop("Not enough significant genes to build a network.")
}

# ── Network construction ──────────────────────────────────────────────────────
# Note: STRINGdb might require internet access to download data if not cached
string_db <- STRINGdb$new(version="11.5", species=STRING_SPECIES, 
                         score_threshold=STRING_SCORE, input_directory=STRING_CACHE_DIR)

mapped <- string_db$map(as.data.frame(sig), "Symbol", removeUnmappedRows = TRUE)
if (is.null(mapped) || nrow(mapped) < 5) {
    message("Warning: Few STRING mappings found. Attempting to use top 100 DEGs by p-value instead.")
    sig_expanded <- meta |> arrange(rem_pvalue) |> head(100) |> select(Symbol, pooled_log2FC, Confidence_Level)
    mapped <- string_db$map(as.data.frame(sig_expanded), "Symbol", removeUnmappedRows = TRUE)
}

if (is.null(mapped) || nrow(mapped) < 2) stop("Insufficient STRING mappings")

interactions <- string_db$get_interactions(mapped$STRING_id)
if (is.null(interactions) || nrow(interactions) == 0) stop("No interactions found")

g <- graph_from_data_frame(
  d        = interactions[, c("from", "to", "combined_score")],
  directed = FALSE,
  vertices = mapped[, c("STRING_id", "Symbol", "pooled_log2FC", "Confidence_Level")]
)
V(g)$name <- V(g)$Symbol
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# Metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)

# Main component
comps <- components(g)
g_main <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))

# Export metrics
metrics_df <- data.frame(
  Symbol      = V(g_main)$name,
  Degree      = V(g_main)$degree,
  Betweenness = V(g_main)$betweenness,
  Confidence  = V(g_main)$Confidence_Level,
  Log2FC      = V(g_main)$pooled_log2FC
) |> arrange(desc(Degree))
write.csv(metrics_df, paste0(out_csv, "network_metrics.csv"), row.names = FALSE)

# Plot
hub_syms <- metrics_df |> filter(Degree >= quantile(Degree, 0.9)) |> pull(Symbol)

p <- ggraph(g_main, layout = "fr") +
  geom_edge_link(color = "grey80", alpha = 0.2, show.legend = FALSE) +
  geom_node_point(aes(color = pooled_log2FC, size = degree)) +
  geom_node_label(aes(label = ifelse(name %in% hub_syms, name, NA)),
                  size = 3, repel = TRUE, fontface = "bold", label.padding = unit(0.1, "lines")) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "COVID-19 Robust PPI Network",
       subtitle = "Sized by Degree, Colored by Pooled Log2FC",
       size = "Degree", color = "Log2FC") +
  theme_graph()

ggsave(paste0(out_fig, "PPI_network.png"), p, width = 12, height = 10, dpi = 300)

write.csv(metrics_df |> head(HUB_N), paste0(out_csv, "hub_genes.csv"), row.names = FALSE)

message("Network analysis complete.")
