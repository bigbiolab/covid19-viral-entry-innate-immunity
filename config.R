# config.R - Project configuration for meta-analysis

# Output directories
OUT_DIR <- "outputs/"
DESEQ_DIR <- paste0(OUT_DIR, "DESeq2/")
META_DIR <- paste0(OUT_DIR, "meta_analysis/")
FIGURES_DIR <- paste0(OUT_DIR, "figures/")

# Cohorts
COHORTS <- c(
  "GSE201325","GSE245922", "GSE275240"
)

# Tissue-specific cohorts (based on docs/Group1_Viral Entry & Host-Pathogen Interection.csv)
BLOOD_COHORTS <- c("GSE245922")
LUNG_COHORTS  <- c("GSE201325", "GSE207923", "GSE275240", "GSE255647", "GSE244484")

# Thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 1
ROBUST_N_STUDIES <- 3

# Network Parameters
STRING_SPECIES <- 9606
STRING_SCORE <- 400
STRING_CACHE_DIR <- "string_cache"
HUB_N <- 20
