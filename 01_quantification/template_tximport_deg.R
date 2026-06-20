# ==============================================================================
# RNA-seq Analysis Pipeline: tximport, Annotations, and DESeq2
# Template Adapter for Multi-Study Metadata Mapping
# ==============================================================================

# 1. INSTALLATION & ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
# Install required packages if missing
# pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)

# Define Core Directory Variables
STUDY_ID    <- "YOUR_GSE_ID"          # e.g., "GSE201325" or "GSE207923"
SALMON_DIR  <- "PATH/TO/SALMON_OUT"    # Path containing your sample folders
OUTPUT_DIR  <- "outputs"               # Base output directory

# Create defensive output directory tree
dir.create(file.path(OUTPUT_DIR, "metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "counts_data/raw_counts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "counts_data/tpm_counts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "PCA/plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "PCA/data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "DESeq2"), recursive = TRUE, showWarnings = FALSE)


# 2. LOCATE QUANTIFICATION FILES
# ------------------------------------------------------------------------------
# Dynamically retrieve folder names from the study directory
samples <- list.dirs(file.path(SALMON_DIR, STUDY_ID), recursive = FALSE, full.names = FALSE)
if(length(samples) == 0) stop("No sample directories found! Double-check STUDY_ID or SALMON_DIR.")

# Build complete path references to each quant.sf file
quant_files <- file.path(SALMON_DIR, STUDY_ID, samples, 'quant.sf')
names(quant_files) <- samples

# Verify files exist locally before processing
if(!all(file.exists(quant_files))) {
  stop("Missing quant.sf files detected for some directories!")
}


# 3. METADATA & EXPERIMENTAL DESIGN DESIGNATION
# ------------------------------------------------------------------------------
# Construct the experimental mapping layout directly inside data.frame()
col_data <- data.frame(
  row.names = samples,
  sample    = samples,
  
  # OPTION A: Single Factor Experiment (e.g., control vs treated)
  condition = factor(
    c(rep("CONTROL_LABEL", EACH_N), rep("TREATMENT_LABEL", EACH_N)), 
    levels = c("CONTROL_LABEL", "TREATMENT_LABEL") # Establishes baseline
  )
  
  # OPTION B: Multi-Factor Experiment (Uncomment and adjust if using Strain x Time)
  # strain = factor(c(rep("Mock", N1), rep("Variant_A", N2)), levels = c("Mock", "Variant_A")),
  # time   = factor(c(rep("6hpi", M1), rep("12hpi", M2)), levels = c("6hpi", "12hpi"))
)

# Optional: Add interaction groupings or batch tracking vectors if applicable
# col_data$group <- factor(paste0(col_data$strain, "_", col_data$time))
# col_data$batch <- factor(rep(c("batch1", "batch2"), length.out = nrow(col_data)))

# Export metadata backup
write.csv(col_data, file.path(OUTPUT_DIR, "metadata", paste0(STUDY_ID, "_metadata.csv")), row.names = FALSE)


# 4. TRANSCRIPT-TO-GENE MAPPING & SUMMARIZATION
# ------------------------------------------------------------------------------
# Query the Ensembl database for Transcript IDs and Gene Symbols
tx2gene <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86, 
  keys = keys(EnsDb.Hsapiens.v86),
  columns = c('TXID', 'SYMBOL')
)
tx2gene <- dplyr::select(tx2gene, -GENEID) # Remove default internal Ensembl Gene ID

# Compile transcript counts into gene-level abundances via tximport
txi <- tximport(
  files = quant_files, 
  type = 'salmon', 
  tx2gene = tx2gene, 
  ignoreTxVersion = TRUE
)


# 5. COUNT EXPORTS WITH EXPLICIT GENE IDENTIFIERS
# ------------------------------------------------------------------------------
# Export Raw Counts Matrix (with Gene names safely kept as a column)
raw_counts_df <- data.frame(SYMBOL = rownames(txi$counts), txi$counts)
write.csv(raw_counts_df, file.path(OUTPUT_DIR, "counts_data/raw_counts", paste0(STUDY_ID, "_raw_counts.csv")), row.names = FALSE)

# Export Normalized TPM Matrix
tpm_counts_df <- data.frame(SYMBOL = rownames(txi$abundance), txi$abundance)
write.csv(tpm_counts_df, file.path(OUTPUT_DIR, "counts_data/tpm_counts", paste0(STUDY_ID, "_tpm_counts.csv")), row.names = FALSE)


# 6. QUALITY CONTROL & PRINCIPAL COMPONENT ANALYSIS
# ------------------------------------------------------------------------------
# Initialize the primary DESeq2 Dataset Object
# NOTE: Update design to ~group or ~batch + condition depending on your setup
dds <- DESeqDataSetFromTximport(txi = txi, colData = col_data, design = ~condition)

# Perform variance-stabilizing Regularized Logarithm transformation
rlog_dds <- rlog(dds)

# Inject tracking features back into the rlog object for safe plotting calls
colData(rlog_dds) <- colData(dds)

# Generate and save the PCA Plot using the major design factor column
plotPCA(rlog_dds, intgroup = "condition") # Swap with "group" or c("strain","time") if multi-factor
ggsave(file.path(OUTPUT_DIR, "PCA/plot", paste0(STUDY_ID, "_PCA.png")), width = 7, height = 5)

# Export raw PCA coordinates for custom ggplot modifications later
pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, file.path(OUTPUT_DIR, "PCA/data", paste0(STUDY_ID, "_pca_coordinates.csv")), row.names = FALSE)


# 7. DIFFERENTIAL EXPRESSION & GENOMIC ANNOTATION
# ------------------------------------------------------------------------------
# Run core Differential Gene Expression Analysis
dds <- DESeq(dds)

# Extract contrast results and shift to a standard data frame format
res_df <- as.data.frame(results(dds)) # Use contrast argument here if dealing with multi-group designs
res_df$SYMBOL <- rownames(res_df)

# Fetch functional annotation descriptors from the library database
annotations <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86, 
  keys = res_df$SYMBOL,
  keytype = "SYMBOL",
  columns = c("GENENAME", "GENEBIOTYPE")
)
annotations <- annotations[!duplicated(annotations$SYMBOL), ] # Deduplicate map

# Merge annotations back into the DESeq2 statistics output
annotated_res <- merge(res_df, annotations, by = "SYMBOL", all.x = TRUE) %>%
  dplyr::relocate(SYMBOL, GENENAME, GENEBIOTYPE)

# Export finalized fully annotated table
write.csv(annotated_res, file.path(OUTPUT_DIR, "DESeq2", paste0(STUDY_ID, "_deseq2_results.csv")), row.names = FALSE)
print(paste("Pipeline completed successfully for study:", STUDY_ID))