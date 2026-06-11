# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: May 2026
# Description:
#   Imports transcript-level quantifications from Salmon
#   and summarizes to gene-level counts for DESeq2. 

# Dataset  : GSE275240 — iPSC-derived lung cells (Alveolar vs Airway)
# Cell type: Alveolar (Alv), Airway (Air)
# Condition: WT (control), x484, x1371
# Replicates: 3 per group (total 18 samples)
# BioProject: PRJDB15620

# Install Bioconductor Packages 
pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)

# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/salmon_out/GSE275240', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE275240', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create Metadata (col_data)
# Dataset  : GSE275240 — iPSC-derived lung cells (Alveolar vs Airway)
# Cell type: Alveolar (Alv), Airway (Air)
# Condition: WT (control), x484, x1371
# Replicates: 3 per group (total 18 samples)
condition_map <- c(
  "DRR456169"="WT",   "DRR456170"="x484",  "DRR456171"="x1371",  # 1-Alv
  "DRR456172"="WT",   "DRR456173"="x484",  "DRR456174"="x1371",  # 2-Alv
  "DRR456175"="WT",   "DRR456176"="x484",  "DRR456177"="x1371",  # 3-Alv
  "DRR456178"="WT",   "DRR456179"="x484",  "DRR456180"="x1371",  # 1-Air
  "DRR456181"="WT",   "DRR456182"="x484",  "DRR456183"="x1371",  # 2-Air
  "DRR456184"="WT",   "DRR456185"="x484",  "DRR456186"="x1371"   # 3-Air
)
celltype_map <- c(
  "DRR456169"="Alv", "DRR456170"="Alv", "DRR456171"="Alv",
  "DRR456172"="Alv", "DRR456173"="Alv", "DRR456174"="Alv",
  "DRR456175"="Alv", "DRR456176"="Alv", "DRR456177"="Alv",
  "DRR456178"="Air", "DRR456179"="Air", "DRR456180"="Air",
  "DRR456181"="Air", "DRR456182"="Air", "DRR456183"="Air",
  "DRR456184"="Air", "DRR456185"="Air", "DRR456186"="Air"
)
gsm_map <- c(
  "DRR456169" = "GSM8473291", "DRR456170" = "GSM8473292", "DRR456171" = "GSM8473293",
  "DRR456172" = "GSM8473294", "DRR456173" = "GSM8473295", "DRR456174" = "GSM8473296",
  "DRR456175" = "GSM8473297", "DRR456176" = "GSM8473298", "DRR456177" = "GSM8473299",
  "DRR456178" = "GSM8473300", "DRR456179" = "GSM8473301", "DRR456180" = "GSM8473302",
  "DRR456181" = "GSM8473303", "DRR456182" = "GSM8473304", "DRR456183" = "GSM8473305",
  "DRR456184" = "GSM8473306", "DRR456185" = "GSM8473307", "DRR456186" = "GSM8473308"
)
# Create the data frame with row names AND a explicit sample column
col_data <- data.frame(
  row.names = samples,
  sample    = samples,
  gsm       = gsm_map[samples],
  cell_line = "iPSC",
  celltype  = factor(celltype_map[samples], levels = c("Alv", "Air")),
  condition = factor(condition_map[samples], levels = c("WT", "x484", "x1371")))  # WT = reference level for all comparisons


# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE275240_metadata.csv", row.names = FALSE)

# Get the mapping from transcript IDs to gene symbols 
# What are the columns in the database?
columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)

# Get the TXID and SYMBOL columns for all entries in database
tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                 keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c('TXID', 'SYMBOL'))

# check tx2gene 
head(tx2gene)

# Remove the gene ID column
tx2gene <- dplyr::select(tx2gene, -GENEID)
head(tx2gene)

# Compile the tximport counts object and make DESeq dataset
# Get tximport counts object
txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# class of txi 
class(txi)

# explore raw counts 
txi$counts

# explore normalizec counts 
txi$abundance

# raw counts 
raw_counts <- txi$counts
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE275240_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE275240_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))


# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~celltype + condition)


# Principal Component Analysis
rlog_dds <- rlog(dds)

# PCA Plot
plotPCA(rlog_dds, intgroup = c("condition", "celltype"))
ggsave("outputs/PCA/plot/GSE275240_PCA.png")

## PCA data
pca_data <- plotPCA(rlog_dds, intgroup = c("condition", "celltype"), returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE275240_data.csv", row.names = FALSE)


# Differential Gene Expression Analysis
dds <- DESeq(dds)


# Extract Results for Each Comparison
# Get the results and immediately convert to a standard dataframe

# Extract DESeq2 results and add prefix to stat columns
get_results_df <- function(dds, contrast, prefix) {
  res    <- results(dds, contrast = contrast)
  res_df <- as.data.frame(res)
  res_df$SYMBOL <- rownames(res_df)
  
  # Add prefix to stat columns so both comparisons can sit side by side
  # baseMean is shared (same model), so kept only once
  res_df <- res_df %>%
    dplyr::select(SYMBOL, baseMean,
                  log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    dplyr::rename(
      !!paste0(prefix, "_log2FC")  := log2FoldChange,
      !!paste0(prefix, "_lfcSE")   := lfcSE,
      !!paste0(prefix, "_stat")    := stat,
      !!paste0(prefix, "_pvalue")  := pvalue,
      !!paste0(prefix, "_padj")    := padj
    )
  return(res_df)
}

res_x484  <- get_results_df(dds,
                            contrast = c("condition", "x484",  "WT"),
                            prefix   = "x484_vs_WT")

res_x1371 <- get_results_df(dds,
                            contrast = c("condition", "x1371", "WT"),
                            prefix   = "x1371_vs_WT") %>%
  dplyr::select(-baseMean)  # baseMean already in res_x484, drop duplicate


# ── 12. Merge Both Comparisons into One Table ─────────────────────────────────
combined_res <- full_join(res_x484, res_x1371, by = "SYMBOL")


# ── 13. Annotate ──────────────────────────────────────────────────────────────
annotations <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys    = combined_res$SYMBOL,
  keytype = "SYMBOL",
  columns = c("GENENAME", "GENEBIOTYPE")
)
annotations <- annotations[!duplicated(annotations$SYMBOL), ]

annotated_res <- merge(combined_res, annotations, by = "SYMBOL", all.x = TRUE) %>%
  dplyr::relocate(SYMBOL, GENENAME, GENEBIOTYPE, baseMean)


# ── 14. Save Final Results ────────────────────────────────────────────────────
write.csv(annotated_res,
          "outputs/DESeq2/GSE275240_deseq2_results.csv",
          row.names = FALSE)



# Extract Results for Each Comparison
# Get the results and immediately convert to a standard dataframe
resdf  <- results(dds)
res_df <- as.data.frame(resdf)

# Rescue the row names (which contain your Gene Symbols/IDs) into a column
res_df$SYMBOL <- rownames(res_df)

# Fetch gene annotations (Full Description, Gene Biotype) from EnsDb
annotations <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                     keys = res_df$SYMBOL,
                                     keytype = "SYMBOL",
                                     columns = c("GENENAME", "GENEBIOTYPE"))

# Remove any accidental duplicate rows from the annotation mapping
annotations <- annotations[!duplicated(annotations$SYMBOL), ]

# Merge annotations into your DESeq2 results data frame
annotated_res <- merge(res_df, annotations, by = "SYMBOL", all.x = TRUE)

# Clean up the column layout (Move identifiers to the front)
annotated_res <- annotated_res %>%
  dplyr::relocate(SYMBOL, GENENAME, GENEBIOTYPE)

# Save the final annotated dataset safely!
write.csv(annotated_res, "outputs/DESeq2/GSE275240_deseq2_results.csv", row.names = FALSE)



