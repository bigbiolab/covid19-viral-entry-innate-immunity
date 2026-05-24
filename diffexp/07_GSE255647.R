# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: May 2026

# Description:
#   Imports transcript-level quantifications from Salmon
#   and summarizes to gene-level counts for DESeq2. 
#     GSE255647: SARS-CoV-1 and SARS-CoV-2 infection in Calu-3/2B4 cells
#     Conditions: Mock (control), SARS-CoV-1, SARS-CoV-2
#     Time points: 12 hpi, 24 hpi, 48 hpi
#     Replicates: 3 per condition per time point (total 27 samples)
#     Cell line  : Calu-3/2B4 (bronchial epithelial cells)
#     BioProject: PRJNA1075838

# Install Bioconductor Packages 
pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)


# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/salmon_out/GSE255647', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE255647', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create Metadata (col_data)
# GSE255647: SARS-CoV-1 and SARS-CoV-2 infection in Calu-3/2B4 cells
# Conditions: Mock (control), SARS-CoV-1, SARS-CoV-2
# Time points: 12 hpi, 24 hpi, 48 hpi
# Replicates: 3 per condition per time point (total 27 samples)
condition_map <- c(
  # Mock
  "SRR27961772" = "Mock", "SRR27961773" = "Mock", "SRR27961774" = "Mock",
  "SRR27961775" = "Mock", "SRR27961776" = "Mock", "SRR27961777" = "Mock",
  "SRR27961778" = "Mock", "SRR27961779" = "Mock", "SRR27961780" = "Mock",
  # SARS-CoV-2
  "SRR27961781" = "SARS2", "SRR27961782" = "SARS2", "SRR27961783" = "SARS2",
  "SRR27961784" = "SARS2", "SRR27961785" = "SARS2", "SRR27961786" = "SARS2",
  "SRR27961787" = "SARS2", "SRR27961788" = "SARS2", "SRR27961789" = "SARS2",
  # SARS-CoV-1
  "SRR27961790" = "SARS1", "SRR27961791" = "SARS1", "SRR27961792" = "SARS1",
  "SRR27961793" = "SARS1", "SRR27961794" = "SARS1", "SRR27961795" = "SARS1",
  "SRR27961796" = "SARS1", "SRR27961797" = "SARS1", "SRR27961798" = "SARS1"
)
time_map <- c(
  "SRR27961772" = "48hpi", "SRR27961773" = "48hpi", "SRR27961774" = "48hpi",
  "SRR27961775" = "24hpi", "SRR27961776" = "24hpi", "SRR27961777" = "24hpi",
  "SRR27961778" = "12hpi", "SRR27961779" = "12hpi", "SRR27961780" = "12hpi",
  "SRR27961781" = "48hpi", "SRR27961782" = "48hpi", "SRR27961783" = "48hpi",
  "SRR27961784" = "24hpi", "SRR27961785" = "24hpi", "SRR27961786" = "24hpi",
  "SRR27961787" = "12hpi", "SRR27961788" = "12hpi", "SRR27961789" = "12hpi",
  "SRR27961790" = "48hpi", "SRR27961791" = "48hpi", "SRR27961792" = "48hpi",
  "SRR27961793" = "24hpi", "SRR27961794" = "24hpi", "SRR27961795" = "24hpi",
  "SRR27961796" = "12hpi", "SRR27961797" = "12hpi", "SRR27961798" = "12hpi"
)
gsm_map <- c(
  "SRR27961772" = "GSM8076559", "SRR27961773" = "GSM8076558",
  "SRR27961774" = "GSM8076557", "SRR27961775" = "GSM8076556",
  "SRR27961776" = "GSM8076555", "SRR27961777" = "GSM8076554",
  "SRR27961778" = "GSM8076553", "SRR27961779" = "GSM8076552",
  "SRR27961780" = "GSM8076551", "SRR27961781" = "GSM8076550",
  "SRR27961782" = "GSM8076549", "SRR27961783" = "GSM8076548",
  "SRR27961784" = "GSM8076547", "SRR27961785" = "GSM8076546",
  "SRR27961786" = "GSM8076545", "SRR27961787" = "GSM8076544",
  "SRR27961788" = "GSM8076543", "SRR27961789" = "GSM8076542",
  "SRR27961790" = "GSM8076541", "SRR27961791" = "GSM8076540",
  "SRR27961792" = "GSM8076539", "SRR27961793" = "GSM8076538",
  "SRR27961794" = "GSM8076537", "SRR27961795" = "GSM8076536",
  "SRR27961796" = "GSM8076535", "SRR27961797" = "GSM8076534",
  "SRR27961798" = "GSM8076533"
)
# Create the data frame with row names AND a explicit sample column
col_data <- data.frame(
   row.names = samples,
   sample    = samples,
   gsm       = gsm_map[samples],
   cell_line = "Calu-3/2B4",
   condition = factor(condition_map[samples],
               levels = c("Mock", "SARS1", "SARS2")),   # Mock = reference level
   timepoint = factor(time_map[samples],
              levels = c("12hpi", "24hpi", "48hpi")))

# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE255647_metadata.csv", row.names = FALSE)

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
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE255647_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE255647_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~timepoint + condition)

# Principal Component Analysis 
rlog_dds <- rlog(dds)

# PCA Plot 
plotPCA(rlog_dds, intgroup = "condition")
ggsave("outputs/PCA/plot/GSE255647_PCA.png")

# PCA data 
pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE255647_data.csv", row.names = F)

# Differential Gene Expression Analysis 
dds <- DESeq(dds)

# Get the results and immediately convert to a standard dataframe
resdf <- results(dds)
res_df <- as.data.frame(resdf)

# Rescue the row names (which contain your Gene Symbols/IDs) into a column
res_df$SYMBOL <- rownames(res_df)

# Get the results and immediately convert to a standard dataframe
# Get results for all three contrasts
add_meta <- function(dds, num, denom) {
  res <- results(dds, contrast = c("condition", num, denom))
  df  <- as.data.frame(res)
  df$SYMBOL   <- rownames(df)
  df$contrast <- paste0(num, "_vs_", denom)
  df
}

res_df <- bind_rows(
  add_meta(dds, "SARS2", "Mock"),
  add_meta(dds, "SARS1", "Mock"),
  add_meta(dds, "SARS2", "SARS1")) |> 
  dplyr::relocate(SYMBOL, contrast)

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
write.csv(annotated_res, "outputs/DESeq2/GSE255647_deseq2_results.csv", row.names = FALSE)
