# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: May 2026
# Description:
#   Imports transcript-level quantifications from Salmon
#   and summarizes to gene-level counts for DESeq2. 

# Install Bioconductor Packages 
pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)


# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/salmon_out/GSE217504', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE217504', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create the data frame with row names AND a explicit sample column
# GSE217504: SARS-CoV-2 infection in Caco-2 cells
# Conditions: mock (control), infected (SARS-CoV-2)
# Time points: mock → 4h, 12h, 48h | infected → 0h, 1h, 2h, 4h, 7h, 12h, 24h, 48h
# Replicates: 3 per condition per time point (total 33 samples)
col_data <- data.frame(
  row.names = samples,
  condition = factor(c(rep("mock", 9),       # SRR22223232–240
                       rep("infected", 24)),   # SRR22223241–264
              levels = c("mock", "infected")),
  timepoint = factor(c(
              # mock (9 samples)
              "48h", "48h", "48h", "12h", "12h", "12h", "4h",  "4h",  "4h",
              # infected (24 samples)
              "48h", "48h", "48h", "24h", "24h", "24h", "12h", "12h", "12h", "7h",  "7h",  "7h", 
              "4h",  "4h",  "4h", "2h",  "2h",  "2h", "1h",  "1h",  "1h", "0h",  "0h",  "0h"), 
              levels = c("0h", "1h", "2h", "4h", "7h", "12h", "24h", "48h")))

# condition as factor (mock = reference)
col_data$condition <- factor(col_data$condition,
                             levels = c("mock", "infected"))

# timepoint as factor (0h = reference)
col_data$timepoint <- factor(col_data$timepoint,
                             levels = c("0h", "1h", "2h", "4h", "7h", "12h", "24h", "48h"))


# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE217504_metadata.csv", row.names = FALSE)

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
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE217504_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE217504_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition)

# Principal Component Analysis 
rlog_dds <- rlog(dds)

# PCA Plot 
plotPCA(rlog_dds)
ggsave("outputs/PCA/plot/GSE217504_PCA.png")

# PCA data 
pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE217504_data.csv", row.names = F)

# Differential Gene Expression Analysis 
dds <- DESeq(dds)

# Get the results and immediately convert to a standard dataframe
resdf <- results(dds)
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
write.csv(annotated_res, "outputs/DESeq2/GSE217504_deseq2_results.csv", row.names = FALSE)
