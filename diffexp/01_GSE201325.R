# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: May 2026

# Description:
#  Imports transcript-level quantification from Salmon
#  and summarizes to gene-level counts for DESeq2. 

#    Dataset: GSE201325 — SARS-CoV-2 spike protein treatment in Calu-3 cells
#    Tissue : Lung
#    Condition: control (control plasmid), treated (SARS-CoV-2 spike protein 100nM)
#    Replicates: 3 per group (total 6 samples)
#    BioProject: PRJNA830876

# Install Bioconductor Packages 
pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)


# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/salmon_out/GSE201325', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE201325', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create Metadata (col_data)
# Dataset: GSE201325 — SARS-CoV-2 spike protein treatment in Calu-3 cells
# Condition: control (control plasmid), treated (SARS-CoV-2 spike protein 100nM)
# Replicates: 3 per group (total 6 samples)
condition_map <- c(
  "SRR18889440" = "treated",     # spike protein 100nM
  "SRR18889441" = "treated",
  "SRR18889442" = "control",     # control plasmid
  "SRR18889443" = "treated",
  "SRR18889444" = "control",
  "SRR18889445" = "control"
)
gsm_map <- c(
  "SRR18889442" = "GSM6058589",
  "SRR18889444" = "GSM6058588",
  "SRR18889445" = "GSM6058587",
  "SRR18889440" = "GSM6058592",
  "SRR18889441" = "GSM6058590",
  "SRR18889443" = "GSM6058591"
)
# Create the data frame with row names AND a explicit sample column
col_data <- data.frame(
  row.names = samples,
  sample    = samples,
  gsm       = gsm_map[samples],
  tissue    = "lung",
  cell_line = "Calu-3",
  condition = factor(condition_map[samples],
                     levels = c("control", "treated")))  # control = reference


# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE201325_metadata.csv", row.names = FALSE)

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
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE201325_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE201325_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition)

# Principal Component Analysis 
rlog_dds <- rlog(dds)

# PCA Plot 
plotPCA(rlog_dds, intgroup = "condition")
ggsave("outputs/PCA/plot/GSE201325_PCA.png")

# PCA data 
pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE201325_data.csv", row.names = F)

# Differential Gene Expression Analysis 
dds <- DESeq(dds)

# Get the results and immediately convert to a standard dataframe
resdf <- results(dds, contrast = c("condition", "treated", "control"))
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
write.csv(annotated_res, "outputs/DESeq2/GSE201325_deseq2_results.csv", row.names = FALSE)
