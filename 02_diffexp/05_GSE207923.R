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
samples <- list.dirs('outputs/salmon_out/GSE207923', recursive = FALSE, full.names = FALSE)
samples

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE207923', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create the data frame with row names AND a explicit sample column
# GSE207923: SARS-CoV-2 infection in NHBE cells
# Conditions: Mock, USA-WA1/2020 (WA1), New York 1-PV08001/2020 (NY1)
# Time points: 6hpi, 12hpi, 24hpi
# Note: Some GSMs have multiple SRRs (split runs)
infection_map <- c(
  "SRR20079510"="Mock", "SRR20079511"="Mock", "SRR20079512"="Mock",
  "SRR20079513"="Mock", "SRR20079514"="Mock", "SRR20079515"="Mock",
  "SRR20079516"="Mock", "SRR20079517"="Mock", "SRR20079518"="Mock",
  "SRR20079519"="Mock", "SRR20079520"="Mock", "SRR20079521"="Mock",
  "SRR20079522"="WA1",  "SRR20079523"="WA1",  "SRR20079524"="WA1",
  "SRR20079525"="WA1",  "SRR20079526"="WA1",  "SRR20079527"="WA1",
  "SRR20079528"="WA1",  "SRR20079529"="WA1",  "SRR20079530"="WA1",
  "SRR20079531"="WA1",  "SRR20079532"="NY1",  "SRR20079533"="NY1",
  "SRR20079534"="NY1",  "SRR20079535"="NY1",  "SRR20079536"="NY1",
  "SRR20079537"="NY1",  "SRR20079538"="NY1",  "SRR20079539"="NY1",
  "SRR20079540"="NY1",  "SRR20079541"="NY1",  "SRR20079542"="NY1",
  "SRR20079543"="NY1"
)
time_map <- c(
  "SRR20079510"="24hpi", "SRR20079511"="24hpi", "SRR20079512"="24hpi",
  "SRR20079513"="24hpi", "SRR20079514"="12hpi", "SRR20079515"="12hpi",
  "SRR20079516"="12hpi", "SRR20079517"="12hpi", "SRR20079518"="6hpi",
  "SRR20079519"="6hpi",  "SRR20079520"="6hpi",  "SRR20079521"="6hpi",
  "SRR20079522"="24hpi", "SRR20079523"="24hpi", "SRR20079524"="24hpi",
  "SRR20079525"="12hpi", "SRR20079526"="12hpi", "SRR20079527"="12hpi",
  "SRR20079528"="6hpi",  "SRR20079529"="6hpi",  "SRR20079530"="6hpi",
  "SRR20079531"="6hpi",  "SRR20079532"="24hpi", "SRR20079533"="24hpi",
  "SRR20079534"="24hpi", "SRR20079535"="24hpi", "SRR20079536"="12hpi",
  "SRR20079537"="12hpi", "SRR20079538"="12hpi", "SRR20079539"="12hpi",
  "SRR20079540"="6hpi",  "SRR20079541"="6hpi",  "SRR20079542"="6hpi",
  "SRR20079543"="6hpi"
)
gsm_map <- c(
  "SRR20079510"="GSM6323035", "SRR20079511"="GSM6323034",
  "SRR20079512"="GSM6323033", "SRR20079513"="GSM6323033",
  "SRR20079514"="GSM6323032", "SRR20079515"="GSM6323032",
  "SRR20079516"="GSM6323031", "SRR20079517"="GSM6323030",
  "SRR20079518"="GSM6323029", "SRR20079519"="GSM6323028",
  "SRR20079520"="GSM6323027", "SRR20079521"="GSM6323027",
  "SRR20079522"="GSM6323026", "SRR20079523"="GSM6323025",
  "SRR20079524"="GSM6323024", "SRR20079525"="GSM6323023",
  "SRR20079526"="GSM6323022", "SRR20079527"="GSM6323021",
  "SRR20079528"="GSM6323020", "SRR20079529"="GSM6323019",
  "SRR20079530"="GSM6323018", "SRR20079531"="GSM6323018",
  "SRR20079532"="GSM6323017", "SRR20079533"="GSM6323016",
  "SRR20079534"="GSM6323015", "SRR20079535"="GSM6323015",
  "SRR20079536"="GSM6323014", "SRR20079537"="GSM6323013",
  "SRR20079538"="GSM6323012", "SRR20079539"="GSM6323012",
  "SRR20079540"="GSM6323011", "SRR20079541"="GSM6323010",
  "SRR20079542"="GSM6323009", "SRR20079543"="GSM6323009"
)

col_data <- data.frame(
  row.names = samples,
  sample    = samples,
  GSM       = gsm_map[samples],
  condition = factor(infection_map[samples], levels = c("Mock", "WA1", "NY1")),
  timepoint = factor(time_map[samples],      levels = c("6hpi", "12hpi", "24hpi"))
)

# infection as factor (Mock = reference)
col_data$infection <- factor(col_data$condition,
                             levels = c("Mock", "WA1", "NY1"))

# timepoint as factor (6hpi = reference)
col_data$timepoint <- factor(col_data$timepoint,
                             levels = c("6hpi", "12hpi", "24hpi"))

# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE207923_metadata.csv", row.names = FALSE)

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
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE207923_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE207923_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~timepoint + condition)

dds_collapsed <- collapseReplicates(dds,
                                    groupby = col_data$GSM,
                                    run     = col_data$sample)

# Principal Component Analysis 
rlog_dds <- rlog(dds_collapsed)

# PCA Plot 
plotPCA(rlog_dds)
ggsave("outputs/PCA/plot/GSE207923_PCA.png")

# PCA data 
pca_data <- plotPCA(rlog_dds, intgroup = "condition", returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE207923_data.csv", row.names = F)

# Differential Gene Expression Analysis 
dds_collapsed <- DESeq(dds_collapsed)

# Get the results and immediately convert to a standard dataframe
resdf <- results(dds_collapsed)
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
write.csv(annotated_res, "outputs/DESeq2/GSE207923_deseq2_results.csv", row.names = FALSE)
