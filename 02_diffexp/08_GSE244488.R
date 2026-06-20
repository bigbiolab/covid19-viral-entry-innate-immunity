# RNA-seq Analysis in R: tximport and Gene-level Summarization
# Author: Md. Jubayer Hossain
# Affiliation: DeepBio Limited | CHIRAL Bangladesh
# Date: May 2026

# Description:
#   Imports transcript-level quantifications from Salmon
#   and summarizes to gene-level counts for DESeq2. 

#   Dataset   : GSE244488 (GEO) / PRJNA1023258, PRJNA1023260, PRJNA1023261 (SRA)
#   Cell line : A549-ACE2 (human alveolar basal epithelial cells;
#               stably overexpressing ACE2 receptor)
#   Conditions:
#     - mock      : uninfected control
#     - bystander : cells exposed to infected culture but not directly infected
#     - infected  : SARS-CoV-2 infection at MOI 1 during 24 h
#
#   Sub-projects:
#     SRP464218 (PRJNA1023258) — 3 mock + 3 bystander + 3 infected
#                                (standard library)
#     SRP464220 (PRJNA1023261) — 6 mock + 6 bystander + 6 infected
#                                (total CoV-depleted library)
#     SRP464222 (PRJNA1023260) — 6 mock + 6 bystander + 6 infected
#                                (total rRNA-depleted library)
#
#   Design consideration:
#     Three sub-projects differ in library preparation type and sequencing depth;
#     a batch covariate (~ condition + SRA_Study) should be included in the
#     DESeq2 design to correct for study-level variation.
#
#   Note: ChIP-seq arm (SRP464232 / SRR26261235–36; anti-p65 and IgG control
#         under TNF-alpha treatment) is excluded from this RNA-seq analysis.

# Install Bioconductor Packages 
pak::pkg_install(c("tidyverse", "tximport", "DESeq2", "EnsDb.Hsapiens.v86"))

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)


# Get the quant files and metadata
# Collect the sample quant files
samples <- list.dirs('outputs/salmon_out/GSE244488', recursive = FALSE, full.names = FALSE)
samples
# Exclude ChIP-seq samples
samples <- samples[!samples %in% c("SRR26261235", "SRR26261236")]

# check quant files 
quant_files <- file.path('outputs/salmon_out/GSE244488', samples, 'quant.sf')
quant_files

# sample names 
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
# all should be TRUE
file.exists(quant_files)  

# Create Metadata (col_data)
# GSE244488: A549-ACE2 human lung epithelial cells infected with SARS-CoV-2
# Conditions:
#     - mock      : uninfected control
#     - bystander : cells exposed to infected culture but not directly infected
#     - infected  : SARS-CoV-2 infection at MOI 1 during 24 h
# Replicates:
#     - SRP464218 : 3 mock + 3 bystander + 3 infected (standard library)
#     - SRP464220 : 6 mock + 6 bystander + 6 infected (total CoV-depleted library)
#     - SRP464222 : 6 mock + 6 bystander + 6 infected (total rRNA-depleted library)
#     - Total RNA-seq samples : 45 (15 mock, 15 bystander, 15 infected)
#     - Note: 2 ChIP-seq samples (SRP464232) excluded from this analysis

condition_map <- c(
  # SRP464218 — standard library 
  "SRR26260925" = "mock",      "SRR26260926" = "mock",      "SRR26260927" = "mock",     
  "SRR26260922" = "bystander", "SRR26260923" = "bystander", "SRR26260924" = "bystander", 
  "SRR26260919" = "infected",  "SRR26260920" = "infected",  "SRR26260921" = "infected",  
  
  # SRP464220 — total CoV-depleted library
  "SRR26260951" = "mock",      "SRR26260949" = "mock",      "SRR26260948" = "mock",     
  "SRR26260947" = "mock",      "SRR26260946" = "mock",      "SRR26260945" = "mock",     
  "SRR26260944" = "bystander", "SRR26260943" = "bystander", "SRR26260942" = "bystander", 
  "SRR26260941" = "bystander", "SRR26260940" = "bystander", "SRR26260939" = "bystander", 
  "SRR26260938" = "infected",  "SRR26260937" = "infected",  "SRR26260936" = "infected",
  "SRR26260935" = "infected",  "SRR26260934" = "infected",  "SRR26260950" = "infected",
  
  # SRP464222 — total rRNA-depleted library
  "SRR26261005" = "mock",      "SRR26261004" = "mock",      "SRR26261003" = "mock",     
  "SRR26261002" = "mock",      "SRR26261001" = "mock",      "SRR26261000" = "mock",     
  "SRR26260999" = "bystander", "SRR26260998" = "bystander", "SRR26260997" = "bystander", 
  "SRR26260996" = "bystander", "SRR26260995" = "bystander", "SRR26260994" = "bystander", 
  "SRR26260993" = "infected",  "SRR26260992" = "infected",  "SRR26260991" = "infected", 
  "SRR26260990" = "infected",  "SRR26260989" = "infected",  "SRR26261006" = "infected"   
)
gsm_map <- c(
  "SRR26260919" = "GSM7817929", "SRR26260920" = "GSM7817928", "SRR26260921" = "GSM7817927", 
  "SRR26260922" = "GSM7817926", "SRR26260923" = "GSM7817925", "SRR26260924" = "GSM7817924",
  "SRR26260925" = "GSM7817923", "SRR26260926" = "GSM7817922", "SRR26260927" = "GSM7817921",
  
  "SRR26260934" = "GSM7817946", "SRR26260935" = "GSM7817945", "SRR26260936" = "GSM7817944", 
  "SRR26260937" = "GSM7817943", "SRR26260938" = "GSM7817942", "SRR26260939" = "GSM7817941",
  "SRR26260940" = "GSM7817940", "SRR26260941" = "GSM7817939", "SRR26260942" = "GSM7817938", 
  "SRR26260943" = "GSM7817937", "SRR26260944" = "GSM7817936", "SRR26260945" = "GSM7817935", 
  "SRR26260946" = "GSM7817934", "SRR26260947" = "GSM7817933", "SRR26260948" = "GSM7817932", 
  "SRR26260949" = "GSM7817931", "SRR26260950" = "GSM7817947", "SRR26260951" = "GSM7817930",
  
  "SRR26260989" = "GSM7817964", "SRR26260990" = "GSM7817963", "SRR26260991" = "GSM7817962", 
  "SRR26260992" = "GSM7817961", "SRR26260993" = "GSM7817960", "SRR26260994" = "GSM7817959",
  "SRR26260995" = "GSM7817958", "SRR26260996" = "GSM7817957", "SRR26260997" = "GSM7817956", 
  "SRR26260998" = "GSM7817955", "SRR26260999" = "GSM7817954", "SRR26261000" = "GSM7817953", 
  "SRR26261001" = "GSM7817952", "SRR26261002" = "GSM7817951", "SRR26261003" = "GSM7817950", 
  "SRR26261004" = "GSM7817949", "SRR26261005" = "GSM7817948", "SRR26261006" = "GSM7817965"
)
study_map <- c(
  "SRR26260919" = "SRP464218", "SRR26260920" = "SRP464218", "SRR26260921" = "SRP464218", 
  "SRR26260922" = "SRP464218", "SRR26260923" = "SRP464218", "SRR26260924" = "SRP464218",
  "SRR26260925" = "SRP464218", "SRR26260926" = "SRP464218", "SRR26260927" = "SRP464218",
  
  "SRR26260934" = "SRP464220", "SRR26260935" = "SRP464220", "SRR26260936" = "SRP464220", 
  "SRR26260937" = "SRP464220", "SRR26260938" = "SRP464220", "SRR26260939" = "SRP464220",
  "SRR26260940" = "SRP464220", "SRR26260941" = "SRP464220", "SRR26260942" = "SRP464220", 
  "SRR26260943" = "SRP464220", "SRR26260944" = "SRP464220", "SRR26260945" = "SRP464220", 
  "SRR26260946" = "SRP464220", "SRR26260947" = "SRP464220", "SRR26260948" = "SRP464220", 
  "SRR26260949" = "SRP464220", "SRR26260950" = "SRP464220", "SRR26260951" = "SRP464220",
  
  "SRR26260989" = "SRP464222", "SRR26260990" = "SRP464222", "SRR26260991" = "SRP464222", 
  "SRR26260992" = "SRP464222", "SRR26260993" = "SRP464222", "SRR26260994" = "SRP464222",
  "SRR26260995" = "SRP464222", "SRR26260996" = "SRP464222", "SRR26260997" = "SRP464222", 
  "SRR26260998" = "SRP464222", "SRR26260999" = "SRP464222", "SRR26261000" = "SRP464222", 
  "SRR26261001" = "SRP464222", "SRR26261002" = "SRP464222", "SRR26261003" = "SRP464222", 
  "SRR26261004" = "SRP464222", "SRR26261005" = "SRP464222", "SRR26261006" = "SRP464222"
)

# Create the data frame with row names AND a explicit sample column
col_data <- data.frame(
  row.names = samples,
  sample    = samples,
  gsm       = gsm_map[samples],
  cell_line = "A549-ACE2",
  condition = factor(condition_map[samples],
                     levels = c("mock", "bystander", "infected")),  # mock = reference
  SRA_Study = factor(study_map[samples],
                     levels = c("SRP464218", "SRP464220", "SRP464222")),
  stringsAsFactors = FALSE
)

# Verify: 15 mock, 15 bystander, 15 infected per condition
table(col_data$condition, col_data$SRA_Study)

# Export metadata for later use 
write.csv(col_data, "outputs/metadata/GSE244488_metadata.csv", row.names = FALSE)

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

# explore normalize counts 
txi$abundance

# raw counts 
raw_counts <- txi$counts
write.csv(raw_counts, "outputs/counts_data/raw_counts/GSE244488_raw_counts.csv", row.names = FALSE)

# TPM 
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/counts_data/tpm_counts/GSE244488_tpm_counts.csv", row.names = FALSE)


# This must return TRUE before you proceed
all(colnames(txi) == rownames(col_data))

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition+SRA_Study)

# Principal Component Analysis 
rlog_dds <- rlog(dds)

# PCA Plot 
plotPCA(rlog_dds, intgroup = "condition")
ggsave("outputs/PCA/plot/GSE244488_PCA.png")

# PCA data 
pca_data <- plotPCA(rlog_dds, intgroup = c("condition", "SRA_Study"), returnData = TRUE)
write.csv(pca_data, "outputs/PCA/data/GSE244488_data.csv", row.names = F)

# Differential Gene Expression Analysis 
dds <- DESeq(dds)

# Get the results and immediately convert to a standard dataframe
# Contrast 1: infected vs mock
res_inf_vs_mock <- results(dds,
                           contrast = c("condition", "infected", "mock"),
                           alpha    = 0.05)

# Contrast 2: bystander vs mock
res_bys_vs_mock <- results(dds,
                           contrast = c("condition", "bystander", "mock"),
                           alpha    = 0.05)

# Convert to data frame
res_df_inf <- as.data.frame(res_inf_vs_mock)
res_df_bys <- as.data.frame(res_bys_vs_mock)

# Rescue the row names (which contain your Gene Symbols/IDs) into a column
res_df_inf$SYMBOL <- rownames(res_df_inf)
res_df_bys$SYMBOL <- rownames(res_df_bys)

# Fetch gene annotations (Full Description, Gene Biotype) from EnsDb
annotations_inf <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                     keys = res_df_inf$SYMBOL,
                                     keytype = "SYMBOL",
                                     columns = c("GENENAME", "GENEBIOTYPE"))
annotations_bys <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                     keys    = res_df_bys$SYMBOL,
                                     keytype = "SYMBOL",
                                     columns = c("GENENAME", "GENEBIOTYPE"))
# Remove any accidental duplicate rows from the annotation mapping
annotations_inf <- annotations_inf[!duplicated(annotations_inf$SYMBOL), ]
annotations_bys <- annotations_bys[!duplicated(annotations_bys$SYMBOL), ]

# Merge annotations into your DESeq2 results data frame
annotated_inf <- merge(res_df_inf, annotations_inf, by = "SYMBOL", all.x = TRUE)
annotated_bys <- merge(res_df_bys, annotations_bys, by = "SYMBOL", all.x = TRUE)

# Clean up the column layout (Move identifiers to the front)
annotated_inf <- annotated_inf %>%
  dplyr::relocate(SYMBOL, GENENAME, GENEBIOTYPE)
annotated_bys <- annotated_bys %>%
  dplyr::relocate(SYMBOL, GENENAME, GENEBIOTYPE)

# Save the final annotated dataset safely!
write.csv(annotated_inf, "outputs/DESeq2/GSE244488_infected_vs_mock_deseq2_results.csv", row.names = FALSE)
write.csv(annotated_bys, "outputs/DESeq2/GSE244488_bystander_vs_mock_deseq2_results.csv", row.names = FALSE)
