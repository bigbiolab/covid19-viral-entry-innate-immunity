# covid19-viral-entry-innate-immunity

COVID-19 viral entry and innate host–pathogen interactions — a cross-study
RNA-seq meta-analysis of the human transcriptomic response to SARS-CoV-2
infection, comparing **control vs COVID/infected** conditions.

## Overview

This project integrates publicly available GEO RNA-seq cohorts (Group 1:
*Viral Entry & Innate Host–Pathogen Interaction*) into a single random-effects
meta-analysis. Scope is restricted to **control vs COVID** contrasts — SARS-CoV-2
infected cells or COVID-19 patient samples versus uninfected/healthy controls.
Each cohort is quantified with Salmon, summarized to gene level with `tximport`,
and tested for differential expression with DESeq2. Per-study effect sizes
(log2 fold changes and standard errors) are pooled across studies to identify
**robust, reproducible** differentially expressed genes (DEGs), then carried
into network, functional-enrichment, and immune-deconvolution analyses.

### Cohorts (control vs COVID)

| GSE | Model / cells | Comparison | Tissue |
|-----|---------------|------------|--------|
| GSE207923 | NHBE | SARS-CoV-2 infected vs mock (6/12/24 hpi) | Lung |
| GSE217504 | Caco-2 | SARS-CoV-2 infected vs mock (time course) | Gut |
| GSE244488 | A549-ACE2 | SARS-CoV-2 infected (S+) vs mock | Lung |
| GSE245922 | Primary monocytes | COVID-19 vs healthy | Blood |
| GSE255647 | Calu-3/2B4 | SARS-CoV-1/-2 infected vs mock (time course) | Lung |

Mechanistic single-arm cohorts without an infection-vs-control contrast are
**excluded** from this analysis: GSE201325 (Spike protein treatment vs control),
GSE211851 (nsp13 vs vector), and GSE275240 (ACE2 knockout vs wild-type).

## Study Rationale

SARS-CoV-2 infection triggers a host transcriptomic response spanning viral
entry, innate sensing, and interferon signaling. Individual studies are powered
to detect their own effect, but small sample sizes and differences in cell
model, platform, and infection time point make it hard to know **which host
responses are genuinely conserved features of SARS-CoV-2 infection versus
model-specific artefacts.**

A meta-analysis restricted to matched **control-vs-COVID** contrasts across lung,
gut, and blood models lets the shared infection signal rise above per-study
noise. Pooling effect sizes with a random-effects model explicitly accounts for
between-study heterogeneity (different cells, isolates, and time courses) and
yields a confidence-ranked set of DEGs more likely to generalize than any single
study's hit list. That robust core is the right substrate for downstream network
hub detection, pathway enrichment, and immune cell-type deconvolution.

## Research Questions

1. **Conserved response.** Which host genes are robustly and reproducibly
   differentially expressed between SARS-CoV-2 infected and control samples
   across independent RNA-seq cohorts, after accounting for between-study
   heterogeneity?
2. **Tissue specificity.** How do robust infection responses differ between
   lung, gut, and blood/monocyte compartments?
3. **Pathways & networks.** Which biological pathways and protein–protein
   interaction hubs are enriched among the robust DEGs, and do innate immune /
   interferon programs dominate?
4. **Immune landscape.** Do immune cell-type signatures (interferon response,
   monocyte/macrophage, NK, T-cell) shift consistently across cohorts as
   inferred by ssGSEA deconvolution?

## Hypotheses

- **H1.** A conserved core of DEGs — enriched for type I/III interferon and
  NF-κB innate-immune signaling — is reproducibly dysregulated between infected
  and control samples despite differences in cell model, isolate, and time
  point.
- **H2.** Robust DEGs organize into a densely connected PPI network whose hub
  genes are innate-immunity and antiviral-response regulators.
- **H3.** Innate/interferon-associated immune signatures are consistently
  up-shifted in COVID/infected versus control samples, with monocyte/macrophage
  programs most pronounced in the blood cohort (GSE245922).

## Pipeline

| Stage | Directory | Output |
|-------|-----------|--------|
| 1. Quantification | `01_quantification/` | Salmon quant → gene-level counts (`tximport`) |
| 2. Differential expression | `02_diffexp/` | Per-cohort DESeq2 results (infected vs control) |
| 3. Meta-analysis | `03_meta_analysis/` | Random-effects pooled DEGs (`metafor`) |
| 4. Network analysis | `04_network_analysis/` | STRING PPI network, hub genes |
| 5. Enrichment | `05_enrichment_analysis/` | ORA + GSEA (clusterProfiler) |
| 6. Immune analysis | `06_immune_analysis/` | ssGSEA immune-signature deconvolution |

Project-wide parameters (cohort list, tissue groupings, significance thresholds,
network settings) are defined in `config.R`. Robust DEGs are defined as genes
significant (BH-adjusted meta p < 0.05) in at least `ROBUST_N_STUDIES` (3)
cohorts, with |pooled log2FC| ≥ 1 marking high-confidence hits.
