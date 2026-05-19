# Gene2Protein

A comprehensive guide for analyzing the CYP2E1 gene, its transcripts, proteins, and genomic mappings using R and IGV.

## Gene Overview

| Property | Value |
|----------|-------|
| **Gene Name** | CYP2E1 |
| **Ensembl ID** | ENSG00000130649 |
| **Chromosome** | 10 (minus strand) |
| **Genome Build** | GRCh38 / hg38 |
| **Primary Tissue** | Liver |
| **Function** | Breaks down toxins and alcohol (cytochrome P450 enzyme) |

---

## Table of Contents

- [Requirements](#requirements)
- [File Overview](#file-overview)
- [How to Run](#how-to-run)
- [Task Details](#task-details)
  - [Task 1: IGV View](#task-1--igv-view)
  - [Task 2: List All Transcripts](#task-2--list-all-transcripts)
  - [Task 3: List All Proteins](#task-3--list-all-proteins)
  - [Task 4: Protein-to-Genome Map](#task-4--protein-to-genome-map)
- [Key Notes](#key-notes)
- [References](#references)

---

## Requirements

### Software

- **R** (version ≥ 4.0)
- **IGV** — Download from [Broad Institute](https://software.broadinstitute.org/software/igv/download)
  -  Install locally, NOT on Google Drive

### R Packages

All scripts install their own dependencies automatically on first run. Core packages required:

```r
BiocManager::install(c(
  "ensembldb",
  "EnsDb.Hsapiens.v86",   # Ensembl v86, GRCh38 human annotation
  "AnnotationFilter",
  "IRanges",
  "GenomicRanges"
))
