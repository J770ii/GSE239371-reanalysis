# GSE239371-reanalysis

Batch-aware reanalysis of maternal microbiome-associated DNA methylation 
changes in mouse hippocampus.

## Overview

This repository contains the R scripts used for the reanalysis of the 
GSE239371 dataset. This work was submitted to *Basic and Clinical 
Neuroscience* as:

> Shegofte, J. (2026). Maternal microbiome programs developmental-stage-
> specific DNA methylation in mouse hippocampus: a batch-aware reanalysis. 
> *Basic and Clinical Neuroscience*.

## Files

- `GSE239371_reanalysis.R` — Complete R pipeline including:
  - Data loading and quality control
  - Sentrix batch assessment
  - Dam differential methylation analysis
  - Offspring differential methylation analysis (batch-controlled)
  - CpG annotation with TxDb.Mmusculus.UCSC.mm10.knownGene
  - Gene Ontology enrichment (custom hypergeometric test)
  - Figure generation (PCA, Volcano, Heatmap, GO barplot)

## Requirements

- R version 4.3.0 or higher
- Bioconductor packages:
  - limma (>= 3.56.0)
  - TxDb.Mmusculus.UCSC.mm10.knownGene
  - org.Mm.eg.db
  - GO.db
  - GenomicRanges
  - GEOquery
- CRAN packages:
  - ggplot2
  - ggrepel
  - pheatmap
  - R.utils

## Data source

- **GEO accession:** GSE239371
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239371
- **Platform:** Illumina Infinium MouseMethylation-12 v1.0 BeadChip (GPL32685)

## Contact

Jalal Shegofte — [your-email@example.com]
