# HER2 Amplification Analysis in TCGA-BRCA Pan-Cancer Atlas
This script is a comprehensive gene expression analysis identifying HER2-specific transcriptional programs and prognostic signatures in breast invasive carcinoma using TCGA Pan-Cancer Atlas data from cBioPortal <https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018>.

## Overview

This R Markdown script performs differential expression analysis between ERBB2-amplified and normal breast tumors from the TCGA-BRCA PanCanAtlas cohort (~1,000 primary tumors). Key analyses include:

1. **Multi-omics integration**: RNA-seq, CNA, clinical data alignment
2. **HER2 stratification**: CNA-defined amplification status
3. **DESeq2 DE analysis**: 403 significant genes (padj<0.05)
4. **Pathway enrichment**: Immune activation, ECM remodeling signatures
5. **Prognostic modeling**: 17-gene elastic net Cox signature (p<0.0001)

## Workflow Sections

### 1. Data Loading
```
Data Loading
├── untar- Extract TCGA-BRCA tarball (optional)
├── read-RNAseq: Load RSEM RNA-seq counts (genes × samples)
├── read-patientData: Clinical data w/ survival (skip 4 header rows)
└── read-cna: GISTIC copy number aberrations
```

### 2. Sample Matching & ERBB2 Metadata
```
Sample Matching and ERBB2 Metadata
├── match-ids: Standardize TCGA barcodes (12-char format), find 1068 common patients
├── ERBB2-level: Extract ERBB2 CNA row, define amplification (>0 log2 ratio)
│               Output: 738 Normal, 330 Amplified (31%)
```

### 3. DESeq2 Differential Expression
```
DESeq2 Analysis
├── dds-her2: Filter low-count genes, create DESeqDataSet(~ERBB2_amp)
├── deseq-normalize: Median-of-ratios normalization + Wald testing
│                   Output: res_her2 (403 DE genes padj<0.05)
│                   Top: CSN3 (log2FC=-5.91, padj=2.2e-29)
└── top10_fc: Rank by |log2FC| (casein family dominates)
```

### 4. Visualization
```
Visualization
├── plotPCA-her2: VST-transformed PCA (PC1=20% separates HER2 status)
└── clustering: pheatmap top-20 DE genes (row-scaled, HER2 annotation)
```

### 5. Pathway Enrichment Analysis
```
Pathway Analysis (clusterProfiler)
├── enrichment-1: GO-BP ("antimicrobial humoral response", p=8.9e-08)
├── enrichment-2: KEGG (complement/coagulation cascades)
├── enrichment-3: Reactome (ECM organization, R-HSA-1474244)
└── enrichment-4: Similarity treeplots
```

### 6. Prognostic Cox Regression
```
Survival Analysis
├── Survival matching: Map OS time/status to HER2 cohort
├── glmnet: Elastic net Cox (α=0.5, 5-fold CV)
│           Output: 17 prognostic genes at optimal λ
└── ggsurvplot: Risk tertiles (p<0.0001 survival separation)
```

## Prerequisites

```r
# Core analysis packages
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "ReactomePA"))
install.packages(c("pheatmap", "glmnet", "survival", "survminer", "ggpubr", "dplyr", "enrichplot"))
```

## Data Requirements

Download TCGA-BRCA PanCanAtlas from cBioPortal: <https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018>

```
brca_tcga_pan_can_atlas_2018/
├── data_mrna_seq_v2_rsem.txt     # RNA-seq RSEM counts
├── data_clinical_patient.txt     # Clinical + survival
└── data_cna.txt                 # GISTIC CNA
```
