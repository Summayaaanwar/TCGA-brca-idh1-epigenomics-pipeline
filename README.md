# TCGA-BRCA IDH1 Epigenomics Pipeline (Matched Tumor–Normal)

This repository contains an R pipeline to analyze **TCGA-BRCA matched tumor–normal samples** by integrating:
- **RNA-seq expression** (IDH1 stratification)
- **Illumina 450K DNA methylation** (QC + DMPs + DMRs)
- **Promoter enrichment (mCSEA)**
- **Circos/chord plots and gene–CpG network visualization**
---

## What the script does (overview)
1. Loads pre-prepared TCGA objects from `.RData`:
   - expression: `expr_data`
   - methylation: `brca_met`
2. Standardizes TCGA barcodes and identifies Tumor/Normal using sample-type codes.
3. Builds matched Tumor–Normal pairs.
4. Maps **IDH1** (HGNC symbol) to Ensembl using **biomaRt** and defines IDH1 strata:
   - Median-based High/Low
   - Quantile groups (e.g., 0.85 / 0.15)
   - ΔIDH1 (Tumor − Normal) grouping
   - k-means clustering (Tumor and Δ)
5. Expression preprocessing using **edgeR** and **DESeq2**.
6. Methylation QC:
   - removes NA probes
   - removes chrX/chrY probes
   - removes SNP probes (`SNPs.137CommonSingle`)
   - removes cross-reactive probes (Chen et al. list)
7. Differential methylation:
   - Beta → M-values
   - **limma** modeling (including Type × IDH1 interaction)
   - exports full DMP tables
8. DMR calling using **DMRcate** (`lambda=1000`, `C=2`, `pcutoff=0.05`)
9. Promoter enrichment using **mCSEA** (`mCSEATest`, promoters)
10. Visualization: (plots)

---

> Important: The current script uses absolute Desktop paths. For reproducibility, update paths to use `data/`, `resources/`, and `results/`.

---

## Requirements
R (recommended ≥ 4.2)

Packages used in this script include:
- Bioconductor/analysis: `TCGAbiolinks`, `SummarizedExperiment`, `minfi`, `limma`, `edgeR`, `DESeq2`,
  `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `DMRcate`, `mCSEA`
- Data handling: `dplyr`, `tidyr`, `stringr`, `reshape2`
- Visualization: `ggplot2`, `ggpubr`, `ggrepel`, `patchwork`, `ComplexHeatmap`, `pheatmap`,
  `circlize`, `RColorBrewer`
- Networks: `igraph`, `ggraph`, `visNetwork`
- Other: `meta`, `curl`, `sesame`

---

## Setup
### Install packages
```r
install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks","SummarizedExperiment","minfi","limma","edgeR","DESeq2",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19","DMRcate","mCSEA"
))

install.packages(c(
  "dplyr","tidyr","stringr","reshape2",
  "ggplot2","ggpubr","ggrepel","patchwork",
  "ComplexHeatmap","pheatmap",
  "circlize","RColorBrewer",
  "igraph","ggraph","visNetwork",
  "meta","curl","sesame"
))

