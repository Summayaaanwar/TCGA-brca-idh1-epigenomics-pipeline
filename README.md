# TCGA-BRCA IDH1 Epigenomics Pipeline (Matched Tumor–Normal)

This repository contains the R code used for the **TCGA-BRCA computational analyses** in our study. The pipeline integrates:

- **RNA-seq expression** (IDH1 stratification)
- **Illumina HumanMethylation450 (450K) DNA methylation** (QC + DMPs + DMRs)
- **Promoter enrichment** (mCSEA)
- **Circos/chord plots** and **gene–CpG network / hub-gene** visualization
- 
## Data source (TCGA / GDC)
TCGA-BRCA data are accessed through the **Genomic Data Commons (GDC)** using `TCGAbiolinks` (`GDCquery`, `GDCdownload`, `GDCprepare`).  
TCGA data and large derived matrices are **not redistributed** in this repository; users should download TCGA data locally and run the pipeline as described below.

---
## Scripts

### 1) TCGA IDH1 Script pipeline

**`TCGA_BRCA_IDH1_pipeline.R`**  
End-to-end TCGA pipeline including:
- TCGA data retrieval/preparation (RNA-seq + 450K methylation)
- Barcode harmonization and matched Tumor/Normal pairing
- IDH1 mapping (HGNC → Ensembl via `biomaRt`) and IDH1 stratification
- Expression preprocessing (edgeR / DESeq2 as used in the script)
- Methylation QC:
  - removes NA probes
  - removes chrX/chrY probes
  - removes SNP probes (e.g., `SNPs.137CommonSingle`)
  - removes cross-reactive probes (Chen et al. list)
- Differential methylation:
  - Beta → M-values
  - `limma` modeling (including Type × IDH1 structure used in the script)
  - exports DMP tables
- DMR calling using `DMRcate` (parameters as implemented in the script)
- Promoter enrichment using `mCSEA`
- Figures and summaries: IDH1 expression plots, global methylation plots, Manhattan plot, heatmaps, circos/chord plots, and gene–CpG network/hub plots

### 2) Workflow diagram script
**`IDH1 Work Flow.R`**  
Generates the **IDH1 workflow figure** (schematic overview). This script is for figure generation and does not perform wet-lab analyses.

---

## Important note about file paths
The current scripts use **absolute Desktop paths**. For reproducibility, update paths to use a consistent local structure such as:
- `data/` (local input/output objects)
- `resources/` (probe list files)
- `results/` (tables/figures)

---

## Requirements
R (recommended ≥ 4.2) and the packages used in the scripts, including:

**Bioconductor/analysis:** `TCGAbiolinks`, `SummarizedExperiment`, `minfi`, `limma`, `edgeR`, `DESeq2`,  
`IlluminaHumanMethylation450kanno.ilmn12.hg19`, `DMRcate`, `mCSEA`  
**Data handling:** `dplyr`, `tidyr`, `stringr`, `reshape2`  
**Visualization:** `ggplot2`, `ggpubr`, `ggrepel`, `patchwork`, `ComplexHeatmap`, `pheatmap`, `circlize`, `RColorBrewer`  
**Networks:** `igraph`, `ggraph`, `visNetwork`  
**Other:** `meta`, `curl`, `sesame`

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


