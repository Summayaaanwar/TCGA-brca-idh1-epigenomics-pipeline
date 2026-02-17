#TCGA IDH1 Meth Analysis

setwd("~/Desktop/R/TCGA") 

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "limma", "minfi", 
                       "dplyr", "ggplot2", "ComplexHeatmap"))

# Load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(minfi)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)

#### creating barcode for downloading paired data#########

# Query metadata for all methylation samples (no download yet)
query_all <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value"
)

# Extract metadata
metadata <- getResults(query_all)

# Extract the first 3 fields of barcode for patient-level ID
metadata$patient_id <- substr(metadata$cases, 1, 12)

# Split by sample type
tumor_samples  <- metadata[metadata$sample_type == "Primary Tumor", ]
normal_samples <- metadata[metadata$sample_type == "Solid Tissue Normal", ]

# Find patient IDs with both tumor and normal samples
common_patients <- intersect(tumor_samples$patient_id, normal_samples$patient_id)

# Subset to 50 matched patients (or fewer if not enough)
matched_patients <- head(common_patients, 60)

# Get corresponding barcodes
tumor_barcodes  <- tumor_samples$cases[tumor_samples$patient_id %in% matched_patients]
normal_barcodes <- normal_samples$cases[normal_samples$patient_id %in% matched_patients]

# Combine tumor and normal barcodes
selected_barcodes <- c(tumor_barcodes, normal_barcodes)
# Check how many times each patient appears
table(substr(selected_barcodes, 1, 12))

# Function to select one tumor and one normal per patient
select_one_pair <- function(tumor_df, normal_df, patient_ids) {
  selected_tumors  <- do.call(rbind, lapply(patient_ids, function(pid) tumor_df[tumor_df$patient_id == pid, ][1, ]))
  selected_normals <- do.call(rbind, lapply(patient_ids, function(pid) normal_df[normal_df$patient_id == pid, ][1, ]))
  rbind(selected_tumors, selected_normals)
}

# Apply to your matched patients
final_df <- select_one_pair(tumor_samples, normal_samples, matched_patients)

# Final selected barcodes
selected_barcodes <- final_df$cases
length(selected_barcodes)  # Should be exactly 120

write.table(selected_barcodes, 
            file = "selected_TCGA_barcodes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
rm(final_df,metadata,normal_samples,tumor_samples,query_all)

#### GDC query with selected barcodes### perform it after verifying data for transcriptome.
query_meth <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  barcode = selected_barcodes
)

GDCdownload(query_meth, files.per.chunk = 10, directory ="~/Desktop/R/TCGA/meth/")
brca_met <- GDCprepare(query_meth, summarizedExperiment = TRUE, directory = "~/Desktop/R/TCGA/meth/")
save(brca_met, file = "~/Desktop/R/TCGA/meth_BRCA_matched_120.RData")

table(colData(brca_met)$sample_type)

beta <- assay(brca_met)
write.csv(as.data.frame(beta), "meth-TCGA.csv")

metadata <- colData(brca_met)
metadata_df <- as.data.frame(metadata)

# Convert list columns to character
metadata_df[] <- lapply(metadata_df, function(col) {
  if (is.list(col)) {
    sapply(col, function(x) paste(as.character(x), collapse = "; "))
  } else {
    col
  }
})

# Write to CSV
write.csv(metadata_df, "rowData-meth.csv", row.names = FALSE)


row_anno <- rowData(brca_met)
write.csv(as.data.frame(row_anno), "meth-probe_annotations.csv")


##### .... verifying if the barcode samples also have transcriptome data...##3

# Check expression availability
# Get just the first 15 characters for each barcode
short_barcodes <- substr(selected_barcodes, 1, 15)

# Remove duplicates
short_barcodes <- unique(short_barcodes)

# Try querying again
query_expr_check <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = short_barcodes
)

# Check how many matched
matched_barcodes <- unique(query_expr_check$results[[1]]$cases)
length(matched_barcodes)

query_expr <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = matched_barcodes
)

GDCdownload(query_expr, files.per.chunk = 10, directory ="~/Desktop/R/TCGA/expr/")
expr_data <- GDCprepare(query_expr, summarizedExperiment = TRUE, directory = "~/Desktop/R/TCGA/expr/")
save(expr_data, file = "~/Desktop/R/TCGA/Exp_BRCA_matched_120.RData")
# Extract expression matrix
expr_matrix <- assay(expr_data)

# Check if gene IDs are Ensembl (e.g., "ENSG00000138413") or gene symbols
# If Ensembl IDs, map them to symbols:
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_ids <- rownames(expr_matrix)

# Remove version numbers
gene_ids_clean <- gsub("\\..*", "", gene_ids)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Attach gene symbols as rownames or metadata
names(gene_symbols) <- rownames(expr_matrix)

idh1_ensembl <- mapIds(org.Hs.eg.db, keys = "IDH1", column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
idh1_row <- which(gene_ids_clean == idh1_ensembl)
idh1_expr <- expr_matrix[idh1_row, ]
idh1_expr <- as.numeric(idh1_expr)
names(idh1_expr) <- colnames(expr_matrix)

#====.......Tumor data only.....####
tumor_samples <- colData(expr_data)$sample_type == "Primary Tumor"
expr_matrix_tumor <- expr_matrix[, tumor_samples]

#...identify IDH1 row for tumor only......###

idh1_ensembl <- mapIds(org.Hs.eg.db, keys = "IDH1", 
                       column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
idh1_index <- which(gene_ids_clean == idh1_ensembl)
# Extract IDH1 expression in tumors
tumor_idh1_expr_values <- expr_matrix_tumor[idh1_index, ]
# strtify be median value
median_idh1 <- median(tumor_idh1_expr_values, na.rm = TRUE)
idh1_status <- ifelse(tumor_idh1_expr_values >= median_idh1, "High", "Low")

tumor_coldata <- colData(expr_data)[tumor_samples, ]
tumor_coldata$IDH1_group <- idh1_status
dim(tumor_coldata)

### phenotype information for tumor samples
# Get tumor-only sample barcodes
tumor_barcodes <- colnames(expr_matrix)[tumor_samples]

# Extract IDH1 groups
idh1_status <- ifelse(expr_matrix[idh1_index, tumor_samples] >= median_idh1, "High", "Low")

# Make phenotype dataframe
pheno_df <- data.frame(
  Sample = tumor_barcodes,
  IDH1_group = idh1_status
)
rownames(pheno_df) <- pheno_df$Sample





# Download and prepare the data
GDCdownload(query_met, files.per.chunk = 10, directory = "~/Desktop/R/TCGA/meth/")

brca_met <- GDCprepare(query_met, summarizedExperiment = TRUE, directory = "~/Desktop/R/TCGA/meth/")

# save summarized experiment object
save(brca_met, file = "~/Desktop/R/TCGA/meth_BRCA_SummarizedExperiment.RData")
rm(brca_met)
# load summarized experiment object 
load("meth_BRCA_SummarizedExperiment.RData")
brca_met
# Extract beta values
beta <- assay(brca_met)

# Get sample metadata
metadata <- colData(brca_met)

a <- as.data.frame(assay(brca_met))
write.csv(a, "meth-TCGA.csv")
b <- rowData(brca_met)
write.csv(b, "rowData-meth-TCGA.csv")

# Step 2: Query mutation data for IDH1
query_mut <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "MuTect2"
)

mRNA_df <- read.csv("TCGA-transcriptome.csv", row.names = 1)
rowData <- read.csv("rowData-TCGA-transcriptome.csv")

GDCdownload(query_mut)
maf <- GDCprepare(query_mut)

setwd("~/Desktop/R/TCGA") 
# Load required packages
library(TCGAbiolinks)    # TCGA data access
library(SummarizedExperiment) 
library(minfi)           # Methylation QC
library(limma)           # Differential analysis
library(DESeq2)          # RNA-seq analysis
library(ggplot2)         # Visualization
library(ggpubr)          # Publication-ready plots
library(ComplexHeatmap)  # Heatmaps
library(sesame)          # Methylation probe filters



#.....Step 1: Load data and match sample IDs....

# Load expression and methylation data
load("~/Desktop/R/TCGA/Exp_BRCA_matched_120.RData")   # expression_data
load("~/Desktop/R/TCGA/meth_BRCA_matched_120.RData")  # methylation_data

# Standardize sample names to 15-character short TCGA barcodes
shorten_tcga_barcode <- function(barcodes) {
  substr(barcodes, 1, 15)
}

colnames(expr_data) <- shorten_tcga_barcode(colnames(expr_data))
colnames(brca_met) <- shorten_tcga_barcode(colnames(brca_met))

# Find common samples
common_samples <- intersect(colnames(expr_data), colnames(brca_met))
expr_data <- expr_data[, common_samples]
brca_met <- brca_met[, common_samples]
dim(expr_data)
dim(brca_met)

# Extract tumor/normal labels using barcode positions 14-15
get_sample_type <- function(barcode) substr(barcode, 14, 15)

sample_types <- sapply(common_samples, get_sample_type)
tumor_samples <- common_samples[sample_types == "01"]
normal_samples <- common_samples[sample_types == "11"]

# Match tumor-normal pairs by patient ID (first 12 characters)
get_patient_id <- function(barcode) substr(barcode, 1, 12)

tumor_patients <- get_patient_id(tumor_samples)
normal_patients <- get_patient_id(normal_samples)
matched_ids <- intersect(tumor_patients, normal_patients)


# Create a dataframe of matched tumor/normal pairs
matched_pairs <- data.frame(
  PatientID = matched_ids,
  Tumor = sapply(matched_ids, function(id) tumor_samples[get_patient_id(tumor_samples) == id][1]),
  Normal = sapply(matched_ids, function(id) normal_samples[get_patient_id(normal_samples) == id][1]),
  stringsAsFactors = FALSE
)
dim(matched_pairs)
#...... Maping Gene IDs...............#
# Install and load biomaRt if needed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

# Connect to Ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Map gene symbol "IDH1" to Ensembl gene ID
idh1_mapping <- getBM(
  filters = "hgnc_symbol",
  values = "IDH1",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart
)

idh1_ensembl <- idh1_mapping$ensembl_gene_id[1]

# Remove version numbers
# Extract IDH1 expression (clean ENSEMBL IDs)
rownames(expr_data) <- sub("\\..*", "", rownames(expr_data))


# Check if Ensembl ID exists in expression data
if (!(idh1_ensembl %in% rownames(expr_data_qc))) {
  stop("IDH1 Ensembl ID not found in expression_data")
}


# Filter out lowly expressed genes
expr_mat <- assay(expr_data)
dim(expr_mat)

library(edgeR)

# Calculate CPM (Counts per Million)
cpm_mat <- cpm(expr_mat)

# Keep genes expressed at CPM > 1 in at least 10% of samples
keep_genes <- rowSums(cpm_mat > 1) >= ncol(cpm_mat) * 0.10
expr_mat_filtered <- expr_mat[keep_genes, ]

# Sample-wise total expression
total_expr <- colSums(expr_mat_filtered)

# Identify low-quality samples (e.g., bottom 1%)
low_expr_samples <- names(total_expr[total_expr < quantile(total_expr, 0.01)])

# Remove them from the expression matrix
expr_mat_filtered <- expr_mat_filtered[, !(colnames(expr_mat_filtered) %in% low_expr_samples)]
# Keep only high-quality genes and samples in orignal summarized experimet
expr_data_qc <- expr_data[keep_genes, !(colnames(expr_data) %in% low_expr_samples)]
dim(expr_data_qc)

# Get IDH1 expression for tumors and normals
matched_pairs$Tumor %in% colnames(expr_data_qc)
valid_samples <- colnames(expr_data_qc) # Current barcodes in filtered expression data

# Keep only pairs where both tumor and normal are in expr_data_qc
matched_pairs <- matched_pairs[
  matched_pairs$Tumor %in% valid_samples &
    matched_pairs$Normal %in% valid_samples,
]
 dim(matched_pairs) ## should 100 samples 
 
idh1_tumor <- assay(expr_data_qc)[idh1_ensembl, matched_pairs$Tumor]  # Correct if barcodes match
idh1_normal <- assay(expr_data_qc)[idh1_ensembl, matched_pairs$Normal]
a <- as.data.frame(idh1_tumor)
b<- as.data.frame(idh1_normal)
a$IDH1 <- ifelse(a$idh1_tumor > median(a$idh1_tumor), "High", "Low")
b$IDH1 <- ifelse(b$idh1_normal > median(b$idh1_normal), "High", "Low")

#quantile to distribute data (Tumor)
q_high <- quantile(a$idh1_tumor, 0.85, na.rm = TRUE)
q_low  <- quantile(a$idh1_tumor, 0.15, na.rm = TRUE)

a$IDH_mean <- ifelse(a$idh1_tumor >= q_high, "High",
                     ifelse(a$idh1_tumor <= q_low, "Low", "Medium"))
#quantile to distribute data (Nrmal)
q_high_n <- quantile(b$idh1_normal, 0.85, na.rm = TRUE)
q_low_n  <- quantile(b$idh1_normal, 0.15, na.rm = TRUE)

b$IDH_mean <- ifelse(b$idh1_normal >= q_high_n, "High",
                     ifelse(b$idh1_normal <= q_low_n, "Low", "Medium"))
# Calculate ΔIDH1 (Tumor - Normal)
delta_idh1 <- idh1_tumor - idh1_normal

# ΔIDH1 = Tumor - Normal
matched_pairs$delta_IDH1 <- idh1_tumor - idh1_normal
rm(delta_idh1)
library(ggpubr)

idh1_exp <- assay(expr_data_qc)[idh1_ensembl,]# Include comma for all samples
clinical <- colData(expr_data_qc)
idh1_exp <- as.data.frame(idh1_exp)
# Group samples
groups <- ifelse(clinical$sample_type == "Primary Tumor", "Tumor", "Normal")

#  Boxplot for Tumor vs Normal
tiff("IDH1_expression_TCGA.tiff", units= "in", width= 10, height = 4, res = 300)
ggplot(idh1_exp, aes(x = groups, y = idh1_exp, fill = groups)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "IDH1 expression: Tumor vs Normal",
       x = "Sample Type",
       y = "Expression") +
  stat_compare_means(method = "wilcox.test") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"))
dev.off()


# Stratify based on median ΔIDH1
median_delta <- median(matched_pairs$delta_IDH1, na.rm = TRUE)

matched_pairs$IDH1_group <- ifelse(
  matched_pairs$delta_IDH1 > median_delta,
  "High_dIDH1",
  "Low_dIDH1"
)
dim(matched_pairs)

#quantile to distribute data (Delata IDH1)
q_high_dIDH1 <- quantile(matched_pairs$delta_IDH1, 0.85, na.rm = TRUE)
q_low_dIDH1  <- quantile(matched_pairs$delta_IDH1, 0.15, na.rm = TRUE)

matched_pairs$Quat_dIDH1 <- ifelse(matched_pairs$delta_IDH1 >= q_high_dIDH1, "High",
                     ifelse(matched_pairs$delta_IDH1 <= q_low_dIDH1, "Low", "Medium"))
str(matched_pairs$delta_IDH1)


kmeans_result <- kmeans(matrix(matched_pairs$delta_IDH1, ncol = 1), centers = 3) # You can adjust 'centers' based on the desired number of clusters
clusters <- factor(kmeans_result$cluster)
# Step 3: Create a data frame with original ΔIDH1 values and clusters
cluster_df <- data.frame(
  Sample = row.names(matched_pairs),  # If delta_IDH1 is named
  delta_IDH1 = matched_pairs$delta_IDH1,
  Cluster = clusters
)

# Step 4: Assign proper labels (High, Medium, Low)
# First, determine mean ΔIDH1 per cluster
cluster_means <- tapply(cluster_df$delta_IDH1, cluster_df$Cluster, mean)

# Sort clusters by increasing mean
sorted_clusters <- names(sort(cluster_means))  # e.g., "2", "3", "1"

# Assign labels
cluster_labels <- setNames(c("dIDH1_Low", "dIDH1_Medium", "dIDH1_High"), sorted_clusters)

# Apply labels to Cluster column
cluster_df$dIDH1_Group <- cluster_labels[as.character(cluster_df$Cluster)]


##.....Prepare group info for methylation analysis.....
###.....
# Prepare metadata for methylation analysis
group_info <- data.frame(
  Sample = c(matched_pairs$Normal, matched_pairs$Tumor),
  IDH1 =c(b$IDH1,a$IDH1),
  IDH1_Quant = c(b$IDH_mean,a$IDH_mean),
  Type = rep(c("Normal", "Tumor"), each = nrow(matched_pairs)),
  PatientID = rep(matched_pairs$PatientID, 2),
  IDH1_group = rep(matched_pairs$IDH1_group, 2),
  dIDH1_Quant = rep(matched_pairs$Quat_dIDH1, 2),
  Cluster = rep(cluster_df$dIDH1_Group,2),
  stringsAsFactors = FALSE
)

dim(group_info)
#write.csv(group_info, "TCGA_IDH1_Group_info.csv")

#.......... ✅ Step 5: QC and Differential Methylation Analysis
#..............................................................

# Load annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(annotation)
# Extract methylation matrix
meth_mat <- assay(brca_met)
head(meth_mat)
dim(meth_mat)
# Keep only probes present in annotation
common_probes <- intersect(rownames(meth_mat), rownames(annotation))
meth_mat <- meth_mat[common_probes, ]
annotation <- annotation[common_probes, ]

# 1. Remove probes with NA
meth_mat <- meth_mat[complete.cases(meth_mat), ]
# 2. Remove chrX and chrY probes

sex_chr <- annotation[rownames(meth_mat), "chr"] %in% c("chrX", "chrY")
meth_mat <- meth_mat[!sex_chr, ]

# 3. Remove SNP-related probes
data("SNPs.137CommonSingle")  # From minfiData
snp_probes <- as.character(SNPs.137CommonSingle$name)
meth_mat <- meth_mat[!rownames(meth_mat) %in% snp_probes, ]

# 4. Remove cross-reactive probes (if you have a list)
crs.reac <- read.csv("/Users/Saeed_1/Desktop/R/cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]
meth_mat <- meth_mat[!rownames(meth_mat) %in% crs.reac, ]

#.....Prepare meta data for stratified methylation analysis
..........
# Prepare metadata: Tumor samples only
dim(matched_pairs) 
dim(group_info)
common <- intersect(colnames(meth_mat), group_info$Sample)
meth_mat <- meth_mat[,common]
dim(meth_mat)

all(colnames(meth_mat) %in% group_info$Sample)  # Should return TRUE
length(colnames(meth_mat)) == length(group_info$Sample)  # Should be TRUE

#barcodes <- readLines("~/Desktop/R/TCGA/selected_TCGA_barcodes.txt")
# Make sure the matched tumor and normal samples exist in methylation matrix
valid_pairs_meth <- matched_pairs$Tumor %in% colnames(meth_mat) & matched_pairs$Normal %in% colnames(meth_mat)
paired_samples <- matched_pairs[valid_pairs_meth, ]


# Subset methylation matrix to just tumor and normal samples
tumor_mat <- meth_mat[, paired_samples$Tumor]
normal_mat <- meth_mat[, paired_samples$Normal]

# Ensure columns are aligned for subtraction
tumor_mat <- tumor_mat[, match(paired_samples$Tumor, colnames(tumor_mat))]
normal_mat <- normal_mat[, match(paired_samples$Normal, colnames(normal_mat))]

# Compute Δmethylation matrix (Tumor - Normal)
delta_meth <- tumor_mat - normal_mat
colnames(delta_meth) <- paired_samples$Tumor  # Label with tumor sample ID

rm(valid_samples, valid_pairs_meth)

#........box plot meth tumor vs control / IDH1 high and low...

library(dplyr)
library(ggplot2)
library(ggpubr)

# Step 1: Calculate mean methylation for each sample
sample_means <- colMeans(meth_mat, na.rm = TRUE)

# Step 2: Convert to data frame
df_means <- data.frame(
  Sample = names(sample_means),
  Mean_Methylation = sample_means
)

# Step 3: Merge with metadata
#group_info$Sample <- rownames(group_info)  # ensure Sample column exists in meta1
plot_data <- left_join(df_means, group_info, by = "Sample")

# Check structure
head(plot_data)

# Step 4: Boxplot for Tumor vs Normal
ggplot(plot_data, aes(x = Type, y = Mean_Methylation, fill = IDH1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "Global Methylation: Tumor vs Normal",
       x = "Sample Type",
       y = "Mean Beta Value") +
  stat_compare_means(method = "wilcox.test") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"))

# Recomended test for more that two categories is Kruskal. test
tiff("Glob_meth_in_IDH1_groups.tiff", units="in", width=12, height=4, res=300)
ggplot(plot_data, aes(x = Type, y = Mean_Methylation, fill = IDH1_Quant)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "Global Methylation by IDH1 Expression",
       x = "Sample Type",
       y = "Mean Beta Value") +
  stat_compare_means(method = "kruskal.test") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"))

dev.off()

#...Meth in dIDH1 quantile groups
# Recomended test for more that two categories is Kruskal. test
tiff("Glob_meth_in_IDH1_groups.tiff", units="in", width=12, height=4, res=300)
ggplot(plot_data, aes(x = Type, y = Mean_Methylation, fill = dIDH1_Quant)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(title = "Global Methylation by IDH1 Expression",
       x = "Sample Type",
       y = "Mean M-Value") +
  stat_compare_means(method = "kruskal.test") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"))

dev.off()


#..............DMR Pipeline........
# Ensure alignment of metadata and methylation data
dim(group_info)
meta1 <- group_info[match(colnames(meth_mat), group_info$Sample), ]
dim(meta1)

write.csv(meta1,"clinical_IDH1_tcga.csv")

all(colnames(meth_mat) == meta1$Sample)  # Should return TRUE
length(colnames(meth_mat)) == length(meta1$Sample)  # Should be TRUE

meta1$IDH1_group <- factor(meta1$IDH1_group)
meta1$Type <- factor(meta1$Type)
meta1$IDH1 <- factor(meta1$IDH1)
meta1$IDH1_Quant <- factor(meta1$IDH1_Quant)
#design <- model.matrix(~ Type + IDH1_group, data = meta1) #This models the effect of IDH1_group after adjusting for SampleType.
#design <- model.matrix(~ Type * IDH1_Quant, data = meta1) # If you want to test interaction (e.g., does the effect of IDH1 differ between tumor and normal)
design <- model.matrix(~ Type * IDH1, data = meta1)         
 #This expands to:
                #(Intercept): Baseline (Control + IDH1High)
                #TypeTumor: Tumor effect (vs. Control)
                #IDH1Low: IDH1Low effect (vs. IDH1High)
                #TypeTumor:IDH1Low: Interaction (tumor-specific IDH1Low effect)

#Column	                             Meaning
#(Intercept)	                       Baseline group (probably Normal + IDH1_QuantHigh, depending on ref levels)
#TypeTumor	                         Effect of being Tumor vs Normal
#IDH1_QuantLow	                     Effect of Low IDH1 vs reference (likely High)
#IDH1_QuantMedium	                   Effect of Medium IDH1 vs reference (likely High)
#TypeTumor:IDH1_QuantLow	           Interaction: Tumor-specific effect in IDH1 Low
#TypeTumor:IDH1_QuantMedium	         Interaction: Tumor-specific effect in IDH1 Medium
library(limma)

# 2. Define contrasts
colnames(design)
colnames(design)[4] <- "TypeTumor.IDH1Low"
contrasts <- makeContrasts(
  Tumor_vs_Control_IDH1High = TypeTumor,
  Tumor_vs_Control_IDH1Low = TypeTumor + `TypeTumor.IDH1Low`,
 IDH1Low_vs_High_InTumor = IDH1Low + `TypeTumor.IDH1Low`,
  Interaction = `TypeTumor.IDH1Low`,
  levels = design
)
# Quatile based IDH1 contrast
#colnames(design) <- make.names(colnames(design))
#contrasts <- makeContrasts(TvC= TypeTumor, 
  #    Low_vs_High_in_Tumor = TypeTumor + IDH1_QuantLow + TypeTumor.IDH1_QuantLow, # Tumor: IDH1 Low vs High
   #   Medium_vs_High_in_Tumor= TypeTumor + IDH1_QuantMedium + TypeTumor.IDH1_QuantMedium, # Tumor: IDH1 Medium vs High
    #  Low_vs_High_in_Normal = IDH1_QuantLow,  # Normal: IDH1 Low vs High
     # Medium_vs_High_in_Normal = IDH1_QuantMedium,  # Normal: IDH1 Medium vs High
      #  levels = design
     # )

Mval1 <- log2(meth_mat / (1 - meth_mat))
rownames(Mval1)
library(limma)

# 3. Fit model and apply contrasts
fit_meth <- lmFit(Mval1, design)
fit_meth2 <- contrasts.fit(fit_meth, contrasts)
fit_meth2 <- eBayes(fit_meth2)


# 4. Extract results ( from first contrast based on IDH1 mean)
results_TvC_High <- topTable(fit_meth2, coef = "Tumor_vs_Control_IDH1High", number = Inf)
results_TvC_Low <- topTable(fit_meth2, coef = "Tumor_vs_Control_IDH1Low", number = Inf)
results_IDH1effect <- topTable(fit_meth2, coef = "IDH1Low_vs_High_InTumor", number = Inf)
results_interaction <- topTable(fit_meth2, coef = "Interaction", number = Inf)

# 4. Extract results (IDH1 Quant contrast)
#fit_meth2$coefficients
#results_TvC <- topTable(fit_meth2, coef = "TvC", number = Inf)
#results_Low_Vs_High_Tumor <- topTable(fit_meth2, coef = "Low_vs_High_in_Tumor", number = Inf)
#results_Med_Vs_High_Tumor <- topTable(fit_meth2, coef = "Medium_vs_High_in_Tumor", number = Inf)


head(annotation)
colnames(annotation)

ann1 <- annotation[match(rownames(meth_mat),annotation$Name),
                   c(1:3,19,24, 26)]

DMPs <- topTable(fit_meth2, num=Inf, coef="IDH1Low_vs_High_InTumor", genelist=ann1) #coef = 3
head(DMPs)
write.csv(DMPs, file = "IDH1Low_vs_High_InTumor_results.csv", row.names = TRUE)
rm(DMP_TC_High)
DMP_TC_int <- topTable(fit_meth2, num=Inf, coef="Interaction", genelist=ann1) #coef = 4

##3 plot differntial methylations (cpG)
par(mfrow=c(2,5))
par(mar=c(2,2,2,2))
sapply(rownames(DMPs)[1:10], function(cpg){
  plotCpg(meth_mat, cpg=cpg, pheno=meta1$Type, ylab = "Beta values")
})
dev.off()

sapply(rownames(DMPs)[1:10], function(cpg) {
  plotCpg(meth_mat, cpg = cpg, pheno = meta1$Type, ylab = "Beta values")
  title(main = paste0("\n", DMPs[cpg, "UCSC_RefGene_Name"]), cex.main = 0.8)
})
 # improving and coloring 

library(ggplot2)
library(patchwork)

plots <- lapply(rownames(DMPs)[1:10], function(cpg) {
  gene <- ifelse(
    is.na(DMPs[cpg, "UCSC_RefGene_Name"]) || DMPs[cpg, "UCSC_RefGene_Name"] == "",
    "intergenic/ unannotated",
    DMPs[cpg, "UCSC_RefGene_Name"]
  )
  
  df <- data.frame(
    Beta = meth_mat[cpg, ],
    Group = factor(meta1$Type, levels = c("Normal", "Tumor"))
  )
  
  ggplot(df, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, shape = 21, size = 1.5, aes(color = Group), alpha = 0.7) +
    scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
    scale_color_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
    labs(
      title = paste0(cpg, "\n", gene),
      y = "Beta values",
      x = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      legend.position = "none"
    )
})

# Combine into a 2-row, 5-column layout
wrap_plots(plots, nrow = 2)

ggsave("/Users/Saeed_1/Desktop/R/TCGA/IDH1_top 10_Cpg.tiff", width = 12, height = 4, dpi = 600)


plotCpg(meth_mat, cpg="cg05268564", pheno=meta1$Type, ylab = "Beta values")
# making a volcano plot

library(ggplot2)
library(reshape2)  # For melting the data
library(ggpubr)



#....
# Use in DMRcate
library(DMRcate)
# Assume you have annotation data (from minfi or ExperimentHub)

annot1 <- cpg.annotate(
  object = Mval1,
  datatype = "array",
  what = "M",
  analysis.type = "differential",
  design = design,
  contrasts = TRUE,
  cont.matrix = contrasts,
 # coef = "Interaction",
  coef = "IDH1Low_vs_High_InTumor",  # adjust based on your contrast
  #coef = "High_vs_Low_dIDH1",
  genome = "hg19",
  arraytype = "450K",  # or "EPIC"
 fdr = 0.05
 )

dmr <- dmrcate(annot1, lambda = 1000, C = 2, pcutoff = 0.05)

resultsRanges1 <- extractRanges(dmr, genome = "hg19")
dmr.IDH1 <- data.frame(resultsRanges1)
write.csv(as.data.frame(resultsRanges1), "~/Desktop/R/TCGA/DMRs_BC_50_samples.csv")

# Plot first DMR
groups <- c(High="magenta", Low="forestgreen")
cols <- groups[as.character(meta1$IDH1)]
cols
names(cols) <- meta1$Sample  # assuming this column matches colnames(Mval)


DMR.plot(ranges=resultsRanges1, dmr=33575, CpGs=meth_mat, what="M",
          arraytype = "450K", phen.col=meta1$IDH1_group, genome="hg19")

....
#.... Manhattan plot
.....
library(ggplot2)
library(dplyr)
library(ggrepel)

library(stringr)  

# 1. Prepare autosome-only data
dmp_plot <- DMPs %>%
  # Standardize chromosome names
  mutate(
    chr = ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))
  ) %>%
  # Filter to only chr1-22
  filter(chr %in% paste0("chr", 1:22)) %>%
  # Create plotting variables
  mutate(
    log10p = -log10(pmax(P.Value, 1e-20)), # Handle zero p-values
    abs_logFC = abs(logFC),
    chr = factor(chr, levels = paste0("chr", 1:22)) # Proper ordering
  )

# 2. Create plot with improved labeling
tiff("Manhattan_DMPs_IDH.tiff", units="in", width=25, height=14, res=600)
ggplot(dmp_plot, aes(x = pos, y = log10p)) +
  
  # Chromosome background
  geom_rect(
    data = . %>% group_by(chr) %>% summarize(xmin = min(pos), xmax = max(pos)),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
        fill = as.numeric(chr) %% 2 == 0),
    alpha = 0.2, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("grey90", "white"), guide = "none") +
  
  # Data points
  geom_point(aes(color = Relation_to_Island, size = abs_logFC), alpha = 0.7) +
  scale_color_brewer(palette = "Set1", na.value = "grey50", name = "CpG Context") +
  scale_size_continuous(range = c(1, 3), name = "|Methylation Δβ|") +
  
  # Significance thresholds
  geom_hline(yintercept = -log10(0.05/nrow(dmp_plot)), 
             linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dotted", color = "blue") +
  
  # Top hits with improved labeling
  {
    top_hits <- dmp_plot %>% 
      arrange(P.Value) %>% 
      head(10) %>%
      mutate(label = ifelse(is.na(UCSC_RefGene_Name), 
                            probe, 
                            paste0(probe, " (", UCSC_RefGene_Name, ")")))
    
    list(
      geom_point(
        data = top_hits,
        aes(x = pos, y = log10p),
        color = "gold", size = 4, shape = 1
      ),
  
        geom_text_repel(
          data = top_hits,
          #aes(label = label),
          aes(label = ifelse(nchar(as.character(label)) > 13, 
                             stringr::str_wrap(label, width = 10), 
                             label)),  # Wrap text at 10 chars
          size = 4.3,
          angle = 15,                   # 45-degree angled labels
          #nudge_y = 0.1,                # Slight upward nudge
          direction = "both",           # Allows labels to move in all directions
          segment.angle = 20,           # Angled connecting lines
          segment.curvature = -0.1,     # Gentle curves
          segment.ncp = 3,              # Smooth curve points
          segment.inflect = TRUE,       # Better pathfinding
          box.padding = 1,              # More space around labels
          max.overlaps = Inf,           # Force all labels to show
          min.segment.length = 0,       # Always show connecting lines
          force = 3,                    # Stronger repulsion force
          point.padding = 0.5,          # Space from points
          hjust = 0,                    # Left-align angled text
          vjust = 0.4,                  # Center vertically
         bg.color = "white",           # Opaque background
          bg.r = 0.1,                   # Background corner radius
          position = position_jitter(width = 0.2, height = 0.2)  # Random starting position
      )
    )
  } +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+  # Add top margin
  
  # Faceting and themes
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  labs(
    x = "Genomic Position",
    y = "-log10(p-value)",
    title = "Genome-Wide Differential Methylation (Autosomes Only)",
    caption = paste("Top hits labeled with probe ID and gene symbol\n",
                    "Red line: Bonferroni threshold (p =", format(0.05/nrow(dmp_plot), scientific = TRUE), ")")
  ) +
theme_minimal(base_size = 20) +  # Base font size increased to 14
  theme(
    strip.text.x = element_text(
      angle = 70,          # 45-degree angle
      hjust = 0.5,         # Horizontal centering
      vjust = 0.5,         # Vertical centering
      size = 16,           # Adjust size as needed
      face = "bold"        # Optional: bold font
    ),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),  # Larger axis titles
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 16, face = "bold"),  # Chromosome labels
    panel.grid.major.x = element_blank(), legend.position = "bottom",
    plot.caption = element_text(size = 15, color = "grey40")
  )

dev.off()

#.....heatmap
#.....

str(DMPs)
# Choose a lenient p-value threshold for exploration
sig_cpgs <- DMPs[DMPs$adj.P.Val < 0.001, ]


# Confirm how many CpGs passed
nrow(sig_cpgs)

top_cpg_ids <- head(rownames(sig_cpgs),71)  # or 100
heatmap_mat <- meth_mat[top_cpg_ids, , drop = FALSE]
library(pheatmap)

# Optionally scale for better visualization
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# Annotate samples by SampleType and IDH1_group
sample_anno <- meta1[, c("IDH1", "Type")]

head(sample_anno)
# Only keep annotations for samples in heatmap
annotation_col <- sample_anno[colnames(heatmap_mat), c("Type", "IDH1")]

# Plot heatmap
tiff("Heatmap_topDMPs_IDH.tiff", units="in", width=10, height=6, res=600)
pheatmap(
  heatmap_mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  scale = "none",
  clustering_method = "ward.D2",
  main = "Top Differentially Methylated CpGs"
)

dev.off()
#.........Functional analysis.....
#.............................
library(dplyr)

# For a vector or dataframe column
DMPs %>% distinct(UCSC_RefGene_Group) %>% pull()

# Add a RegionGroup column
DMPs$RegionGroup <- case_when(
  grepl("TSS200|TSS1500|5'UTR", DMPs$UCSC_RefGene_Group, ignore.case = TRUE) ~ "Promoter",
  grepl("1stExon|Body", DMPs$UCSC_RefGene_Group, ignore.case = TRUE) ~ "Gene Body",
  grepl("3'UTR", DMPs$UCSC_RefGene_Group, ignore.case = TRUE) ~ "3'UTR",
  TRUE ~ "Other"
)

#...Plot distribution
library(ggplot2)

ggplot(DMPs, aes(x = RegionGroup)) +
  geom_bar(fill = "darkcyan") +
  theme_minimal() +
  labs(title = "Distribution of DMPs Across Functional Regions",
       x = "Functional Region",
       y = "Number of DMPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"))
#...Split by direction
DMPs$Direction <- ifelse(DMPs$logFC > 0, "Hypermethylated", "Hypomethylated")
tiff("DMP-Dist.tiff", units="in", width=10, height=4, res=600)

ggplot(DMPs, aes(x = RegionGroup, fill = Direction)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Hypermethylated" = "firebrick", "Hypomethylated" = "steelblue")) +
  theme_minimal() +
  labs(title = "IDH1-mediated altrations in DMP distribution among Functional Region",
       x = "Functional Region",
       y = "Number of DMPs",
       fill = "Methylation Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # removes major grid lines
        panel.grid.minor = element_blank(),  # removes minor grid lines
        panel.border = element_blank(),      # remove panel border
        axis.line = element_line(),          # adds axis lines (x and y)
        axis.text = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),)
dev.off()
.....
#.....Enrichment analysis...
.....

head(annotation)
ann_df <- data.frame(
  probe = rownames(annotation),
  UCSC_RefGene_Group = annotation$UCSC_RefGene_Group,
  stringsAsFactors = FALSE
)

# Define functional region categories for background
ann_df$RegionGroup <- dplyr::case_when(
  grepl("TSS200|TSS1500|5'UTR", ann_df$UCSC_RefGene_Group, ignore.case = TRUE) ~ "Promoter",
  grepl("1stExon|Body", ann_df$UCSC_RefGene_Group, ignore.case = TRUE) ~ "Gene Body",
  grepl("3'UTR", ann_df$UCSC_RefGene_Group, ignore.case = TRUE) ~ "3'UTR",
  TRUE ~ "Other"
)

DMPs$probe <- rownames(DMPs)

# Annotate DMPs with RegionGroup
DMPs_annotated <- merge(DMPs, ann_df[, c("probe", "RegionGroup")], by = "probe")
head(DMPs_annotated)

# Total number of probes per region (background)
background_counts <- table(ann_df$RegionGroup)

# DMPs per region (foreground)
# Replace the incorrect RegionGroup (from DMPs) with the correct one from ann_df
DMPs_annotated$RegionGroup <- DMPs_annotated$RegionGroup.y

DMP_counts <- table(DMPs_annotated$RegionGroup)

# Combine into matrix
enrichment_table <- data.frame(
  Region = names(background_counts),
  DMPs = as.numeric(DMP_counts[names(background_counts)]),
  Background = as.numeric(background_counts)
)

# Replace NAs (regions not present in DMPs) with 0
enrichment_table$DMPs[is.na(enrichment_table$DMPs)] <- 0

# Perform Fisher's Exact Test for each region

enrichment_results <- apply(enrichment_table, 1, function(row) {
  dmp <- as.numeric(row["DMPs"])
  bg <- as.numeric(row["Background"])
  not_dmp <- sum(enrichment_table$DMPs) - dmp
  not_bg <- sum(enrichment_table$Background) - bg
  
  fisher_mat <- matrix(c(dmp, bg, not_dmp, not_bg), nrow = 2)
  fisher_test <- fisher.test(fisher_mat)
  
  return(c(OR = fisher_test$estimate,
           p.value = fisher_test$p.value))
})

# Format results
enrichment_results_df <- cbind(
  enrichment_table,
  t(enrichment_results)
)

# Adjust p-values
enrichment_results_df$FDR <- p.adjust(enrichment_results_df$p.value, method = "BH")

enrichment_results_df

#Plot ODD's ratio This will tell you which regions are overrepresented among DMPs compared to the full array. For example:

#OR > 1: enrichment
#OR < 1: depletion
#FDR < 0.05: statistically significant 
colnames(enrichment_results_df)
colnames(enrichment_results_df)[colnames(enrichment_results_df) == "OR.odds ratio"] <- "OR"


library(ggplot2)

ggplot(enrichment_results_df, aes(x = Region, y = OR)) +
  geom_bar(stat = "identity", fill = "purple") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Enrichment of DMPs in Functional Regions",
       y = "Odds Ratio (DMP enrichment)",
       x = "Region Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"))

# Add stars based on FDR
enrichment_results_df$signif <- cut(enrichment_results_df$FDR,
                                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                    labels = c("***", "**", "*", "")
)

ggplot(enrichment_results_df, aes(x = Region, y = OR)) +
  geom_bar(stat = "identity", fill = "purple") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = signif), vjust = -0.5, size = 5) +
  theme_minimal() +
  labs(title = "Enrichment of DMPs in Functional Regions",
       y = "Odds Ratio (DMP enrichment)",
       x = "Region Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"))

#.... Forest Plot ODDs ration

library(dplyr)
# How many DMPs are present at FDR < 0.05?
sum(DMPs_annotated$adj.P.Val < 0.05, na.rm = TRUE)

# 1. Label significant vs. not
DMPs_annotated <- DMPs_annotated %>%
  mutate(sig = ifelse(adj.P.Val < 0.05, "DMP", "nonDMP"))

# 2. Count significant & non-significant in each RegionGroup
region_table <- DMPs_annotated %>%
  count(RegionGroup, sig) %>%
  tidyr::pivot_wider(names_from = sig, values_from = n, values_fill = 0)

# 3. Add totals across all regions
total_DMPs <- sum(region_table$DMP, na.rm = TRUE)
total_nonDMPs <- sum(region_table$nonDMP, na.rm = TRUE)

# Apply Fisher test per row
get_fisher_stats <- function(dmp, non_dmp, total_dmp, total_non_dmp) {
  m <- matrix(c(dmp, total_dmp - dmp, non_dmp, total_non_dmp - non_dmp),
              nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  list(OR = as.numeric(ft$estimate),
       CI_lower = ft$conf.int[1],
       CI_upper = ft$conf.int[2],
       p.value = ft$p.value)
}

# Apply to each region row
region_results <- region_table %>%
  rowwise() %>%
  mutate(stats = list(get_fisher_stats(DMP, nonDMP, total_DMPs, total_nonDMPs))) %>%
  tidyr::unnest_wider(stats) %>%
  mutate(label = sprintf("%.2f (%.2f–%.2f)", OR, CI_lower, CI_upper),
         RegionGroup = factor(RegionGroup, levels = unique(RegionGroup)))

install.packages("meta")
library(meta)

# Prepare your data (as in region_results)
region_df <- as.data.frame(region_table)


# Create the complementary group: everything *outside* the region (i.e., background)
region_df$nonRegion_DMP <- sum(region_df$DMP) - region_df$DMP
region_df$nonRegion_nonDMP <- sum(region_df$nonDMP) - region_df$nonDMP

# Run metabin (treat each region as a "study")
meta_obj <- metabin(
  event.e = DMP, 
  n.e = DMP + nonDMP, 
  event.c = nonRegion_DMP, 
  n.c = nonRegion_DMP + nonRegion_nonDMP,
  studlab = RegionGroup,
  data = region_df,
  sm = "OR", method = "Inverse"
)

# Plot forest
tiff("Forest-region-meth.tiff", units="in", width=10, height=4, res=600)
forest(meta_obj,
       comb.fixed = FALSE,      # Don't show combined effect
       comb.random = FALSE,     # Keep per-region effect only
       col.square = "black",
       col.diamond = "red",
       print.I2 = TRUE,
       print.tau2 = TRUE,
       print.Q = FALSE,
       rightcols = c("effect", "ci"),
       leftcols = c("studlab", "event.e", "event.c"),
       leftlabs = c("Region", "Sig DMPs", "Other DMPs"),
       text.random = "Overall",
       xlab = "Odds Ratio",
       #smlab = "Enrichment of DMPs by Region"
       fontsize = 10,
       backtransf = TRUE
)

dev.off()

#.........Circose plot..........
......................
install.packages("circlize")
install.packages("RColorBrewer")

library(circlize)
library(RColorBrewer)
library(dplyr)
#...Prepare data (circos)
str(DMPs_annotated)


circos_df <- DMPs_annotated %>%
  dplyr::select(probe, chr, pos, logFC, RegionGroup, Direction, sig) %>%
  filter(!is.na(chr), !is.na(pos)) %>%
  mutate(
    chr = as.character(chr),
    chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),  # ensure 'chr' prefix
    start = pos,
    end = pos + 1,
    RegionGroup = as.factor(RegionGroup),
    Direction = factor(Direction, levels = c("Hypomethylated", "Hypermethylated"))
  )
circos_df <- circos_df %>%
  mutate(y = 1)  # constant height, used for plotting

#... Deep seek solutions
# Input data
df <- circos_df
library(circlize)

df <- df[order(-abs(df$logFC)), ]  # Sort by absolute logFC
# 1. Subset top hyper/hypo probes to avoid freezing
N <- 100  # Reduced from 500 to test
dim(df)
df_top <- df[1:100, ]  # Take top 100
dim(df_top)

df_hyper <- df_top[df_top$logFC > 0, ]
df_hypo <- df_top[df_top$logFC < 0, ]


cat("Hypermethylated:", nrow(df_hyper), 
    "\nHypomethylated:", nrow(df_hypo),
    "\nTotal:", nrow(df_top)) 
df_top <- rbind(
  df_hyper[order(-df_hyper$logFC)[1:min(N, nrow(df_hyper))], ],
  df_hypo[order(df_hypo$logFC)[1:min(N, nrow(df_hypo))], ]
)

dim(df_top)

# Prepare index for links
df_top <- df_top %>%
  group_by(RegionGroup) %>%
  mutate(region_index = row_number()) %>%
  ungroup() %>%
  group_by(chr) %>%
  mutate(chrom_index = row_number()) %>%
  ungroup()

# 2. Define sectors
region_sectors <- aggregate(region_index ~ RegionGroup, df_top, max)
chrom_sectors <- aggregate(chrom_index ~ chr, df_top, max)
colnames(region_sectors) <- c("sector", "end")
colnames(chrom_sectors) <- c("sector", "end")

region_sectors$start <- 1
chrom_sectors$start <- 1

sectors_df <- rbind(region_sectors[, c("sector", "start", "end")],
                    chrom_sectors[, c("sector", "start", "end")])

# 1. Convert ALL sector-related columns to CHARACTER (removes factor issues)
df_top$RegionGroup <- as.character(df_top$RegionGroup)
df_top$chr <- as.character(df_top$chr)
sectors_df$sector <- as.character(sectors_df$sector)

# 2. Verify there are no NA sectors
sum(is.na(df_top$RegionGroup)) # Should be 0
sum(is.na(df_top$chr)) # Should be 0

png("chord_plot100LFC.png", width = 3500, height = 3500, res = 300)
# 1. Set ultra-tight layout parameters
circos.clear()
circos.par(
  cell.padding = c(0, 0, 0, 0),  # Zero padding all sides
  track.margin = c(0, 0),         # No additional margin
  gap.degree = 0.9,                 # Minimal gap between sectors
  start.degree = 90               # Standard start position
)

# 2. Ensure all sectors have minimum width
min_width <- 4  # Minimum sector width
sectors_df$end <- pmax(sectors_df$end, sectors_df$start + min_width)

# 3. Initialize with safe parameters
circos.initialize(
  factors = sectors_df$sector,
  xlim = sectors_df[, c("start", "end")]
)

# 4. Add basic track to verify sectors
circos.trackPlotRegion(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ycenter,
      CELL_META$sector.index,
      cex = 1, #size
      font = 2, #bold
      facing = "clockwise", # "inside", # # Curves text to follow circle
      niceFacing = TRUE
    )
  }
)


# Compute and draw links within defined limits
# 1. Compute RegionGroup band: 25% to75%
region_band <- sectors_df %>%
  dplyr::filter(sector %in% unique(df_top$RegionGroup)) %>%
  dplyr::mutate(
    span = end - start,
    mid_start = start + 0.25 * span,
    mid_end   = start + 0.75 * span
  ) %>%
  dplyr::select(sector, mid_start, mid_end)

# 2. Compute Chromosome band: 10% to 90%
chrom_band <- sectors_df %>%
  dplyr::filter(sector %in% unique(df_top$chr)) %>%
  dplyr::mutate(
    span = end - start,
    band_start = start + 0.10 * span,
    band_end   = start + 0.90 * span
  ) %>%
  dplyr::select(sector, band_start, band_end)

# 3. Loop and draw links
for (reg in unique(df_top$RegionGroup)) {
  reg_links <- df_top %>% filter(RegionGroup == reg)
  n_links <- nrow(reg_links)
  if (n_links == 0) next
  
  reg_band <- region_band %>% filter(sector == reg)
  if (nrow(reg_band) == 0) next
  
  # Divide 40–60% span into N equal slots
  reg_positions <- seq(reg_band$mid_start, reg_band$mid_end, length.out = n_links + 1)
  
  for (i in seq_len(n_links)) {
    chr <- reg_links$chr[i]
    chrom_index <- reg_links$chrom_index[i]
    if (!(chr %in% chrom_band$sector)) next
    
    chr_band <- chrom_band %>% filter(sector == chr)
    
    # Clamp chromosome point to [5%, 95%] and keep 1 unit width
    chrom_start <- max(chrom_index, chr_band$band_start)
    chrom_end   <- min(chrom_index + 1, chr_band$band_end)
    
    # Region point: one slot between reg_positions[i] and reg_positions[i+1]
    point1 <- c(reg_positions[i], reg_positions[i + 1])
    point2 <- c(chrom_start, chrom_end)
    
    # Draw link
    circos.link(
      sector.index1 = reg,
      point1 = point1,
      sector.index2 = chr,
      point2 = point2,
      col = ifelse(reg_links$logFC[i] > 0, "#FF0000DD", "#0000FFDD"),
      border = NA,
      lwd = 0.8,
      rou1 = 0.8,
      rou2 = 0.8
    )
  }
}

# Add legend
legend(
  x = 0.45,                   # Right side outside the plot
  y = 1,                    # Vertical center
  legend = c("Hypermethylated", "Hypomethylated"),
  fill = c("#E41A1C", "#0000FFDD"),  # Match your link colors
  border = NA,
  bty = "n",                  # No legend box
  xpd = TRUE,                 # Allow plotting outside main area
  cex = 1.2,                  # Font size
  title = "Top 100 differential Methylated probes",
  title.adj = 0.25
)

# For PDF/PNG output, adjust coordinates:
# legend(x = 0.85, y = 1.1, ...)  # Typical position for saved plots
dev.off()

#................................................
 #...... Methyaltion frequency Circose plot....
#................................................
#.................

library(dplyr)
library(circlize)

# Prepare data
freq_table <- circos_df %>%
  filter(RegionGroup != "Gene Body") %>%
  mutate(
    chr = factor(paste0("chr", gsub("chr", "", chr))), 
    RegionGroup = ifelse(RegionGroup == "Promoter", "Promoter", "Non-Promoter")
  ) %>%
  dplyr::count(chr, RegionGroup, Direction, name = "Freq")

chromosomes <- levels(freq_table$chr)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 2,
  track.margin = c(0.02, 0.02), # Slightly larger margin
  cell.padding = c(0, 0, 0, 0)
)

# 1. Simulate ideogram first
png("IDH1_media_Meth_freq_Pro.png", width = 3500, height = 3500, res = 300)

circos.initialize(
  factors = chromosomes,
  xlim = cbind(rep(0, length(chromosomes)), rep(1, length(chromosomes)))
)

# 2. Plot Promoter track (outside ideogram)

{region_type <- "Promoter"
  data_subset <- freq_table %>%
    filter(RegionGroup == region_type) %>%
    mutate(chr = factor(chr, levels = chromosomes))
  ymax <- max(data_subset$Freq, na.rm = TRUE) * 1.3
  track_index <- NULL
  
  circos.trackPlotRegion(
    factors = data_subset$chr,
    ylim = c(0, ymax),
    track.height = 0.15,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      current_data <- data_subset %>% filter(as.character(chr) == chr)
      
      circos.lines(c(0, 1), c(0, 0), col = "black", lwd = 0.5)
      
      if (nrow(current_data) > 0) {
        if ("Hypermethylated" %in% current_data$Direction) {
          val <- current_data$Freq[current_data$Direction == "Hypermethylated"]
          circos.rect(0.2, 0, 0.4, val, col = "#e41a1c", border = NA)
        }
        if ("Hypomethylated" %in% current_data$Direction) {
          val <- current_data$Freq[current_data$Direction == "Hypomethylated"]
          circos.rect(0.6, 0, 0.8, val, col = "#377eb8", border = NA)
        }
      }
      
      circos.yaxis(
        side = "left",
        at = c(0, round(ymax / 2), round(ymax)),
        labels.cex = 0.3,
        tick.length = 0.02,
        labels.niceFacing = TRUE
      )
      
      track_index <<- CELL_META$track.index
    }
  )
  
  circos.text(
    x = 0.5,
    y = ymax * 0.9,
    labels = region_type,
    sector.index = chromosomes[1],
    track.index = track_index,
    cex = 0.8,
    font = 2,
    niceFacing = TRUE
  )
}

# 3. Plot Non-Promoter track (inside ideogram)
{region_type <- "Non-Promoter"
  data_subset <- freq_table %>%
    filter(RegionGroup == region_type) %>%
    mutate(chr = factor(chr, levels = chromosomes))
  
  ymax <- max(data_subset$Freq, na.rm = TRUE) * 1.3
  track_index <- NULL
  
  circos.trackPlotRegion(
    factors = data_subset$chr,
    ylim = c(0, ymax),
    track.height = 0.15,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      current_data <- data_subset %>% filter(as.character(chr) == chr)
      
      circos.lines(c(0, 1), c(0, 0), col = "black", lwd = 0.5)
      
      if (nrow(current_data) > 0) {
        if ("Hypermethylated" %in% current_data$Direction) {
          val <- current_data$Freq[current_data$Direction == "Hypermethylated"]
          circos.rect(0.2, 0, 0.4, val, col = "#e41a1c", border = NA)
        }
        if ("Hypomethylated" %in% current_data$Direction) {
          val <- current_data$Freq[current_data$Direction == "Hypomethylated"]
          circos.rect(0.6, 0, 0.8, val, col = "#377eb8", border = NA)
        }
      }
      
      circos.yaxis(
        side = "left",
        at = c(0, round(ymax / 2), round(ymax)),
        labels.cex = 0.3,
        tick.length = 0.02,
        labels.niceFacing = TRUE
      )
      
      track_index <<- CELL_META$track.index
    }
  )
  
  circos.text(
    x = 0.5,
    y = ymax * 0.9,
    labels = region_type,
    sector.index = chromosomes[1],
    track.index = track_index,
    cex = 0.8,
    font = 2,
    niceFacing = TRUE
  )
} 



#4. ...ideogram

circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    circos.rect(0, 0, 1, 1, col = "#F5F5F5", border = "black")
    circos.text(0.5, 0.5, chr, cex = 0.8, font = 2,niceFacing = TRUE)
  }
)

#legend(
# x = 0.35,                   # Right side outside the plot
# y = 1,                    # Vertical center
#legend = c("Hypermethylated", "Hypomethylated"),
# fill = c("#e41a1c", "#377eb8"),  # Match your link colors
#  border = NA,
# bty = "n",                  # No legend box
# xpd = TRUE,                 # Allow plotting outside main area
#cex = 1.2,                  # Font size
# title = "Genome-wide methylation in regulatory regions",
#  title.adj = 0.25
#)

# Create legend object
lgd <- Legend(
  labels = c("Hypermethylated", "Hypomethylated"),
  legend_gp = gpar(fill = c("#e41a1c", "#377eb8")),
  title = "Genome-Wide Methylation\nin Regulatory Regions",
  title_gp = gpar(fontface = "bold",fontsize = 12 ),
  # label_gp = gpar(fontsize = 10),
  title_position = "topcenter"
)

# Draw in center
circos.clear()  # Make sure to call this AFTER your circos plot is fully drawn
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")

# For PDF/PNG output, adjust coordinates:
# legend(x = 0.85, y = 1.1, ...)  # Typical position for saved plots
dev.off()

circos.clear()

### mCSEA for enrichment

...........

# Read data
DMPs <- read.csv("IDH1Low_vs_High_InTumor_results.csv", header=TRUE)

# Create t-stat vector with correct names
t_stat_vector <- DMPs$t
names(t_stat_vector) <- DMPs$X  # Use CpG IDs from column X

# Read and prepare phenotype
meta1 <- read.csv("clinical_IDH1_tcga.csv", header=TRUE)
rownames(meta1) <- meta1$Sample

# Find common samples and CpGs
common_samples <- intersect(colnames(meth_mat), rownames(meta1))
common_cpgs <- intersect(names(t_stat_vector), rownames(meth_mat))

# Subset data
meth_mat <- meth_mat[common_cpgs, common_samples]
t_stat_vector <- t_stat_vector[common_cpgs]

# Create phenotype as a DATA FRAME (not vector)
pheno <- meta1[common_samples, "IDH1", drop = FALSE]  # Keep as data frame
pheno$IDH1 <- factor(pheno$IDH1)  # Convert to factor within the data frame

# Check the structure
str(pheno)
table(pheno$IDH1)

install.packages("curl")
BiocManager::install("mCSEA")
library(curl)
library(mCSEA)

# Run mCSEA
res_mCSEA <- mCSEATest(
  t_stat_vector,
  meth_mat,
  pheno,
  minCpGs = 5,
  regionsTypes = "promoters",
  platform = "450k"
)

# Check results
summary(res_mCSEA)

head(res_mCSEA$promoters) #The leading edge CpGs are the real drivers of the ES; these can be considered the most important CpGs with the largest logFC. 
#Kolmogorov-Smirnov (KS) Test The KS test is a non-parametric statistical test used to compare two distributions, often employed in gene set enrichment analysis (GSEA).

# Extract CpGs from the leading edge of IDH1
leading_cpgs <- unlist(strsplit(res_mCSEA$promoters["IDH1", "leadingEdge"], ", "))
# Check all regions for IDH1

prom <- as.data.frame(res_mCSEA$promoters)
prom$Gene <- rownames(prom)

str(prom)

write.csv(prom, "prom_IDH1_tcga.csv")

mCSEAPlot(
  res_mCSEA,
  regionType = "promoters",
  genes = FALSE,
  dmrName = "IDH1",
transcriptAnnotation = "symbol",
  makePDF = TRUE) # Save as pdf

#.... Top N bar plot of NES
library(dplyr)
library(ggplot2)

# Get top 20 most significant
top_prom <- prom %>%
 # mutate(Gene = rownames(prom)) %>%  # If rownames are genes; skip if Gene column already exists
  arrange(padj) %>%
  dplyr::slice(1:20)

ggplot(top_prom, aes(x = reorder(Gene, NES), y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Top 20 Enriched Promoter Sets",
       x = "Promoter Set",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.title = element_text(, size= 12, face = "bold"), 
        axis.text.y = element_text(size = 10, color = "black", family = "arial")
        # Change size or other properties here
  )
ggsave("T20_enr_pro.png", width = 6, height = 5, unit= "in", dpi = 300, bg= "white")

#.... Volcano Plot (NES vs -log10 FDR)
ggplot(prom, aes(x = NES, y = -log10(padj))) +
  geom_point(alpha = 0.4, color = "darkblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = "mCSEA Results: Promoter Enrichment",
       x = "Normalized Enrichment Score (NES)",
       y = "-log10 Adjusted P-value") +
  theme_minimal()

#... . Dot Plot: NES vs Size
ggplot(prom, aes(x = NES, y = size)) +
  geom_point(aes(color = -log10(padj)), alpha = 0.7) +
  scale_color_gradient(low = "gray", high = "red") +
  labs(title = "mCSEA Enrichment: Promoter Sets",
       x = "Normalized Enrichment Score",
       y = "Region Size (Number of CpGs)",
       color = "-log10 FDR") +
  theme_minimal()

#...... Gene NEtwork analysis and leading CpG clustering
#.....................................................
#.....1. Prepare Gene CpG maping table

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)

# Build the data cleanly
gene_cpg_df <- do.call(rbind, lapply(seq_len(nrow(prom)), function(i) {
  gene <- as.character(prom$Gene[i])
  cpgs <- unlist(strsplit(gsub("\\s+", "", prom$leadingEdge[i]), ","))
  data.frame(Gene = rep(gene, length(cpgs)), CpG = cpgs, stringsAsFactors = FALSE)
}))

# Convert columns to character explicitly (in case of lingering list types)
gene_cpg_df$Gene <- as.character(gene_cpg_df$Gene)
gene_cpg_df$CpG <- as.character(gene_cpg_df$CpG)

# Confirm column type
str(gene_cpg_df)

gene_cpg_df <- as.data.frame(gene_cpg_df)  # Just to be sure


# Step 2: Check how many genes each CpG links to
cpg_counts <- gene_cpg_df %>%
  dplyr::count(CpG, name = "GeneCount")

shared_cpgs <- cpg_counts %>% # calculate CpG with 3 or more connections
  filter(GeneCount >= 3)

gene_cpg_shared <- gene_cpg_df %>%
  filter(CpG %in% shared_cpgs$CpG)

# ...Build and Plot the Network

# Create edge list
edges <- gene_cpg_shared  # or full gene_cpg_df

# Create node list
nodes <- data.frame(
  name = unique(c(edges$Gene, edges$CpG)),
  type = ifelse(unique(c(edges$Gene, edges$CpG)) %in% edges$Gene, "Gene", "CpG")
)

#...... vis netwoork interactive map
install.packages("visNetwork")
library(visNetwork)


# Prepare nodes data frame with required columns
vis_nodes <- data.frame(
  id = unique(c(edges$Gene, edges$CpG)),
  label = unique(c(edges$Gene, edges$CpG)),
  group = ifelse(unique(c(edges$Gene, edges$CpG)) %in% edges$Gene, "Gene", "CpG"),
  stringsAsFactors = FALSE
)

# Prepare edges data frame with required columns
vis_edges <- data.frame(
  from = edges$Gene,
  to = edges$CpG,
  stringsAsFactors = FALSE
)

# Create the interactive network
visNetwork(vis_nodes, vis_edges) %>%
  visOptions(highlightNearest = TRUE, 
             nodesIdSelection = TRUE) %>%
  visGroups(groupname = "Gene", color = "steelblue") %>%
  visGroups(groupname = "CpG", color = "tomato") %>%
  visLegend() %>%
  visPhysics(solver = "forceAtlas2Based",  # Better layout for large networks
             forceAtlas2Based = list(gravitationalConstant = -50))

#....... using plot and igraph

library(igraph)

# Create igraph object. Simple network
g <- graph_from_data_frame(d = vis_edges, vertices = vis_nodes, directed = FALSE)

# Assign visual properties
V(g)$color <- ifelse(V(g)$group == "Gene", "steelblue", "tomato")
V(g)$label.color <- "black"
V(g)$size <- 5  # smaller nodes
V(g)$label.cex <- 0.5  # smaller text
V(g)$frame.color <- "gray90"

# Improved layout: Kamada-Kawai
set.seed(42)
layout_clean <- layout_with_kk(g)

# Open PNG with higher resolution
png("gene_cpg_network_cleaner.png", width = 1800, height = 1600, res = 200)

# Plot cleanly
plot(g,
     layout = layout_clean,
     vertex.shape = "circle",
     edge.color = rgb(0.5, 0.5, 0.5, 0.3),  # soft gray with transparency
     edge.width = 1,
     main = "Gene–CpG Network"
)

# Add legend
legend("bottomleft",
       legend = c("Gene", "CpG"),
       col = c("steelblue", "tomato"),
       pch = 21,
       pt.bg = c("steelblue", "tomato"),
       pt.cex = 2,
       bty = "n"
)

dev.off()

#.....hub genes highlight and gene CpG cluster

library(igraph)

# Create graph
g <- graph_from_data_frame(d = vis_edges, vertices = vis_nodes, directed = FALSE)

# Compute degree (number of connections per node)
deg <- degree(g)

a <- as.data.frame(deg)
write.csv(a, "connections_genes_cpg_tcga.csv")

# Identify top 10 hub genes
hub_genes <- names(sort(deg[V(g)$group == "Gene"], decreasing = TRUE))[1:10]

# Assign default node properties
V(g)$color <- ifelse(V(g)$group == "Gene", "steelblue", "tomato")
V(g)$size <- 7
V(g)$label.cex <- 0.7
V(g)$frame.color <- "gray90"
V(g)$label.color <- "black"

# Highlight hub genes
V(g)$color[V(g)$name %in% hub_genes] <- "gold"
V(g)$size[V(g)$name %in% hub_genes] <- 18
V(g)$label.cex[V(g)$name %in% hub_genes] <- 1
V(g)$label.font <- ifelse(V(g)$name %in% hub_genes, 2, 1)  # bold font for hubs

# Layout and plot
set.seed(123)
layout_clean <- layout_with_kk(g)
png("gene_cpg_network_hubs1.png", width = 3000, height = 3000, res = 200)

plot(g,
     layout = layout_clean,
     vertex.shape = "circle",
     edge.color = rgb(0.5, 0.5, 0.5, 0.3),
     edge.width = 1.5,
     main = "Gene–CpG Network with Highlighted Hub Genes"
)

# Add legend
legend("bottomleft",
       legend = c("Hub Genes", "Other Genes", "CpGs"),
       col = c("gold", "steelblue", "tomato"),
       pch = 21,
       pt.bg = c("gold", "steelblue", "tomato"),
       pt.cex = 2,
       bty = "n"
)

dev.off()

#... separate hub genes and label nodes
#.....
library(igraph)

# Build the full graph
g_full <- graph_from_data_frame(d = vis_edges, vertices = vis_nodes, directed = FALSE)

# Identify hub genes (top 10 by degree)
is_gene <- V(g_full)$group == "Gene"
gene_degrees <- degree(g_full)[is_gene]

# Step 1: Get hub gene vertices by name
hub_genes <- names(sort(degree(g_full)[V(g_full)$group == "Gene"], decreasing = TRUE))[1:10]

# Step 2: Get all neighbors of hub genes (CpGs), convert to names
hub_neighbors <- unlist(neighborhood(g_full, 1, nodes = hub_genes))
hub_neighbor_names <- V(g_full)$name[hub_neighbors]

# Step 3: Combine and get unique node names
hub_nodes <- unique(c(hub_genes, hub_neighbor_names))

# Step 4: Subgraph using valid vertex names
g_sub <- induced_subgraph(g_full, vids = V(g_full)[name %in% hub_nodes])


# Set visual properties
V(g_sub)$shape <- "circle"
V(g_sub)$size <- ifelse(V(g_sub)$name %in% hub_genes, 18, 8)
V(g_sub)$color <- ifelse(V(g_sub)$name %in% hub_genes, "gold", "steelblue")
V(g_sub)$label <- V(g_sub)$name
V(g_sub)$label.cex <- 0.7
V(g_sub)$label.font <- 2
V(g_sub)$label.color <- "black"
V(g_sub)$frame.color <- "gray80"


# Identify CpGs in the subgraph
is_cpg <- V(g_sub)$group == "CpG"

# Count how many hub genes each CpG is connected to
cpg_neighbors <- sapply(V(g_sub)$name[is_cpg], function(cpg) {
  connected_genes <- neighbors(g_sub, cpg)$name
  sum(connected_genes %in% hub_genes)
})

# Store into a named vector for easy lookup
cpg_counts <- setNames(cpg_neighbors, V(g_sub)$name[is_cpg])

# Assign colors:
V(g_sub)$color <- ifelse(
  V(g_sub)$name %in% hub_genes, "gold",
  ifelse(
    V(g_sub)$name %in% names(cpg_counts) & cpg_counts[V(g_sub)$name] > 3,
    "bisque",
    "steelblue"
  )
)


# Layout: arrange each hub gene cluster separately
layout_list <- lapply(hub_genes, function(hub) {
  nodes <- c(hub, neighbors(g_sub, hub)$name)
  sub <- induced_subgraph(g_sub, nodes)
  layout <- layout_in_circle(sub)
  layout <- layout + matrix(runif(length(nodes)*2, -0.5, 0.5), ncol = 2)  # jitter to avoid exact overlap
  rownames(layout) <- nodes
  layout
})

# Combine layouts with spacing between clusters
cluster_spacing <- 30
full_layout <- matrix(NA, nrow = vcount(g_sub), ncol = 2)
rownames(full_layout) <- V(g_sub)$name
x_offset <- 0

for (lay in layout_list) {
  lay[,1] <- lay[,1] + x_offset
  full_layout[rownames(lay), ] <- lay
  x_offset <- x_offset + cluster_spacing
}

# Plot
png("hub_gene_clusters_labeled1.png", width = 4500, height = 4000, res = 300)

plot(g_sub,
     layout = full_layout,
     vertex.shape = V(g_sub)$shape,
     vertex.color = V(g_sub)$color,
     vertex.size = V(g_sub)$size,
     vertex.label = V(g_sub)$label,
     vertex.label.cex = V(g_sub)$label.cex,
     vertex.label.font = V(g_sub)$label.font,
     vertex.label.color = V(g_sub)$label.color,
     edge.color = "gray60",
     edge.width = 1.2,
     main = "Hub Gene Clusters with Connected CpGs"
)

dev.off()

