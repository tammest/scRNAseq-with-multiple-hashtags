```r
# ============================================================
# 
#
# Create a Seurat object from Cell Ranger output, separate RNA,
# antibody hashtag, and lipid hashtag count matrices, perform
# RNA QC, filtering, normalization, variable feature selection,
# and scaling.
# ============================================================

library(Seurat)
library(tidyverse)

# -----------------------------
# Paths
# -----------------------------

input_h5 <- "data/filtered_feature_bc_matrix.h5"

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load Cell Ranger data
# -----------------------------

kody.sc <- Read10X_h5(filename = input_h5)

# The Cell Ranger H5 contains RNA and hashtag/feature barcode data.
# Adjust indices if your object structure is different.
kody.rna <- as.matrix(kody.sc[[1]])
kody.htos <- as.matrix(kody.sc[[2]])

# First 3 hashtags are antibody hashtags.
# Rows 4-9 are lipid hashtags.
kody.htos.antibody <- kody.htos[1:3, ]
kody.htos.lipid <- kody.htos[4:9, ]

# Keep only shared cell barcodes between RNA and hashtag matrices.
joint.bcs <- intersect(colnames(kody.rna), colnames(kody.htos))

kody.rna <- kody.rna[, joint.bcs]
kody.htos.antibody <- kody.htos.antibody[, joint.bcs]
kody.htos.lipid <- kody.htos.lipid[, joint.bcs]

# Check hashtag names
print(rownames(kody.htos.antibody))
print(rownames(kody.htos.lipid))

# -----------------------------
# Create Seurat object
# -----------------------------

kody.hashtag <- CreateSeuratObject(counts = kody.rna)

# -----------------------------
# RNA quality control
# -----------------------------

counts_data <- GetAssayData(
  kody.hashtag,
  assay = "RNA",
  layer = "counts"
)

# Mouse mitochondrial genes usually start with "mt-"
mt_genes <- rownames(counts_data)[grepl("^mt-", rownames(counts_data))]

kody.hashtag <- PercentageFeatureSet(
  kody.hashtag,
  features = mt_genes,
  assay = "RNA",
  col.name = "percent.mt"
)

# Save QC plot
qc_plot <- VlnPlot(
  kody.hashtag,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(
  filename = file.path(output_dir, "qc_violin_plot.png"),
  plot = qc_plot,
  width = 10,
  height = 4,
  dpi = 300
)

# Filter cells
kody.hashtag <- subset(
  kody.hashtag,
  subset = nFeature_RNA > 350 &
    nFeature_RNA < 4950 &
    percent.mt < 10
)

# -----------------------------
# RNA normalization and scaling
# -----------------------------

kody.hashtag <- NormalizeData(
  kody.hashtag,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

rna_data <- GetAssayData(
  kody.hashtag,
  assay = "RNA",
  slot = "data"
)

num_zeros_normalized <- sum(rna_data == 0)
print(paste("Number of zeros after normalization:", num_zeros_normalized))

kody.hashtag <- FindVariableFeatures(
  kody.hashtag,
  nfeatures = 2000
)

kody.hashtag <- ScaleData(kody.hashtag)

rna_data_scaled <- GetAssayData(
  kody.hashtag,
  assay = "RNA",
  slot = "scale.data"
)

num_zeros_scaled <- sum(rna_data_scaled == 0)
print(paste("Number of zeros after scaling:", num_zeros_scaled))

# -----------------------------
# Subset HTO matrices to filtered cells
# -----------------------------

filtered_cells <- colnames(kody.hashtag)

common_cells_antibody <- intersect(filtered_cells, colnames(kody.htos.antibody))
common_cells_lipid <- intersect(filtered_cells, colnames(kody.htos.lipid))

print(paste("Antibody HTOs match filtered cells:", length(common_cells_antibody) == length(filtered_cells)))
print(paste("Lipid HTOs match filtered cells:", length(common_cells_lipid) == length(filtered_cells)))

kody.htos.antibody_subset <- kody.htos.antibody[, filtered_cells]
kody.htos.lipid_subset <- kody.htos.lipid[, filtered_cells]

# -----------------------------
# Save outputs
# -----------------------------

saveRDS(
  kody.hashtag,
  file = file.path(output_dir, "kody_hashtag_qc.rds")
)

saveRDS(
  kody.htos.antibody_subset,
  file = file.path(output_dir, "hto_antibody_counts_filtered.rds")
)

saveRDS(
  kody.htos.lipid_subset,
  file = file.path(output_dir, "hto_lipid_counts_filtered.rds")
)

write.csv(
  kody.hashtag@meta.data,
  file = file.path(output_dir, "metadata_qc_filtered.csv")
)

print("Saved QC-filtered Seurat object and filtered HTO matrices.")
```
