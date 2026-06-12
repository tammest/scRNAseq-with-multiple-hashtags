# ============================================================
# Add antibody and lipid hashtag assays to the QC-filtered
# Seurat object, normalize each HTO assay using CLR, and run
# hashtag demultiplexing separately for antibody and lipid tags.
# ============================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsci)

# -----------------------------
# Paths
# -----------------------------

input_dir <- "outputs"
output_dir <- "outputs"

seurat_rds <- file.path(input_dir, "kody_hashtag_qc.rds")
hto_antibody_rds <- file.path(input_dir, "hto_antibody_counts_filtered.rds")
hto_lipid_rds <- file.path(input_dir, "hto_lipid_counts_filtered.rds")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load QC-filtered object and HTO matrices
# -----------------------------

kody.hashtag <- readRDS(seurat_rds)
kody.htos.antibody_subset <- readRDS(hto_antibody_rds)
kody.htos.lipid_subset <- readRDS(hto_lipid_rds)

# Check that HTO matrices match the filtered Seurat object
filtered_cells <- colnames(kody.hashtag)

if (!all(filtered_cells %in% colnames(kody.htos.antibody_subset))) {
  stop("Not all filtered cells are present in the antibody HTO matrix.")
}

if (!all(filtered_cells %in% colnames(kody.htos.lipid_subset))) {
  stop("Not all filtered cells are present in the lipid HTO matrix.")
}

# Reorder HTO matrices to match Seurat object cell order
kody.htos.antibody_subset <- kody.htos.antibody_subset[, filtered_cells]
kody.htos.lipid_subset <- kody.htos.lipid_subset[, filtered_cells]

print("Antibody HTO matrix dimensions:")
print(dim(kody.htos.antibody_subset))

print("Lipid HTO matrix dimensions:")
print(dim(kody.htos.lipid_subset))

print("Antibody HTO names:")
print(rownames(kody.htos.antibody_subset))

print("Lipid HTO names:")
print(rownames(kody.htos.lipid_subset))

# -----------------------------
# Add HTO assays
# -----------------------------

kody.hashtag[["HTOantibody"]] <- CreateAssayObject(
  counts = kody.htos.antibody_subset
)

kody.hashtag[["HTOlipid"]] <- CreateAssayObject(
  counts = kody.htos.lipid_subset
)

# -----------------------------
# Normalize HTO assays
# -----------------------------

kody.hashtag <- NormalizeData(
  kody.hashtag,
  assay = "HTOantibody",
  normalization.method = "CLR"
)

kody.hashtag <- NormalizeData(
  kody.hashtag,
  assay = "HTOlipid",
  normalization.method = "CLR"
)

# -----------------------------
# Demultiplex antibody hashtags
# -----------------------------

kody.hashtag <- HTODemux(
  kody.hashtag,
  assay = "HTOantibody",
  positive.quantile = 0.99
)

# Save antibody classification before running lipid demultiplexing,
# because HTODemux writes to hash.ID and classification columns.
hash_id_antibody <- kody.hashtag@meta.data$hash.ID
hto_antibody_classification <- kody.hashtag@meta.data$HTOantibody_classification
hto_antibody_classification_global <- kody.hashtag@meta.data$HTOantibody_classification.global

# -----------------------------
# Demultiplex lipid hashtags
# -----------------------------

kody.hashtag <- HTODemux(
  kody.hashtag,
  assay = "HTOlipid",
  positive.quantile = 0.99
)

# Store lipid classification after lipid demultiplexing
hash_id_lipid <- kody.hashtag@meta.data$hash.ID
hto_lipid_classification <- kody.hashtag@meta.data$HTOlipid_classification
hto_lipid_classification_global <- kody.hashtag@meta.data$HTOlipid_classification.global

# Add saved antibody and lipid labels to metadata with clear names
kody.hashtag@meta.data$hash.ID.antibody <- hash_id_antibody
kody.hashtag@meta.data$hash.ID.lipid <- hash_id_lipid

kody.hashtag@meta.data$HTOantibody_classification_saved <- hto_antibody_classification
kody.hashtag@meta.data$HTOantibody_classification_global_saved <- hto_antibody_classification_global

kody.hashtag@meta.data$HTOlipid_classification_saved <- hto_lipid_classification
kody.hashtag@meta.data$HTOlipid_classification_global_saved <- hto_lipid_classification_global

# -----------------------------
# Summary tables
# -----------------------------

antibody_table <- table(kody.hashtag@meta.data$HTOantibody_classification_saved)
lipid_table <- table(kody.hashtag@meta.data$HTOlipid_classification_saved)

print("Antibody HTO classification counts:")
print(antibody_table)

print("Lipid HTO classification counts:")
print(lipid_table)

write.csv(
  as.data.frame(antibody_table),
  file = file.path(output_dir, "antibody_hto_classification_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(lipid_table),
  file = file.path(output_dir, "lipid_hto_classification_counts.csv"),
  row.names = FALSE
)

# -----------------------------
# UMAP/DimPlot placeholders
# -----------------------------
# These plots will only work after UMAP has been run.
# They are kept here as optional checks if the object already has embeddings.

if ("umap" %in% names(kody.hashtag@reductions)) {
  p_antibody <- DimPlot(
    kody.hashtag,
    label = FALSE,
    group.by = "HTOantibody_classification_saved"
  ) + scale_color_igv()

  p_lipid <- DimPlot(
    kody.hashtag,
    label = FALSE,
    group.by = "HTOlipid_classification_saved"
  ) + scale_color_igv()

  ggsave(
    filename = file.path(output_dir, "umap_antibody_hto_classification.png"),
    plot = p_antibody,
    width = 7,
    height = 6,
    dpi = 300
  )

  ggsave(
    filename = file.path(output_dir, "umap_lipid_hto_classification.png"),
    plot = p_lipid,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# -----------------------------
# Save demultiplexed object
# -----------------------------

saveRDS(
  kody.hashtag,
  file = file.path(output_dir, "kody_hashtag_demuxed.rds")
)

write.csv(
  kody.hashtag@meta.data,
  file = file.path(output_dir, "metadata_demuxed.csv")
)

print("Saved demultiplexed Seurat object and metadata.")
