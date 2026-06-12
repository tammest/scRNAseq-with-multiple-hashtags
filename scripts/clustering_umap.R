# ============================================================
# Run PCA, select PCs, perform clustering, and generate UMAP
# embeddings for the full Seurat object and for phenotype-only
# cells.
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
figure_dir <- file.path(output_dir, "figures")

input_rds <- file.path(input_dir, "kody_hashtag_phenotypes.rds")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper function: select PCs
# -----------------------------

select_pcs <- function(seurat_obj) {
  stdv <- seurat_obj[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)

  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(
    which((percent_stdv[1:(length(percent_stdv) - 1)] -
             percent_stdv[2:length(percent_stdv)]) > 0.1),
    decreasing = TRUE
  )[1] + 1

  min_pc <- min(co1, co2, na.rm = TRUE)

  if (is.infinite(min_pc) || is.na(min_pc)) {
    min_pc <- 14
  }

  return(min_pc)
}

# -----------------------------
# Load object
# -----------------------------

kody.hashtag <- readRDS(input_rds)

# -----------------------------
# PCA, clustering, and UMAP on full object
# -----------------------------

kody.hashtag <- RunPCA(
  kody.hashtag,
  assay = "RNA",
  features = VariableFeatures(kody.hashtag),
  npcs = 50
)

min_pc <- select_pcs(kody.hashtag)

print(paste("Number of PCs selected for full object:", min_pc))

set.seed(1007)

kody.hashtag <- FindNeighbors(
  object = kody.hashtag,
  dims = 1:min_pc
)

kody.hashtag <- FindClusters(
  kody.hashtag,
  resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1, 2, 3, 4)
)

kody.hashtag <- RunUMAP(
  kody.hashtag,
  dims = 1:min_pc,
  min.dist = 0.5
)

# Save full object plots

p_full_clusters <- DimPlot(
  kody.hashtag,
  group.by = "RNA_snn_res.0.4",
  label = TRUE,
  label.box = TRUE
) +
  scale_color_igv()

ggsave(
  filename = file.path(figure_dir, "full_object_umap_clusters_res_0_4.png"),
  plot = p_full_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

p_full_phenotype <- DimPlot(
  kody.hashtag,
  group.by = "phenotype",
  label = FALSE
) +
  scale_color_manual(
    values = c(
      "Phenotype 1" = "blue",
      "Phenotype 2" = "gray",
      "Phenotype 3" = "red",
      "Other" = "lightgray"
    )
  )

ggsave(
  filename = file.path(figure_dir, "full_object_umap_by_phenotype.png"),
  plot = p_full_phenotype,
  width = 7,
  height = 6,
  dpi = 300
)

# -----------------------------
# Subset to phenotype-assigned cells
# -----------------------------

justphenotypes <- subset(
  kody.hashtag,
  subset = phenotype %in% c("Phenotype 1", "Phenotype 2", "Phenotype 3")
)

# Re-run RNA preprocessing on phenotype-only object

justphenotypes <- NormalizeData(
  justphenotypes,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

justphenotypes <- FindVariableFeatures(
  justphenotypes,
  nfeatures = 2000
)

justphenotypes <- ScaleData(justphenotypes)

# PCA, clustering, and UMAP on phenotype-only object

justphenotypes <- RunPCA(
  justphenotypes,
  assay = "RNA",
  features = VariableFeatures(justphenotypes),
  npcs = 50
)

min_pc_phenotype <- select_pcs(justphenotypes)

print(paste("Number of PCs selected for phenotype-only object:", min_pc_phenotype))

set.seed(1007)

justphenotypes <- FindNeighbors(
  object = justphenotypes,
  dims = 1:min_pc_phenotype
)

justphenotypes <- FindClusters(
  justphenotypes,
  resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1, 2, 3, 4)
)

justphenotypes <- RunUMAP(
  justphenotypes,
  dims = 1:min_pc_phenotype,
  min.dist = 0.5
)

Idents(justphenotypes) <- "RNA_snn_res.0.4"

# -----------------------------
# Save phenotype-only plots
# -----------------------------

p_pheno <- DimPlot(
  justphenotypes,
  label = FALSE,
  group.by = "phenotype"
) +
  scale_color_manual(
    values = c(
      "Phenotype 1" = "blue",
      "Phenotype 2" = "gray",
      "Phenotype 3" = "red"
    )
  )

ggsave(
  filename = file.path(figure_dir, "phenotype_only_umap_by_phenotype.png"),
  plot = p_pheno,
  width = 7,
  height = 6,
  dpi = 300
)

p_pheno_clusters <- DimPlot(
  justphenotypes,
  label = TRUE,
  label.box = TRUE,
  group.by = "RNA_snn_res.0.4"
) +
  scale_color_igv()

ggsave(
  filename = file.path(figure_dir, "phenotype_only_umap_clusters_res_0_4.png"),
  plot = p_pheno_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

# Optional hashtag classification plots

if ("HTOantibody_classification_saved" %in% colnames(justphenotypes@meta.data)) {
  p_antibody <- DimPlot(
    justphenotypes,
    label = FALSE,
    group.by = "HTOantibody_classification_saved"
  ) +
    scale_color_igv()

  ggsave(
    filename = file.path(figure_dir, "phenotype_only_umap_antibody_hashtag.png"),
    plot = p_antibody,
    width = 7,
    height = 6,
    dpi = 300
  )
}

if ("HTOlipid_classification_saved" %in% colnames(justphenotypes@meta.data)) {
  p_lipid <- DimPlot(
    justphenotypes,
    label = FALSE,
    group.by = "HTOlipid_classification_saved"
  ) +
    scale_color_igv()

  ggsave(
    filename = file.path(figure_dir, "phenotype_only_umap_lipid_hashtag.png"),
    plot = p_lipid,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# -----------------------------
# Save objects
# -----------------------------

saveRDS(
  kody.hashtag,
  file = file.path(output_dir, "kody_hashtag_clustered.rds")
)

saveRDS(
  justphenotypes,
  file = file.path(output_dir, "justphenotypes_clustered.rds")
)

write.csv(
  kody.hashtag@meta.data,
  file = file.path(output_dir, "metadata_full_clustered.csv")
)

write.csv(
  justphenotypes@meta.data,
  file = file.path(output_dir, "metadata_phenotype_only_clustered.csv")
)

print("Saved clustered full object and phenotype-only object.")
