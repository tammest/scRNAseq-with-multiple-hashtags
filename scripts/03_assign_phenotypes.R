# ============================================================
# 03_assign_phenotypes.R
#
# Assign phenotype labels using combined lipid and antibody
# hashtag classifications.
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

input_rds <- file.path(input_dir, "kody_hashtag_demuxed.rds")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load demultiplexed object
# -----------------------------

kody.hashtag <- readRDS(input_rds)

# -----------------------------
# Define classification columns
# -----------------------------
# These columns were created in 02_demultiplex_lipid_antibody_hashtags.R

lipid_col <- "HTOlipid_classification_saved"
antibody_col <- "HTOantibody_classification_saved"

if (!lipid_col %in% colnames(kody.hashtag@meta.data)) {
  stop(paste("Missing metadata column:", lipid_col))
}

if (!antibody_col %in% colnames(kody.hashtag@meta.data)) {
  stop(paste("Missing metadata column:", antibody_col))
}

# -----------------------------
# Assign phenotypes
# -----------------------------
# Phenotype definitions:
#
# Phenotype 1:
#   Lipid hashtags: CMO301, CMO302
#   Antibody hashtag: TotalseqB301
#
# Phenotype 2:
#   Lipid hashtags: CMO303, CMO304
#   Antibody hashtag: TotalseqB302
#
# Phenotype 3:
#   Lipid hashtags: CMO305, CMO306
#   Antibody hashtag: TotalseqB303
#
# Negative lipid calls are allowed if the antibody hashtag is
# clearly assigned to the corresponding phenotype.

meta <- kody.hashtag@meta.data

meta$phenotype <- case_when(
  meta[[lipid_col]] %in% c("CMO301", "CMO302", "CMO301_CMO302", "Negative") &
    meta[[antibody_col]] %in% c("TotalseqB301", "Negative") ~ "Phenotype 1",

  meta[[lipid_col]] %in% c("CMO303", "CMO304", "CMO303_CMO304", "Negative") &
    meta[[antibody_col]] %in% c("TotalseqB302", "Negative") ~ "Phenotype 2",

  meta[[lipid_col]] %in% c("CMO305", "CMO306", "CMO305_CMO306", "Negative") &
    meta[[antibody_col]] %in% c("TotalseqB303", "Negative") ~ "Phenotype 3",

  TRUE ~ "Other"
)

kody.hashtag@meta.data <- meta

# -----------------------------
# Summary tables
# -----------------------------

phenotype_counts <- table(kody.hashtag@meta.data$phenotype)
print("Phenotype counts:")
print(phenotype_counts)

phenotype_by_lipid <- table(
  kody.hashtag@meta.data$phenotype,
  kody.hashtag@meta.data[[lipid_col]]
)

phenotype_by_antibody <- table(
  kody.hashtag@meta.data$phenotype,
  kody.hashtag@meta.data[[antibody_col]]
)

write.csv(
  as.data.frame(phenotype_counts),
  file = file.path(output_dir, "phenotype_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(phenotype_by_lipid),
  file = file.path(output_dir, "phenotype_by_lipid_hashtag.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(phenotype_by_antibody),
  file = file.path(output_dir, "phenotype_by_antibody_hashtag.csv"),
  row.names = FALSE
)

# -----------------------------
# Optional UMAP plots
# -----------------------------
# These will only be generated if UMAP has already been run.

if ("umap" %in% names(kody.hashtag@reductions)) {
  p1 <- DimPlot(
    kody.hashtag,
    label = FALSE,
    group.by = "phenotype"
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
    filename = file.path(output_dir, "umap_by_phenotype.png"),
    plot = p1,
    width = 7,
    height = 6,
    dpi = 300
  )

  p2 <- DimPlot(
    kody.hashtag,
    label = FALSE,
    group.by = lipid_col
  ) +
    scale_color_igv()

  ggsave(
    filename = file.path(output_dir, "umap_by_lipid_hashtag.png"),
    plot = p2,
    width = 7,
    height = 6,
    dpi = 300
  )

  p3 <- DimPlot(
    kody.hashtag,
    label = FALSE,
    group.by = antibody_col
  ) +
    scale_color_igv()

  ggsave(
    filename = file.path(output_dir, "umap_by_antibody_hashtag.png"),
    plot = p3,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# -----------------------------
# Save object
# -----------------------------

saveRDS(
  kody.hashtag,
  file = file.path(output_dir, "kody_hashtag_phenotypes.rds")
)

write.csv(
  kody.hashtag@meta.data,
  file = file.path(output_dir, "metadata_with_phenotypes.csv")
)

print("Saved Seurat object with phenotype labels.")
