# ============================================================
# Marker visualization and differential expression analysis
# for phenotype-assigned and reclustered Seurat objects.
# ============================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(ggsci)
library(Nebulosa)

# -----------------------------
# Paths
# -----------------------------

input_dir <- "outputs"
output_dir <- "outputs"
figure_dir <- file.path(output_dir, "figures")
marker_dir <- file.path(output_dir, "markers")

input_rds <- file.path(input_dir, "justphenotypes_clustered.rds")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load object
# -----------------------------

justphenotypes <- readRDS(input_rds)

cluster_col <- "RNA_snn_res.0.4"

if (!cluster_col %in% colnames(justphenotypes@meta.data)) {
  stop(paste("Missing clustering column:", cluster_col))
}

Idents(justphenotypes) <- cluster_col

# -----------------------------
# Marker gene visualization
# -----------------------------

genes_to_plot <- c(
  "Itga5", "Itga6", "Krt10", "Krt14", "Pdgfra", "Pecam1",
  "Acta2", "Lyve1", "Pdpn", "Ptprc", "Cd3e", "Trac",
  "Trbc1", "Trdc", "Cd19", "Ly6g", "Fcer1a", "Siglecf",
  "Itgam", "Itgax", "Ly6c1", "Adgre1", "Cd68", "Cd163",
  "Mrc1", "H2-Ab1", "H2-Eb1", "Cd80", "Cd86", "Csf1r",
  "Cdh5", "Vwf"
)

genes_present <- genes_to_plot[genes_to_plot %in% rownames(justphenotypes)]

if (length(genes_present) > 0) {
  dot_plot <- DotPlot(
    justphenotypes,
    features = genes_present
  ) +
    coord_flip() +
    scale_color_gradientn(colors = c("grey85", brewer.pal(7, "Reds"))) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0))

  ggsave(
    filename = file.path(figure_dir, "marker_genes_dotplot.png"),
    plot = dot_plot,
    width = 10,
    height = 8,
    dpi = 300
  )

  violin_plot <- VlnPlot(
    justphenotypes,
    features = genes_present,
    group.by = cluster_col,
    stack = TRUE,
    flip = TRUE
  )

  ggsave(
    filename = file.path(figure_dir, "marker_genes_violinplot.png"),
    plot = violin_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# -----------------------------
# Optional endothelial marker density plot
# -----------------------------

endo_genes <- c(
  "B2m", "Epas1", "Hspa1a", "Tmsb4x", "Vwf", "Rgcc",
  "Igfbp7", "Sparc", "Vim", "Timp3", "Ifitm3", "Tmsb10",
  "Zfp36", "Socs3", "Pecam1", "A2m", "Sparcl1", "Itm2b",
  "Junb", "Myh9", "Ifi27", "Aqp1", "Col4a2", "Tm4sf1",
  "Hspa1b", "Gnai2", "Col4a1", "Pcdh17", "Mgp", "Heg1"
)

endo_genes_present <- endo_genes[endo_genes %in% rownames(justphenotypes)]

if (length(endo_genes_present) > 0) {
  density_plot <- plot_density(
    justphenotypes,
    features = endo_genes_present
  )

  ggsave(
    filename = file.path(figure_dir, "endothelial_marker_densityplot.png"),
    plot = density_plot,
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
}

# -----------------------------
# Find markers across clusters
# -----------------------------

all_genes_justphenotypes <- FindAllMarkers(
  justphenotypes,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  all_genes_justphenotypes,
  file = file.path(marker_dir, "all_genes_justphenotypes_by_cluster.csv"),
  row.names = FALSE
)

# -----------------------------
# Find markers within each phenotype
# -----------------------------

get_all_genes_by_phenotype <- function(seurat_obj, cluster_column) {
  phenotypes <- unique(seurat_obj$phenotype)
  phenotypes <- phenotypes[!is.na(phenotypes)]

  results <- list()

  for (pheno in phenotypes) {
    phenotype_obj <- subset(seurat_obj, subset = phenotype == pheno)

    if (ncol(phenotype_obj) < 10) {
      warning(paste("Skipping", pheno, "- too few cells."))
      next
    }

    Idents(phenotype_obj) <- cluster_column

    markers <- FindAllMarkers(
      phenotype_obj,
      only.pos = FALSE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )

    results[[pheno]] <- markers
  }

  return(results)
}

all_genes_by_phenotype <- get_all_genes_by_phenotype(
  justphenotypes,
  cluster_col
)

for (pheno in names(all_genes_by_phenotype)) {
  safe_name <- gsub(" ", "_", pheno)

  write.csv(
    all_genes_by_phenotype[[pheno]],
    file = file.path(marker_dir, paste0(safe_name, "_markers_by_cluster.csv")),
    row.names = FALSE
  )
}

# -----------------------------
# Pairwise phenotype comparisons
# -----------------------------

compare_genes_between_phenotypes <- function(seurat_obj) {
  Idents(seurat_obj) <- "phenotype"

  phenotypes <- unique(seurat_obj$phenotype)
  phenotypes <- phenotypes[!is.na(phenotypes)]

  results <- list()

  for (i in seq_len(length(phenotypes) - 1)) {
    for (j in (i + 1):length(phenotypes)) {
      phenotype_1 <- phenotypes[i]
      phenotype_2 <- phenotypes[j]

      comparison_name <- paste(
        gsub(" ", "_", phenotype_1),
        "vs",
        gsub(" ", "_", phenotype_2),
        sep = "_"
      )

      markers <- FindMarkers(
        seurat_obj,
        ident.1 = phenotype_1,
        ident.2 = phenotype_2,
        only.pos = FALSE,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )

      markers$gene <- rownames(markers)
      results[[comparison_name]] <- markers
    }
  }

  return(results)
}

all_genes_between_phenotypes <- compare_genes_between_phenotypes(justphenotypes)

# Save pairwise phenotype comparisons to CSV and Excel workbook

wb <- createWorkbook()

for (comparison_name in names(all_genes_between_phenotypes)) {
  comparison_result <- all_genes_between_phenotypes[[comparison_name]]

  write.csv(
    comparison_result,
    file = file.path(marker_dir, paste0(comparison_name, ".csv")),
    row.names = FALSE
  )

  sheet_name <- substr(comparison_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, comparison_result)
}

saveWorkbook(
  wb,
  file = file.path(marker_dir, "all_genes_between_phenotypes.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Compare specific phenotype-cluster pairs
# -----------------------------

compare_phenotype_clusters <- function(
    seurat_obj,
    phenotype_1,
    cluster_1,
    phenotype_2,
    cluster_2,
    cluster_column
) {
  subset_obj <- subset(
    seurat_obj,
    subset =
      (phenotype == phenotype_1 & .data[[cluster_column]] == cluster_1) |
      (phenotype == phenotype_2 & .data[[cluster_column]] == cluster_2)
  )

  subset_obj$comparison_group <- case_when(
    subset_obj$phenotype == phenotype_1 &
      subset_obj@meta.data[[cluster_column]] == cluster_1 ~ paste0(
        gsub(" ", "_", phenotype_1), "_cluster_", cluster_1
      ),

    subset_obj$phenotype == phenotype_2 &
      subset_obj@meta.data[[cluster_column]] == cluster_2 ~ paste0(
        gsub(" ", "_", phenotype_2), "_cluster_", cluster_2
      ),

    TRUE ~ NA_character_
  )

  subset_obj <- subset(subset_obj, subset = !is.na(comparison_group))

  Idents(subset_obj) <- "comparison_group"

  ident_1 <- paste0(gsub(" ", "_", phenotype_1), "_cluster_", cluster_1)
  ident_2 <- paste0(gsub(" ", "_", phenotype_2), "_cluster_", cluster_2)

  markers <- FindMarkers(
    subset_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  markers$gene <- rownames(markers)

  return(markers)
}

# Example comparison from your original code:
# Phenotype 1 cluster 4 vs Phenotype 2 cluster 4

specific_comparison <- compare_phenotype_clusters(
  seurat_obj = justphenotypes,
  phenotype_1 = "Phenotype 1",
  cluster_1 = "4",
  phenotype_2 = "Phenotype 2",
  cluster_2 = "4",
  cluster_column = cluster_col
)

write.csv(
  specific_comparison,
  file = file.path(marker_dir, "Phenotype_1_cluster_4_vs_Phenotype_2_cluster_4.csv"),
  row.names = FALSE
)

wb_specific <- createWorkbook()
addWorksheet(wb_specific, "P1_c4_vs_P2_c4")
writeData(wb_specific, "P1_c4_vs_P2_c4", specific_comparison)

saveWorkbook(
  wb_specific,
  file = file.path(marker_dir, "phenotype_cluster_comparison.xlsx"),
  overwrite = TRUE
)

print("Marker analysis completed.")
