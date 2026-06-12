# ============================================================
# GO enrichment analysis for phenotype differential expression
# results generated in 05_marker_analysis.R.
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(GOplot)

# -----------------------------
# Paths
# -----------------------------

input_dir <- "outputs/markers"
output_dir <- "outputs/go_enrichment"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper function: run GO enrichment
# -----------------------------

run_go_enrichment <- function(marker_file, comparison_name, ont = "BP") {
  markers <- read.csv(marker_file)

  if (!"gene" %in% colnames(markers)) {
    stop(paste("Missing gene column in:", marker_file))
  }

  if (!"p_val_adj" %in% colnames(markers)) {
    stop(paste("Missing p_val_adj column in:", marker_file))
  }

  significant_genes <- markers %>%
    filter(p_val_adj < 0.05) %>%
    pull(gene) %>%
    unique()

  if (length(significant_genes) == 0) {
    warning(paste("No significant genes found for:", comparison_name))
    return(NULL)
  }

  go_results <- enrichGO(
    gene = significant_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  go_df <- as.data.frame(go_results)

  write.csv(
    go_df,
    file = file.path(output_dir, paste0(comparison_name, "_GO_", ont, ".csv")),
    row.names = FALSE
  )

  return(go_results)
}

# -----------------------------
# Helper function: prepare top GO gene table
# -----------------------------

prepare_top_go_gene_table <- function(marker_file, go_results, comparison_name, top_n_terms = 15, top_n_genes = 10) {
  markers <- read.csv(marker_file)

  if (is.null(go_results)) {
    return(NULL)
  }

  go_df <- as.data.frame(go_results)

  if (nrow(go_df) == 0) {
    warning(paste("No GO terms found for:", comparison_name))
    return(NULL)
  }

  top_go_terms <- go_df %>%
    arrange(p.adjust) %>%
    head(top_n_terms)

  gene_list <- unlist(strsplit(top_go_terms$geneID, "/"))

  gene_fc <- markers %>%
    filter(gene %in% gene_list) %>%
    select(gene, avg_log2FC)

  gene_go_mapping <- data.frame(
    gene = gene_list,
    GO_term = rep(
      top_go_terms$Description,
      sapply(strsplit(top_go_terms$geneID, "/"), length)
    )
  )

  top_genes <- merge(gene_fc, gene_go_mapping, by = "gene")

  top_genes <- top_genes %>%
    arrange(desc(avg_log2FC)) %>%
    head(top_n_genes)

  write.csv(
    top_go_terms,
    file = file.path(output_dir, paste0(comparison_name, "_top_GO_terms.csv")),
    row.names = FALSE
  )

  write.csv(
    top_genes,
    file = file.path(output_dir, paste0(comparison_name, "_top_GO_genes.csv")),
    row.names = FALSE
  )

  return(list(
    top_go_terms = top_go_terms,
    top_genes = top_genes
  ))
}

# -----------------------------
# Run GO enrichment on phenotype comparisons
# -----------------------------

marker_files <- list.files(
  input_dir,
  pattern = "Phenotype_.*_vs_Phenotype_.*\\.csv$",
  full.names = TRUE
)

if (length(marker_files) == 0) {
  stop("No phenotype comparison marker files found in outputs/markers.")
}

go_results_list <- list()
top_go_tables <- list()

for (marker_file in marker_files) {
  comparison_name <- tools::file_path_sans_ext(basename(marker_file))

  print(paste("Running GO enrichment for:", comparison_name))

  go_results <- run_go_enrichment(
    marker_file = marker_file,
    comparison_name = comparison_name,
    ont = "BP"
  )

  go_results_list[[comparison_name]] <- go_results

  top_go_tables[[comparison_name]] <- prepare_top_go_gene_table(
    marker_file = marker_file,
    go_results = go_results,
    comparison_name = comparison_name,
    top_n_terms = 15,
    top_n_genes = 10
  )
}

# -----------------------------
# Save all GO results to workbook
# -----------------------------

wb <- createWorkbook()

for (comparison_name in names(go_results_list)) {
  go_obj <- go_results_list[[comparison_name]]

  if (is.null(go_obj)) {
    next
  }

  go_df <- as.data.frame(go_obj)

  if (nrow(go_df) == 0) {
    next
  }

  sheet_name <- substr(comparison_name, 1, 31)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, go_df)
}

saveWorkbook(
  wb,
  file = file.path(output_dir, "GO_enrichment_results.xlsx"),
  overwrite = TRUE
)

print("GO enrichment analysis completed.")
