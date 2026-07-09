# ============================================================
# Enrichment helper functions for leptin cell-culture analyses
# ============================================================

suppressPackageStartupMessages({
  library(enrichR)
  library(dplyr)
  library(writexl)
})

# ------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------

clean_gene_vector <- function(genes) {
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  unique(genes)
}

run_enrichr_safe <- function(genes, dbs, min_genes = 1) {
  genes <- clean_gene_vector(genes)
  if (length(genes) < min_genes) {
    warning("Gene set has fewer than min_genes; returning empty tables.")
    empty <- setNames(vector("list", length(dbs)), dbs)
    empty <- lapply(empty, function(x) data.frame())
    return(empty)
  }
  enrichR::enrichr(genes, dbs)
}

add_label_to_enrichr_result <- function(enrich_result, label_col, label_value) {
  lapply(enrich_result, function(df) {
    if (!is.null(df) && nrow(df) > 0) {
      df[[label_col]] <- label_value
    }
    df
  })
}

combine_one_database <- function(enrich_results,
                                 db,
                                 label_col = NULL,
                                 label_values = NULL,
                                 p_col = "P.value",
                                 p_cutoff = NULL,
                                 padj_col = "Adjusted.P.value",
                                 padj_cutoff = NULL) {
  dfs <- Map(function(x, label_value) {
    df <- x[[db]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!is.null(label_col)) df[[label_col]] <- label_value
    df
  }, enrich_results, label_values)

  combined <- dplyr::bind_rows(dfs)

  if (nrow(combined) == 0) return(combined)

  if (!is.null(p_cutoff) && p_col %in% colnames(combined)) {
    combined <- combined[combined[[p_col]] < p_cutoff, , drop = FALSE]
  }

  if (!is.null(padj_cutoff) && padj_col %in% colnames(combined)) {
    combined <- combined[combined[[padj_col]] < padj_cutoff, , drop = FALSE]
  }

  combined
}

# ------------------------------------------------------------
# A) WT k-means cluster enrichment
# ------------------------------------------------------------

write_cluster_enrichment_csvs <- function(enrich_result,
                                          cluster_file_label,
                                          outdir,
                                          suffix = "wt_time_effect") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  file_map <- c(
    GO_Biological_Process_2025 = "BP",
    MSigDB_Hallmark_2020 = "Hallmark",
    KEGG_2021_Human = "kegg",
    Reactome_Pathways_2024 = "reac"
  )

  for (db in intersect(names(file_map), names(enrich_result))) {
    outfile <- file.path(
      outdir,
      paste0(cluster_file_label, "_", file_map[[db]], "_", suffix, ".csv")
    )
    write.csv(enrich_result[[db]], file = outfile, row.names = FALSE)
  }
}

run_cluster_enrichment <- function(clu_df,
                                   gene_col = "GeneID",
                                   cluster_col = "Cluster",
                                   cluster_order = c("cluster_3", "cluster_4", "cluster_2", "cluster_1", "cluster_5"),
                                   cluster_labels = c("C3", "C4", "C2", "C1", "C5"),
                                   cluster_file_labels = c("cl3", "cl4", "cl2", "cl1", "cl5"),
                                   dbs = c("GO_Biological_Process_2025", "MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_Pathways_2024"),
                                   outdir = "commom_leptin_response/wt_time_effect",
                                   write_csvs = TRUE,
                                   csv_suffix = "wt_time_effect",
                                   write_hallmark_xlsx = TRUE,
                                   hallmark_xlsx = "commom_leptin_response/manuscript/SupplementaryFiles/SF3_hallmark_enrichment_WT_clusters.xlsx") {

  stopifnot(length(cluster_order) == length(cluster_labels))
  stopifnot(length(cluster_order) == length(cluster_file_labels))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  enrich_results <- setNames(vector("list", length(cluster_order)), cluster_file_labels)

  for (i in seq_along(cluster_order)) {
    genes <- clu_df[clu_df[[cluster_col]] == cluster_order[[i]], gene_col]
    enr <- run_enrichr_safe(genes, dbs)
    enr <- add_label_to_enrichr_result(enr, "cluster", cluster_labels[[i]])
    enrich_results[[cluster_file_labels[[i]]]] <- enr

    if (write_csvs) {
      write_cluster_enrichment_csvs(
        enrich_result = enr,
        cluster_file_label = cluster_file_labels[[i]],
        outdir = outdir,
        suffix = csv_suffix
      )
    }
  }

  # Combined tables matching your filtering choices.
  bp <- combine_one_database(
    enrich_results = enrich_results,
    db = "GO_Biological_Process_2025",
    label_col = NULL,
    label_values = cluster_labels,
    padj_cutoff = 0.05,
    p_cutoff = 0.01
  )

  hallmark <- combine_one_database(
    enrich_results = enrich_results,
    db = "MSigDB_Hallmark_2020",
    label_col = NULL,
    label_values = cluster_labels,
    padj_cutoff = 0.05
  )

  kegg <- combine_one_database(
    enrich_results = enrich_results,
    db = "KEGG_2021_Human",
    label_col = NULL,
    label_values = cluster_labels,
    padj_cutoff = 0.05
  )

  reactome <- combine_one_database(
    enrich_results = enrich_results,
    db = "Reactome_Pathways_2024",
    label_col = NULL,
    label_values = cluster_labels,
    padj_cutoff = 0.05
  )

  if (write_hallmark_xlsx && "MSigDB_Hallmark_2020" %in% dbs) {
    dir.create(dirname(hallmark_xlsx), recursive = TRUE, showWarnings = FALSE)

    # Match your supplementary workbook style: cl1, cl2, cl3, cl4, cl5.
    hallmark_list <- list(
      cl1 = enrich_results[["cl1"]][["MSigDB_Hallmark_2020"]],
      cl2 = enrich_results[["cl2"]][["MSigDB_Hallmark_2020"]],
      cl3 = enrich_results[["cl3"]][["MSigDB_Hallmark_2020"]],
      cl4 = enrich_results[["cl4"]][["MSigDB_Hallmark_2020"]],
      cl5 = enrich_results[["cl5"]][["MSigDB_Hallmark_2020"]]
    )
    writexl::write_xlsx(hallmark_list, hallmark_xlsx)
  }

  list(
    raw = enrich_results,
    bp = bp,
    hallmark = hallmark,
    kegg = kegg,
    reactome = reactome
  )
}

# Optional plotting helper for the WT cluster Hallmark dotplot.
plot_cluster_hallmark_dotplot <- function(hallmark,
                                          outfile,
                                          remove_row_indices = NULL,
                                          remove_terms = NULL,
                                          width = 6,
                                          height = 7) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(forcats)
    library(stringr)
    library(grid)
  })

  df <- hallmark %>%
    mutate(
      term = stringr::str_wrap(Term, width = 40),
      overlap_num = as.numeric(sub("/.*", "", Overlap)),
      overlap_den = as.numeric(sub(".*/", "", Overlap)),
      min_p = ave(Adjusted.P.value, term, FUN = function(x) min(x, na.rm = TRUE))
    ) %>%
    mutate(term = forcats::fct_reorder(term, min_p, .desc = FALSE))

  if (!is.null(remove_row_indices)) {
    terms_to_remove <- df$Term[remove_row_indices]
    df <- df[!df$Term %in% terms_to_remove, , drop = FALSE]
  }

  if (!is.null(remove_terms)) {
    df <- df[!df$Term %in% remove_terms, , drop = FALSE]
  }

  pdf(file = outfile, width = width, height = height)
  p <- ggplot(df, aes(x = cluster, y = term)) +
    geom_point(aes(size = overlap_num, color = Adjusted.P.value), alpha = 0.9) +
    scale_size(range = c(2, 12), name = "Overlap (# genes)") +
    scale_color_distiller(
      palette = "Reds",
      name = "Adj. p-value",
      direction = -1,
      limits = c(0, 0.04)
    ) +
    guides(color = guide_colorbar(barheight = 4)) +
    labs(x = "", y = "", title = "Cluster enrichment") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 8), size = 12),
      axis.text.y  = element_text(face = "bold", size = 12, angle = 0),
      legend.key.size    = unit(0.4, "lines"),
      legend.spacing.y   = unit(0.2, "lines"),
      legend.box.spacing = unit(0.2, "lines")
    )
  print(p)
  dev.off()

  invisible(df)
}

# ------------------------------------------------------------
# B) Up/down DEG enrichment for any comparison
# ------------------------------------------------------------

run_up_down_enrichment <- function(up_genes,
                                   down_genes,
                                   dbs,
                                   outdir,
                                   comparison_label,
                                   up_file = NULL,
                                   down_file = NULL,
                                   p_cutoff = 0.05,
                                   write_raw_xlsx = TRUE) {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  up_enrich <- run_enrichr_safe(up_genes, dbs)
  down_enrich <- run_enrichr_safe(down_genes, dbs)

  up_enrich <- add_label_to_enrichr_result(up_enrich, "regulation", "UP")
  down_enrich <- add_label_to_enrichr_result(down_enrich, "regulation", "Down")

  if (write_raw_xlsx) {
    if (is.null(up_file)) {
      up_file <- file.path(outdir, paste0("UP_enrichment_", comparison_label, "_comp.xlsx"))
    }
    if (is.null(down_file)) {
      down_file <- file.path(outdir, paste0("Down_enrichment_", comparison_label, "_comp.xlsx"))
    }
    writexl::write_xlsx(up_enrich, up_file)
    writexl::write_xlsx(down_enrich, down_file)
  }

  combined_by_db <- setNames(vector("list", length(dbs)), dbs)
  for (db in dbs) {
    combined <- dplyr::bind_rows(up_enrich[[db]], down_enrich[[db]])
    if (nrow(combined) > 0 && "P.value" %in% colnames(combined)) {
      combined <- combined[combined$P.value < p_cutoff, , drop = FALSE]
    }
    combined$comparison <- comparison_label
    combined_by_db[[db]] <- combined
  }

  list(
    up = up_enrich,
    down = down_enrich,
    combined_by_db = combined_by_db
  )
}

run_up_down_enrichment_set <- function(gene_sets,
                                       dbs = c("GO_Biological_Process_2025"),
                                       outdir,
                                       p_cutoff = 0.05,
                                       write_raw_xlsx = TRUE) {

  # gene_sets must be a named list, for example:
  # list(
  #   `0min` = list(up = Up_Stat3_vsWT_0$id, down = Down_Stat3_vsWT_0$id),
  #   `30min` = list(up = Up_Stat3_vsWT_30$id, down = Down_Stat3_vsWT_30$id),
  #   `180min` = list(up = Up_Stat3_vsWT_180$id, down = Down_Stat3_vsWT_180$id)
  # )

  results <- setNames(vector("list", length(gene_sets)), names(gene_sets))

  for (comparison_label in names(gene_sets)) {
    results[[comparison_label]] <- run_up_down_enrichment(
      up_genes = gene_sets[[comparison_label]]$up,
      down_genes = gene_sets[[comparison_label]]$down,
      dbs = dbs,
      outdir = outdir,
      comparison_label = comparison_label,
      p_cutoff = p_cutoff,
      write_raw_xlsx = write_raw_xlsx
    )
  }

  combined_by_db <- setNames(vector("list", length(dbs)), dbs)
  for (db in dbs) {
    combined_by_db[[db]] <- dplyr::bind_rows(
      lapply(results, function(x) x$combined_by_db[[db]])
    )
  }

  list(
    raw = results,
    combined_by_db = combined_by_db
  )
}
