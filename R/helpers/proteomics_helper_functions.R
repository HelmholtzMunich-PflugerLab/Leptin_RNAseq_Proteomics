############################################################
# Script: proteomics_helper_functions.R
# Project: Temporal organization of leptin signaling reveals distinct
#          regulatory roles for STAT3 and ETV6
#
# Purpose:
#   Provide helper functions for the proteomics analysis workflow.
#
# Description:
#   This script contains reusable functions for data preparation,
#   PCA plotting, limma result extraction, heatmap generation,
#   and proteome-transcriptome comparison plots used by
#   proteomics_main_clean.R.
#
# Usage:
#   Source this file from 1_proteomics_main.R before running
#   the main analysis sections.
############################################################
#--------------------------------------------------------------
# Create sample metadata from the metadata sheet.
#
# Expected columns:
#   SampleID  : sample identifier matching abundance columns
#   Condition : sample condition, e.g. WT_minusLep, S15_plusLep
#
# Output:
#   data frame with genotype and leptin factors added.
#--------------------------------------------------------------
make_sample_info <- function(metadata) {
  sample_info <- metadata %>%
    dplyr::mutate(
      genotype = dplyr::case_when(
        stringr::str_detect(Condition, "^WT_")  ~ "WT",
        stringr::str_detect(Condition, "^S15_") ~ "Stat3Δ",
        TRUE ~ NA_character_
      ),
      leptin = dplyr::case_when(
        stringr::str_detect(Condition, "minusLep$") ~ "0min",
        stringr::str_detect(Condition, "plusLep$")  ~ "180min",
        TRUE ~ NA_character_
      )
    )

  if (any(is.na(sample_info$genotype)) || any(is.na(sample_info$leptin))) {
    stop("Could not parse genotype/leptin for some samples. Check Condition names.")
  }

  sample_info$genotype <- factor(sample_info$genotype, levels = c("WT", "Stat3Δ"))
  sample_info$leptin   <- factor(sample_info$leptin,   levels = c("0min", "180min"))

  sample_info
}

#--------------------------------------------------------------
# Prepare log2 protein abundance matrix.
#
# The original analysis used log2(normalized abundance + 1), followed by
# filtering for proteins quantified in at least 50% of samples and in at
# least two samples per genotype/leptin group.
#--------------------------------------------------------------
make_log2_protein_matrix <- function(protein_tbl, sample_info,
                                     accession_col = "Accession",
                                     remove_col_index = 2,
                                     offset = 1,
                                     min_prop = 0.5,
                                     min_per_group = 2) {
  protein_tbl <- protein_tbl[, -remove_col_index]

  E <- protein_tbl %>%
    tibble::column_to_rownames(accession_col) %>%
    dplyr::select(dplyr::all_of(sample_info$SampleID)) %>%
    as.matrix()

  stopifnot(identical(colnames(E), sample_info$SampleID))

  E <- log2(E + offset)

  keep <- rowMeans(!is.na(E)) >= min_prop
  E <- E[keep, , drop = FALSE]

  grp <- interaction(sample_info$genotype, sample_info$leptin, drop = TRUE)
  keep2 <- apply(E, 1, function(x) {
    all(tapply(!is.na(x), grp, sum) >= min_per_group)
  })

  E[keep2, , drop = FALSE]
}

#--------------------------------------------------------------
# Make PCA plot from a protein abundance matrix.
#--------------------------------------------------------------
plot_proteomics_pca <- function(E, sample_info, title,
                                color_var = "genotype",
                                shape_var = "leptin",
                                point_size = 3,
                                color_values = NULL,
                                shape_values = c("0min" = 16, "180min" = 17),
                                fixed_point_color = NULL) {
  X_pca <- t(E)
  pca_res <- prcomp(X_pca, center = TRUE, scale. = TRUE)

  pca_df <- data.frame(
    Sample   = rownames(X_pca),
    PC1      = pca_res$x[, 1],
    PC2      = pca_res$x[, 2],
    genotype = sample_info$genotype,
    leptin   = sample_info$leptin
  )

  if (is.null(fixed_point_color)) {
    p <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]], shape = .data[[shape_var]])
    ) +
      ggplot2::geom_point(size = point_size)
  } else {
    p <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, shape = .data[[shape_var]])
    ) +
      ggplot2::geom_point(size = point_size, color = fixed_point_color)
  }

  p <- p +
    ggplot2::scale_shape_manual(values = shape_values) +
    ggplot2::labs(
      title = title,
      x = paste0("PC1: ", round(100 * summary(pca_res)$importance[2, 1], 1), "% variance"),
      y = paste0("PC2: ", round(100 * summary(pca_res)$importance[2, 2], 1), "% variance")
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()

  if (!is.null(color_values) && is.null(fixed_point_color)) {
    p <- p + ggplot2::scale_color_manual(values = color_values)
  }

  p
}

#--------------------------------------------------------------
# Fit limma model for genotype, leptin, and genotype-by-leptin effects.
#
# Model:
#   abundance ~ genotype * leptin
#
# Contrasts reproduce the original proteomics analysis.
#--------------------------------------------------------------
run_proteomics_limma <- function(E, sample_info) {
  stopifnot(identical(colnames(E), sample_info$SampleID))

  sample_info$genotype <- factor(sample_info$genotype, levels = c("WT", "Stat3Δ"))
  sample_info$leptin   <- factor(sample_info$leptin,   levels = c("0min", "180min"))

  design <- stats::model.matrix(~ genotype * leptin, data = sample_info)
  colnames(design) <- make.names(colnames(design))

  fit <- limma::lmFit(E, design)
  fit <- limma::eBayes(fit)

  cmat <- limma::makeContrasts(
    WT_leptin_180_vs_0        = leptin180min,
    Stat3DEL_leptin_180_vs_0  = leptin180min + genotypeStat3Δ.leptin180min,
    Baseline_Stat3DEL_vs_WT   = genotypeStat3Δ,
    At180_Stat3DEL_vs_WT      = genotypeStat3Δ + genotypeStat3Δ.leptin180min,
    Interaction_DiffInDiff    = genotypeStat3Δ.leptin180min,
    levels = design
  )

  fit2 <- limma::contrasts.fit(fit, cmat)
  fit2 <- limma::eBayes(fit2)

  list(fit = fit, fit2 = fit2, design = design, contrasts = cmat)
}

#--------------------------------------------------------------
# Extract limma topTable and add protein identifiers.
#--------------------------------------------------------------
extract_dap_table <- function(fit2, coef, protein_annot) {
  limma::topTable(fit2, coef = coef, number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("Accession") %>%
    dplyr::left_join(
      protein_annot %>% dplyr::select(dplyr::any_of(c("Accession", "Gene Symbol", "Description"))),
      by = "Accession"
    )
}

#--------------------------------------------------------------
# Prepare signed -log10(P) comparison between proteomics and RNA-seq.
#
# This reproduces the original scatterplot metric:
#   signed score = -log10(P value) * sign(log2FC)
#--------------------------------------------------------------
prepare_signed_comparison <- function(dap_table, deg_table,
                                      deg_p_overrides = NULL,
                                      threshold = 1.30103) {
  daps <- dap_table

  # Preserve the original positional renaming logic after topTable/add_ids.
  colnames(daps)[2] <- "DAP_log2fc"
  colnames(daps)[5] <- "DAP_pvalue"
  colnames(daps)[6] <- "DAP_adj_pvalue"
  colnames(daps)[8] <- "id"

  daps$id <- stringr::str_to_title(daps$id)
  daps <- daps %>%
    dplyr::mutate(
      logDAP_pvalue  = -log10(DAP_pvalue),
      sign_DAPpvalue = logDAP_pvalue * sign(DAP_log2fc)
    )

  degs <- deg_table
  colnames(degs)[3] <- "DEG_log2fc"
  colnames(degs)[6] <- "DEG_pvalue"
  colnames(degs)[7] <- "DEG_adj_pvalue"

  if (!is.null(deg_p_overrides)) {
    for (i in seq_along(deg_p_overrides)) {
      degs[i, 6] <- deg_p_overrides[[i]]
    }
  }

  degs <- degs %>%
    dplyr::mutate(
      logDEG_pvalue  = -log10(DEG_pvalue),
      sign_DEGpvalue = logDEG_pvalue * sign(DEG_log2fc)
    )

  com_data <- dplyr::inner_join(daps, degs, by = "id") %>%
    tibble::column_to_rownames("id")

  com_data$color <- ifelse(
    com_data$sign_DEGpvalue > threshold & com_data$sign_DAPpvalue > threshold,
    "firebrick3",
    ifelse(
      com_data$sign_DEGpvalue < -threshold & com_data$sign_DAPpvalue < -threshold,
      "royalblue3",
      "grey"
    )
  )

  com_data
}

#--------------------------------------------------------------
# Plot signed -log10(P) proteome/transcriptome concordance.
#--------------------------------------------------------------
plot_signed_scatter <- function(com_data) {
  ggplot2::ggplot(com_data, ggplot2::aes(x = sign_DEGpvalue, y = sign_DAPpvalue)) +
    ggplot2::geom_point(data = subset(com_data, color == "grey"), color = "grey", size = 1.2) +
    ggplot2::geom_point(data = subset(com_data, color == "firebrick3"), color = "firebrick3", size = 1.6) +
    ggplot2::geom_point(data = subset(com_data, color == "royalblue3"), color = "royalblue3", size = 1.6) +
    ggplot2::geom_hline(yintercept = c(1.30103, -1.30103), linetype = "dashed", color = "black", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = c(1.30103, -1.30103), linetype = "dashed", color = "black", linewidth = 0.25) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "−log10(P) × sign(FC) transcriptome",
      y = "−log10(P) × sign(FC) proteome"
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 10),
      axis.text  = ggplot2::element_text(size = 10)
    )
}

#--------------------------------------------------------------
# Run signed scatterplot analysis, save PDF and Excel output.
#--------------------------------------------------------------
run_signed_scatter_analysis <- function(dap_table, deg_path, deg_sheet,
                                        pdf_file, xlsx_file,
                                        deg_p_overrides = NULL,
                                        annotate_x = -30,
                                        annotate_y = 6) {
  deg_table <- readxl::read_xlsx(deg_path, sheet = deg_sheet)

  com_data <- prepare_signed_comparison(
    dap_table = dap_table,
    deg_table = deg_table,
    deg_p_overrides = deg_p_overrides
  )

  p <- plot_signed_scatter(com_data)
  res <- stats::cor.test(com_data$sign_DEGpvalue, com_data$sign_DAPpvalue, method = "spearman")

  grDevices::pdf(file = pdf_file, width = 4, height = 4)
  print(
    p +
      ggplot2::theme_classic() +
      ggplot2::labs(
        x = "-log10(Pval) × sign(FC) transcriptome",
        y = "-log10(Pval) × sign(FC) proteome"
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(face = "plain", size = 10),
        axis.text  = ggplot2::element_text(face = "plain")
      ) +
      ggplot2::annotate(
        "text", x = annotate_x, y = annotate_y,
        label = paste("rho= ", round(res$estimate, 4)),
        size = 3, color = "grey40"
      )
  )
  grDevices::dev.off()

  writexl::write_xlsx(tibble::rownames_to_column(com_data, var = "id"), xlsx_file)

  list(data = com_data, plot = p, cor_test = res)
}

#--------------------------------------------------------------
# Build a row-scaled heatmap matrix from significant DAPs.
#--------------------------------------------------------------
make_proteomics_heatmap_matrix <- function(E, protein_annot, dap_table,
                                           sample_order,
                                           fdr_cutoff = 0.05,
                                           lfc_cutoff = 1) {
  sig <- dap_table %>%
    dplyr::filter(adj.P.Val < fdr_cutoff) %>%
    dplyr::filter(logFC >= lfc_cutoff | logFC <= -lfc_cutoff)

  log2_data_ht <- as.data.frame(E) %>%
    tibble::rownames_to_column("Accession") %>%
    dplyr::left_join(
      protein_annot %>% dplyr::select(dplyr::any_of(c("Accession", "Gene Symbol", "Description"))),
      by = "Accession"
    )

  ht_data <- log2_data_ht %>%
    dplyr::filter(`Gene Symbol` %in% sig$`Gene Symbol`) %>%
    dplyr::filter(!is.na(`Gene Symbol`))

  m <- ht_data %>%
    dplyr::select(-dplyr::any_of(c("Description"))) %>%
    tibble::column_to_rownames("Gene Symbol") %>%
    dplyr::select(-dplyr::any_of("Accession")) %>%
    as.matrix()

  m <- t(scale(t(m)))
  m[, sample_order, drop = FALSE]
}

#--------------------------------------------------------------
# Create ComplexHeatmap annotation for proteomics samples.
#--------------------------------------------------------------
make_proteomics_col_annotation <- function(sample_info, sample_order,
                                           include_genotype = TRUE) {
  ann <- data.frame(
    genotype = sample_info$genotype,
    samples  = sample_info$SampleID,
    leptin   = sample_info$leptin
  )
  ann <- ann[match(sample_order, ann$samples), , drop = FALSE]
  ann$leptin <- factor(ann$leptin, levels = c("0min", "180min"))
  ann$genotype <- factor(ann$genotype, levels = c("WT", "Stat3Δ"))

  if (include_genotype) {
    ComplexHeatmap::HeatmapAnnotation(
      genotype = ann$genotype,
      leptin = ann$leptin,
      col = list(
        genotype = c(WT = "#46464E", Stat3Δ = "#1F7A8C"),
        leptin = c(`0min` = "grey75", `180min` = "#5BC0D0")
      ),
      annotation_name_side = "left"
    )
  } else {
    ComplexHeatmap::HeatmapAnnotation(
      leptin = ann$leptin,
      col = list(leptin = c("0min" = "#A4A4AE", "180min" = "#46464E")),
      annotation_legend_param = list(leptin = list(nrow = 1, title = "leptin"))
    )
  }
}

#--------------------------------------------------------------
# Draw and save proteomics heatmap.
#--------------------------------------------------------------
save_proteomics_heatmap <- function(m, colAnn, pdf_file, title,
                                    width = 5, height = 7,
                                    color_limits = c(-3, 0, 3),
                                    show_row_names = TRUE,
                                    show_column_names = FALSE,
                                    row_fontsize = 10,
                                    heatmap_width_mm = 3,
                                    heatmap_height_mm = 3,
                                    column_split = NULL) {
  grDevices::pdf(file = pdf_file, width = width, height = height)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      m,
      name = "Abundance",
      col = circlize::colorRamp2(color_limits, c("royalblue3", "white", "firebrick3")),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      column_split = column_split,
      top_annotation = colAnn,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      show_row_dend = FALSE,
      row_names_gp = grid::gpar(fontsize = row_fontsize),
      width = ncol(as.matrix(m)) * grid::unit(heatmap_width_mm, "mm"),
      height = nrow(as.matrix(m)) * grid::unit(heatmap_height_mm, "mm"),
      column_title = title,
      heatmap_legend_param = list(legend_direction = "horizontal")
    ),
    heatmap_legend_side = "bottom"
  )
  grDevices::dev.off()
}
