# ============================================================
# DEG helper functions for DESeq2 leptin time-course analyses
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(writexl)
  library(EnhancedVolcano)
})

# -----------------------------
# Small utilities
# -----------------------------

safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  substr(x, 1, 31)
}

check_coef_exists <- function(dds, coef_names) {
  coef_names <- unique(coef_names[!is.na(coef_names) & nzchar(coef_names)])
  missing <- setdiff(coef_names, resultsNames(dds))
  if (length(missing) > 0) {
    stop(
      "Missing coefficient(s) in resultsNames(dds):\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

coef_names_from_spec <- function(spec) {
  if (!is.null(spec$name)) {
    return(spec$name)
  }

  if (!is.null(spec$contrast)) {
    # For contrast = c("time", "180min", "30min"), these are factor levels,
    # not coefficient names, so do not check against resultsNames(dds).
    if (is.character(spec$contrast) && length(spec$contrast) == 3) {
      return(character(0))
    }

    # For contrast = list(c("coef1", "coef2"), c("coef3")), these are coefficient names.
    if (is.list(spec$contrast)) {
      return(unlist(spec$contrast, use.names = FALSE))
    }
  }

  character(0)
}

# -----------------------------
# Extract one DESeq2 result
# -----------------------------

get_deseq_result <- function(dds, spec) {
  if (!is.null(spec$name)) {
    return(results(dds, name = spec$name))
  }

  if (!is.null(spec$contrast)) {
    return(results(dds, contrast = spec$contrast))
  }

  stop("Each comparison spec must contain either `name` or `contrast`.", call. = FALSE)
}

# -----------------------------
# Process one DEG result table
# -----------------------------

process_deg_result <- function(res,
                               comparison,
                               padj_cutoff = 0.05,
                               lfc_cutoff = 1) {

  res_df <- as.data.frame(res)
  res_df$id <- rownames(res_df)

  res_df <- res_df %>%
    dplyr::select(id, dplyr::everything()) %>%
    dplyr::mutate(comparison = comparison) %>%
    dplyr::arrange(padj, dplyr::desc(log2FoldChange))

  sig_df <- res_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff)

  up_df <- sig_df %>%
    dplyr::filter(log2FoldChange >= lfc_cutoff)

  down_df <- sig_df %>%
    dplyr::filter(log2FoldChange <= -lfc_cutoff)

  list(
    all = res_df,
    sig = sig_df,
    up = up_df,
    down = down_df
  )
}

# -----------------------------
# Optional custom colours for volcano plots
# -----------------------------

make_volcano_keyvals <- function(res,
                                 padj_cutoff = 0.05,
                                 lfc_cutoff = 1,
                                 up_colour = "#B56336",
                                 down_colour = "#E0B9A0",
                                 ns_colour = "grey70",
                                 na_colour = "grey90",
                                 up_label = "Up",
                                 down_label = "Down",
                                 ns_label = "NS") {

  keyvals <- ifelse(
    !is.na(res$padj) & res$padj < padj_cutoff & res$log2FoldChange >= lfc_cutoff,
    up_colour,
    ifelse(
      !is.na(res$padj) & res$padj < padj_cutoff & res$log2FoldChange <= -lfc_cutoff,
      down_colour,
      ns_colour
    )
  )

  keyvals[is.na(keyvals)] <- na_colour

  names(keyvals)[keyvals == up_colour] <- up_label
  names(keyvals)[keyvals == down_colour] <- down_label
  names(keyvals)[keyvals %in% c(ns_colour, na_colour)] <- ns_label

  keyvals
}

# -----------------------------
# Basic volcano plot function
# -----------------------------

plot_deg_volcano <- function(res,
                             output_file,
                             title,
                             xlim = c(-5, 5),
                             padj_cutoff = 0.05,
                             lfc_cutoff = 1,
                             width = 6,
                             height = 8,
                             custom_colours = FALSE,
                             up_label = "Up",
                             down_label = "Down",
                             lab_size = 4,
                             title_size = 12,
                             legend_size = 10,
                             axis_size = 12) {

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  volcano_args <- list(
    toptable = res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    labSize = lab_size,
    ylab = "-log10 adjusted P value",
    title = title,
    subtitle = "",
    titleLabSize = title_size,
    subtitleLabSize = 10,
    legendLabSize = legend_size,
    axisLabSize = axis_size,
    xlim = xlim
  )

  if (isTRUE(custom_colours)) {
    volcano_args$colCustom <- make_volcano_keyvals(
      res,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff,
      up_label = up_label,
      down_label = down_label
    )
  }

  pdf(file = output_file, width = width, height = height)
  print(do.call(EnhancedVolcano, volcano_args))
  dev.off()
}

# -----------------------------
# Run and export a set of DEG comparisons
# -----------------------------

run_deg_set <- function(dds,
                        comparisons,
                        outdir,
                        xlsx_file,
                        padj_cutoff = 0.05,
                        lfc_cutoff = 1,
                        make_volcano = TRUE,
                        custom_volcano_colours = FALSE) {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Check coefficient names where possible.
  needed <- unique(unlist(lapply(comparisons, coef_names_from_spec), use.names = FALSE))
  check_coef_exists(dds, needed)

  result_objects <- list()
  processed <- list()

  for (comparison_name in names(comparisons)) {
    spec <- comparisons[[comparison_name]]

    res <- get_deseq_result(dds, spec)
    result_objects[[comparison_name]] <- res

    processed[[comparison_name]] <- process_deg_result(
      res = res,
      comparison = comparison_name,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    )

    if (isTRUE(make_volcano)) {
      plot_deg_volcano(
        res = res,
        output_file = file.path(outdir, paste0("Volcano_", comparison_name, ".pdf")),
        title = if (is.null(spec$title)) comparison_name else spec$title,
        xlim = if (is.null(spec$xlim)) c(-5, 5) else spec$xlim,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        custom_colours = custom_volcano_colours,
        up_label = if (is.null(spec$up_label)) "Up" else spec$up_label,
        down_label = if (is.null(spec$down_label)) "Down" else spec$down_label,
        lab_size = ifelse(isTRUE(custom_volcano_colours), 8, 4),
        title_size = ifelse(isTRUE(custom_volcano_colours), 18, 12),
        legend_size = ifelse(isTRUE(custom_volcano_colours), 14, 10),
        axis_size = ifelse(isTRUE(custom_volcano_colours), 18, 12)
      )
    }
  }

  combined_all <- dplyr::bind_rows(lapply(processed, `[[`, "all"))
  combined_sig <- dplyr::bind_rows(lapply(processed, `[[`, "sig"))
  combined_up <- dplyr::bind_rows(lapply(processed, `[[`, "up"))
  combined_down <- dplyr::bind_rows(lapply(processed, `[[`, "down"))

  export_list <- list()

  for (comparison_name in names(processed)) {
    short <- gsub("[^A-Za-z0-9_]+", "_", comparison_name)
    export_list[[safe_sheet_name(paste0("all_", short))]] <- processed[[comparison_name]]$all
    export_list[[safe_sheet_name(paste0("sig_", short))]] <- processed[[comparison_name]]$sig
    export_list[[safe_sheet_name(paste0("up_", short))]] <- processed[[comparison_name]]$up
    export_list[[safe_sheet_name(paste0("down_", short))]] <- processed[[comparison_name]]$down
  }

  export_list[["combined_all"]] <- combined_all
  export_list[["combined_sig"]] <- combined_sig
  export_list[["combined_up"]] <- combined_up
  export_list[["combined_down"]] <- combined_down

  output_xlsx <- file.path(outdir, xlsx_file)
  writexl::write_xlsx(export_list, output_xlsx)

  invisible(list(
    results = result_objects,
    tables = processed,
    combined_all = combined_all,
    combined_sig = combined_sig,
    combined_up = combined_up,
    combined_down = combined_down,
    output_xlsx = output_xlsx
  ))
}

# ============================================================
# Comparison specification helpers
# ============================================================

# 1) WT within-genotype time effects
make_wt_time_specs <- function() {
  list(
    WT_30v0 = list(
      name = "time_30min_vs_0min",
      title = "WT stimulated for 30 min vs unstimulated",
      xlim = c(-3, 7),
      up_label = "Up at 30 min",
      down_label = "Up at 0 min"
    ),
    WT_180v0 = list(
      name = "time_180min_vs_0min",
      title = "WT stimulated for 180 min vs unstimulated",
      xlim = c(-5, 7),
      up_label = "Up at 180 min",
      down_label = "Up at 0 min"
    ),
    WT_180v30 = list(
      contrast = c("time", "180min", "30min"),
      title = "WT stimulated for 180 min vs 30 min",
      xlim = c(-10, 10),
      up_label = "Up at 180 min",
      down_label = "Up at 30 min"
    )
  )
}

# 2) KO within-genotype time effects: KO 30v0, KO 180v0, KO 180v30
make_ko_time_specs <- function(ko_label,
                               interaction_prefix,
                               coef_30v0 = "time_30min_vs_0min",
                               coef_180v0 = "time_180min_vs_0min") {

  int_30 <- paste0(interaction_prefix, ".time30min")
  int_180 <- paste0(interaction_prefix, ".time180min")

  specs <- list(
    `30v0` = list(
      contrast = list(c(coef_30v0, int_30)),
      title = paste0(ko_label, " stimulated for 30 min vs unstimulated"),
      xlim = c(-3, 7),
      up_label = "Up at 30 min",
      down_label = "Up at 0 min"
    ),
    `180v0` = list(
      contrast = list(c(coef_180v0, int_180)),
      title = paste0(ko_label, " stimulated for 180 min vs unstimulated"),
      xlim = c(-5, 7),
      up_label = "Up at 180 min",
      down_label = "Up at 0 min"
    ),
    `180v30` = list(
      contrast = list(
        c(coef_180v0, int_180),
        c(coef_30v0, int_30)
      ),
      title = paste0(ko_label, " stimulated for 180 min vs 30 min"),
      xlim = c(-10, 10),
      up_label = "Up at 180 min",
      down_label = "Up at 30 min"
    )
  )

  stats::setNames(specs, paste0(ko_label, "_", names(specs)))
}

# 3) Interaction / delta-response effects relative to WT
# 30v0 = KO 30v0 response minus WT 30v0 response
# 180v0 = KO 180v0 response minus WT 180v0 response
# 180v30 = (KO 180v30 response) minus (WT 180v30 response)
make_interaction_specs <- function(ko_label,
                                   interaction_prefix) {

  int_30 <- paste0(interaction_prefix, ".time30min")
  int_180 <- paste0(interaction_prefix, ".time180min")

  specs <- list(
    `Interac_30v0` = list(
      name = int_30,
      title = paste0("Difference in 0 to 30 min response between ", ko_label, " and WT"),
      xlim = c(-7, 5),
      up_label = paste0("Higher response in ", ko_label),
      down_label = "Higher response in WT"
    ),
    `Interac_180v0` = list(
      name = int_180,
      title = paste0("Difference in 0 to 180 min response between ", ko_label, " and WT"),
      xlim = c(-11, 10),
      up_label = paste0("Higher response in ", ko_label),
      down_label = "Higher response in WT"
    ),
    `Interac_180v30` = list(
      contrast = list(c(int_180), c(int_30)),
      title = paste0("Difference in 30 to 180 min response between ", ko_label, " and WT"),
      xlim = c(-10, 10),
      up_label = paste0("Higher response in ", ko_label),
      down_label = "Higher response in WT"
    )
  )

  stats::setNames(specs, paste0(ko_label, "_", names(specs)))
}

# 4) Genotype difference versus WT at each time point
make_vs_wt_specs <- function(ko_label,
                             type_coef,
                             interaction_prefix) {

  int_30 <- paste0(interaction_prefix, ".time30min")
  int_180 <- paste0(interaction_prefix, ".time180min")

  specs <- list(
    `vsWT_0` = list(
      name = type_coef,
      title = paste0(ko_label, " vs WT: 0 min"),
      xlim = c(-15, 12),
      up_label = paste0("Higher in ", ko_label),
      down_label = "Higher in WT"
    ),
    `vsWT_30` = list(
      contrast = list(c(type_coef, int_30)),
      title = paste0(ko_label, " vs WT: 30 min"),
      xlim = c(-15, 12),
      up_label = paste0("Higher in ", ko_label),
      down_label = "Higher in WT"
    ),
    `vsWT_180` = list(
      contrast = list(c(type_coef, int_180)),
      title = paste0(ko_label, " vs WT: 180 min"),
      xlim = c(-15, 12),
      up_label = paste0("Higher in ", ko_label),
      down_label = "Higher in WT"
    )
  )

  stats::setNames(specs, paste0(ko_label, "_", names(specs)))
}

# ============================================================
# Example calls
# ============================================================

# WT time effects
# wt_time <- run_deg_set(
#   dds = dds,
#   comparisons = make_wt_time_specs(),
#   outdir = "commom_leptin_response/wt_time_effect",
#   xlsx_file = "degs_wt_time_effect.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE,
#   custom_volcano_colours = FALSE
# )

# Stat3DEL within-genotype time effects
# stat3_time <- run_deg_set(
#   dds = dds,
#   comparisons = make_ko_time_specs("Stat3DEL", "typeStat3_DEL"),
#   outdir = "commom_leptin_response/stat3_del_time_effect",
#   xlsx_file = "degs_stat3del_time_effect.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )

# Etv6DEL within-genotype time effects
# etv6_time <- run_deg_set(
#   dds = dds,
#   comparisons = make_ko_time_specs("Etv6DEL", "typeEtv6_DEL"),
#   outdir = "commom_leptin_response/Etv6_del_time_effect",
#   xlsx_file = "degs_etv6del_time_effect.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )

# Stat3DEL delta-response interaction effects relative to WT
# stat3_interaction <- run_deg_set(
#   dds = dds,
#   comparisons = make_interaction_specs("Stat3DEL", "typeStat3_DEL"),
#   outdir = "commom_leptin_response/stat3_interaction_effect",
#   xlsx_file = "degs_stat3_interaction_effect.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )

# Etv6DEL delta-response interaction effects relative to WT
# etv6_interaction <- run_deg_set(
#   dds = dds,
#   comparisons = make_interaction_specs("Etv6DEL", "typeEtv6_DEL"),
#   outdir = "commom_leptin_response/Etv6_interaction_effect",
#   xlsx_file = "degs_etv6_interaction_effect.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )

# Stat3DEL vs WT at each time point
# stat3_vs_wt <- run_deg_set(
#   dds = dds,
#   comparisons = make_vs_wt_specs("Stat3DEL", "type_Stat3_DEL_vs_WT", "typeStat3_DEL"),
#   outdir = "commom_leptin_response/Stat3_vsWT",
#   xlsx_file = "degs_stat3_vsWT.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )

# Etv6DEL vs WT at each time point
# etv6_vs_wt <- run_deg_set(
#   dds = dds,
#   comparisons = make_vs_wt_specs("Etv6DEL", "type_Etv6_DEL_vs_WT", "typeEtv6_DEL"),
#   outdir = "commom_leptin_response/Etv6_vsWT",
#   xlsx_file = "degs_etv6_vsWT.xlsx",
#   padj_cutoff = 0.05,
#   lfc_cutoff = 1,
#   make_volcano = TRUE
# )
