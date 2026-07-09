############################################################
# Script: 1_proteomics_main.R
# Project: Temporal organization of leptin signaling reveals distinct
#          regulatory roles for STAT3 and ETV6
#
# Purpose:
#   Analyse normalized proteomics abundance data from the leptin
#   stimulation time course and Stat3 deletion experiment.
#
# Description:
#   This script imports normalized proteomics abundance data, prepares
#   sample metadata, performs PCA, fits limma models for leptin-induced
#   and Stat3-dependent protein abundance changes, compares proteomic
#   responses with RNA-seq-derived signed p-value statistics, and
#   generates proteomics heatmaps and summary plots.
#
# Required upstream analyses:
#   Run the 1_RNAseq_degs_heatmap_pathways.R analysis scripts before
#   this script. This workflow expects
#   RNA-seq-derived differential expression results, pathway results,
#   and activity/correlation objects used for proteome-transcriptome
#   comparisons.
#
# Helper script:
#   proteomics_helper_functions.R
############################################################

# ============================================================
# 1. Load packages and helper functions
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(limma)
library(readxl)
library(writexl)
library(ComplexHeatmap)
library(circlize)
library(grid)

source("R/helpers/proteomics_helper_functions.R")

# ============================================================
# 2. File paths and analysis settings
# ============================================================

proteomics_file <- "Data/Proteomics_abundance.xlsx"
out_dir <- "proteomic_dta_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Sample order used in the original heatmaps.
wt_sample_order <- c(
  "F3_WT_minusLep", "F13_WT_minusLep", "F14_WT_minusLep", "F15_WT_minusLep",
  "F1_WT_plusLep", "F6_WT_plusLep", "F10_WT_plusLep"
)

all_sample_order <- c(
  "F3_WT_minusLep", "F13_WT_minusLep", "F14_WT_minusLep", "F15_WT_minusLep",
  "F1_WT_plusLep", "F6_WT_plusLep", "F10_WT_plusLep",
  "F2_S15_minusLep", "F4_S15_minusLep", "F8_S15_minusLep", "F9_S15_minusLep",
  "F7_S15_plusLep", "F11_S15_plusLep", "F12_S15_plusLep", "F16_S15_plusLep"
)

# ============================================================
# 3. Load proteomics data and metadata
# ============================================================

protein_raw <- readxl::read_xlsx(proteomics_file, sheet = 3)
metadata <- readxl::read_xlsx(proteomics_file, sheet = 4)

sample_info <- make_sample_info(metadata)

E <- make_log2_protein_matrix(
  protein_tbl = protein_raw,
  sample_info = sample_info,
  accession_col = "Accession",
  remove_col_index = 2,
  offset = 1,
  min_prop = 0.5,
  min_per_group = 2
)

protein_annot <- protein_raw

# ============================================================
# 4. PCA analysis
# ============================================================

p_pca_all <- plot_proteomics_pca(
  E = E,
  sample_info = sample_info,
  title = "Proteomics PCA",
  color_var = "genotype",
  shape_var = "leptin",
  point_size = 3,
  color_values = c(WT = "#46464E", Stat3Δ = "#1F7A8C")
)

ggplot2::ggsave(
  filename = file.path(out_dir, "proteomic_pca.pdf"),
  plot = p_pca_all,
  width = 6,
  height = 6
)

wt_idx <- sample_info$genotype == "WT"

p_pca_wt <- plot_proteomics_pca(
  E = E[, wt_idx, drop = FALSE],
  sample_info = sample_info[wt_idx, , drop = FALSE],
  title = "WT:Proteomics PCA",
  color_var = "leptin",
  shape_var = "leptin",
  point_size = 5,
  fixed_point_color = "#46464E"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "proteomic_pca_WTonly_1.pdf"),
  plot = p_pca_wt,
  width = 6,
  height = 6
)

# ============================================================
# 5. Differential abundance analysis with limma
# ============================================================

limma_out <- run_proteomics_limma(E, sample_info)
fit2 <- limma_out$fit2

res_wt_lep <- extract_dap_table(fit2, "WT_leptin_180_vs_0", protein_annot)
res_ko_lep <- extract_dap_table(fit2, "Stat3DEL_leptin_180_vs_0", protein_annot)
res_base <- extract_dap_table(fit2, "Baseline_Stat3DEL_vs_WT", protein_annot)
res_180 <- extract_dap_table(fit2, "At180_Stat3DEL_vs_WT", protein_annot)
res_interaction <- extract_dap_table(fit2, "Interaction_DiffInDiff", protein_annot)

# ============================================================
# 6. Transcriptome-proteome signed p-value scatterplots
# ============================================================

scatter_wt <- run_signed_scatter_analysis(
  dap_table = res_wt_lep,
  deg_path = "transcriptomics_data_analysis/wt_time_effect/degs_wt_time_effect_test.xlsx",#
  deg_sheet = 5,
  pdf_file = file.path(out_dir, "DEG_WT180vsWT0_vs_DAP_WT180vsWT0_SIGNED_pvalue.pdf"),
  xlsx_file = file.path(out_dir, "DEG_WT180vsWT0_vs_DAP_WT180vsWT0_SIGNED_pvalue.xlsx"),
  deg_p_overrides = c(10^(-276))
)

scatter_base <- run_signed_scatter_analysis(
  dap_table = res_base,
  deg_path = "transcriptomics_data_analysis/Stat3_vsWT/degs_stat3_vsWT_test.xlsx",#
  deg_sheet = 1,
  pdf_file = file.path(out_dir, "DEG_Stat3vsWT0_vs_DAP_Stat3vsWT0_SIGNED_pvalue.pdf"),
  xlsx_file = file.path(out_dir, "DEG_Stat3vsWT0_vs_DAP_Stat3vsWT0_SIGNED_pvalue.xlsx"),
  deg_p_overrides = rep(10^(-267), 3)
)

scatter_180 <- run_signed_scatter_analysis(
  dap_table = res_180,
  deg_path = "transcriptomics_data_analysis/Stat3_vsWT/degs_stat3_vsWT_test.xlsx",#
  deg_sheet = 9,
  pdf_file = file.path(out_dir, "DEG_Stat3vsWT180_vs_DAP_Stat3vsWT180_SIGNED_pvalue.pdf"),
  xlsx_file = file.path(out_dir, "DEG_Stat3vsWT180_vs_DAP_Stat3vsWT180_SIGNED_pvalue.xlsx"),
  deg_p_overrides = rep(10^(-302), 6)
)

scatter_interaction <- run_signed_scatter_analysis(
  dap_table = res_interaction,
  deg_path = "transcriptomics_data_analysis/stat3_interaction_effect/degs_stat3_interaction_effect.xlsx",
  deg_sheet = 5,
  pdf_file = file.path(out_dir, "DEG_Interaction_vs_DAP_Interaction_SIGNED_pvalue.pdf"),
  xlsx_file = file.path(out_dir, "DEG_Interaction_vs_DAP_Interaction_SIGNED_pvalue.xlsx")
)

# ============================================================
# 7. Proteomics heatmaps
# ============================================================

# WT leptin-responsive proteins.
m_wt_lep <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_wt_lep,
  sample_order = wt_sample_order,
  fdr_cutoff = 0.05,
  lfc_cutoff = 1
)

ann_wt <- make_proteomics_col_annotation(
  sample_info = sample_info,
  sample_order = wt_sample_order,
  include_genotype = FALSE
)

save_proteomics_heatmap(
  m = m_wt_lep,
  colAnn = ann_wt,
  pdf_file = file.path(out_dir, "heatmap_daps_WTleptin.pdf"),
  title = "WT:Leptin stimulation effect",
  width = 4,
  height = 4,
  color_limits = c(-2, 0, 2),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 10,
  heatmap_width_mm = 3,
  heatmap_height_mm = 3
)

# Common annotation for all genotype/time heatmaps.
ann_all <- make_proteomics_col_annotation(
  sample_info = sample_info,
  sample_order = all_sample_order,
  include_genotype = TRUE
)

ann_for_split <- data.frame(
  samples = sample_info$SampleID,
  leptin = sample_info$leptin
)
ann_for_split <- ann_for_split[match(all_sample_order, ann_for_split$samples), , drop = FALSE]
ann_for_split$leptin <- factor(ann_for_split$leptin, levels = c("0min", "180min"))

# Stat3 deletion × leptin interaction proteins.
m_interaction <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_interaction,
  sample_order = all_sample_order,
  fdr_cutoff = 0.05,
  lfc_cutoff = 1
)

save_proteomics_heatmap(
  m = m_interaction,
  colAnn = ann_all,
  pdf_file = file.path(out_dir, "Heatmap_proteomic_Interaction_effect.pdf"),
  title = "Stat3 Del X time effect proteins",
  width = 5,
  height = 6,
  color_limits = c(-3, 0, 3),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 12,
  column_split = ann_for_split$leptin
)

# Stat3Δ leptin response proteins.
m_ko_lep <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_ko_lep,
  sample_order = all_sample_order,
  fdr_cutoff = 0.05,
  lfc_cutoff = 1
)

save_proteomics_heatmap(
  m = m_ko_lep,
  colAnn = ann_all,
  pdf_file = file.path(out_dir, "Heatmap_proteomic_Stat3del_180_0min_effect.pdf"),
  title = "Stat3 180vs0min protein",
  width = 5,
  height = 7,
  color_limits = c(-3, 0, 3),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 12,
  column_split = ann_for_split$leptin
)

# Stat3Δ vs WT baseline proteins.
m_base <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_base,
  sample_order = all_sample_order,
  fdr_cutoff = 0.05,
  lfc_cutoff = 1
)

save_proteomics_heatmap(
  m = m_base,
  colAnn = ann_all,
  pdf_file = file.path(out_dir, "Heatmap_proteomic_Stat3delvsWT_0min_effect_withnames.pdf"),
  title = "Stat3 vs WT 0min protein",
  width = 5,
  height = 22,
  color_limits = c(-3, 0, 3),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 10,
  heatmap_width_mm = 3,
  heatmap_height_mm = 2.5,
  column_split = ann_for_split$leptin
)

# Stat3Δ vs WT at 180 min proteins.
m_180 <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_180,
  sample_order = all_sample_order,
  fdr_cutoff = 0.05,
  lfc_cutoff = 1
)

save_proteomics_heatmap(
  m = m_180,
  colAnn = ann_all,
  pdf_file = file.path(out_dir, "Heatmap_proteomic_Stat3delvsWT_180min_effect_withnames.pdf"),
  title = "Stat3 vs WT 180min protein",
  width = 5,
  height = 12,
  color_limits = c(-3, 0, 3),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 10,
  heatmap_width_mm = 3,
  heatmap_height_mm = 3,
  column_split = ann_for_split$leptin
)

# Combined Stat3Δ vs WT proteins from 0 and 180 min contrasts.
combined_symbols <- unique(c(
  res_base %>% dplyr::filter(adj.P.Val < 0.05, logFC >= 1 | logFC <= -1) %>% dplyr::pull(`Gene Symbol`),
  res_180  %>% dplyr::filter(adj.P.Val < 0.05, logFC >= 1 | logFC <= -1) %>% dplyr::pull(`Gene Symbol`)
))

res_combined <- res_base %>%
  dplyr::filter(`Gene Symbol` %in% combined_symbols)

m_combined <- make_proteomics_heatmap_matrix(
  E = E,
  protein_annot = protein_annot,
  dap_table = res_combined,
  sample_order = all_sample_order,
  fdr_cutoff = 1,
  lfc_cutoff = 0
)

save_proteomics_heatmap(
  m = m_combined,
  colAnn = ann_all,
  pdf_file = file.path(out_dir, "Heatmap_proteomic_Stat3delvsWT_180_0_min_combined_effect_withnames.pdf"),
  title = "Stat3 vs WT protein",
  width = 5,
  height = 24,
  color_limits = c(-3, 0, 3),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_fontsize = 10,
  heatmap_width_mm = 2.5,
  heatmap_height_mm = 2.5,
  column_split = ann_for_split$leptin
)
