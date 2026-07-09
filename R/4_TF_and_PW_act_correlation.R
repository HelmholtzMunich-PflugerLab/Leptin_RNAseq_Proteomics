############################################################
# Script: 4_TF_and_PW_act_correlation.R
# Project: Temporal organization of leptin signaling reveals distinct
#          regulatory roles for STAT3 and ETV6
#
# Purpose:
#   Correlate WT transcription factor and pathway activities with
#   WT time-course gene-expression modules.
#
# Description:
#   This script compares regulatory activity dynamics with expression
#   modules derived from WT leptin-responsive genes. For each WT gene
#   cluster, a module score is calculated as the mean row-scaled,
#   batch-corrected VST expression across genes in that cluster.
#   Spearman correlations are then computed between these module scores
#   and inferred TF/pathway activities across WT samples.
#
# Required input objects:
#   vsd                  : VST-transformed DESeq2 object.
#   mat_bc               : Batch-corrected, non-scaled VST expression matrix.
#   clu_df               : WT DEG k-means cluster assignment data frame.
#   M_tf                 : TF activity matrix from DoRothEA/decoupleR.
#   M_pw                 : Pathway activity matrix from PROGENy/decoupleR.
#   features_filtered_tf : TFs selected from WT TF activity dynamics.
#   features_filtered_pw : Pathways selected from WT pathway activity dynamics.
#   clu_df_tf71          : TF-to-cluster assignment from WT TF heatmap.
#
# Required metadata columns in colData(vsd):
#   type, time, batch, samples
#
# Main outputs:
#   cor_mat_tf_show : TF activity vs WT DEG module correlation matrix.
#   cor_mat_pw_show : Pathway activity vs WT DEG module correlation matrix.
#   Heatmaps showing TF/pathway correlations with WT expression modules.
#
# Notes:
#   - Run 1_RNAseq_degs_heatmap_pathways.R and 2_Tfs_Pathways_activity_dynamics.R
#     before this script
#   - Correlations are calculated across WT samples only.
#   - Spearman correlation is used because it captures monotonic
#     relationships without assuming linearity.
#   - Expression module scores are based on batch-corrected VST values.
############################################################


# ============================================================
# 1. Load required packages
# ============================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)


# ============================================================
# 2. Check required input objects
# ============================================================

stopifnot(exists("vsd"))
stopifnot(exists("mat_bc"))
stopifnot(exists("clu_df"))
stopifnot(exists("M_tf"))
stopifnot(exists("M_pw"))
stopifnot(exists("features_filtered_tf"))
stopifnot(exists("features_filtered_pw"))
stopifnot(exists("clu_df_tf71"))


# ============================================================
# 3. Prepare sample metadata
# ============================================================

# Extract sample metadata from the VST object and enforce the
# same factor levels used in the main activity-dynamics script.
meta <- as.data.frame(colData(vsd))

stopifnot(all(c("type", "time", "batch", "samples") %in% colnames(meta)))

meta$time  <- factor(meta$time, levels = c("0min", "30min", "180min"))
meta$type  <- factor(meta$type, levels = c("WT", "Stat3_DEL", "Etv6_DEL"))
meta$batch <- factor(meta$batch)


# Select WT samples and order them chronologically.
wt_cols <- meta %>%
  mutate(time = factor(time, levels = c("0min", "30min", "180min"))) %>%
  filter(type == "WT") %>%
  arrange(time) %>%
  pull(samples)


# ============================================================
# 4. Calculate WT gene-expression module scores
# ============================================================

# Subset batch-corrected VST expression to ordered WT samples.
# mat_bc should contain batch-corrected, non-scaled VST values
# generated in the upstream PCA/DEG analysis script.
expr_wt_bc <- mat_bc[, wt_cols, drop = FALSE]

# Row-scale genes so that module scores reflect temporal expression
# patterns rather than absolute expression magnitude.
Z_wt <- t(scale(t(expr_wt_bc)))
Z_wt[is.na(Z_wt)] <- 0


# Split WT leptin-responsive genes into k-means expression clusters.
cluster_sets <- split(clu_df$GeneID, clu_df$Cluster)


# Calculate one per-sample module score for each expression cluster.
# Each score is the mean z-scaled expression of genes in that cluster.
module_scores <- sapply(cluster_sets, function(gs) {
  gs <- intersect(gs, rownames(Z_wt))
  
  # Ignore clusters with too few genes represented in the expression matrix.
  if (length(gs) < 5) return(rep(NA_real_, ncol(Z_wt)))
  
  colMeans(Z_wt[gs, , drop = FALSE])
})


# Rows are WT samples; columns are expression modules.
MS <- as.data.frame(module_scores)
rownames(MS) <- colnames(Z_wt)


# ============================================================
# 5. Align TF/pathway activities to WT module scores
# ============================================================

# Ensure that TF and pathway activity matrices are aligned to the
# same WT samples used for module scoring.
stopifnot(all(rownames(MS) %in% colnames(M_tf)))
stopifnot(all(rownames(MS) %in% colnames(M_pw)))

TF_wt <- t(M_tf[, rownames(MS), drop = FALSE])  # samples x TFs
PW_wt <- t(M_pw[, rownames(MS), drop = FALSE])  # samples x pathways


# ============================================================
# 6. Correlate module scores with TF/pathway activities
# ============================================================

# For each expression module, calculate Spearman correlation with
# every inferred TF activity across WT samples.
cor_mat_tf <- sapply(MS, function(ms) {
  apply(TF_wt, 2, function(x) {
    cor(ms, x, method = "spearman", use = "pair")
  })
})


# For each expression module, calculate Spearman correlation with
# every inferred pathway activity across WT samples.
cor_mat_pw <- sapply(MS, function(ms) {
  apply(PW_wt, 2, function(x) {
    cor(ms, x, method = "spearman", use = "pair")
  })
})


# Keep only TFs/pathways selected in the upstream activity-dynamics analysis.
cor_mat_tf_show <- cor_mat_tf[
  rownames(cor_mat_tf) %in% features_filtered_tf,
  ,
  drop = FALSE
]

cor_mat_pw_show <- cor_mat_pw[
  rownames(cor_mat_pw) %in% features_filtered_pw,
  ,
  drop = FALSE
]


# ============================================================
# 7. Order rows and columns for heatmap display
# ============================================================

# Keep expression modules in their original cluster order.
ord_mod <- colnames(cor_mat_tf_show)

# Rank TFs/pathways by strongest absolute module correlation.
ord_tf <- order(apply(abs(cor_mat_tf_show), 1, max), decreasing = TRUE)
ord_pw <- order(apply(abs(cor_mat_pw_show), 1, max), decreasing = TRUE)


# ============================================================
# 8. Plot pathway activity correlations with expression modules
# ============================================================

Heatmap(
  cor_mat_pw_show[ord_pw, ord_mod, drop = FALSE],
  name = "ρ",
  col = colorRamp2(c(-1, 0, 1), c("royalblue3", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Pathway activities correlated with WT expression modules"
)


# ============================================================
# 9. Plot TF activity correlations with expression modules
# ============================================================

# Use TF heatmap clusters as row splits for easier biological
# interpretation of the TF correlation heatmap.
tf_split_map <- setNames(clu_df_tf71$Cluster, clu_df_tf71$GeneID)

row_split_tf <- factor(
  tf_split_map[rownames(cor_mat_tf_show)],
  levels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
)


Heatmap(
  cor_mat_tf_show[ord_tf, ord_mod, drop = FALSE],
  name = "corr",
  col = colorRamp2(c(-1, 0, 1), c("royalblue3", "white", "firebrick3")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = row_split_tf[ord_tf],
  row_names_gp = gpar(fontsize = 10),
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title = "TF activities correlated with WT expression modules"
)
