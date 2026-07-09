
############################################################
# Script: 2_Tfs_Pathways_activity_dynamics.R
# Project: Temporal organization of leptin signaling reveals distinct
#          regulatory roles for STAT3 and ETV6
#
# Purpose:
#   This script estimates transcription factor (TF) and pathway activity
#   dynamics from VST-normalized RNA-seq expression data.
#
# Main analyses:
#   1. Prepare VST expression matrix and sample metadata from a DESeq2 object.
#   2. Estimate pathway activity using PROGENy through decoupleR.
#   3. Estimate TF activity using mouse DoRothEA regulons through decoupleR.
#   4. Model WT leptin time-course activity changes with limma.
#   5. Generate WT TF/pathway activity heatmaps and WT pathway line plots.
#      TF line plots are generated in the downstream script WT_71_Tf_lineplot.R.
#   6. Test genotype-by-time interaction effects for Stat3_DEL and Etv6_DEL.
#
# Required input object:
#   dds: DESeq2 object containing raw counts and sample metadata.
#
# Required metadata columns in colData(dds):
#   type    : genotype/cell type; expected values include WT, Stat3_DEL, Etv6_DEL
#   time    : leptin stimulation time; expected values include 0min, 30min, 180min
#   batch   : batch covariate used in limma models
#   samples : sample names used for plotting/ordering
#
# Important notes:
#   - Use this script after running 1_RNAseq_degs_heatmap_pathways.R
#   - Gene identifiers are assumed to be mouse gene symbols.
#   - The statistical calculations and plotting choices are kept consistent
#     with the original analysis so that the same inputs should give the
#     same outputs.
#   - Downstream scripts used with this workflow:
#       3_WT_71_Tf_lineplot.R
#       4_TF_and_PW_act_correlation.R
############################################################

# ============================================================
# Load required packages
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(decoupleR)
library(dorothea)
library(limma)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(stringr)


# ============================================================
# Helper functions
# ============================================================

#--------------------------------------------------------------
# Convert decoupleR output into a feature-by-sample matrix.
#
# Input:
#   tbl  : decoupleR result table.
#   stat : Activity statistic to extract (default = "norm_wmean").
#
# Output:
#   Matrix with features (TFs/pathways) as rows and samples as
#   columns, suitable for downstream limma analysis.
#--------------------------------------------------------------

to_mat <- function(tbl, stat = "norm_wmean") {
  tbl |>
    dplyr::filter(statistic == stat) |>
    dplyr::select(source, condition, score) |>
    tidyr::pivot_wider(names_from = condition, values_from = score) |>
    tibble::column_to_rownames("source") |>
    as.matrix()
}

#--------------------------------------------------------------
# Fit a limma model for the WT leptin time course.
#
# The model adjusts for batch effects and estimates differential
# activity between:
#   - 30 min vs 0 min
#   - 180 min vs 0 min
#   - 180 min vs 30 min
#
# Input:
#   M    : Feature activity matrix (rows = features, columns = samples).
#   meta : Sample metadata.
#
# Output:
#   Empirical Bayes-moderated limma model.
#--------------------------------------------------------------

run_wt_limma <- function(M, meta) {
  wt <- meta$type == "WT"
  
  design_wt <- model.matrix(~ 0 + time + batch, data = meta[wt, ])
  
  contr_wt <- limma::makeContrasts(
    `t30-t0`   = time30min - time0min,
    `t180-t0`  = time180min - time0min,
    `t180-t30` = time180min - time30min,
    levels = design_wt
  )
  
  fit <- limma::lmFit(M[, wt, drop = FALSE], design_wt)
  fit <- limma::contrasts.fit(fit, contr_wt)
  fit <- limma::eBayes(fit)
  
  fit
}

#--------------------------------------------------------------
# Select the most responsive features from the WT analysis.
#
# Features are ranked by limma B-statistic for two specified
# contrasts and merged into a unique feature list.
#
# Input:
#   fit     : limma model.
#   coef_a  : First contrast.
#   coef_b  : Second contrast.
#   n       : Number of top features per contrast.
#
# Output:
#   Character vector of unique feature names.
#--------------------------------------------------------------

get_top_features <- function(fit, coef_a = "t30-t0", coef_b = "t180-t0", n = 50) {
  top_a <- rownames(limma::topTable(fit, coef = coef_a, sort.by = "B", number = n))
  top_b <- rownames(limma::topTable(fit, coef = coef_b, sort.by = "B", number = n))
  
  unique(c(top_a, top_b))
}



#--------------------------------------------------------------
# Remove features with consistently low activity in WT samples.
#
# Mean activity is calculated separately for each WT time point.
# Features with mean activity below the specified cutoff at every
# time point are excluded from downstream visualisation.
#
# Typical cutoffs:
#   TF activity      : 1
#   Pathway activity : 0
#
# Input:
#   M        : Feature activity matrix.
#   meta     : Sample metadata.
#   features : Candidate feature list.
#   cutoff   : Activity threshold.
#
# Output:
#   Filtered feature vector.
#--------------------------------------------------------------


filter_low_activity <- function(M, meta, features, cutoff) {
  wt_flag <- meta$type == "WT"
  
  wt_means <- M[, wt_flag, drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column("feature") |>
    tidyr::pivot_longer(-feature, names_to = "sample", values_to = "score") |>
    dplyr::mutate(time = meta[sample, "time"]) |>
    dplyr::group_by(feature, time) |>
    dplyr::summarise(mean = mean(score, na.rm = TRUE), .groups = "drop")
  
  low_features <- wt_means |>
    dplyr::group_by(feature) |>
    dplyr::summarise(all_low = all(mean < cutoff, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(all_low) |>
    dplyr::pull(feature)
  
  setdiff(features, low_features)
}

#--------------------------------------------------------------
# Plot WT activity dynamics for selected features.
#
# Generates a time-course line plot of mean feature activity
# across WT samples. Feature labels are placed at the final
# time point using ggrepel, and optional custom colours can be
# supplied for individual features.
#
# Input:
#   M            : Feature activity matrix.
#   features     : Features to plot.
#   title        : Plot title.
#   meta         : Sample metadata.
#   wt_flag      : Logical vector indicating WT samples.
#   time_levels  : Ordered time points.
#   colors       : Optional named vector of colours.
#   show_legend  : Display legend if TRUE.
#
# Output:
#   ggplot object.
#--------------------------------------------------------------

plot_time_lines <- function(M, features, title, meta, wt_flag,
                            time_levels = c("0min","30min","180min"),
                            colors = NULL, show_legend = FALSE) {
  
  df <- M[features, wt_flag, drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column("feature") |>
    tidyr::pivot_longer(-feature, names_to = "sample", values_to = "score") |>
    dplyr::mutate(time = factor(meta[sample, "time"], levels = time_levels, ordered = TRUE)) |>
    dplyr::group_by(feature, time) |>
    dplyr::summarise(mean = mean(score), .groups = "drop")
  
  last_t <- tail(time_levels, 1)
  lab_df <- df |> dplyr::filter(time == last_t) |> dplyr::mutate(label = feature)
  
  p <- ggplot(df, aes(x = time, y = mean, group = feature, color = feature)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = label),
      nudge_x = 0.25, direction = "y", hjust = 0,
      segment.color = NA, size = 3, show.legend = FALSE
    ) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(5.5, 40, 5.5, 5.5)
    ) +
    labs(x = "Time (min)", y = "Activity (norm_wmean)", title = title)
  
  if (!is.null(colors)) {
    # keep only colors for features present
    use_cols <- colors[names(colors) %in% unique(df$feature)]
    p <- p + scale_color_manual(values = use_cols)
  }
  p
}
#--------------------------------------------------------------
# Fit genotype × time interaction models.
#
# The model tests whether temporal activity responses differ
# between WT and a specified knockout genotype while adjusting
# for batch effects.
#
# Interaction coefficients represent genotype-specific changes
# relative to WT at:
#   - 30 min
#   - 180 min
#
# Input:
#   M        : Feature activity matrix.
#   meta     : Sample metadata.
#   genotype : Knockout genotype name (must match meta$type).
#   prefix   : Prefix used to name output contrasts.
#
# Output:
#   Empirical Bayes-moderated interaction model.
#--------------------------------------------------------------

run_interaction_limma <- function(M, meta, genotype, prefix) {
  
  meta$type <- relevel(factor(meta$type), "WT")
  meta$time <- relevel(factor(meta$time), "0min")
  
  design_all <- model.matrix(~ batch + type * time, data = meta)
  colnames(design_all) <- make.names(colnames(design_all))
  
  fit_all <- limma::lmFit(M, design_all)
  
  contrast_30  <- paste0("type", genotype, ".time30min")
  contrast_180 <- paste0("type", genotype, ".time180min")
  
  stopifnot(contrast_30 %in% colnames(design_all))
  stopifnot(contrast_180 %in% colnames(design_all))
  
  cont_int <- limma::makeContrasts(
    contrasts = c(
      contrast_30,
      contrast_180
    ),
    levels = design_all
  )
  
  colnames(cont_int) <- c(
    paste0(prefix, "30"),
    paste0(prefix, "180")
  )
  
  fit_int <- limma::contrasts.fit(fit_all, cont_int)
  fit_int <- limma::eBayes(fit_int)
  
  fit_int
}

#--------------------------------------------------------------
# Extract interaction-effect estimates for plotting.
#
# Retrieves limma coefficients corresponding to the 30- and
# 180-minute interaction contrasts for selected features and
# reshapes them into long format.
#
# Input:
#   fit      : limma interaction model.
#   feats    : Features to extract.
#   coef30   : 30-minute interaction coefficient.
#   coef180  : 180-minute interaction coefficient.
#
# Output:
#   Long-format data frame containing activity differences.
#--------------------------------------------------------------

get_delta <- function(fit, feats, coef30, coef180) {
  feats <- unique(feats)
  feats <- intersect(feats, rownames(fit$coefficients))
  
  co <- fit$coefficients[feats, c(coef30, coef180), drop = FALSE]
  
  data.frame(
    feature = rep(rownames(co), times = 2),
    time = factor(rep(c("30", "180"), each = nrow(co)), levels = c("30", "180")),
    delta = c(co[, coef30], co[, coef180])
  )
}

#--------------------------------------------------------------
# Prepare data for dumbbell plots.
#
# Combines interaction-effect estimates with statistical
# significance from limma and ranks features by the largest
# absolute interaction effect across time points.
#
# Input:
#   fit      : limma interaction model.
#   feats    : Features to display.
#   coef30   : 30-minute interaction coefficient.
#   coef180  : 180-minute interaction coefficient.
#
# Output:
#   Data frame containing:
#     - interaction effects
#     - adjusted P values
#     - raw P values
#     - feature ranking
#--------------------------------------------------------------

make_dumbbell_data <- function(fit, feats, coef30, coef180) {
  delta <- get_delta(fit, feats, coef30, coef180)
  
  db <- delta |>
    tidyr::pivot_wider(names_from = time, values_from = delta) |>
    dplyr::mutate(max_abs = pmax(abs(`30`), abs(`180`))) |>
    dplyr::arrange(dplyr::desc(max_abs))
  
  tt30 <- limma::topTable(fit, coef = coef30, number = Inf) |>
    tibble::rownames_to_column("feature")
  
  tt180 <- limma::topTable(fit, coef = coef180, number = Inf) |>
    tibble::rownames_to_column("feature")
  
  annot <- tt30 |>
    dplyr::select(feature, padj30 = adj.P.Val, pval30 = P.Value) |>
    dplyr::left_join(
      tt180 |> dplyr::select(feature, padj180 = adj.P.Val, pval180 = P.Value),
      by = "feature"
    )
  
  db |>
    dplyr::left_join(annot, by = "feature")
}

#--------------------------------------------------------------
# Plot genotype-specific interaction effects as a dumbbell plot.
#
# Each feature is represented by:
#    blue point  : 30-minute interaction effect
#    red point   : 180-minute interaction effect
#    connecting arrow indicating temporal change
#
# Features can either be displayed in a predefined order or
# ranked automatically by the largest interaction effect.
#
# Input:
#   db            : Data frame generated by make_dumbbell_data().
#   title         : Plot title.
#   feature_order : Optional custom feature ordering.
#
# Output:
#   ggplot object.
#--------------------------------------------------------------

plot_dumbbell <- function(db, title, feature_order = NULL) {
  if (!is.null(feature_order)) {
    db$feature <- factor(db$feature, levels = feature_order)
    y_aes <- aes(y = feature)
  } else {
    y_aes <- aes(y = reorder(feature, max_abs))
  }
  
  ggplot(db, y_aes) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(
      aes(x = `30`, xend = `180`, yend = feature),
      color = "grey50",
      linewidth = 0.6,
      arrow = arrow(length = unit(2.5, "mm"), type = "closed")
    ) +
    geom_point(aes(x = `30`), color = "#1f77b4", size = 2.2) +
    geom_point(aes(x = `180`), color = "#d62728", size = 2.2) +
    scale_x_continuous(
      "Δactivity (interaction)",
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(
      y = NULL,
      title = title,
      subtitle = "30 min (blue) → 180 min (red)"
    ) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank())
}
#--------------------------------------------------------------




# ============================================================
# Prepare VST expression matrix and metadata
# ============================================================

# Check that the expected DESeq2 object is available before starting.
stopifnot(exists("dds"))

vsd <- vst(dds, blind=FALSE)

expr <- assay(vsd)                           # genes x samples 
meta <- as.data.frame(colData(vsd))          #meta should have type, time, batch, samples
stopifnot(all(c("type", "time", "batch", "samples") %in% colnames(meta)))
meta$time  <- factor(meta$time, levels=c("0min","30min","180min"))
meta$type  <- factor(meta$type, levels=c("WT","Stat3_DEL","Etv6_DEL"))
meta$batch <- factor(meta$batch)
stopifnot(identical(colnames(expr), rownames(meta)))


# ============================================================
#  Estimate PROGENy pathway activity
# ============================================================

# Retrieve mouse PROGENy model. top = 500 uses the top 500 target genes per pathway.
res_progeny <- decoupleR::get_progeny(organism = "mouse", top = 500)

# save this object and it call be called again
#saveRDS(res_progeny,"commom_leptin_response/res_progeny.rds")
#res_progeny <- readRDS("commom_leptin_response/res_progeny.rds")

# estimate pathway activity
act_pw <- run_wmean(expr, res_progeny, .source="source", .target="target", .mor="weight",
                    minsize=5, times=1000, seed=1)



# ============================================================
# Estimate DoRothEA TF activity
# ============================================================

data(dorothea_mm)   # mouse regulons (tf, target, mor, confidence)

# Keep A–C (higher stringency)

reg_mouse <- dorothea_mm %>%
  filter(confidence %in% c("A","B","C")) %>%
  transmute(
    tf        = toupper(tf),      # decoupleR is fine with either; uppercase is conventional
    target    = target,           # mouse SYMBOLs
    mor       = mor,              # +1 / -1
    confidence = confidence
  )

# saveRDS(reg_mouse,"commom_leptin_response/reg_mouse.rds")
# reg_mouse <- readRDS("commom_leptin_response/reg_mouse.rds")

# estimate TF activities
act_tf <- run_wmean(expr, reg_mouse, .source="tf", .target="target", .mor="mor",
                    minsize=5, times=1000, seed=1)


#Convert decoupleR output into a feature-by-sample matrix.
M_pw <- to_mat(act_pw, "norm_wmean")
M_tf <- to_mat(act_tf, "norm_wmean")

# align with meta
M_pw <- M_pw[, rownames(meta)]
M_tf <- M_tf[, rownames(meta)]




# ============================================================
# WT TF activity dynamics
# ============================================================

fit_tf_wt <- run_wt_limma(M_tf, meta)

tf_wt_inf_30min <- limma::topTable(fit_tf_wt, coef = "t30-t0", sort.by = "B", number = Inf)
tf_wt_inf_180min <- limma::topTable(fit_tf_wt, coef = "t180-t0", sort.by = "B", number = Inf)

#write.csv(tf_wt_inf_30min, "transcriptomics_data_analysis/tf_wt_inf_30min.csv")
# write.csv(tf_wt_inf_180min, "transcriptomics_data_analysis/tf_wt_inf_180min.csv")


tf_top50 <- get_top_features(fit_tf_wt, n = 50)

features_filtered_tf <- filter_low_activity(
  M = M_tf,
  meta = meta,
  features = tf_top50,
  cutoff = 1
)

length(features_filtered_tf)
# Generate heatmap of selected WT TF activity dynamics.

wt_cols <- meta %>%
  mutate(time = factor(time, levels = c("0min","30min","180min"))) %>%
  filter(type == "WT") %>%
  arrange(time) %>%
  pull(samples)

heatmap_wt_tf <- M_tf[, wt_cols, drop = FALSE]
heatmap_wt_tf <- heatmap_wt_tf[features_filtered_tf,]
heatmap_wt_tf <- t(scale(t(heatmap_wt_tf)))


meta_ord <- meta[wt_cols,]



colAnn <- ComplexHeatmap::HeatmapAnnotation(
  time = meta_ord$time,
  #Exp       = meta_ord$batch,
  col = list(
    time = c("0min"="#A4A4AE", "30min"="#6E6E78", "180min"="#46464E")
    #batch      = c(I="purple", II="magenta")
  ),
  annotation_legend_param = list(
    time = list(
      nrow = 1,                  
      title = "time"
      # optionally: labels_gp = gpar(fontsize = 10)
    )
  )
)


set.seed(123)
tf_71_heatmap <- Heatmap(heatmap_wt_tf,
                         cluster_columns = F,
                         column_title="TF activity dynamics in WT samples",
                         col = circlize::colorRamp2(c(-2,0,2), c("royalblue3","white", "firebrick3")),####
                         top_annotation = colAnn,
                         heatmap_legend_param = list(title = "scaled\nactivity",legend_direction = "horizontal"),
                         row_names_gp = gpar(fontsize = 10),show_column_names = F,
                         column_names_gp = gpar(fontsize = 8),km=4,
                         width = ncol(as.matrix(heatmap_wt_tf))*unit(2.5, "mm"), 
                         height = nrow(as.matrix(heatmap_wt_tf))*unit(3, "mm"))

tf_71_heatmap <- draw(tf_71_heatmap,
                      heatmap_legend_side = "bottom",
                      annotation_legend_side = "bottom")


tf_71_heatmap


clusterlist_tf = row_order(tf_71_heatmap)

clu_df_tf71 <- lapply(names(clusterlist_tf), function(i){
  out <- data.frame(GeneID = rownames(heatmap_wt_tf)[clusterlist_tf[[i]]],
                    Cluster = paste0("cluster_", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  do.call(rbind, .)

#NOTE: run 3_WT_71_Tf_lineplot.R for Tf lineplot

# ============================================================
# WT pathway activity dynamics
# ============================================================

fit_pw_wt <- run_wt_limma(M_pw, meta)


pw_wt_inf_30min <- limma::topTable(fit_pw_wt, coef = "t30-t0", sort.by = "B", number = Inf)
pw_wt_inf_180min <- limma::topTable(fit_pw_wt, coef = "t180-t0", sort.by = "B", number = Inf)

# write.csv(pw_wt_inf_30min, "transcriptomics_data_analysis/pw_wt_inf_30min.csv") 
# write.csv(pw_wt_inf_180min, "transcriptomics_data_analysis/pw_wt_inf_180min.csv")


pw_top <- get_top_features(fit_pw_wt, n = Inf)

features_filtered_pw <- filter_low_activity(
  M = M_pw,
  meta = meta,
  features = pw_top,
  cutoff = 0
)

#heatmap for pathways
heatmap_wt_pw <- M_pw[, wt_cols, drop = FALSE]
heatmap_wt_pw <- heatmap_wt_pw[features_filtered_pw,]
heatmap_wt_pw <- t(scale(t(heatmap_wt_pw)))

draw(Heatmap(heatmap_wt_pw,
             cluster_columns = F,
             column_title="Pathway activity dynamics in WT samples",
             col = circlize::colorRamp2(c(-2,0,2), c("royalblue3","white", "firebrick3")),####
             top_annotation = colAnn,
             heatmap_legend_param = list(title = "scaled\nactivity",legend_direction = "horizontal"),
             row_names_gp = gpar(fontsize = 10),show_column_names = F,
             column_names_gp = gpar(fontsize = 8)),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")


pal_pw <- c(
  "JAK-STAT" = "#1f7a8c",  # STAT3 teal (matches STATs block)
  "NFkB"     = "#d62728",  # NFKB1 red (same as NF-κB)
  "TNFa"     = "#e6550d",  # RELA orange (canonical TNF→NF-κB axis)
  "MAPK"     = "#fdae6b",  # ELK4/AP-1 orange (MAPK→AP-1)
  "TGFb"     = "#6a51a3",  # AHR purple (keeps signaling-family purples)
  "Hypoxia"  = "#6b6ecf",  # E2F7 slate (cool, distinct from STAT teals)
  "Estrogen" = "#CC79A7",  # classical magenta for nuclear receptors
  "Androgen" = "#a05a2c",  # CEBPB brown (steroid-hormone family tone)
  "VEGF"     = "#2ca9bc",  # STAT1 teal (vascular/PI3K–MAPK crosstalk; cool tone)
  "WNT"      = "#7F7F7F"   # neutral grey (morphogen set apart)
)


# Ensure we only plot pathways present in the pathway activity matrix.
features_filtered_pw <- intersect(features_filtered_pw , rownames(M_pw))
wt_flag <- meta$type == "WT"

#line plot for pathways
p <- plot_time_lines(
  M = M_pw,
  features = features_filtered_pw,
  title = "Pathway activity (WT)",
  meta = meta,
  wt_flag = wt_flag,
  colors = pal_pw,            # <-- apply your palette
  show_legend = F         # show legend if you prefer
)

p



#NOTE:  run "4_TF_and_PW_act_correlation.R" for TF and pathway activity correlation(Fig 3G,J) with WT time course degs clusters(of Fig2C)


# ============================================================
# Stat3_DEL interaction effect
# ============================================================

fit_tf_stat3_int <- run_interaction_limma(
  M = M_tf,
  meta = meta,
  genotype = "Stat3_DEL",
  prefix = "S"
)

tt30_tf_stat3 <- limma::topTable(fit_tf_stat3_int, coef = "S30", number = Inf) |>
  tibble::rownames_to_column("feature")

tt180_tf_stat3 <- limma::topTable(fit_tf_stat3_int, coef = "S180", number = Inf) |>
  tibble::rownames_to_column("feature")


# write.csv(tt30_tf_stat3, "transcriptomics_data_analysis/tt30_tf_stat3.csv")
# write.csv(tt180_tf_stat3, "transcriptomics_data_analysis/tt180_tf_stat3.csv")


fit_pw_stat3_int <- run_interaction_limma(
  M = M_pw,
  meta = meta,
  genotype = "Stat3_DEL",
  prefix = "S"
)
tt30_pw_stat3 <- limma::topTable(fit_pw_stat3_int, coef = "S30", number = Inf) |>
  tibble::rownames_to_column("feature")

tt180_pw_stat3 <- limma::topTable(fit_pw_stat3_int, coef = "S180", number = Inf) |>
  tibble::rownames_to_column("feature")

# write.csv(tt30_pw_stat3 , "transcriptomics_data_analysis/tt30_pw_stat3 .csv")
# write.csv(tt180_pw_stat3, "transcriptomics_data_analysis/tt180_pw_stat3.csv")



tt30_sig_tf_stat3 <- tt30_tf_stat3 |> dplyr::filter(adj.P.Val < 0.05)
tt180_sig_tf_stat3 <- tt180_tf_stat3 |> dplyr::filter(adj.P.Val < 0.05)

tfs_stat3del_dumbbellplot <- unique(c(
  intersect(clu_df_tf71$GeneID, tt30_sig_tf_stat3$feature),
  intersect(clu_df_tf71$GeneID, tt180_sig_tf_stat3$feature)
))

db_tf_stat3 <- make_dumbbell_data(
  fit = fit_tf_stat3_int,
  feats = tfs_stat3del_dumbbellplot,
  coef30 = "S30",
  coef180 = "S180"
)

p_tf_stat3 <- plot_dumbbell(
  db = db_tf_stat3,
  title = "Stat3Δ impairment of TF activity (DEL − WT)",
  feature_order = c(
    "PKNOX1","MAFG","KLF6","MEF2A","FOXK2","ZFP217","E2F7","EGR1","CREM","CREB1",
    "RARA","PPARG","ESR1","ELK4","FLI1","ETV4","ETS1","CEBPD","CEBPB","CEBPA",
    "EPAS1","ARNT","AHR","JUNB","JUN","FOSL2","ATF4","RELA","REL","NFKB2","NFKB1",
    "STAT6","STAT5B","STAT5A","STAT4","STAT3","STAT1"
  )
)

p_tf_stat3


# Stat3_DEL pathway interaction

features_pw_stat3 <- c(
  "NFkB", "TNFa", "Androgen", "MAPK", "JAK-STAT",
  "Estrogen", "WNT"
)

db_pw_stat3 <- make_dumbbell_data(
  fit = fit_pw_stat3_int,
  feats = features_pw_stat3,
  coef30 = "S30",
  coef180 = "S180"
)

p_pw_stat3 <- plot_dumbbell(
  db = db_pw_stat3,
  title = "Stat3Δ impairment of pathway activity (DEL − WT)"
)

p_pw_stat3

# ============================================================
# Etv6_DEL interaction effect
# ============================================================

fit_tf_etv6_int <- run_interaction_limma(
  M = M_tf,
  meta = meta,
  genotype = "Etv6_DEL",
  prefix = "E"
)

tt30_tf_etv6 <- limma::topTable(fit_tf_etv6_int, coef = "E30", number = Inf) |>
  tibble::rownames_to_column("feature")

tt180_tf_etv6 <- limma::topTable(fit_tf_etv6_int, coef = "E180", number = Inf) |>
  tibble::rownames_to_column("feature")
# write.csv(tt30_tf_etv6, "transcriptomics_data_analysis/tt30_tf_etv6.csv")
# write.csv(tt180_tf_etv6, "transcriptomics_data_analysis/tt180_tf_etv6.csv")


fit_pw_etv6_int <- run_interaction_limma(
  M = M_pw,
  meta = meta,
  genotype = "Etv6_DEL",
  prefix = "E"
)

tt30_pw_etv6 <- limma::topTable(fit_pw_etv6_int, coef = "E30", number = Inf) |>
  tibble::rownames_to_column("feature")

tt180_pw_etv6 <- limma::topTable(fit_pw_etv6_int, coef = "E180", number = Inf) |>
  tibble::rownames_to_column("feature")
# write.csv(tt30_pw_etv6, "transcriptomics_data_analysis/tt30_pw_etv6.csv")
# write.csv(tt180_pw_etv6, "transcriptomics_data_analysis/tt180_pw_etv6.csv")



feature_tf_etv6 <- c(
  "BATF", "E2F7", "ELK4", "CEBPB", "AHR", "FOXL2", "FOXK2",
  "RELA", "NFKB2", "NFKB1", "STAT6", "STAT5B", "STAT5A",
  "STAT4", "STAT3", "STAT1"
)

db_tf_etv6 <- make_dumbbell_data(
  fit = fit_tf_etv6_int,
  feats = feature_tf_etv6,
  coef30 = "E30",
  coef180 = "E180"
)

p_tf_etv6 <- plot_dumbbell(
  db = db_tf_etv6,
  title = "Etv6Δ impairment of TF activity (DEL − WT)",
  feature_order = feature_tf_etv6
)

p_tf_etv6

# Etv6_DEL pathway interaction

db_pw_etv6 <- make_dumbbell_data(
  fit = fit_pw_etv6_int,
  feats = features_filtered_pw,
  coef30 = "E30",
  coef180 = "E180"
)

p_pw_etv6 <- plot_dumbbell(
  db = db_pw_etv6,
  title = "Etv6Δ impairment of pathway activity (DEL − WT)"
)

p_pw_etv6
