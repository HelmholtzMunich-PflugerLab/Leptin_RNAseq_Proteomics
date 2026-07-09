#****************************************************************************
# Script: 1_RNAseq_degs_heatmap_pathways.R
# Project: Leptin stimulation response in WT, Stat3_DEL, and Etv6_DEL cells
#
# Description:
# This script performs RNA-seq differential expression analysis for a leptin
# stimulation time-course experiment comparing WT, Stat3_DEL, and Etv6_DEL
# samples across 0 min, 30 min, and 180 min time points.
#
# The workflow includes:
#   1. Loading count matrix and sample metadata.
#   2. Creating a DESeq2 dataset with batch, genotype, time, and interaction terms.
#   3. Performing size-factor normalization and variance-stabilizing transformation.
#   4. Assessing sample-level structure using PCA plots.
#   5. Extracting differentially expressed genes for:
#        - WT time effects,
#        - Stat3_DEL and Etv6_DEL time effects,
#        - genotype-vs-WT contrasts at each time point,
#        - genotype-specific interaction effects relative to WT.
#   6. Generating a WT leptin-response heatmap using significant time-responsive genes.
#   7. Applying k-means clustering to WT time-responsive genes.
#   8. Calculating and plotting cluster centroid trajectories.
#   9. Performing pathway enrichment analysis for WT clusters and genotype-vs-WT DEGs.
#  10. Plotting selected genes of interest across conditions.
#
# Input files:
#   - Data/Counts_WT_Stat3del_Etv6del.txt
#   - Data/colData_WT_Stat3del_Etv6del.txt
#   - R/deg_analysis_functions.R
#   - R/enrichment_helpers.R
#
# Main outputs:
#   - Normalized count matrix
#   - PCA plots
#   - DEG tables and volcano plots
#   - WT time-course heatmap with k-means clusters
#   - Cluster membership table
#   - Cluster centroid plot
#   - Enrichment result tables and enrichment plots
#   - Selected gene expression plots
#
# Notes:
#   - Differential expression testing is performed using the original count data.
#   - Batch-corrected VST values are used only for visualization and clustering.
#   - The DESeq2 model includes batch and genotype-by-time interaction effects.
#****************************************************************************



# ============================================================
# Load required packages
# ============================================================
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(vsn)
library(genefilter)
library(ggrepel)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)


# ============================================================
# load raw data and run DEseq2
# ============================================================

##Read the raw count matrix
counts <- read.table("Data/Counts_WT_Stat3del_Etv6del.txt", header=T)

##Read sample metadata.
colData <- read.table("Data/colData_WT_Stat3del_Etv6del.txt", header=T)


##Convert metadata variables to factors
colData$batch <- factor(colData$batch)

##Set WT as the reference genotype
colData$type <- factor(colData$type, levels = c("WT","Stat3_DEL","Etv6_DEL"))

##Set 0min as the reference timepoint
colData$time <- factor(colData$time, levels = c("0min","30min","180min"))


## matches count columns to colData rownames, so sample IDs are assigned as rownames here.
rownames(colData)<- colData$samples
rownames(colData)
## check that count-matrix columns and colData rownames are aligned before constructing dds.
colnames(counts) == rownames(colData)


#Build the DESeq2 object using a design with batch, genotype, time, and genotype-by-time interaction.

dds = DESeqDataSetFromMatrix(countData=counts,
                             colData=colData, 
                             design = ~batch + type*time)

## Filter out very low-count genes before differential expression testing.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
## Estimate library-size normalization factors.
dds<- estimateSizeFactors(dds)
sizeFactors(dds)

##Set the reference to be compared
dds$time = relevel(dds$time,"0min")
dds$time


##Fit the negative-binomial GLM and estimate dispersions and coefficients.
dds = DESeq(dds)


## Extract normalized counts for export/inspection; these are not VST-transformed values.
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)


## Save the normalized count matrix to the output folder.
#dir.create("transcriptomics_data_analysis") # create folder
#write.table(dt, file="transcriptomics_data_analysis/norm_matrix.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)
#------------------------------------------------------------#


# ============================================================
#Data quality assessment using PCA 
# ============================================================


## Variance-stabilizing transformation for PCA, heatmaps, and expression visualization.
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)


dds$time <- factor(dds$time, levels = c("0min","30min","180min"))

## PCA using all samples to inspect major structure across time, genotype, and batch.
pcaData <- plotPCA(vsd, intgroup=c("time","type","batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(PC1, PC2, color=time, shape=type, label=batch)) +
  geom_point(size=2) +
  geom_text_repel(size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()+ labs(title = "PCA plot")


##PCA for WT  samples

##Subset metadata to WT samples for WT-only PCA.
meta0 <- as.data.frame(colData(vsd)) %>% filter(type == "WT")
meta0$time  <- factor(meta0$time, levels = c("0min","30min","180min"))
meta0$batch <- factor(meta0$batch)

vsd2 <- vsd[, rownames(meta0)]
pcaData <- plotPCA(vsd2, intgroup = "time", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(PC1, PC2, color = time, shape = time)) +
  geom_point(size = 4, color = "#46464E") +
  scale_shape_manual(values = c(`0min`=16, `30min`=15, `180min`=17)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  labs(title = "WT: leptin stimulation (time effect)")

## PCA for Stat3del and WT  samples
## batch I only

meta0 <- as.data.frame(colData(vsd))
## Restrict to batch I for the Stat3_DEL versus WT PCA comparison.
meta0 <- meta0 %>%  dplyr::filter(batch == "I")
meta0$time  <- factor(meta0$time,  levels=c("0min","30min","180min"))
meta0$type  <- factor(meta0$type,  levels=c("WT","Stat3_DEL"))   # drop Etv6
meta0$batch <- factor(meta0$batch)

keep <- meta0$type %in% c("WT","Stat3_DEL")
meta  <- droplevels(meta0[keep, ])
vsd2  <- vsd[, rownames(meta)]

pcaData <- plotPCA(vsd2, intgroup=c("time","type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=type, shape=time)) +
  geom_point(size=6) +
  scale_shape_manual(values = c('0min'= 16, '30min' = 15, '180min'=17)) +
  scale_color_manual(values = c(WT = "#46464E", Stat3_DEL = "#1F7A8C"))+
   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()+ 
  theme(title = element_text(size=18))+
  labs(title = "Stat3 DEL: PCA showing deletion and\nLeptin stimulation effect")

## PCA for Etv6del and WT  samples
## batch II only

meta0 <- as.data.frame(colData(vsd))
##  Restrict to batch II for the Etv6_DEL versus WT PCA comparison.
meta0 <- meta0 %>%  dplyr::filter(batch == "II")
meta0$time  <- factor(meta0$time,  levels=c("0min","30min","180min"))
meta0$type  <- factor(meta0$type,  levels=c("WT","Etv6_DEL"))   # drop Etv6
meta0$batch <- factor(meta0$batch)

keep <- meta0$type %in% c("WT","Etv6_DEL")
meta  <- droplevels(meta0[keep, ])
vsd2  <- vsd[, rownames(meta)]

pcaData <- plotPCA(vsd2, intgroup=c("time","type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=type, shape=time)) +
  geom_point(size=6) +
  scale_shape_manual(values = c('0min'= 16, '30min' = 15, '180min'=17)) +
  scale_color_manual(values = c(WT = "#46464E", Etv6_DEL = "#6B58B5"))+
  #geom_text_repel(size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()+ 
  theme(title = element_text(size=18))+
  labs(title = "Etv6 DEL: PCA showing deletion and\nLeptin stimulation effect")

#------------------------------------------------------------#


# ============================================================
#*Extracting DEGs
# ============================================================


## Load custom DEG helper functions used to generate contrast specifications,
#tables, xlsx files, and volcano plots.
source("R/helpers/deg_analysis_functions.R")

#--------WT leptin-response contrasts across time.-----------#
wt_time <- run_deg_set(
  dds = dds,
  comparisons = make_wt_time_specs(),
  outdir = "transcriptomics_data_analysis/wt_time_effect",
  xlsx_file = "degs_wt_time_effect.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE,
  custom_volcano_colours = FALSE
)


#-----------Within-Stat3_DEL leptin-response contrasts across time.-----#
stat3_time <- run_deg_set(
  dds = dds,
  comparisons = make_ko_time_specs("Stat3DEL", "typeStat3_DEL"),
  outdir = "transcriptomics_data_analysis/stat3_del_time_effect",
  xlsx_file = "degs_stat3del_time_effect.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)


#------------Within-Etv6_DEL leptin-response contrasts across time.---#

etv6_time <- run_deg_set(
  dds = dds,
  comparisons = make_ko_time_specs("Etv6DEL", "typeEtv6_DEL"),
  outdir = "transcriptomics_data_analysis/Etv6_del_time_effect",
  xlsx_file = "degs_etv6del_time_effect.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)



#---------Stat3_DEL interaction contrasts test whether the time response differs from WT.---------#

stat3_interaction <- run_deg_set(
  dds = dds,
  comparisons = make_interaction_specs("Stat3DEL", "typeStat3_DEL"),
  outdir = "transcriptomics_data_analysis/stat3_interaction_effect",
  xlsx_file = "degs_stat3_interaction_effect.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)


#---------Etv6_DEL interaction contrasts test whether the time response differs from WT.---------#

etv6_interaction <- run_deg_set(
  dds = dds,
  comparisons = make_interaction_specs("Etv6DEL", "typeEtv6_DEL"),
  outdir = "transcriptomics_data_analysis/Etv6_interaction_effect",
  xlsx_file = "degs_etv6_interaction_effect.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)

# Stat3DEL vs WT at each time point
#---------Stat3_DEL versus WT contrasts at each matched timepoint.---------------#

stat3_vs_wt <- run_deg_set(
  dds = dds,
  comparisons = make_vs_wt_specs("Stat3DEL", "type_Stat3_DEL_vs_WT", "typeStat3_DEL"),
  outdir = "transcriptomics_data_analysis/Stat3_vsWT",
  xlsx_file = "degs_stat3_vsWT.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)


#---------  Etv6_DEL versus WT contrasts at each matched timepoint. ---------#

etv6_vs_wt <- run_deg_set(
  dds = dds,
  comparisons = make_vs_wt_specs("Etv6DEL", "type_Etv6_DEL_vs_WT", "typeEtv6_DEL"),
  outdir = "transcriptomics_data_analysis/Etv6_vsWT",
  xlsx_file = "degs_etv6_vsWT.xlsx",
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  make_volcano = TRUE
)
#------------------------------------------------------------#



# ============================================================
# Heatmap with K-means clustering for WT time course
# ============================================================

## Collect WT time-effect DEG sets that will define genes of interest for the heatmap.
Up_wt_30 <- wt_time$tables$WT_30v0$up
Up_wt_180 <- wt_time$tables$WT_180v0$up
Up_wt_180v30 <- wt_time$tables$WT_180v30$up

Down_wt_30 <- wt_time$tables$WT_30v0$down
Down_wt_180 <- wt_time$tables$WT_180v0$down
Down_wt_180v30 <- wt_time$tables$WT_180v30$down

## Combine up- and down-regulated WT time-effect genes across all WT contrasts.
goi <- c(Up_wt_30$id,Up_wt_180$id, Up_wt_180v30$id,
         Down_wt_30$id,Down_wt_180$id,Down_wt_180v30$id)

length(goi)
#[1]1557

## Remove duplicated genes 
goi <- unique(goi)
length(goi)
# 990



library(ComplexHeatmap)
library(circlize)

## Extract VST expression values for the selected WT time-response genes.
mat0 <- assay(vsd)[intersect(goi, rownames(vsd)), , drop = FALSE]

## Keep WT sample metadata for the WT time-course heatmap.
meta <- as.data.frame(colData(vsd)) %>% filter(type %in% "WT")
meta$time  <- factor(meta$time,  levels = c("0min", "30min", "180min"))
meta$batch <- factor(meta$batch)


## Subset expression matrix columns to WT samples.
mat0 <- mat0[,colnames(mat0) %in%rownames(meta)]

# Ensure identical column order across everything
stopifnot(identical(colnames(mat0), rownames(meta)))

# 1) Batch-correction FOR VISUALIZATION only

design_vis <- model.matrix(~ time, data = meta)

## Batch-correct VST values for heatmap, not for DESeq2 testing.
mat_bc <- limma::removeBatchEffect(mat0, batch = meta$batch, design = design_vis)

# 2) Row z-score (shape-based clustering)

mat_z <- t(scale(t(mat_bc)))


# 3) Order columns only for display 

ord_cols <- order(meta$time, meta$batch)
mat_z    <- mat_z[, ord_cols]
meta_ord <- meta[ord_cols,]

# 4) Heatmap

## Build top annotation showing timepoint colours above the heatmap.
colAnn <- ComplexHeatmap::HeatmapAnnotation(
  time = meta_ord$time,
  col = list(
    time = c("0min"="#A4A4AE", "30min"="#6E6E78", "180min"="#46464E")
    ),
  
  annotation_legend_param = list(
    time = list(
      nrow = 1,                  
      title = "time"
      
    )
  )
)

## Make k-means cluster assignment reproducible.
set.seed(41)
## Draw the WT time-course heatmap with k-means row splitting.
ht_km1<- Heatmap(
  mat_z,
  name = "Exp",
  col = circlize::colorRamp2(c(-3,0,3.5), c("royalblue3","white", "firebrick3")),
  top_annotation = colAnn,
  show_row_names = FALSE,
  show_column_names = F,
  cluster_rows = F,          # we use km
  cluster_columns = F,       # fixed time order
  km=5,
  column_split = meta_ord$time,
  column_title = "Leptin stimulation time efffect",
  heatmap_legend_param = list(
    title = "Exp",
    legend_direction = "horizontal"
  )
)

ht_km1= draw(ht_km1,heatmap_legend_side = "bottom",
             annotation_legend_side = "bottom")

ht_km1


##  Extract row order for each heatmap k-means cluster after drawing the heatmap.
clusterlist = row_order(ht_km1)

## Convert heatmap cluster membership into a gene-to-cluster dataframe.
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(mat_z)[clusterlist[[i]]],
                    Cluster = paste0("cluster_", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  do.call(rbind, .)




#--------------Plotting centriod: cluster mean expression plot-----------------#


## Per-gene mean at each timepoint (use batch-corrected, non scaled values for means)

mat_bc_ord <- mat_bc[, ord_cols]

## Calculate each gene's mean expression at each timepoint.
tmeans <- t(apply(mat_bc_ord, 1, function(v) tapply(v, meta_ord$time, mean)))
tmeans <- tmeans[, c("0min","30min","180min")] 


## Average timepoint means across genes within each k-means cluster to obtain cluster centroids.
centroids <- tmeans |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneID") |>
  inner_join(clu_df, by = "GeneID") |>
  group_by(Cluster) |>
  summarise(`0min` = mean(`0min`),
            `30min`= mean(`30min`),
            `180min`= mean(`180min`), .groups="drop")


# compute changes relative to baseline (WT_0) and gap between earlyup and late up

head(centroids)

## Compute centroid changes relative to 0min and the 30min-to-180min gap.
centroids <- centroids %>%
  mutate(
    rel30  = `30min` - `0min`,
    rel180 = `180min` - `0min`,
    gap    = `30min` - `180min`
  )

# ---- Label clusters and set thresholds ----
up_th   <- 0.15   # "meaningfully up" (VST units)
down_th <- -0.15  # "meaningfully down"
gap_th  <- 0.25   # separates early vs late vs sustained

## Assign a biologically interpretable response label based on centroid trajectory thresholds.
label_shape <- function(rel30, rel180, gap, up_th, down_th, gap_th){
  if (rel30 >= up_th & rel180 >= up_th & abs(gap) <= gap_th) return("SustainedUp")
  if (rel30 >= up_th & rel180 >= up_th & gap >  gap_th)      return("EarlyUp")
  if (rel30 >= up_th & rel180 <  up_th & gap >  gap_th/2)    return("UpThenDown")
  if (rel180 >= up_th & (rel180 - rel30) > gap_th)           return("LateUP")
  if (rel180 <= down_th & rel30  > down_th)                  return("LateDown")
  if (rel30  <= down_th & rel180 >= up_th)                   return("DownThenUp")
  if (rel30  <= down_th & rel180 <= down_th)                 return("Down")
  return("Mixed/Other")
}

## Apply trajectory labels to all k-means clusters.
centroids$label <- mapply(label_shape,
                          centroids$rel30, centroids$rel180, centroids$gap,
                          MoreArgs = list(up_th=up_th, down_th=down_th, gap_th=gap_th))


## Count the number of genes in each cluster for plot labels.
size_map <- clu_df %>%
  count(Cluster, name = "size") 


## Compute centroid changes relative to 0min and the 30min-to-180min gap.
centroids <- centroids %>%
  left_join(size_map, by = "Cluster") %>%
  mutate(
    # C1, C2 … from "cluster_1" etc.
    Cid = sub("^cluster_", "C", Cluster),
    cluster_name = paste0(Cid, " – ", label, " (n=", size, ")")
  )




## Convert centroid table from wide to long format for ggplot.
centroids_long <- pivot_longer(
  centroids, cols = c(`0min`, `30min`, `180min`),
  names_to = "time", values_to = "meanVST"
)
centroids_long$time <- factor(centroids_long$time, levels = c("0min","30min","180min"))


## Create text labels combining cluster ID and trajectory label.
centroids_long$clabel <- paste0(centroids_long$Cid, " – ", centroids_long$label)

# take the last timepoint for labeling
lab_df <- subset(centroids_long, time == "180min")




## Fix factor order for response labels
resp_levels <- c("EarlyUp","SustainedUp","LateUP","LateDown",
                 "UpThenDown","DownThenUp","Down","Mixed/Other")
centroids_long$label <- factor(centroids_long$label, levels = resp_levels)


## Manual colours for centroid response classes.
pal <- c(
  EarlyUp        = "#299DB2",  # teal
  #SustainedUp  = "#1F7A8C",  # deep teal
  LateUP        = "#6B58B5",  # purple
  LateDown     = "#D62728",  # red
  UpThenDown   = "#FF7F0E",  # orange
  # DownThenUp   = "#2CA02C",  # green
  # Down         = "#46464E",  # charcoal
  "Mixed/Other"= "#9E9E9E"   # grey
)


## Plot cluster centroid trajectories over leptin stimulation time.
ggplot(centroids_long, aes(time, meanVST, group = Cluster, color = label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_text_repel(data = lab_df, aes(label = clabel),
                  nudge_x = 0.15, direction = "y",
                  segment.size = 0.2, size = 5.5, show.legend = FALSE) +
  theme_classic() +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.2))) +
  scale_color_manual(name = "Response type", values = pal, drop = FALSE) +
  labs(x = "", y = "Mean VST expression values", title = "Cluster centroids")+
  theme(axis.text = element_text(size=14), axis.title  = element_text(size=14),
        title =element_text(size=16), legend.text =element_text(size=14),
        legend.title =element_text(size=16))
#-----------------------------------------------------------#






# ============================================================
#*Enrichments 
# ============================================================

#-----------for each km-cluster from heatmap-------------
## Load custom enrichment helper functions.
source("R/helpers/enrichment_helpers.R")

## Run pathway enrichment separately for each heatmap-derived k-means cluster.
cluster_enrich <- run_cluster_enrichment(
  clu_df = clu_df,
  gene_col = "GeneID",
  cluster_col = "Cluster",
  cluster_order = c("cluster_3", "cluster_4", "cluster_2", "cluster_1", "cluster_5"),
  cluster_labels = c("C3", "C4", "C2", "C1", "C5"),
  cluster_file_labels = c("cl3", "cl4", "cl2", "cl1", "cl5"),
  dbs = c(
    "GO_Biological_Process_2025",
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "Reactome_Pathways_2024"
  ),
  outdir = "transcriptomics_data_analysis/wt_time_effect",
  write_csvs = TRUE,
  write_hallmark_xlsx = TRUE,
  hallmark_xlsx = "transcriptomics_data_analysis/wt_time_effect/SF3_hallmark_enrichment_WT_clusters.xlsx"
)

## Combine selected Hallmark enrichment sheets into one dataframe for plotting.
hallmark <- rbind(readxl::read_xlsx("SF3_hallmark_enrichment_WT_clusters.xlsx", sheet=1),
             readxl::read_xlsx("SF3_hallmark_enrichment_WT_clusters.xlsx", sheet=2),
             readxl::read_xlsx("SF3_hallmark_enrichment_WT_clusters.xlsx", sheet=3),
             readxl::read_xlsx("SF3_hallmark_enrichment_WT_clusters.xlsx", sheet=4),
             readxl::read_xlsx("SF3_hallmark_enrichment_WT_clusters.xlsx", sheet=5))

## Keep significantly enriched Hallmark terms by adjusted p-value.
hallmark <- hallmark[hallmark$Adjusted.P.value<0.05,]

## Prepare enrichment results for dot-plot visualization.
df <- hallmark %>%
  mutate(
    term = str_wrap(Term, width = 40),
    overlap_num = as.numeric(sub("/.*", "", Overlap)),     # numerator
    overlap_den = as.numeric(sub(".*/", "", Overlap)),     # denominator 
    min_p = ave(Adjusted.P.value, term, FUN = function(x) min(x, na.rm = TRUE))
  ) %>%
  mutate(term = fct_reorder(term, min_p, .desc = FALSE))

##  Selected pathway terms of interest to display.

poi <- c("Interferon Gamma Response",
         "Epithelial Mesenchymal Transition",
         "IL-2/STAT5 Signaling",
         "Inflammatory Response",
         "IL-6/JAK/STAT3 Signaling",
         "TNF-alpha Signaling via NF-kB",
         "Hypoxia", 
         "Interferon Alpha Response",
         "Cholesterol Homeostasis",
         "Apoptosis", 
         "TGF-beta Signaling", 
         "p53 Pathway",
         "Complement")

length(unique(poi))
## Restrict dot plot to selected terms of interest.
df <- df[df$Term %in% poi, ]


## Dot plot showing enrichment strength and overlap size by cluster.
ggplot(df, aes(x = cluster, y = term)) +
  geom_point(aes(size = overlap_num, color = Adjusted.P.value), alpha = 0.9) +
  scale_size(range = c(2, 12), name = "Overlap (# genes)") +
  scale_color_distiller(palette = "Reds", name = "Adj. p-value",
                        direction = -1, limits = c(0, 0.04)) +
  guides(
    color = guide_colorbar(barheight = 4)  
  ) +
  labs(x = "", y = "", title = "Cluster enrichment") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8), size = 12),
    axis.text.y  = element_text(face = "bold",size=12, angle = 0),
    legend.key.size    = unit(0.4, "lines"),
    legend.spacing.y   = unit(0.2, "lines"),
    legend.box.spacing = unit(0.2, "lines")
  )


#-----------Stat3Δ vs WT enrichment using GO BP-------------------

##source("R/helpers/enrichment_helpers.R")

## Build up/down gene lists for Stat3_DEL versus WT at each timepoint.

stat3_vswt_gene_sets <- list(
  `0min` = list(
    up = stat3_vs_wt$tables$Stat3DEL_vsWT_0$up$id,
    down = stat3_vs_wt$tables$Stat3DEL_vsWT_0$down$id
  ),
  `30min` = list(
    up = stat3_vs_wt$tables$Stat3DEL_vsWT_30$up$id,
    down = stat3_vs_wt$tables$Stat3DEL_vsWT_30$down$id
  ),
  `180min` = list(
    up = stat3_vs_wt$tables$Stat3DEL_vsWT_180$up$id,
    down = stat3_vs_wt$tables$Stat3DEL_vsWT_180$down$id
  )
)

##  Run GO BP enrichment for Stat3_DEL versus WT up/down gene sets.

stat3_vswt_enrich <- run_up_down_enrichment_set(
  gene_sets = stat3_vswt_gene_sets,
  dbs = c("GO_Biological_Process_2025"),
  outdir = "transcriptomics_data_analysis/Stat3_vsWT",
  p_cutoff = 0.05,
  write_raw_xlsx = TRUE
)


#---------Etv6Δ vs WT enrichment using GO BP----------------

## Build up/down gene lists for Etv6_DEL versus WT at each timepoint.
etv6_vswt_gene_sets <- list(
  `0min` = list(
    up = etv6_vs_wt$tables$Etv6DEL_vsWT_0$up$id,
    down = etv6_vs_wt$tables$Etv6DEL_vsWT_0$down$id
  ),
  `30min` = list(
    up = etv6_vs_wt$tables$Etv6DEL_vsWT_30$up$id,
    down = etv6_vs_wt$tables$Etv6DEL_vsWT_30$down$id
  ),
  `180min` = list(
    up = etv6_vs_wt$tables$Etv6DEL_vsWT_180$up$id,
    down = etv6_vs_wt$tables$Etv6DEL_vsWT_180$down$id
  )
)


## Run GO BP enrichment for Etv6_DEL versus WT up/down gene sets.

etv6_vswt_enrich <- run_up_down_enrichment_set(
  gene_sets = etv6_vswt_gene_sets,
  dbs = c("GO_Biological_Process_2025"),
  outdir = "transcriptomics_data_analysis/Etv6_vsWT",
  p_cutoff = 0.05,
  write_raw_xlsx = TRUE
)
#-----------------------------------------------------------#



# ============================================================
#Plotting selected genes
# ============================================================



#---------------------------------WT time course--------------------------------
## Start selected-gene plotting using WT batch-corrected VST expression.
dt_wt<- mat_bc
colnames(dt_wt)
dt_wt <- t(dt_wt)
dt_wt <- as.data.frame(dt_wt)         

meta$condition <- meta$time
dt_wt$condition <- meta$condition
unique(dt_wt$condition)  

## Set x-axis order for WT time-course barplots.
dt_wt$condition <- factor(dt_wt$condition,c("0min","30min","180min"))

library(Hmisc)

#EXAMPLE: barplot for Socs3 in the WT time course
wt_colours <- c("#6E6E78", "#5A5A63", "#46464E")

p <- ggplot(dt_wt, aes(x=condition, y=Socs3, fill=condition))+
  stat_summary(fun=mean, geom="bar", width = 0.55)+ 
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar", width=.2)+
  geom_point(position = position_jitter(width = 0.01),size = 0.5, color = "black") +
  scale_fill_manual(values=wt_colours)+
  labs(title="Socs3",x ="", y = "Mean expression values" )+
  theme_bw()

p <- p +
  coord_cartesian(ylim = c(5, 15)) +             
  theme(
    plot.title  = element_text(size=16, face="bold", hjust=0.5),
    axis.title.x = element_text(size=10, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.x  = element_text(face="bold", size=12, angle=0),
    axis.text.y  = element_text(face="bold", size=10),
    plot.margin  = margin(4, 4, 4, 4, unit="mm")
  ) +
  theme(
    panel.border = element_rect(colour="black", linewidth=0.25),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black", linewidth=0),
    legend.position = "none"
  )
p
#------------------------ ETv6 del time course --------------------------------
# barplots for selectd genes for etv6del and WT cell types

meta0 <- as.data.frame(colData(vsd))
## Restrict  to batch II.
meta0 <- meta0 %>% dplyr::filter(batch == "II")
meta0$time  <- factor(meta0$time,  levels = c("0min","30min","180min"))
meta0$type  <- factor(meta0$type,  levels = c("WT","Etv6_DEL"))
meta0$batch <- factor(meta0$batch)

keep <- meta0$type %in% c("WT","Etv6_DEL")
meta  <- droplevels(meta0[keep, ])
vsd2  <- vsd[, rownames(meta)]
expr  <- assay(vsd2)
## Convert expression matrix to sample-by-gene dataframe for plotting.
dt_etv6 <- as.data.frame(t(expr))

## Metadata subset used to define plotting condition labels for Etv6_DEL/WT samples.
etv6_coldata <- colData[25:48,]

dt_etv6$condition <- etv6_coldata$condition

##  Set x-axis order for Etv6_DEL versus WT condition groups.
dt_etv6$condition <- factor(dt_etv6$condition,c("WT_0min","WT_30min","WT_180min","Etv6_DEL_0min","Etv6_DEL_30min","Etv6_DEL_180min"))

## Example selected-gene barplot for Nox4 in Etv6_DEL versus WT.

etv6_wt_cols <- c("#6E6E78", "#5A5A63", "#46464E", "#B9A7DE", "#8D76D0", "#6B58B5")

p <- ggplot(dt_etv6, aes(x=condition, y=Nox4, fill=condition))+
  stat_summary(fun=mean, geom="bar", width = 0.55)+ 
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar", width=.2)+
  geom_point(position = position_jitter(width = 0.01),size = 0.5, color = "black")+
  scale_fill_manual(values=etv6_wt_cols)+
  labs(title="Nox4",x ="Time(min) ", y = "Mean Expression" )+
  theme_bw()

p + theme(plot.title = element_text( size=14, face="bold", hjust = 0.5),axis.title.x = element_text(size=10, face="bold",),axis.title.y = element_text(size=10, face="bold"),
          axis.text.x = element_text(face="bold",size=8,angle = 0,hjust= 1),
          axis.text.y = element_text(face="bold",size=10)) +
  theme(panel.border =element_rect(colour = "black", linewidth=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",linewidth=0))
#----------------------------------------------------------------------------#





# ============================================================
sessionInfo()
# ============================================================

R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
  [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
  [1] tidyr_1.3.1                 ComplexHeatmap_2.20.0      
[3] openxlsx_4.2.8              limma_3.60.6               
[5] pheatmap_1.0.13             RColorBrewer_1.1-3         
[7] dplyr_1.1.4                 genefilter_1.86.0          
[9] vsn_3.72.0                  EnhancedVolcano_1.22.0     
[11] ggrepel_0.9.6               ggplot2_4.0.0              
[13] DESeq2_1.44.0               SummarizedExperiment_1.34.0
[15] Biobase_2.64.0              MatrixGenerics_1.16.0      
[17] matrixStats_1.5.0           GenomicRanges_1.56.2       
[19] GenomeInfoDb_1.40.1         IRanges_2.38.1             
[21] S4Vectors_0.42.1            BiocGenerics_0.50.0        

loaded via a namespace (and not attached):
  [1] DBI_1.2.3               writexl_1.5.4           readxl_1.4.5           
[4] rlang_1.1.6             magrittr_2.0.4          clue_0.3-66            
[7] GetoptLong_1.0.5        compiler_4.4.1          RSQLite_2.4.3          
[10] png_0.1-8               vctrs_0.6.5             pkgconfig_2.0.3        
[13] shape_1.4.6.1           crayon_1.5.3            fastmap_1.2.0          
[16] magick_2.9.0            XVector_0.44.0          labeling_0.4.3         
[19] utf8_1.2.6              UCSC.utils_1.0.0        preprocessCore_1.66.0  
[22] purrr_1.1.0             bit_4.6.0               zlibbioc_1.50.0        
[25] cachem_1.1.0            jsonlite_2.0.0          blob_1.2.4             
[28] DelayedArray_0.30.1     BiocParallel_1.38.0     parallel_4.4.1         
[31] cluster_2.1.8.1         R6_2.6.1                stringi_1.8.7          
[34] pkgload_1.4.1           cellranger_1.1.0        Rcpp_1.1.0             
[37] iterators_1.0.14        Matrix_1.7-4            splines_4.4.1          
[40] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1      
[43] abind_1.4-8             doParallel_1.0.17       codetools_0.2-20       
[46] affy_1.82.0             lattice_0.22-7          tibble_3.3.0           
[49] withr_3.0.2             KEGGREST_1.44.1         S7_0.2.0               
[52] survival_3.8-3          zip_2.3.3               circlize_0.4.16        
[55] Biostrings_2.72.1       pillar_1.11.1           affyio_1.74.0          
[58] BiocManager_1.30.26     foreach_1.5.2           generics_0.1.4         
[61] scales_1.4.0            xtable_1.8-4            glue_1.8.0             
[64] tools_4.4.1             hexbin_1.28.5           annotate_1.82.0        
[67] locfit_1.5-9.12         XML_3.99-0.19           Cairo_1.6-5            
[70] AnnotationDbi_1.66.0    colorspace_2.1-2        GenomeInfoDbData_1.2.12
[73] cli_3.6.5               S4Arrays_1.4.1          gtable_0.3.6           
[76] digest_0.6.37           SparseArray_1.4.8       rjson_0.2.23           
[79] farver_2.1.2            memoise_2.0.1           lifecycle_1.0.4        
[82] httr_1.4.7              GlobalOptions_0.1.2     statmod_1.5.0          
[85] bit64_4.6.0-1          
