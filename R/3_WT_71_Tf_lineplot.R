############################################################
# Script: 3_WT_71_Tf_lineplot.R
# Project: Temporal organization of leptin signaling reveals distinct
#          regulatory roles for STAT3 and ETV6
#
# Purpose:
#   Generate WT TF activity line plots for the 71 transcription
#   factors selected from the WT TF activity heatmap.
#
# Description:
#   This script uses the TF clusters obtained from
#   2_Tfs_Pathways_activity_dynamics.R and plots mean WT TF activity
#   over the leptin time course for each heatmap cluster.
#
# Required input objects from Tfs_Pathways_activity_dynamics.R:
#   M_tf        : TF activity matrix from decoupleR/DoRothEA.
#   meta        : Sample metadata.
#   wt_flag     : Logical vector identifying WT samples.
#   clu_df_tf71 : Data frame containing TF-to-cluster assignments.
#
# Required helper function:
#   plot_time_lines() from Tfs_Pathways_activity_dynamics.R
#
# Main outputs:
#   p1, p2, p3, p4 : WT TF activity line plots for clusters 1–4.
#
# Notes:
#   - TFs are coloured by broad transcription factor family.
#   - Related TFs within the same family receive related colour shades.
#   - Cluster identity is inherited from the WT TF activity heatmap.
############################################################


#=====================================
#Load required packages
#=====================================
library(dplyr)
library(colorspace)


#=====================================
# Sanity checks for required input objects
#=====================================

stopifnot(exists("M_tf"))
stopifnot(exists("meta"))
stopifnot(exists("wt_flag"))
stopifnot(exists("clu_df_tf71"))
stopifnot(exists("plot_time_lines"))

stopifnot(all(c("GeneID", "Cluster") %in% colnames(clu_df_tf71)))




#-------------------------------------------------------------
# Assign each TF to a broad biological/structural TF family.
# This is used only for colouring the line plots.

tf_family <- c(
  # STAT
  STAT1="STAT", STAT3="STAT", STAT4="STAT", STAT5A="STAT", STAT5B="STAT", STAT6="STAT",
  
  # AP-1 / immediate early (stress module)
  JUN="AP1", JUNB="AP1", FOSL2="AP1", ATF2="AP1",
  
  # NFKB / REL
  REL="NFKB", RELA="NFKB", NFKB1="NFKB", NFKB2="NFKB",
  
  # ETS family + close ETS-like
  ETS1="ETS", ETS2="ETS", FLI1="ETS", ELF1="ETS", ETV4="ETS", GABPA="ETS",
  
  # ELK / SRF (MAPK–SRF program; keep as its own functional module)
  ELK1="SRF_ELK", ELK4="SRF_ELK", SRF="SRF_ELK",
  
  # SMAD / TGFβ
  SMAD3="SMAD", SMAD4="SMAD",
  
  # Hypoxia / bHLH-PAS / AhR axis
  HIF1A="bHLH_PAS", EPAS1="bHLH_PAS", AHR="bHLH_PAS", ARNT="bHLH_PAS",
  
  # Nuclear receptors
  ESR1="NuclearReceptor", ESR2="NuclearReceptor", PPARA="NuclearReceptor",
  PPARG="NuclearReceptor", NR4A1="NuclearReceptor", RARA="NuclearReceptor",
  
  # CEBP
  CEBPA="CEBP", CEBPB="CEBP", CEBPD="CEBP",
  
  # Forkhead
  FOXM1="Forkhead", FOXL2="Forkhead", FOXK2="Forkhead",
  
  # Cell-cycle / proliferation
  E2F7="CellCycle", MYC="CellCycle", TCF7L2="CellCycle", ZFP217="CellCycle",
  
  # KLF / SP / related zinc-finger
  KLF3="KLF_SP", KLF6="KLF_SP", SP1="KLF_SP", NFIC="KLF_SP",
  
  # MEF2
  MEF2A="MEF2", MEF2B="MEF2", MEF2C="MEF2",
  
  # NRF2 / small Maf oxidative stress
  NFE2L2="NRF2_MAF", MAFG="NRF2_MAF", NFE2L1="NRF2_MAF",
  
  # CREB/ATF family (cAMP/ISR-like; keep separate from AP1 for clarity)
  CREB1="CREB_ATF", CREB3="CREB_ATF", CREM="CREB_ATF", ATF4="CREB_ATF",
  
  # GATA
  GATA2="GATA",
  
  # BATF (AP-1 related but often treated separately; you can set BATF="AP1" if you prefer)
  BATF="BATF",
  
  # Other / mixed or cofactor/TF-like
  PRDM1="Other", MEIS1="Other", PKNOX1="Other", OTX2="Other",
  NCOA2="Other", RFX2="Other", TBP="Other", TFAP2A="TFAP2"
)



#-------------------------------------------------------------
# Base colour assigned to each TF family.
# Individual TFs within a family are later given darker/lighter
# shades of the corresponding base colour.


family_base_vibrant <- c(
  STAT            = "#00C2A8",  # turquoise (non-grey)
  AP1             = "#E31A1C",  # vivid red
  NFKB            = "#7A00CC",  # vivid violet
  ETS             = "#00A650",  # vivid green
  SRF_ELK         = "#FF2DAA",  # hot magenta
  SMAD            = "#4B00FF",  # electric indigo
  bHLH_PAS        = "#00B8D9",  # bright teal/cyan
  NuclearReceptor = "#FF8C00",  # vivid orange
  CEBP            = "#FF6A00",  # orange-red
  Forkhead        = "#A3D900",  # bright yellow-green
  CellCycle       = "#0066FF",  # vivid blue
  KLF_SP          = "#00B0FF",  # sky blue (instead of grey)
  MEF2            = "#0052CC",  # strong royal blue (still distinct from STAT)
  NRF2_MAF        = "#2ECC71",  # vivid emerald
  CREB_ATF        = "#FFD000",  # vivid gold
  GATA            = "#FF4D00",  # vivid orange-red (replaces deep rose)
  BATF            = "#FF3D00",  # bright vermillion
  TFAP2           = "#FF6FB1",  # bright pink (not muted)
  Other           = "#00AEEF"   # bright cyan (avoid grey)
)


#-------------------------------------------------------------
# Generate a named colour vector for a set of TFs.
#
# TFs from the same family are assigned related shades of the
# same base colour, which helps preserve biological grouping
# while keeping individual TF lines distinguishable.


make_semantic_palette_vibrant <- function(features, tf_family, family_base,
                                          darken_amt = 0.20, lighten_amt = 0.06) {
  features <- unique(as.character(features))
  features <- features[!is.na(features)]
  
  fam <- unname(tf_family[features])
  fam[is.na(fam)] <- "Other"
  names(fam) <- features
  
  out <- setNames(rep(NA_character_, length(features)), features)
  
  for (f in unique(fam)) {
    tfs <- names(fam)[fam == f]
    base <- family_base[[f]]
    if (is.null(base) || is.na(base)) base <- family_base[["Other"]]
    
    n <- length(tfs)
    if (n == 1) {
      out[tfs] <- base
    } else {
      amounts <- seq(darken_amt, -lighten_amt, length.out = n)
      cols <- vapply(amounts, function(a) {
        if (a >= 0) colorspace::darken(base, amount = a) else colorspace::lighten(base, amount = -a)
      }, character(1))
      out[tfs] <- cols
    }
  }
  
  out[is.na(out)] <- family_base[["Other"]]
  out
}

#-------------------------------------------------------------
# Extract TFs belonging to one heatmap cluster.
get_cluster_features <- function(clu_df, clu) {
get_cluster_features <- function(clu_df, clu) {
  clu_df %>% dplyr::filter(Cluster == clu) %>% dplyr::pull(GeneID)
}

#-------------------------------------------------------------

feat1 <- get_cluster_features(clu_df_tf71, "cluster_1")
feat2 <- get_cluster_features(clu_df_tf71, "cluster_2")
feat3 <- get_cluster_features(clu_df_tf71, "cluster_3")
feat4 <- get_cluster_features(clu_df_tf71, "cluster_4")

# Check that all selected cluster TFs are present in the TF activity matrix.
stopifnot(all(unlist(list(feat1, feat2, feat3, feat4)) %in% rownames(M_tf)))


# Use the vibrant family-level palette for all cluster plots.
family_base <- family_base_vibrant

#-------------------------------------------------------------
# Create cluster-specific colour palettes.
pal1 <- make_semantic_palette_vibrant(feat1, tf_family, family_base)
pal2 <- make_semantic_palette_vibrant(feat2, tf_family, family_base)
pal3 <- make_semantic_palette_vibrant(feat3, tf_family, family_base)
pal4 <- make_semantic_palette_vibrant(feat4, tf_family, family_base)

#-------------------------------------------------------------
# Plot mean WT TF activity over time for each heatmap cluster.
p1 <- plot_time_lines(M_tf, feat1, "TF activity (WT) — cluster 1", meta, wt_flag, colors=pal1)
p2 <- plot_time_lines(M_tf, feat2, "TF activity (WT) — cluster 2", meta, wt_flag, colors=pal2)
p3 <- plot_time_lines(M_tf, feat3, "TF activity (WT) — cluster 3", meta, wt_flag, colors=pal3)
p4 <- plot_time_lines(M_tf, feat4, "TF activity (WT) — cluster 4", meta, wt_flag, colors=pal4)

#-------------------------------------------------------------
