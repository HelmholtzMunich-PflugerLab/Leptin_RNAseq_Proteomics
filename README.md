# RNA-seq and Proteomics Analysis of Leptin Signalling

This repository contains the R scripts used to analyse RNA-seq and proteomics data generated from leptin-stimulated wild-type (WT), Stat3_DEL and Etv6_DEL cell lines.

## Associated manuscript

**Temporal organization of leptin signaling reveals distinct regulatory roles for STAT3 and ETV6**

The repository contains sequential analysis scripts together with reusable helper functions for RNA-seq and proteomics data analysis.


## Repository structure

```text
.
├── Data/
│   ├── Counts_WT_Stat3del_Etv6del.txt
│   ├── colData_WT_Stat3del_Etv6del.txt
│   └── Proteomics_abundance.xlsx
│
├── R/
│   ├── 1_RNAseq_degs_heatmap_pathways.R
│   ├── 2_TFs_Pathways_activity_dynamics.R
│   ├── 3_WT_71_TF_lineplot.R
│   ├── 4_TF_and_PW_act_correlation.R
│   ├── 5_proteomics_main.R
│   │
│   └── helpers/
│       ├── deg_analysis_functions.R
│       ├── enrichment_helpers.R
│       └── proteomics_helper_functions.R
│
├── README.md
└── .gitignore
```

---

## Analysis workflow

Run the scripts in the following order.

| Step | Script                              | Description                                                                                                |
| ---- | ----------------------------------- | -----------------------------------------------------------------------------------------------------------|
| 1    | 1_RNAseq_degs_heatmap_pathways.R   | RNA-seq preprocessing, PCA, differential expression analysis, enrichment analysis and heatmap generation    |
| 2    | 2_TFs_Pathways_activity_dynamics.R | Estimate transcription factor and pathway activities, analyse WT dynamics and genotype-by-time interactions |
| 3    | 3_WT_71_TF_lineplot.R              | Generate WT transcription factor activity line plots                                                        |
| 4    | 4_TF_and_PW_act_correlation.R      | Correlate transcription factor and pathway activities with WT expression modules                            |
| 5    | 5_proteomics_main.R                | Proteomics analysis and proteome–transcriptome comparisons                                                  |

---

## Input data

The repository requires the following input files located in the `Data/` directory.

- RNA-seq count matrix
- RNA-seq sample metadata
- Normalized proteomics abundance table

---

## Helper scripts

The helper scripts located in `R/helpers/` provide reusable functions used by the main analysis scripts and are sourced automatically where required.

These helper scripts do **not** need to be run independently.

---

## Software

Analysis was performed in **R 4.4.1**.

Major packages include:

- DESeq2
- limma
- decoupleR
- dorothea
- fgsea
- clusterProfiler
- enrichR
- ComplexHeatmap
- ggplot2