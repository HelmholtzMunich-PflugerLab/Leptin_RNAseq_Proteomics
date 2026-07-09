project/
│
├── README.md
├── .gitignore
│
├── Data/
│   ├── Counts_WT_Stat3del_Etv6del.txt
│   ├── colData_WT_Stat3del_Etv6del.txt
│   └── Proteomics_abundance.xlsx
│
└── R/
    ├── 1_RNAseq_degs_heatmap_pathways.R
    ├── 2_TFs_Pathways_activity_dynamics.R
    ├── 3_WT_71_TF_lineplot.R
    ├── 4_TF_and_PW_act_correlation.R
    ├── 5_proteomics_main.R
    └── helpers/
        ├── deg_analysis_functions.R
        ├── enrichment_helpers.R
    └── proteomics_helper_functions.R