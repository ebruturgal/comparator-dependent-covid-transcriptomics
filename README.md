# Comparator-Dependent Host Transcriptomic Classifiers for COVID-19

This repository contains one complete R analysis script for comparator-aware benchmarking of host transcriptomic machine learning classifiers for COVID-19.

## Repository Contents

```text
.
|-- README.md
|-- scripts/
|   |-- final_stable_analysis.R
|   `-- make_peerj_tables_figures.R
`-- .gitignore
```

## What The Script Does

The single script:

- builds the GSE282464 development count matrix
- defines three comparator-aware classification tasks
- selects up to 200 variance-ranked genes within the training data for each fold/task
- runs internal 5-fold cross-validation across four machine learning models
- writes manuscript Table 1 with all 12 rows
- trains locked Elastic Net models for external validation
- validates externally in GSE199816 and GSE161731
- writes manuscript Table 2
- writes feature-importance tables and figures
- writes ROC, calibration, heatmap, fold-wise, and decision-curve outputs

The submission-builder script:

- checks that Table 1 has 12 rows
- checks that Table 2 contains AUC, threshold metrics, Brier score, calibration intercept, and calibration slope
- writes Supplementary Table S4 as selected model parameters
- writes Supplementary Table S6 with the corrected external DeLong comparison
- copies `Supplementary_Tables_ALL.xlsx` when available; this workbook contains S1-S6
- rebuilds Supplementary Figure S9 with the corrected AUC, Brier, calibration, and DeLong values
- copies/renames the final figure files into a PeerJ-ready submission folder

## Main Tables Produced

```text
Table1_Internal_CV_Performance.csv
Table1_Internal_CV_Performance_Numeric.csv
Table2_ExternalValidation_Performance.csv
Table2_ExternalValidation_Performance_Detailed.csv
```

Table 1 is checked to contain 12 rows:

```text
3 comparator settings x 4 machine learning models
```

## How To Run

From the repository root:

```bash
Rscript scripts/final_stable_analysis.R
Rscript scripts/make_peerj_tables_figures.R
```

Or in R/RStudio:

```r
source("scripts/final_stable_analysis.R")
source("scripts/make_peerj_tables_figures.R")
```

To build the final tables and figures from a specific results folder:

```bash
SOURCE_DIR=/path/to/results OUTPUT_DIR=peerj_submission_ready Rscript scripts/make_peerj_tables_figures.R
```

On Windows PowerShell:

```powershell
$env:SOURCE_DIR = "C:/path/to/results"
$env:OUTPUT_DIR = "peerj_submission_ready"
Rscript scripts/make_peerj_tables_figures.R
```

## Manuscript Table And Figure Mapping

```text
Table 1  -> Table_1_Internal_CV_Performance.csv
Table 2  -> Table_2_ExternalValidation_Performance.csv
Figure 1 -> Figure_1_Study_Workflow.png
Figure 2 -> Figure_2_CrossTask_AUC_Heatmap.png
Figure 3 -> Figure_3_CrossTask_ModelComparison.png
Figure 4 -> Figure_4_ExternalValidation_ROC_Calibration_Combined.png
Figure 5 -> Figure_5_Top_ElasticNet_Features_Across_Comparator_Settings.png
S1       -> Supplementary_Figure_S1_COVID_vs_Influenza_Heatmap.png
S2       -> Supplementary_Figure_S2_COVID_vs_NonInfluenzaViral_Heatmap.png
S3       -> Supplementary_Figure_S3_COVID_vs_Sepsis_Heatmap.png
S4 Fig   -> Supplementary_Figure_S4_Top_RandomForest_Features.png
S5 Fig   -> Supplementary_Figure_S5_Foldwise_Model_Performance.png
S6 Fig   -> Supplementary_Figure_S6_Summary_Barplots_Model_Performance.png
S7 Fig   -> Supplementary_Figure_S7_Decision_Curve_GSE161731.pdf
S8 Fig   -> Supplementary_Figure_S8_Decision_Curve_GSE199816.pdf
S9       -> Supplementary_Figure_S9_External_ROC_Calibration_Summary.png
S1-S6 Tables -> Supplementary_Tables_ALL.xlsx
S4 Table -> Supplementary_Table_S4_Selected_Model_Parameters.csv
S6 Table -> Supplementary_Table_S6_DeLong.csv
```

## Input Data

Before running, place the GSE282464 count files here:

```text
GSE282464/
`-- untarred_raw/
    |-- <sample>.counts.gz
    |-- <sample>.counts.gz
    `-- ...
```

The script uses GEOquery for external validation files where possible.

## Notes

Large raw data files and generated outputs are ignored by `.gitignore`.
