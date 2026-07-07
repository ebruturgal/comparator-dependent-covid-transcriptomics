# Comparator-Dependent Host Transcriptomic Classifiers for COVID-19

This repository contains the R analysis script used for a comparator-aware benchmark of host transcriptomic machine learning classifiers for COVID-19. The analysis evaluates COVID-19 classification against three comparator settings in GSE282464 and performs external validation using GSE199816 and GSE161731.

## Repository Contents

```text
.
|-- README.md
|-- scripts/
|   |-- final_stable_analysis.R
|   `-- make_table2_external_validation.R
`-- .gitignore
```

## Analysis Overview

The script performs:

- construction of a GSE282464 count matrix from raw sample-level count files
- metadata matching from GEO
- three binary comparator tasks:
  - COVID-19 vs Influenza
  - COVID-19 vs non-influenza viral infections
  - COVID-19 vs Sepsis / septic shock
- edgeR filtering, TMM normalization, and logCPM transformation
- 5-fold internal cross-validation
- variance-based top-100 feature selection within each outer training fold
- model comparison using:
  - Elastic Net
  - GPBoost-LightGBM
  - XGBoost
  - Random Forest
- feature-importance summaries
- locked Elastic Net external validation in:
  - GSE199816 for COVID-19 vs sepsis
  - GSE161731 for COVID-19 vs influenza
- ROC, calibration, DeLong comparison, decision-curve, and paper-ready summary outputs

## Requirements

Install R, then run the script from the repository root. The script installs missing CRAN and Bioconductor packages where possible.

Main R packages:

- data.table
- dplyr
- stringr
- ggplot2
- tidyr
- pROC
- glmnet
- xgboost
- randomForest
- pheatmap
- scales
- gpboost
- GEOquery
- edgeR
- rmda
- AnnotationDbi
- org.Hs.eg.db

## Input Data

The main development dataset is GSE282464. Before running the full script, place the raw count files in this structure:

```text
GSE282464/
`-- untarred_raw/
    |-- <sample>.counts.gz
    |-- <sample>.counts.gz
    `-- ...
```

The script downloads supplementary files for GSE199816 and GSE161731 using `GEOquery::getGEOSuppFiles()` when those sections are executed.

## How to Run

From the repository root:

```r
source("scripts/final_stable_analysis.R")
```

Or from a terminal:

```bash
Rscript scripts/final_stable_analysis.R
```

The script writes CSV, PNG, and PDF outputs to the working directory.

To create the manuscript-ready external validation Table 2 after the main workflow finishes:

```r
source("scripts/make_table2_external_validation.R")
```

or:

```bash
Rscript scripts/make_table2_external_validation.R
```

This writes:

- `Table2_ExternalValidation_Performance.csv`
- `Table2_ExternalValidation_Performance_Detailed.csv`

## Important Notes

- Large raw GEO files and generated analysis outputs are intentionally ignored by `.gitignore`.
- The script is provided as a complete analysis workflow. Some later blocks regenerate external validation figures from prediction CSV files created earlier in the workflow.
- The internal benchmark uses 5-fold cross-validation with variance-based top-100 feature selection. It is not a fully nested, fold-restricted preprocessing pipeline.
- If external GEO downloads fail, manually download the supplementary files from GEO and place them in the expected dataset folders.

## Main Outputs

Examples of generated files include:

- `GSE282464_COVID_vs_Influenza_Revised_Summary.csv`
- `GSE282464_COVID_vs_nonInfluenzaViral_Revised_Summary.csv`
- `GSE282464_COVID_vs_Sepsis_Revised_Summary.csv`
- `Paper_Ready_Summary_Revised.csv`
- `ElasticNet_Importance_Revised_AllTasks.csv`
- `RandomForest_Importance_Revised_AllTasks.csv`
- `GSE199816_ExternalValidation_Predictions.csv`
- `GSE161731_ExternalValidation_Predictions.csv`
- `ExternalValidation_Calibration_Summary.csv`
- `DeLong_ExternalValidation_Comparison.csv`
- `Table2_ExternalValidation_Performance.csv`
- ROC, calibration, decision-curve, and heatmap figures

## Citation

If you use this workflow, cite the relevant GEO datasets and software packages used in the analysis, including edgeR, glmnet, xgboost, randomForest, gpboost, pROC, and GEOquery.

### Build Table 2

This script combines the external validation CSV outputs into one manuscript-ready Table 2.

```bash
Rscript scripts/make_table2_external_validation.R
