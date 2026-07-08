# Comparator-Dependent Host Transcriptomic Classifiers for COVID-19

This repository contains one complete R analysis script for comparator-aware benchmarking of host transcriptomic machine learning classifiers for COVID-19.

## Repository Contents

```text
.
|-- README.md
|-- scripts/
|   `-- final_stable_analysis.R
`-- .gitignore
```

## What The Script Does

The single script:

- builds the GSE282464 development count matrix
- defines three comparator-aware classification tasks
- runs internal 5-fold cross-validation across four machine learning models
- writes manuscript Table 1 with all 12 rows
- trains locked Elastic Net models for external validation
- validates externally in GSE199816 and GSE161731
- writes manuscript Table 2
- writes feature-importance tables and figures
- writes ROC, calibration, heatmap, fold-wise, and decision-curve outputs

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
```

Or in R/RStudio:

```r
source("scripts/final_stable_analysis.R")
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

