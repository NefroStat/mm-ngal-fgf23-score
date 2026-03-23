# Multiple Myeloma Composite Risk Score (NGAL & FGF-23)

[![R](https://img.shields.io/badge/R-4.5.1-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Data: Restricted](https://img.shields.io/badge/Data-Restricted_by_Ethics_Committee-red.svg)]()

## Overview
This repository contains the complete analytical pipeline and R scripts for the clinical study investigating the prognostic value of novel biomarkers in Multiple Myeloma (MM). 

The primary objective of this study was to evaluate the independent prognostic value of baseline serum Neutrophil Gelatinase-Associated Lipocalin (NGAL) and Fibroblast Growth Factor-23 (FGF-23) in patients with multiple myeloma. Furthermore, we aimed to develop a composite risk score to predict overall survival independently of baseline estimated glomerular filtration rate (eGFR) and International Staging System (ISS) stage.

## Data Availability Statement
Due to strict patient privacy regulations and restrictions imposed by the local Ethics Committee, the clinical dataset (even in an anonymized format) **cannot be made publicly available** in this repository. 

The data that support the findings of this study are available from the corresponding author upon reasonable request, subject to ethical approval and data-sharing agreements. 

To run the scripts provided in this repository, researchers must obtain the dataset and place it in the `data/processed/` directory under the name `clean_data.xlsx`.

## Expected Data Structure
For transparency and reproducibility, the analytical scripts expect the input dataset (`clean_data.xlsx`) to contain the following key variables:
* `Study_ID`: Unique patient identifier (e.g., MM_PT_001)
* `тривалість`: Overall survival time in months (numeric)
* `Exitus`: Survival status (1 = deceased, 0 = censored)
* `GFR CKD-EPI`: Baseline eGFR in mL/min/1.73m² (numeric)
* `stage ISS`: International Staging System stage (1, 2, or 3)
* `NGAL`: Serum NGAL concentration in ng/mL (numeric)
* `FGF23`: Serum FGF-23 concentration in pg/mL (numeric)
* *Other routine clinical parameters:* `Hb`, `CREA...44`, `Urea`, `CRP`, `alb`, `Total TPROP`, `М-prot g/l`, `β₂-microhlob`.

## Repository Structure

```text
├── data/
│   └── processed/         # Folder where the requested 'clean_data.xlsx' should be placed
├── scripts/               # R scripts for data processing and statistical analysis
│   ├── 01_data_preparation.R
│   ├── 02_baseline_and_boxplots.R
│   ├── 03_cutoff_analysis.R
│   └── 04_survival_analysis.R
├── results/               # Output directory for tables, Kaplan-Meier plots, and Cox models (auto-generated)
└── README.md
Analysis Pipeline (Scripts)
Run the scripts in the following order to reproduce the results:

01_data_preparation.R:
Imports raw clinical data (if available), performs ethical anonymization, standardizes column names, and converts European comma decimals to standard numeric formats.

02_baseline_and_boxplots.R:
Computes baseline clinical characteristics comparing patients by eGFR (<60 vs ≥60 mL/min/1.73m²). Automatically applies Shapiro-Wilk tests for normality to select appropriate statistical tests (Welch's t-test or Mann-Whitney U test). Generates publication-ready boxplots.

03_cutoff_analysis.R:
Determines the optimal prognostic cutoffs for NGAL and FGF-23 using maximally selected rank statistics (survminer package). Generates log-rank statistical cutoff plots.

04_survival_analysis.R:
Generates Kaplan-Meier survival curves for individual biomarkers and the novel composite risk score. Performs Univariate and Multivariate Cox Proportional Hazards Regression, adjusting for baseline eGFR and ISS stage.

System Requirements & Dependencies
The analysis was performed using R software (version 4.5.1). Ensure you have the following R packages installed before running the scripts: install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "readxl", "readr", "writexl", "survival", "survminer"))
How to Run the Code
Clone this repository to your local machine.

Open the project folder in RStudio.

Obtain the dataset from the corresponding author and save it exactly as data/processed/clean_data.xlsx.

Run the scripts sequentially from 01 to 04. All outputs (plots and tables) will be automatically saved in the newly generated results/ folder.

Citation
If you use this code or analytical pipeline in your research, please cite our forthcoming paper:

(Citation placeholder - will be updated upon publication)

Contact
For any questions regarding the analysis pipeline or data access requests, please open an issue in this repository or contact the corresponding author (v.bardash@1tmolviv.com).