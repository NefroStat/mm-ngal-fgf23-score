# ==============================================================================
# Script: 01_data_preparation.R
# Description: Imports raw clinical data, performs ethical anonymization 
#              (removes patient names), standardizes column names, and ensures 
#              correct data types. Exports a clean dataset for downstream analysis.
#              Formatted according to Open Science / Nature Portfolio standards.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
# Uncomment the following line if packages are not installed yet:
# install.packages(c("readxl", "dplyr", "readr", "writexl"))

library(readxl)
library(dplyr)
library(readr)
library(writexl) # Added for saving human-readable Excel files

# Set up project directory structure for data (if it doesn't exist)
if (!dir.exists("data")) {
  dir.create("data")
}
if (!dir.exists("data/raw")) {
  dir.create("data/raw")
}
if (!dir.exists("data/processed")) {
  dir.create("data/processed")
}
message("Verified standard data directory structure (data/raw, data/processed).")

# ------------------------------------------------------------------------------
# 1. Data Import (Using RELATIVE path for GitHub portability)
# ------------------------------------------------------------------------------
# This assumes your working directory is the main project folder (myeloma-ckd-survival)
raw_file_path <- "data/таблиця 95 пац. з тривалістю хвороби.xlsx"

if (file.exists(raw_file_path)) {
  # Read the Excel file, suppressing name repair warnings
  df <- read_excel(raw_file_path, .name_repair = "unique_quiet")
  message("Success: Raw data successfully imported using relative path.")
} else {
  stop("Error: Raw data file not found at 'data/...'. Please ensure your working directory is set to the project root and the file is in the 'data' folder.")
}

# ------------------------------------------------------------------------------
# 2. Ethical Anonymization (HIPAA/GDPR Compliance)
# ------------------------------------------------------------------------------
# Remove the patient name column ("П.І.П.") and assign unique Study IDs
if ("П.І.П." %in% colnames(df)) {
  df <- df %>%
    select(-`П.І.П.`) %>% # Drop names
    mutate(Study_ID = paste0("MM_PT_", sprintf("%03d", row_number())), .before = 1)
  message("Data anonymized: Patient names removed, unique Study_IDs generated.")
}

# ------------------------------------------------------------------------------
# 3. Technical Standardization (Fixing issues for downstream analysis)
# ------------------------------------------------------------------------------
# Fix the problematic hyphen in FGF-23 (which breaks R formulas)
if ("FGF-23" %in% colnames(df)) {
  colnames(df)[colnames(df) == "FGF-23"] <- "FGF23"
}

# Ensure critical survival and numerical variables are formatted strictly as numeric
# This handles European comma decimals (e.g., 3,14 to 3.14)
numeric_columns <- c("тривалість", "Exitus", "GFR CKD-EPI", "NGAL", "FGF23", 
                     "CREA...44", "Urea", "CRP", "Hb", "Са", "alb", "Total TPROP")

for (col in numeric_columns) {
  if (col %in% colnames(df)) {
    # Replace commas with dots and convert to numeric
    df[[col]] <- as.numeric(gsub(",", ".", df[[col]])) 
  }
}

message("Variable types standardized. Hyphens removed from variable names.")

# ------------------------------------------------------------------------------
# 4. Export Clean Data (CSV for GitHub, XLSX for local review)
# ------------------------------------------------------------------------------
# 1. Save as CSV for R/GitHub compatibility (Machine-readable standard)
processed_path_csv <- "data/processed/clean_data.csv"
write.csv(df, processed_path_csv, row.names = FALSE)

# 2. Save as XLSX for local review (Solves regional delimiter issues in Excel)
processed_path_excel <- "data/processed/clean_data.xlsx"
write_xlsx(df, processed_path_excel)

message(paste("Success: Clean dataset exported to", processed_path_csv, "and", processed_path_excel))

# ------------------------------------------------------------------------------
# 5. Session Information
# ------------------------------------------------------------------------------
# sessionInfo()