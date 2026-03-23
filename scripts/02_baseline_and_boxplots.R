# ==============================================================================
# Script: 02_baseline_and_boxplots.R
# Description: Computes baseline characteristics (Table 1) comparing patients 
#              by eGFR (<60 vs >=60 mL/min). Automatically selects statistical 
#              tests (t-test or Mann-Whitney) based on Shapiro-Wilk normality.
#              Calculates group sizes (n) and generates publication-ready boxplots.
#              Formatted according to Nature Portfolio / High-Impact standards.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(writexl)

if (!dir.exists("results")) {
  dir.create("results")
  message("Created directory: 'results/'")
}

# ------------------------------------------------------------------------------
# 1. Data Preparation & Variable Definitions
# ------------------------------------------------------------------------------
# Create English group labels for eGFR (<60 and >=60)
df$eGFR_group_eng <- ifelse(df$`GFR CKD-EPI` >= 60, "eGFR ≥ 60", 
                            ifelse(df$`GFR CKD-EPI` < 60, "eGFR < 60", NA))

# Factorize to ensure consistent plotting order
df$eGFR_group_eng <- factor(df$eGFR_group_eng, levels = c("eGFR ≥ 60", "eGFR < 60"))

# Calculate N for each group to include in the table headers later
n_ge_60 <- sum(df$eGFR_group_eng == "eGFR ≥ 60", na.rm = TRUE)
n_lt_60 <- sum(df$eGFR_group_eng == "eGFR < 60", na.rm = TRUE)

# Define the continuous variables to analyze (Updated to β₂-microhlob)
cont_vars <- c("Years", "М-prot g/l", "β₂-microhlob", "CREA...44", 
               "Urea", "Са", "Hb", "Total TPROP", "alb", "CRP", "NGAL", "FGF23", "ACR")

# Define variable titles with standard measurement units for plots
var_titles <- c(
  "CREA...44"   = "Creatinine (µmol/L)",
  "CRP"         = "C-Reactive Protein (mg/L)",
  "Hb"          = "Hemoglobin (g/L)",
  "М-prot g/l"  = "M-protein (g/L)",
  "NGAL"        = "NGAL (ng/mL)",
  "Urea"        = "Urea (mmol/L)",
  "Years"       = "Age (years)",
  "β₂-microhlob"= "Beta-2 microglobulin (mg/L)",
  "Са"          = "Calcium (mmol/L)",
  "Total TPROP" = "Total Protein (g/L)",
  "alb"         = "Albumin (g/L)",
  "FGF23"       = "FGF-23 (RU/mL)",
  "ACR"         = "ACR (mg/mmol)"
)

# ------------------------------------------------------------------------------
# 2. Compute Baseline Statistics (Table 1)
# ------------------------------------------------------------------------------
table_results <- data.frame(
  Variable = character(), 
  `eGFR_ge_60` = character(), 
  `eGFR_lt_60` = character(), 
  p_value = numeric(), 
  Test_Used = character(), 
  stringsAsFactors = FALSE
)

for (var in cont_vars) {
  
  # ЗАХИСНИЙ МЕХАНІЗМ: Перевіряємо, чи існує колонка в базі
  if (!var %in% colnames(df)) {
    message(paste("Warning: Column '", var, "' not found in dataset. Skipping...", sep=""))
    next # Пропускаємо цю змінну і йдемо до наступної
  }
  
  temp_data <- df[!is.na(df[[var]]) & !is.na(df$eGFR_group_eng), ]
  g1 <- temp_data[[var]][temp_data$eGFR_group_eng == "eGFR ≥ 60"]
  g2 <- temp_data[[var]][temp_data$eGFR_group_eng == "eGFR < 60"]
  
  format_med_iqr <- function(x) {
    if(length(x) == 0) return("NA")
    paste0(round(median(x), 1), " (", round(quantile(x, 0.25), 1), "-", round(quantile(x, 0.75), 1), ")")
  }
  
  val1 <- format_med_iqr(g1)
  val2 <- format_med_iqr(g2)
  
  shap1 <- ifelse(length(g1) >= 3, shapiro.test(g1)$p.value, 0)
  shap2 <- ifelse(length(g2) >= 3, shapiro.test(g2)$p.value, 0)
  
  if (shap1 > 0.05 & shap2 > 0.05) {
    if(length(g1) > 1 & length(g2) > 1) {
      p_val <- t.test(g1, g2)$p.value
      test_used <- "Welch Two Sample t-test"
    } else { p_val <- NA; test_used <- "NA" }
  } else {
    if(length(g1) > 0 & length(g2) > 0) {
      p_val <- wilcox.test(g1, g2, exact = FALSE)$p.value
      test_used <- "Mann-Whitney U test"
    } else { p_val <- NA; test_used <- "NA" }
  }
  
  table_results <- rbind(table_results, data.frame(
    Variable = var, 
    eGFR_ge_60 = val1, 
    eGFR_lt_60 = val2,
    p_value = p_val, 
    Test_Used = test_used
  ))
}

# Додаємо (n) у назви колонок
names(table_results)[names(table_results) == "eGFR_ge_60"] <- paste0("eGFR ≥ 60 (n=", n_ge_60, ")")
names(table_results)[names(table_results) == "eGFR_lt_60"] <- paste0("eGFR < 60 (n=", n_lt_60, ")")

# ------------------------------------------------------------------------------
# 3. Export Tabular Data
# ------------------------------------------------------------------------------
write.csv(table_results, "results/table1_baseline_statistics.csv", row.names = FALSE)
write_xlsx(table_results, "results/table1_baseline_statistics.xlsx")

message("Success: Baseline statistics saved to 'results/'.")

# ------------------------------------------------------------------------------
# 4. Generate Advanced Boxplots (Excluding Creatinine and Urea)
# ------------------------------------------------------------------------------
significant_vars <- table_results$Variable[table_results$p_value < 0.05 & 
                                             !is.na(table_results$p_value) & 
                                             !(table_results$Variable %in% c("CREA...44", "Urea"))]

if (length(significant_vars) > 0) {
  plot_list <- list()
  for (var in significant_vars) {
    p_val <- table_results$p_value[table_results$Variable == var]
    p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", format(round(p_val, 3), nsmall = 3)))
    
    df_plot <- df[!is.na(df[[var]]) & !is.na(df$eGFR_group_eng), ]
    
    p <- ggplot(df_plot, aes(x = eGFR_group_eng, y = .data[[var]], fill = eGFR_group_eng)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black", linewidth = 0.5) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black", shape = 16) +
      scale_fill_manual(values = c("#2E9FDF", "#E7B800")) + 
      theme_classic(base_size = 12) +
      labs(title = var_titles[var], subtitle = p_label, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#d73027", face = "italic"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.line = element_line(color = "black")
      )
    
    # Застосовуємо логарифмічну шкалу до вкрай асиметричних даних
    if (var %in% c("CRP", "ACR", "β₂-microhlob", "М-prot g/l")) {
      p <- p + scale_y_log10() + labs(y = "(log10 scale)") + 
        theme(axis.title.y = element_text(size = 8, color = "gray40", angle = 90))
    }
    plot_list[[var]] <- p
  }
  
  final_plot <- wrap_plots(plot_list, ncol = 3)
  
  # Зберігаємо графіки
  ggsave("results/significant_biomarkers_advanced.pdf", plot = final_plot, width = 11, height = 8, dpi = 300)
  ggsave("results/significant_biomarkers_advanced.png", plot = final_plot, width = 11, height = 8, dpi = 300)
  
  message("Success: Advanced multi-panel figure saved to 'results/'.")
  print(final_plot)
} else {
  message("Warning: No statistically significant variables found to generate plots.")
}