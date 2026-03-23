# ==============================================================================
# Script: 04_survival_analysis.R
# Description: Generates publication-ready Kaplan-Meier survival curves (Overall, 
#              NGAL, FGF23, Composite Score) and performs Cox Proportional 
#              Hazards Regression. Exports plots and summary tables.
#              (Updated: Added multivariate Cox model with eGFR AND ISS stage adjustment)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(writexl)

# ------------------------------------------------------------------------------
# 1. Load Data & Prepare Risk Groups
# ------------------------------------------------------------------------------
df <- read_excel("data/processed/clean_data.xlsx")

# Розподіляємо пацієнтів на групи
df <- df %>%
  mutate(
    NGAL_Group = ifelse(NGAL > 153.8, "High", "Low"),
    FGF23_Group = ifelse(FGF23 > 36.4, "High", "Low"),
    
    # Бали для фінальної моделі (1 бал за кожен перевищений маркер)
    NGAL_Score = ifelse(NGAL > 153.8, 1, 0),
    FGF23_Score = ifelse(FGF23 > 36.4, 1, 0),
    Raw_Score = NGAL_Score + FGF23_Score,
    
    # ОБ'ЄДНАННЯ ГРУП: Зводимо 0 та 1 бал в одну категорію "Score 0-1"
    Composite_Score = ifelse(Raw_Score <= 1, "Score 0-1", "Score 2")
  )

# Перетворюємо в фактори ("Low" або "Score 0-1" завжди референсна група)
df$NGAL_Group <- factor(df$NGAL_Group, levels = c("Low", "High"))
df$FGF23_Group <- factor(df$FGF23_Group, levels = c("Low", "High"))
df$Composite_Score <- factor(df$Composite_Score, levels = c("Score 0-1", "Score 2"))

# --- ДОДАНО: Перетворюємо стадію ISS на фактор ---
# Це критично важливо, щоб R порівнював 2 і 3 стадію відносно 1-ї (референсної)
df$`stage ISS` <- factor(df$`stage ISS`, levels = c(1, 2, 3))

message("Data categorized. Proceeding to Survival Analysis...")

# ------------------------------------------------------------------------------
# 2. Univariate Cox Regression Models
# ------------------------------------------------------------------------------
extract_cox_results <- function(model, model_name) {
  sum_mod <- summary(model)
  
  clean_vars <- gsub("GroupHigh", " (High vs Low)", rownames(sum_mod$conf.int))
  clean_vars <- gsub("Composite_Score", "", clean_vars)
  
  data.frame(
    Model = model_name,
    Variable = clean_vars,
    Hazard_Ratio = round(sum_mod$conf.int[, 1], 2),
    CI_95_Lower = round(sum_mod$conf.int[, 3], 2),
    CI_95_Upper = round(sum_mod$conf.int[, 4], 2),
    p_value = round(sum_mod$coefficients[, 5], 4),
    row.names = NULL
  )
}

cox_ngal <- coxph(Surv(тривалість, Exitus) ~ NGAL_Group, data = df)
cox_fgf23 <- coxph(Surv(тривалість, Exitus) ~ FGF23_Group, data = df)
cox_score <- coxph(Surv(тривалість, Exitus) ~ Composite_Score, data = df)

cox_results_table <- rbind(
  extract_cox_results(cox_ngal, "NGAL Only"),
  extract_cox_results(cox_fgf23, "FGF23 Only"),
  extract_cox_results(cox_score, "Composite risk score")
)

write_xlsx(cox_results_table, "results/Cox_Regression_Results.xlsx")
message("Success: Univariate Cox regression table saved to 'results/Cox_Regression_Results.xlsx'")

# ------------------------------------------------------------------------------
# 2.5. Multivariate Cox Regression Models (Separate Models)
# ------------------------------------------------------------------------------
message("Running multivariate Cox regression models...")

# Модель 1: Коригування на eGFR (як було в оригіналі)
fit_multi_egfr <- coxph(
  formula = Surv(тривалість, Exitus) ~ Composite_Score + `GFR CKD-EPI`, 
  data = df
)

# Модель 2: Коригування на стадію ISS
fit_multi_iss <- coxph(
  formula = Surv(тривалість, Exitus) ~ Composite_Score + `stage ISS`, 
  data = df
)

# Зберігаємо результати обох моделей у зручний текстовий файл
sink("results/Multivariate_Cox_Summary.txt")

cat("=== MULTIVARIATE COX REGRESSION SUMMARY ===\n\n")

cat("--- MODEL 1: ADJUSTED FOR eGFR ---\n")
print(summary(fit_multi_egfr))
cat("\n-------------------------------------------\n\n")

cat("--- MODEL 2: ADJUSTED FOR ISS STAGE ---\n")
print(summary(fit_multi_iss))
cat("\n===========================================\n")

sink()

message("Success: Multivariate Cox summaries saved to 'results/Multivariate_Cox_Summary.txt'")

# ------------------------------------------------------------------------------
# 3. Kaplan-Meier Curves (Publication Quality)
# ------------------------------------------------------------------------------
message("Generating clean Kaplan-Meier plots...")

km_theme <- theme_survminer()

# --- Графік 1: Загальна виживаність ---
fit_overall <- survfit(Surv(тривалість, Exitus) ~ 1, data = df)
p_overall <- ggsurvplot(
  fit_overall, data = df,
  title = "", 
  xlab = "Follow-up time (months)", ylab = "Survival probability",
  palette = "#2C3E50", conf.int = TRUE, ggtheme = km_theme
)

# --- Графік 2: Виживаність за NGAL ---
fit_ngal <- survfit(Surv(тривалість, Exitus) ~ NGAL_Group, data = df)
p_ngal <- ggsurvplot(
  fit_ngal, data = df,
  title = "", 
  pval = TRUE, risk.table = TRUE,
  risk.table.title = "Number at risk",
  xlab = "Follow-up time (months)", ylab = "Survival probability",
  palette = c("#2E9FDF", "#E74C3C"), 
  legend.title = "", 
  legend.labs = c("NGAL ↓", "NGAL ↑"), 
  ggtheme = km_theme
)

# --- Графік 3: Виживаність за FGF-23 ---
fit_fgf23 <- survfit(Surv(тривалість, Exitus) ~ FGF23_Group, data = df)
p_fgf23 <- ggsurvplot(
  fit_fgf23, data = df,
  title = "", 
  pval = TRUE, risk.table = TRUE,
  risk.table.title = "Number at risk",
  xlab = "Follow-up time (months)", ylab = "Survival probability",
  palette = c("#2E9FDF", "#E74C3C"),
  legend.title = "",
  legend.labs = c("FGF-23 ↓", "FGF-23 ↑"),
  ggtheme = km_theme
)

# --- Графік 4: Фінальна виживаність за Балами ---
fit_score <- survfit(Surv(тривалість, Exitus) ~ Composite_Score, data = df)
p_score <- ggsurvplot(
  fit_score, data = df,
  title = "", 
  pval = TRUE, risk.table = TRUE,
  risk.table.title = "Number at risk",
  xlab = "Follow-up time (months)", ylab = "Survival probability",
  palette = c("#2E9FDF", "#E74C3C"), 
  legend.title = "",
  legend.labs = c("Score: 0-1", "Score: 2"), 
  ggtheme = km_theme
)

# ------------------------------------------------------------------------------
# 3.5. Extract Exact Overall Survival Statistics (Median, 3-year, 5-year)
# ------------------------------------------------------------------------------
message("Extracting precise survival milestones for the manuscript...")

surv_milestones <- summary(fit_overall, times = c(36, 60))

os_rates <- data.frame(
  Time_Months = surv_milestones$time,
  Time_Years = surv_milestones$time / 12,
  Survival_Percent = round(surv_milestones$surv * 100, 1),
  CI_95_Lower = round(surv_milestones$lower * 100, 1),
  CI_95_Upper = round(surv_milestones$upper * 100, 1)
)

sink("results/Overall_Survival_Milestones.txt")
cat("=== OVERALL SURVIVAL SUMMARY FOR THE MANUSCRIPT ===\n\n")
cat("1. MEDIAN SURVIVAL (Months):\n")
print(fit_overall)
cat("\n--------------------------------------------------\n")
cat("2. SURVIVAL RATES (3-year and 5-year):\n")
print(os_rates)
cat("\n==================================================\n")
sink()

message("Success: Overall Survival milestones saved to 'results/Overall_Survival_Milestones.txt'")

# ------------------------------------------------------------------------------
# 4. Save Plots Securely to Files (PDF/PNG)
# ------------------------------------------------------------------------------
save_km_plot <- function(plot_obj, file_prefix) {
  pdf(paste0("results/", file_prefix, ".pdf"), width = 8, height = 7, onefile = FALSE)
  print(plot_obj)
  dev.off()
  
  png(paste0("results/", file_prefix, ".png"), width = 800, height = 700, res = 120)
  print(plot_obj)
  dev.off()
}

save_km_plot(p_overall, "KM_01_Overall_Survival")
save_km_plot(p_ngal, "KM_02_NGAL_Survival")
save_km_plot(p_fgf23, "KM_03_FGF23_Survival")
save_km_plot(p_score, "KM_04_Composite_Score_Survival")

message("Success: All clean Kaplan-Meier survival curves saved to 'results/'!")