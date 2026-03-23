# ==============================================================================
# Script: 03_cutoff_analysis.R
# Description: Determination of optimal cutoffs for NGAL & FGF23.
#              Fixed for Ukrainian Excel compatibility and ggplot2 errors.
# ==============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(writexl) # Для збереження в українському форматі (xlsx)

# 1. Завантаження даних
# Використовуємо .xlsx, щоб уникнути проблем із розбиттям колонок
data_path <- "data/processed/clean_data.xlsx"

if (file.exists(data_path)) {
  df <- read_excel(data_path)
  message("Success: Clean dataset (Excel) loaded.")
} else {
  stop("Error: Clean Excel file not found in 'data/processed/'.")
}

# 2. Розрахунок точок відсічення
res.cut <- surv_cutpoint(
  df, 
  time = "тривалість", 
  event = "Exitus", 
  variables = c("NGAL", "FGF23"), 
  minprop = 0.1 
)

# Отримуємо таблицю результатів
cutoff_summary <- as.data.frame(summary(res.cut))
print(cutoff_summary)

# Зберігаємо результати саме в .xlsx (український формат Excel)
write_xlsx(cutoff_summary, "results/optimal_cutoffs_summary.xlsx")
message("Success: Cutoff summary saved as .xlsx")

# ------------------------------------------------------------------------------
# 3. Побудова графіків (Виправлено помилку fortify)
# ------------------------------------------------------------------------------
message("Generating statistical cutoff plots...")

# Функція для створення графіка, яка ГАРАНТОВАНО працює з ggplot2
make_cutoff_plot <- function(res_cut_obj, var_name, title_name, unit) {
  
  # ВИПРАВЛЕННЯ: Дістаємо дані та ПРИМУСОВО робимо їх простою таблицею
  # В об'єкті res.cut дані лежать у списку під назвою змінної
  raw_data <- as.data.frame(res_cut_obj[[var_name]])
  
  # Отримуємо значення точки відсічення
  cutoff_val <- summary(res_cut_obj)[var_name, "cutpoint"]
  
  ggplot(raw_data, aes(x = cuts, y = stats)) +
    geom_line(color = "#2E9FDF", linewidth = 1) +
    geom_vline(xintercept = cutoff_val, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    annotate("text", x = cutoff_val, y = max(raw_data$stats, na.rm = TRUE), 
             label = paste("Cutoff:", round(cutoff_val, 1)), 
             vjust = -0.5, color = "#E74C3C", fontface = "bold") +
    labs(
      title = title_name,
      x = paste(var_name, unit),
      y = "Log-rank Statistic"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10)
    )
}

# Створюємо графіки без жодних помилок
p_ngal  <- make_cutoff_plot(res.cut, "NGAL", "NGAL Cutoff", "(ng/mL)")
p_fgf23 <- make_cutoff_plot(res.cut, "FGF23", "FGF-23 Cutoff", "(RU/mL)")

# Об'єднуємо та зберігаємо
final_fig <- p_ngal + p_fgf23

ggsave("results/cutoff_evidence_final.png", plot = final_fig, width = 10, height = 5, dpi = 300)
# Зберігаємо також PDF для високої якості друку
pdf("results/cutoff_evidence_final.pdf", width = 10, height = 5)
print(final_fig)
dev.off()

message("Success: Clean plots generated and saved to 'results/'.")
print(final_fig)