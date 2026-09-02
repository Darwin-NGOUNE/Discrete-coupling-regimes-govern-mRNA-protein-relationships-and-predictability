# ==============================================================================
# SCRIPT: generate_baseline_pearson_boxplots_Proc1.R
# PURPOSE: Baseline Model Pearson Correlation Boxplots with Transparent Dots for Proc 1
#
# SPECIFICATIONS (UPDATED):
# 1. Procedure 1 (ComBat Batch-Corrected)
# 2. Baseline Model ONLY (Condition: "a. Baseline")
# 3. Pearson Correlation ONLY with expression(rho[BP]) on Y-axis
# 4. Transparent jittered points overlaid on boxplots inside a clear Gitter (Grid)
# 5. CCl4 formatted with subscript 4 (CCl[4]) in plot titles
# 6. X-axis: 4 DiPa groups side-by-side (DiPa 1 & 2, DiPa 3 & 4, DiPa 5 & 6, DiPa 8)
# 7. Percentages (>=0.8: XX%, >=0.5: YY%) placed AT THE BOTTOM of each DiPa box
# 8. TWO SEPARATE PDF FILES:
#    - PDF 1 (Richtung 1): Train BDL / Test CCl4
#    - PDF 2 (Richtung 2): Train CCl4 / Test BDL
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE OUTPUT DIRECTORY & PATHS
# ------------------------------------------------------------------------------
output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prime_dir   <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/"
dipa_groups <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")

# ------------------------------------------------------------------------------
# STEP 2: HELPER TO EXTRACT PEARSON CORRELATION FOR BASELINE MODEL
# ------------------------------------------------------------------------------
extract_pearson <- function(model_obj) {
  if (is.null(model_obj) || is.character(model_obj)) return(NA_real_)
  val <- model_obj$correlation.pearson
  if (is.null(val)) val <- model_obj$Pearson
  if (is.null(val)) val <- model_obj$pearson
  if (!is.null(val)) return(as.numeric(val))
  return(NA_real_)
}

extract_baseline_pearson <- function(filepath, dipa_label) {
  if (!file.exists(filepath)) return(data.table())
  env <- new.env()
  load(filepath, envir = env)
  vars <- ls(env)
  
  obj_rna <- NULL
  for (v in vars) {
    if (grepl("rna", tolower(v))) obj_rna <- env[[v]]
  }
  if (is.null(obj_rna)) return(data.table())
  
  proteins <- names(obj_rna)
  if (is.null(proteins)) proteins <- paste0("Prot_", seq_along(obj_rna))
  
  rows <- list()
  for (i in seq_along(obj_rna)) {
    x <- obj_rna[[i]]
    node <- if (!is.null(x$prediction.obj)) x$prediction.obj else x
    base_m <- node$baseline.model
    if (!is.null(base_m)) {
      val_pearson <- extract_pearson(base_m)
      if (!is.na(val_pearson)) {
        rows[[length(rows) + 1]] <- data.table(
          Protein    = proteins[i],
          Pearson    = val_pearson,
          DiPa_Group = dipa_label
        )
      }
    }
  }
  return(rbindlist(rows, fill = TRUE))
}

# ------------------------------------------------------------------------------
# STEP 3: LOAD BASELINE DATA FOR RICHTUNG 1 AND RICHTUNG 2
# ------------------------------------------------------------------------------
prime_files_r1 <- c(
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_1_2_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_3_4_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_5_6_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_0_sub.RData")
)

prime_files_r2 <- c(
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_1_2_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_3_4_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_5_6_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_0_sub.RData")
)

list_r1 <- list()
for (i in 1:4) {
  dt <- extract_baseline_pearson(prime_files_r1[i], dipa_groups[i])
  if (nrow(dt) > 0) list_r1[[length(list_r1) + 1]] <- dt
}
dt_r1 <- rbindlist(list_r1, fill = TRUE)
dt_r1$DiPa_Group <- factor(dt_r1$DiPa_Group, levels = dipa_groups)

list_r2 <- list()
for (i in 1:4) {
  dt <- extract_baseline_pearson(prime_files_r2[i], dipa_groups[i])
  if (nrow(dt) > 0) list_r2[[length(list_r2) + 1]] <- dt
}
dt_r2 <- rbindlist(list_r2, fill = TRUE)
dt_r2$DiPa_Group <- factor(dt_r2$DiPa_Group, levels = dipa_groups)

# ------------------------------------------------------------------------------
# STEP 4: PLOTTING FUNCTION WITH GITTER, RHO_BP AND CCL4 SUBSCRIPT
# ------------------------------------------------------------------------------
create_baseline_boxplots <- function(data_dt, plot_title_expr) {
  
  # Calculate summary percentages per DiPa group with true mathematical \u2265 symbol
  summary_dt <- data_dt[!is.na(Pearson), .(
    n_total  = .N,
    pct_high = round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100),
    pct_suff = round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100),
    label    = paste0("\u2265 0.8: ", round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100), "%\n\u2265 0.5: ", round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100), "%")
  ), by = DiPa_Group]
  
  p <- ggplot(data_dt, aes(x = DiPa_Group, y = Pearson)) +
    # 1. Errorbars and Boxplots
    stat_boxplot(geom = "errorbar", width = 0.25, color = "black", linewidth = 0.8) +
    geom_boxplot(fill = "#E5E7E9", color = "black", alpha = 0.75, outlier.shape = NA, width = 0.45, linewidth = 0.8) +
    
    # 2. Transparent Jittered Dots
    geom_jitter(color = "#2980B9", alpha = 0.40, size = 2.2, width = 0.18) +
    
    # 3. Reference Dashed Lines at 0.5 (Red) and 0.8 (Green)
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#E74C3C", linewidth = 0.9) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#27AE60", linewidth = 0.9) +
    
    # 4. Vertical separator lines between DiPa groups stopping STRICTLY at y = -1.0
    geom_segment(data = data.frame(x = c(1.5, 2.5, 3.5), xend = c(1.5, 2.5, 3.5), y = -1.0, yend = 1.0),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "solid", color = "grey80", linewidth = 0.6, inherit.aes = FALSE) +
    
    # 5. PERCENTAGE LABELS AT THE BOTTOM OF EACH BOX (sitting cleanly below the y = -1.0 line)
    geom_text(data = summary_dt, aes(x = DiPa_Group, y = -1.22, label = label),
              size = 5.2, color = "black", fontface = "bold", lineheight = 0.92, inherit.aes = FALSE) +
    
    # 6. Y-Axis Label with bold Pearson correlation, rho_BP and limits
    scale_y_continuous(
      name = expression(bold("Pearson correlation, "*rho[BP])),
      breaks = seq(-1.0, 1.0, by = 0.25),
      limits = c(-1.38, 1.05)
    ) +
    coord_cartesian(ylim = c(-1.38, 1.05), clip = "off") +
    
    # 7. Theme with horizontal grid lines, no vertical grid lines cutting into bottom text
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 22, color = "black", margin = margin(b = 15)),
      plot.margin = margin(t = 20, r = 20, b = 25, l = 20),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(face = "bold", size = 16, color = "black"),
      axis.title.y = element_text(face = "bold", size = 18, color = "black"),
      axis.text.y  = element_text(face = "bold", size = 15, color = "black"),
      axis.ticks   = element_line(color = "black", linewidth = 0.9),
      
      # Clear horizontal grid, clean vertical
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25, linetype = "dashed"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1.3)
    ) +
    labs(title = plot_title_expr)
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 5: SAVE THE TWO SEPARATE PDF FILES WITH FORMATTED CCL4 TITLES
# ------------------------------------------------------------------------------

# PDF 1: Richtung 1 (Train BDL, test CCl4)
title_r1 <- expression(bold("Train BDL, test "*CCl[4]))
pdf_r1   <- paste0(output_dir, "Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf")

cairo_pdf(pdf_r1, width = 10, height = 8)
print(create_baseline_boxplots(dt_r1, title_r1))
dev.off()
cat("PDF 1 successfully generated:", pdf_r1, "\n")

# PDF 2: Richtung 2 (Train CCl4, test BDL)
title_r2 <- expression(bold("Train "*CCl[4]*", test BDL"))
pdf_r2   <- paste0(output_dir, "Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf")

cairo_pdf(pdf_r2, width = 10, height = 8)
print(create_baseline_boxplots(dt_r2, title_r2))
dev.off()
cat("PDF 2 successfully generated:", pdf_r2, "\n")
