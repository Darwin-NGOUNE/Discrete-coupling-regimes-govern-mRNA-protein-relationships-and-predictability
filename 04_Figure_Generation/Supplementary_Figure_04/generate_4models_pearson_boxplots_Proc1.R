# ==============================================================================
# SCRIPT: generate_3models_pearson_boxplots_Proc1.R
# PURPOSE: 3 Models Comparison (Baseline, Mastery 50, Protein) Pearson Boxplots for Proc 1
#
# SPECIFICATIONS (UNIFORM FONT & CLEAR ASCII >= OPERATOR):
# 1. Procedure 1 (ComBat Batch-Corrected)
# 2. 3 Models: Baseline, Mastery 50, Protein
# 3. 3 Shades of Grey (Gris clair, Gris sombre, Gris tres tres sombre)
# 4. Pearson Correlation ONLY with expression(rho[BP]) on Y-axis
# 5. Transparent jittered points overlaid on boxplots inside a clear Gitter (Grid)
# 6. CCl4 formatted with subscript 4 (CCl[4]) in plot titles
# 7. X-axis: 4 DiPa groups side-by-side (DiPa 1 & 2, DiPa 3 & 4, DiPa 5 & 6, DiPa 8)
# 8. For each DiPa group, 3 model boxplots side-by-side (dodged)
# 9. Clean ASCII ">= 0.8: XX%\n>= 0.5: YY%" placed AT THE BOTTOM under each boxplot
# 10. UNIFORM FONT FAMILY ('sans') across all text elements for maximum readability
# 11. TWO SEPARATE PDF FILES:
#    - PDF 1 (Richtung 1): Train BDL / Test CCl4
#    - PDF 2 (Richtung 2): Train CCl4 / Test BDL
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE DIRECTORIES & PATHS
# ------------------------------------------------------------------------------
output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prime_dir     <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/"
m50_dir       <- paste0(prime_dir, "Mastery50/")
prime_rfrf_dir <- paste0(prime_dir, "RF_RF/")

dipa_groups <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")
model_names <- c("Baseline", "Mastery 50", "Protein", "Protein (RF + RF)")

# ------------------------------------------------------------------------------
# STEP 2: HELPER TO EXTRACT PEARSON CORRELATION
# ------------------------------------------------------------------------------
extract_pearson <- function(model_obj) {
  if (is.null(model_obj) || is.character(model_obj)) return(NA_real_)
  val <- model_obj$correlation.pearson
  if (is.null(val)) val <- model_obj$Pearson
  if (is.null(val)) val <- model_obj$pearson
  if (!is.null(val)) return(as.numeric(val))
  return(NA_real_)
}

extract_model_data <- function(filepath, model_type, dipa_label) {
  if (!file.exists(filepath)) return(data.table())
  env <- new.env()
  load(filepath, envir = env)
  vars <- ls(env)
  
  obj_target <- NULL
  for (v in vars) {
    v_low <- tolower(v)
    if (model_type == "Baseline" && grepl("rna", v_low) && !grepl("plus_rna|combined|mastery", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "Protein" && grepl("protein", v_low) && !grepl("plus_rna|combined|mastery", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "Mastery 50" && grepl("protein", v_low) && !grepl("plus_rna|combined", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "Protein (RF + RF)" && grepl("protein", v_low) && !grepl("plus_rna|combined", v_low)) {
      obj_target <- env[[v]]
    }
  }
  if (is.null(obj_target)) return(data.table())
  
  proteins <- names(obj_target)
  if (is.null(proteins)) proteins <- paste0("Prot_", seq_along(obj_target))
  
  rows <- list()
  for (i in seq_along(obj_target)) {
    x <- obj_target[[i]]
    node <- if (!is.null(x$prediction.obj)) x$prediction.obj else x
    
    if (model_type == "Baseline") {
      base_m <- node$baseline.model
      if (!is.null(base_m)) {
        val_p <- extract_pearson(base_m)
        if (!is.na(val_p)) {
          rows[[length(rows) + 1]] <- data.table(
            Protein = proteins[i], Pearson = val_p, Model = "Baseline", DiPa_Group = dipa_label
          )
        }
      }
    } else {
      lasso_m <- node$model
      if (!is.null(lasso_m)) {
        val_p <- extract_pearson(lasso_m)
        if (!is.na(val_p)) {
          rows[[length(rows) + 1]] <- data.table(
            Protein = proteins[i], Pearson = val_p, Model = model_type, DiPa_Group = dipa_label
          )
        }
      }
    }
  }
  return(rbindlist(rows, fill = TRUE))
}

# ------------------------------------------------------------------------------
# STEP 3: LOAD DATA FOR RICHTUNG 1 AND RICHTUNG 2
# ------------------------------------------------------------------------------

# Richtung 1 (Train BDL / Test CCl4)
r1_prime_files <- c(
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_1_2_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_3_4_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_5_6_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_new_Batch_prime_0_sub.RData")
)

r1_m50_files <- c(
  paste0(m50_dir, "Models_mastery50_BDL_CCL4_1_2_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_BDL_CCL4_3_4_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_BDL_CCL4_5_6_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_BDL_CCL4_0_prime_sub.RData")
)

r1_rfrf_files <- c(
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_new_Batch_prime_1_2_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_new_Batch_prime_3_4_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_new_Batch_prime_5_6_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_new_Batch_prime_0_sub.RData")
)

# Richtung 2 (Train CCl4 / Test BDL)
r2_prime_files <- c(
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_1_2_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_3_4_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_5_6_sub.RData"),
  paste0(prime_dir, "Models_rf_lasso_full_testing_Batch_prime_0_sub.RData")
)

r2_m50_files <- c(
  paste0(m50_dir, "Models_mastery50_CCL4_BDL_1_2_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_CCL4_BDL_3_4_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_CCL4_BDL_5_6_prime_sub.RData"),
  paste0(m50_dir, "Models_mastery50_CCL4_BDL_0_prime_sub.RData")
)

r2_rfrf_files <- c(
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_Batch_prime_1_2_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_Batch_prime_3_4_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_Batch_prime_5_6_sub.RData"),
  paste0(prime_rfrf_dir, "Models_rf_rf_rf_full_testing_Batch_prime_0_sub.RData")
)

load_richtung_data <- function(prime_files, m50_files, rfrf_files) {
  dt_list <- list()
  for (i in 1:4) {
    # Baseline
    dt_base <- extract_model_data(prime_files[i], "Baseline", dipa_groups[i])
    if (nrow(dt_base) > 0) dt_list[[length(dt_list) + 1]] <- dt_base
    
    # Mastery 50
    dt_m50  <- extract_model_data(m50_files[i], "Mastery 50", dipa_groups[i])
    if (nrow(dt_m50) > 0) dt_list[[length(dt_list) + 1]] <- dt_m50
    
    # Protein
    dt_prot <- extract_model_data(prime_files[i], "Protein", dipa_groups[i])
    if (nrow(dt_prot) > 0) dt_list[[length(dt_list) + 1]] <- dt_prot
    
    # Protein (RF + RF)
    dt_rfrf <- extract_model_data(rfrf_files[i], "Protein (RF + RF)", dipa_groups[i])
    if (nrow(dt_rfrf) > 0) dt_list[[length(dt_list) + 1]] <- dt_rfrf
  }
  res <- rbindlist(dt_list, fill = TRUE)
  res$DiPa_Group <- factor(res$DiPa_Group, levels = dipa_groups)
  res$Model      <- factor(res$Model, levels = model_names)
  return(res)
}

dt_r1 <- load_richtung_data(r1_prime_files, r1_m50_files, r1_rfrf_files)
dt_r2 <- load_richtung_data(r2_prime_files, r2_m50_files, r2_rfrf_files)

# ------------------------------------------------------------------------------
# STEP 4: PLOTTING FUNCTION WITH UNIFORM FONT & CLEAN OPERATORS
# ------------------------------------------------------------------------------
create_4models_boxplots <- function(data_dt, plot_title_expr) {
  
  summary_dt <- data_dt[!is.na(Pearson), .(
    n_total  = .N,
    pct_high = round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100),
    pct_suff = round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100),
    label    = paste0(">= 0.8: ", round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100), "%\n>= 0.5: ", round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100), "%")
  ), by = .(DiPa_Group, Model)]
  
  palette_fill <- c(
    "Baseline"          = "#E5E7E9", # Light Grey
    "Mastery 50"        = "#85929E", # Medium Grey
    "Protein"           = "#2C3E50", # Dark Grey
    "Protein (RF + RF)" = "#1B2631"  # Very Dark Grey
  )
  
  p <- ggplot(data_dt, aes(x = DiPa_Group, y = Pearson, fill = Model)) +
    stat_boxplot(geom = "errorbar", position = position_dodge(0.8), width = 0.25, color = "#1A252F", linewidth = 0.6) +
    geom_boxplot(position = position_dodge(0.8), color = "#1A252F", alpha = 0.85, outlier.shape = NA, width = 0.62) +
    
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               color = "#154360", alpha = 0.35, size = 1.3) +
    
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#E74C3C", linewidth = 0.85) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#27AE60", linewidth = 0.85) +
    
    geom_text(data = summary_dt, 
              aes(x = DiPa_Group, y = -1.24, label = label, group = Model),
              position = position_dodge(width = 0.8),
              family = "sans", size = 3.2, color = "#1B4F72", fontface = "bold", lineheight = 0.88) +
    
    scale_fill_manual(values = palette_fill) +
    
    scale_y_continuous(
      name = expression(paste("Bravais-Pearson Correlation (", rho[BP], ")")),
      breaks = seq(-1.0, 1.0, by = 0.25),
      limits = c(-1.42, 1.05)
    ) +
    coord_cartesian(ylim = c(-1.36, 1.05), clip = "off") +
    
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 17, margin = margin(b = 15)),
      plot.margin = margin(t = 20, r = 20, b = 50, l = 20),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(family = "sans", face = "bold", size = 13.5, color = "#2C3E50", margin = margin(t = 5)),
      axis.title.y = element_text(family = "sans", face = "bold", size = 14, color = "#2C3E50"),
      axis.text.y  = element_text(family = "sans", size = 12, color = "#2C3E50"),
      
      legend.position = "top",
      legend.title    = element_text(family = "sans", face = "bold", size = 13),
      legend.text     = element_text(family = "sans", size = 12),
      
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25, linetype = "dashed"),
      panel.border     = element_rect(color = "#2C3E50", fill = NA, linewidth = 1.0)
    ) +
    labs(title = plot_title_expr, fill = "Model:")
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 5: SAVE THE TWO SEPARATE PDF FILES WITH FORMATTED CCL4 TITLES
# ------------------------------------------------------------------------------

# PDF 1: Richtung 1 (Train BDL / Test CCl4)
title_r1 <- expression(paste("Train BDL / Test CCl"[4], " : (Subset)"))
pdf_r1   <- paste0(output_dir, "Procedure_1_4Models_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf")

pdf(pdf_r1, width = 16.5, height = 8.5)
print(create_4models_boxplots(dt_r1, title_r1))
dev.off()
cat("4 Models PDF 1 successfully generated:", pdf_r1, "\n")

# PDF 2: Richtung 2 (Train CCl4 / Test BDL)
title_r2 <- expression(paste("Train CCl"[4], " / Test BDL : (Subset)"))
pdf_r2   <- paste0(output_dir, "Procedure_1_4Models_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf")

pdf(pdf_r2, width = 16.5, height = 8.5)
print(create_4models_boxplots(dt_r2, title_r2))
dev.off()
cat("4 Models PDF 2 successfully generated:", pdf_r2, "\n")


