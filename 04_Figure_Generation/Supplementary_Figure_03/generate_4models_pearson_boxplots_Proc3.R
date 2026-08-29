# ==============================================================================
# SCRIPT: generate_4models_pearson_boxplots_Proc3.R
# PURPOSE: 4 Models Comparison (Baseline, Mastery 50, Protein, K9 (RF + RF))
#          Pearson Boxplots for Procedure 3 with EXACT Proc 1 Layout & Styling
#
# SPECIFICATIONS (IDENTICAL TO PROC 1):
# 1. Procedure 3 (Merged Batch-Corrected Overlap)
# 2. 4 Models: Baseline, Mastery 50, Protein, K9 (RF + RF)
# 3. Palette:
#    - Baseline:      #E5E7E9 (Light Grey)
#    - Mastery 50:    #85929E (Medium Grey)
#    - Protein:       #2C3E50 (Dark Grey)
#    - K9 (RF + RF):  #E67E22 (Orange)
# 4. Error bars on boxplots (stat_boxplot errorbar)
# 5. Jitter points in #154360 with alpha = 0.35, size = 1.3
# 6. Dashed threshold lines:
#    - Green line at y = 0.8 (#27AE60)
#    - Red line at y = 0.5 (#E74C3C)
# 7. Y-axis:
#    - Bravais-Pearson Correlation (rho_BP)
#    - breaks: seq(-1.0, 1.0, by = 0.25)
#    - coord_cartesian(ylim = c(-1.36, 1.05), clip = "off")
# 8. Clean text annotations at y = -1.24 in bold #1B4F72
# 9. PDF Dimensions: width = 16.5, height = 8.5
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE DIRECTORIES & PATHS
# ------------------------------------------------------------------------------
output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prime_dir      <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP/"
m50_dir        <- paste0(prime_dir, "Mastery50/")
prime_rfrf_dir <- paste0(prime_dir, "RF_RF/")

dipa_groups <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")
group_codes <- c("1_2", "3_4", "5_6", "0")
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
# STEP 3: LOAD DATA FOR PROCEDURE 3
# ------------------------------------------------------------------------------
prime_files <- sprintf(file.path(prime_dir, "Models_rf_preselection_full_objects_mergedata_blind_batch_prime_over_%s_24.RData"), group_codes)
m50_files   <- sprintf(file.path(m50_dir, "Models_mastery50_mergedata_blind_batch_over_%s_24.RData"), group_codes)
rfrf_files  <- sprintf(file.path(prime_rfrf_dir, "Models_rf_rf_preselection_full_objects_mergedata_blind_batch_prime_over_%s_24.RData"), group_codes)

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

dt_proc3_4m <- rbindlist(dt_list, fill = TRUE)
dt_proc3_4m$DiPa_Group <- factor(dt_proc3_4m$DiPa_Group, levels = dipa_groups)
dt_proc3_4m$Model      <- factor(dt_proc3_4m$Model, levels = model_names)

# ------------------------------------------------------------------------------
# STEP 4: PLOTTING FUNCTION WITH EXACT PROC 1 STYLING & OPERATORS
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
# STEP 5: SAVE PDF & PNG PREVIEWS
# ------------------------------------------------------------------------------
title_proc3 <- expression(paste("BDL + CCl"[4], ": (Subset)"))
pdf_proc3   <- paste0(output_dir, "Procedure_3_4Models_Pearson_Merged_Batch.pdf")

pdf(pdf_proc3, width = 16.5, height = 8.5)
print(create_4models_boxplots(dt_proc3_4m, title_proc3))
dev.off()
cat("Procedure 3 4-Models PDF successfully generated:", pdf_proc3, "\n")

# png_proc3 <- paste0(output_dir, "Procedure_3_4Models_Pearson_Merged_Batch.png")
# png(png_proc3, width = 3300, height = 1700, res = 200)
# print(create_4models_boxplots(dt_proc3_4m, title_proc3))
# dev.off()

cat("\nALL 4-MODEL PROCEDURE 3 BOXPLOTS COMPLETED WITH EXACT PROC 1 STYLING!\n")
