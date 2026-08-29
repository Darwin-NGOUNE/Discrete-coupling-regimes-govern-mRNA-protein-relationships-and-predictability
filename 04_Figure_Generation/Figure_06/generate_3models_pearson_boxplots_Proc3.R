# ==============================================================================
# SCRIPT: generate_3models_pearson_boxplots_Proc3.R
# PURPOSE: 3 Models Comparison (Baseline, Mastery 50, Protein) Pearson Boxplots for Proc 3
#
# SPECIFICATIONS (TITLE UPDATED TO BDL + CCl4 Batch (Overlap-Subset)):
# 1. Procedure 3 (Merged Batch-Corrected 24h Subset)
# 2. 3 Models: Baseline, Mastery 50, Protein
# 3. 3 Shades of Grey (Gris clair, Gris sombre, Gris tres tres sombre)
# 4. Pearson Correlation ONLY with expression(rho[BP]) on Y-axis
# 5. Transparent jittered points overlaid on boxplots inside a clear Gitter (Grid)
# 6. X-axis: 4 DiPa groups side-by-side (DiPa 1 & 2, DiPa 3 & 4, DiPa 5 & 6, DiPa 8)
# 7. For each DiPa group, 3 model boxplots side-by-side (dodged)
# 8. Clean ASCII ">= 0.8: XX%\n>= 0.5: YY%" placed AT THE BOTTOM under each boxplot
# 9. UNIFORM FONT FAMILY ('sans') across all text elements for maximum readability
# 10. OUTPUT PDF: Procedure_3_3Models_Pearson_Merged_Batch.pdf
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 1: DEFINE DIRECTORIES & PATHS
# ------------------------------------------------------------------------------
output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prime_dir   <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP/"
dipa_groups <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")
group_codes <- c("1_2", "3_4", "5_6", "0")
model_names <- c("Baseline", "Mastery 50", "Protein")

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
m50_files   <- sprintf(file.path(prime_dir, "Mastery50/Models_mastery50_mergedata_blind_batch_over_%s_24.RData"), group_codes)

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
}

dt_proc3_3m <- rbindlist(dt_list, fill = TRUE)
dt_proc3_3m$DiPa_Group <- factor(dt_proc3_3m$DiPa_Group, levels = dipa_groups)
dt_proc3_3m$Model      <- factor(dt_proc3_3m$Model, levels = model_names)

# ------------------------------------------------------------------------------
# STEP 4: PLOTTING FUNCTION FOR 3-MODEL COMPARISON
# ------------------------------------------------------------------------------
create_3models_boxplots <- function(data_dt, plot_title_expr) {
  
  summary_dt <- data_dt[!is.na(Pearson), .(
    n_total  = .N,
    pct_high = round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100),
    pct_suff = round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100),
    label    = paste0(">= 0.8: ", round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100), "%\n>= 0.5: ", round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100), "%")
  ), by = .(DiPa_Group, Model)]
  
  grey_palette_fill <- c("Baseline" = "#E5E7E9", "Mastery 50" = "#85929E", "Protein" = "#2C3E50")
  
  p <- ggplot(data_dt, aes(x = DiPa_Group, y = Pearson, fill = Model)) +
    # 1. Boxplots grouped side-by-side
    stat_boxplot(geom = "errorbar", position = position_dodge(0.8), width = 0.25, color = "#1A252F", linewidth = 0.6) +
    geom_boxplot(position = position_dodge(0.8), color = "#1A252F", alpha = 0.85, outlier.shape = NA, width = 0.65) +
    
    # 2. Transparent Jittered Dots Dodged with Boxplots
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               color = "#154360", alpha = 0.35, size = 1.6) +
    
    # 3. Reference Dashed Lines at 0.5 (Red) and 0.8 (Green)
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#E74C3C", linewidth = 0.85) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#27AE60", linewidth = 0.85) +
    
    # 4. PROPERLY DODGED & UNIFORM FONT PERCENTAGE LABELS AT THE BOTTOM
    geom_text(data = summary_dt, 
              aes(x = DiPa_Group, y = -1.24, label = label, group = Model),
              position = position_dodge(width = 0.8),
              family = "sans", size = 3.7, color = "#1B4F72", fontface = "bold", lineheight = 0.88) +
    
    # 5. Fill Colors (3 Shades of Grey)
    scale_fill_manual(values = grey_palette_fill) +
    
    # 6. Y-Axis Label with expression(rho[BP])
    scale_y_continuous(
      name = expression(paste("Bravais-Pearson Correlation (", rho[BP], ")")),
      breaks = seq(-1.0, 1.0, by = 0.25),
      limits = c(-1.42, 1.05)
    ) +
    coord_cartesian(ylim = c(-1.36, 1.05), clip = "off") +
    
    # 7. Theme with UNIFORM SANS FONT & VISIBLE GRID (GITTER)
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
      
      # Clear Grid / Gitter styling
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5, linetype = "solid"),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25, linetype = "dashed"),
      panel.border     = element_rect(color = "#2C3E50", fill = NA, linewidth = 1.0)
    ) +
    labs(title = plot_title_expr, fill = "Model:")
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 5: SAVE THE PDF FILE FOR PROCEDURE 3 3-MODEL COMPARISON
# ------------------------------------------------------------------------------
title_p3_3m  <- expression(paste("BDL + CCl"[4], ": (Subset)"))
pdf_proc3_3m <- paste0(output_dir, "Procedure_3_3Models_Pearson_Merged_Batch.pdf")

pdf(pdf_proc3_3m, width = 16, height = 8)
print(create_3models_boxplots(dt_proc3_3m, title_p3_3m))
dev.off()
cat("Procedure 3 3-Model PDF successfully generated:", pdf_proc3_3m, "\n")
