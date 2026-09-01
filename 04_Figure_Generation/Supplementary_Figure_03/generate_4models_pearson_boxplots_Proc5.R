# ==============================================================================
# SCRIPT: generate_4models_pearson_boxplots_Proc5.R
# PURPOSE: Generate 4-Models Pearson Correlation Boxplots for Procedure 5 (Intra-Cohort LOOCV)
#          Models: Baseline, Mastery 50, RF + LASSO, RF + RF
#          Panels: Intra BDL (Panel A) and Intra CCl4 (Panel B) for Supplementary Figure 3
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

intra_dir      <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung/Intra_new"
m50_dir        <- file.path(intra_dir, "Mastery50")
intra_rfrf_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung/Intra_RF_RF"

dipa_groups <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")
group_codes <- c("1_2", "3_4", "5_6", "0")
model_names <- c("Baseline", "Mastery 50", "RF + LASSO", "RF + RF")

extract_pearson <- function(model_obj) {
  if (is.null(model_obj) || is.character(model_obj)) return(NA_real_)
  val <- model_obj$correlation.pearson
  if (is.null(val)) val <- model_obj$Pearson
  if (is.null(val)) val <- model_obj$pearson
  if (!is.null(val)) return(as.numeric(val))
  
  y_test <- if (!is.null(model_obj$y_test)) model_obj$y_test else if (!is.null(model_obj$y.t)) model_obj$y.t else model_obj$y
  pred   <- if (!is.null(model_obj$prediction)) model_obj$prediction else model_obj$predictions.insample
  if (!is.null(y_test) && !is.null(pred)) {
    len <- min(length(y_test), length(pred))
    if (len >= 3) {
      v1 <- as.numeric(y_test[1:len])
      v2 <- as.numeric(pred[1:len])
      valid <- which(!is.na(v1) & !is.na(v2))
      if (length(valid) >= 3) return(cor(v1[valid], v2[valid]))
    }
  }
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
    } else if (model_type == "RF + LASSO" && grepl("protein", v_low) && !grepl("plus_rna|combined|mastery|rf_rf|rf\\.rf", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "Mastery 50" && grepl("mastery", v_low) && !grepl("plus_rna|combined", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "RF + RF" && grepl("protein", v_low) && (grepl("rf_rf", v_low) || grepl("rf\\.rf", v_low) || grepl("rf", v_low)) && !grepl("plus_rna|combined|mastery|lasso", v_low)) {
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
      base_m <- if (!is.null(node$baseline.model)) node$baseline.model else node
      val_p  <- extract_pearson(base_m)
      if (!is.na(val_p)) {
        rows[[length(rows) + 1]] <- data.table(
          Protein = proteins[i], Pearson = val_p, Model = "Baseline", DiPa_Group = dipa_label
        )
      }
    } else {
      lasso_m <- if (!is.null(node$model)) node$model else node
      val_p   <- extract_pearson(lasso_m)
      if (!is.na(val_p)) {
        rows[[length(rows) + 1]] <- data.table(
          Protein = proteins[i], Pearson = val_p, Model = model_type, DiPa_Group = dipa_label
        )
      }
    }
  }
  return(rbindlist(rows, fill = TRUE))
}

run_proc5_4models <- function(dataset_name, plot_title, out_filename) {
  dt_list <- list()
  
  for (i in 1:4) {
    if (dataset_name == "BDL") {
      f_prime <- file.path(intra_dir, sprintf("Models_rf_preselection_full_objects_CD_new_LCPM_prime_cluster_%s_sub.RData", group_codes[i]))
      f_m50   <- file.path(m50_dir,   sprintf("Models_mastery50_LCPM_cluster_%s_sub.RData", group_codes[i]))
      f_rfrf  <- file.path(intra_rfrf_dir, sprintf("Models_rf_rf_preselection_full_objects_CD_new_LCPM_prime_cluster_%s_sub.RData", group_codes[i]))
    } else {
      f_prime <- file.path(intra_dir, sprintf("Models_rf_preselection_full_objects_CD_new_DTccl4_prime_cluster_%s_sub.RData", group_codes[i]))
      f_m50   <- file.path(m50_dir,   sprintf("Models_mastery50_DTccl4_cluster_%s_sub.RData", group_codes[i]))
      f_rfrf  <- file.path(intra_rfrf_dir, sprintf("Models_rf_rf_preselection_full_objects_CD_new_DTccl4_prime_cluster_%s_sub.RData", group_codes[i]))
    }
    
    dt_base <- extract_model_data(f_prime, "Baseline",   dipa_groups[i])
    dt_m50  <- extract_model_data(f_m50,   "Mastery 50",  dipa_groups[i])
    dt_prot <- extract_model_data(f_prime, "RF + LASSO", dipa_groups[i])
    dt_rfrf <- extract_model_data(f_rfrf,  "RF + RF",     dipa_groups[i])
    
    if (nrow(dt_base) > 0) dt_list[[length(dt_list) + 1]] <- dt_base
    if (nrow(dt_m50) > 0)  dt_list[[length(dt_list) + 1]] <- dt_m50
    if (nrow(dt_prot) > 0) dt_list[[length(dt_list) + 1]] <- dt_prot
    if (nrow(dt_rfrf) > 0) dt_list[[length(dt_list) + 1]] <- dt_rfrf
  }
  
  dt_all <- rbindlist(dt_list, fill = TRUE)
  dt_all[, DiPa_Group := factor(DiPa_Group, levels = dipa_groups)]
  dt_all[, Model := factor(Model, levels = model_names)]
  
  summary_dt <- dt_all[!is.na(Pearson), .(
    pct_high = round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100),
    pct_suff = round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100),
    label    = paste0(">= 0.8: ", round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100), "%\n>= 0.5: ", round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100), "%")
  ), by = .(DiPa_Group, Model)]
  
  palette_fill <- c(
    "Baseline"   = "#E5E7E9",
    "Mastery 50" = "#85929E",
    "RF + LASSO" = "#2C3E50",
    "RF + RF"    = "#1B2631"
  )
  
  p <- ggplot(dt_all, aes(x = DiPa_Group, y = Pearson, fill = Model)) +
    stat_boxplot(geom = "errorbar", position = position_dodge(0.8), width = 0.25, color = "#1A252F", linewidth = 0.6) +
    geom_boxplot(position = position_dodge(0.8), color = "#1A252F", alpha = 0.85, outlier.shape = NA, width = 0.62) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               color = "#154360", alpha = 0.35, size = 1.3) +
    
    geom_segment(data = data.frame(x = c(1.5, 2.5, 3.5), xend = c(1.5, 2.5, 3.5), y = -1.0, yend = 1.0),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "solid", color = "grey80", linewidth = 0.6, inherit.aes = FALSE) +
    
    geom_hline(yintercept = seq(-1.0, 1.0, by = 0.25), color = "grey90", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#E74C3C", linewidth = 0.85) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#27AE60", linewidth = 0.85) +
    
    geom_text(data = summary_dt, 
              aes(x = DiPa_Group, y = -1.22, label = label, group = Model),
              position = position_dodge(width = 0.8),
              family = "sans", size = 3.1, color = "black", fontface = "bold", lineheight = 0.88, inherit.aes = FALSE) +
    
    scale_fill_manual(values = palette_fill) +
    
    scale_y_continuous(
      name = expression(bold("Pearson correlation, "*rho[BP])),
      breaks = seq(-1.0, 1.0, by = 0.25),
      limits = c(-1.42, 1.05)
    ) +
    coord_cartesian(ylim = c(-1.42, 1.05), clip = "off") +
    
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(family = "sans", face = "bold", hjust = 0.5, size = 19, color = "black", margin = margin(b = 12)),
      plot.margin = margin(t = 15, r = 20, b = 32, l = 20),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(family = "sans", face = "bold", size = 15, color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(family = "sans", face = "bold", size = 15, color = "black"),
      axis.text.y  = element_text(family = "sans", face = "bold", size = 12, color = "black"),
      axis.ticks   = element_line(color = "black", linewidth = 0.7),
      legend.position = "top",
      legend.title    = element_text(family = "sans", face = "bold", size = 14, color = "black"),
      legend.text     = element_text(family = "sans", face = "bold", size = 13, color = "black"),
      legend.margin   = margin(b = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.1)
    ) +
    labs(
      title = plot_title,
      fill  = "Model:"
    )
  
  out_path <- file.path(output_dir, out_filename)
  cairo_pdf(out_path, width = 16.5, height = 8.5)
  print(p)
  dev.off()
  return(out_path)
}

# 1. Intra BDL (Panel A)
out_bdl <- run_proc5_4models(
  dataset_name = "BDL",
  plot_title   = expression(bold("Intra BDL")),
  out_filename = "Procedure_5_4Models_Pearson_Intra_BDL.pdf"
)
cat("Procedure 5 Intra BDL PDF successfully generated:", out_bdl, "\n")

# 2. Intra CCl4 (Panel B)
out_ccl4 <- run_proc5_4models(
  dataset_name = "CCL4",
  plot_title   = expression(bold("Intra "*CCl[4])),
  out_filename = "Procedure_5_4Models_Pearson_Intra_CCl4.pdf"
)
cat("Procedure 5 Intra CCl4 PDF successfully generated:", out_ccl4, "\n")
