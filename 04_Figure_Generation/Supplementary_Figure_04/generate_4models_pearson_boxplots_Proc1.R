# ==============================================================================
# SCRIPT: generate_4models_pearson_boxplots_Proc1.R
# PURPOSE: Generate 4-Models Pearson Correlation Boxplots for Supplementary Figure 4
#          Procedure 1 Direction 2 (Train CCl4, test BDL)
#          Models: Baseline, Mastery 50, RF + LASSO, RF + RF
# ==============================================================================

library(data.table)
library(ggplot2)
library(dplyr)

output_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

master_fig_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
if (!dir.exists(master_fig_dir)) dir.create(master_fig_dir, recursive = TRUE)

prime_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis"
m50_dir   <- file.path(prime_dir, "Mastery50")
rf_rf_dir <- file.path(prime_dir, "RF_RF")

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
    } else if (model_type == "Mastery 50" && (grepl("mastery", v_low) || grepl("protein", v_low)) && !grepl("plus_rna|combined", v_low)) {
      obj_target <- env[[v]]
    } else if (model_type == "RF + RF" && (grepl("rf_rf", v_low) || grepl("protein", v_low)) && !grepl("plus_rna|combined|mastery|lasso", v_low)) {
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

run_supp4_4models <- function() {
  dt_list <- list()
  
  for (i in 1:4) {
    f_prime <- file.path(prime_dir, sprintf("Models_rf_lasso_full_testing_Batch_prime_%s_sub.RData", group_codes[i]))
    f_m50   <- file.path(m50_dir,   sprintf("Models_mastery50_CCL4_BDL_%s_prime_sub.RData", group_codes[i]))
    f_rfrf  <- file.path(rf_rf_dir, sprintf("Models_rf_rf_full_testing_Batch_prime_%s_sub.RData", group_codes[i]))
    
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
    label    = paste0("\u2265 0.8: ", round(sum(Pearson >= 0.8, na.rm = TRUE) / .N * 100), "%\n\u2265 0.5: ", round(sum(Pearson >= 0.5, na.rm = TRUE) / .N * 100), "%")
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
              family = "sans", size = 3.35, color = "black", fontface = "bold", lineheight = 0.90, inherit.aes = FALSE) +
    
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
      plot.margin = margin(t = 15, r = 20, b = 25, l = 20),
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
      title = expression(bold("Train "*CCl[4]*", test BDL")),
      fill  = "Model:"
    )
  
  # 1. Output Component PDF
  out_path_comp <- file.path(output_dir, "Procedure_1_4Models_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf")
  cairo_pdf(out_path_comp, width = 16.5, height = 8.5)
  print(p)
  dev.off()
  
  # 2. Output Master Supplementary Figure 4 PDF
  out_path_master <- file.path(master_fig_dir, "Supplementary_Figure_4.pdf")
  cairo_pdf(out_path_master, width = 16.5, height = 8.5)
  print(p)
  dev.off()
  
  cat("Supplementary Figure 4 PDF successfully generated:\n", out_path_comp, "\n", out_path_master, "\n")
}

run_supp4_4models()
