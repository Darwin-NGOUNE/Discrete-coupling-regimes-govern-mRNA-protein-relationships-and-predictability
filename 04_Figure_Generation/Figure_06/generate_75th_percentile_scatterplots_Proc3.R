# ==============================================================================
# SCRIPT: generate_75th_percentile_scatterplots_Proc3.R
# PURPOSE: Generate 75th Percentile Scatterplot PDFs for Procedure 3 (Overlap-Subset Dataset)
#          strictly using 2 MOUSE STATUS GROUPS (Control vs Disease) - Figure 6 Panel B
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

dipa_groups <- c("1_2", "3_4", "5_6", "0")
dipa_labels <- c("DiPa 1 & 2", "DiPa 3 & 4", "DiPa 5 & 6", "DiPa 8")

base_dir  <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot"
prime_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP"
out_dir   <- file.path(base_dir, "Scatterplots_75th_Percentile")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

extract_model_node <- function(item, model_type = "Baseline") {
  if (is.null(item)) return(NULL)
  obj <- if (!is.null(item$prediction.obj)) item$prediction.obj else item
  if (model_type == "Baseline") {
    if (!is.null(obj$baseline.model)) return(obj$baseline.model)
  } else {
    if (!is.null(obj$model)) return(obj$model)
  }
  return(NULL)
}

get_exp <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.null(x$experiments.test)) return(x$experiments.test)
  if (!is.null(x$prediction.obj$experiments.test)) return(x$prediction.obj$experiments.test)
  if (!is.null(x$experiments)) return(x$experiments)
  if (!is.null(x$prediction.obj$experiments)) return(x$prediction.obj$experiments)
  return(NULL)
}

get_prot_id <- function(x) {
  if (is.null(x)) return(NA_character_)
  obj <- if (!is.null(x$prediction.obj)) x$prediction.obj else x
  if (!is.null(obj$baseline.model$protein)) return(as.character(obj$baseline.model$protein))
  if (!is.null(obj$model$protein)) return(as.character(obj$model$protein))
  return(NA_character_)
}

safe_extract_pred_y <- function(raw_item, model_type = "Baseline") {
  if (is.null(raw_item)) return(NULL)
  
  node <- extract_model_node(raw_item, model_type = model_type)
  if (is.null(node) || is.character(node)) return(NULL)
  
  y_test_full  <- if (!is.null(node$y_test)) node$y_test else if (!is.null(node$y.t)) node$y.t else node$y
  pred_full    <- if (!is.null(node$prediction)) node$prediction else node$predictions.insample
  y_train_full <- if (!is.null(node$y_train)) node$y_train else node$y
  exp_full     <- get_exp(raw_item)
  
  if (is.null(y_test_full) || is.null(pred_full)) return(NULL)
  
  len <- min(length(y_test_full), length(pred_full))
  y_t_sub <- as.numeric(y_test_full[1:len])
  p_sub   <- as.numeric(pred_full[1:len])
  exp_vec <- if (!is.null(exp_full)) as.character(exp_full[1:len]) else rep("control", len)
  
  valid_idx <- which(!is.na(y_t_sub) & !is.na(p_sub))
  y_t_sub <- y_t_sub[valid_idx]
  p_sub   <- p_sub[valid_idx]
  exp_vec <- exp_vec[valid_idx]
  
  len <- length(y_t_sub)
  group_aligned <- character(len)
  
  for (j in seq_len(len)) {
    e_val <- exp_vec[j]
    s_lower <- tolower(e_val)
    is_ctrl <- grepl("oil|control|ctrl|month0", s_lower) & !grepl("ccl4|bdl", s_lower)
    group_aligned[j] <- if (is_ctrl) "Control" else "Disease"
  }
  
  res_dt <- data.table(
    y_true = y_t_sub,
    y_pred = p_sub,
    Group  = factor(group_aligned, levels = c("Control", "Disease"))
  )
  attr(res_dt, "y_train") <- as.numeric(y_train_full)
  return(res_dt)
}

calc_full_metrics_dt <- function(dt_xy) {
  if (is.null(dt_xy) || nrow(dt_xy) < 3) return(c(r = NA, rmse = NA, nrmse = NA, r2_test = NA))
  
  y_t  <- dt_xy$y_true
  y_p  <- dt_xy$y_pred
  y_tr <- attr(dt_xy, "y_train")
  
  r_val <- cor(y_t, y_p, use = "complete.obs")
  mse   <- mean((y_t - y_p)^2)
  rmse  <- sqrt(mse)
  
  mean_tr <- if (!is.null(y_tr)) mean(y_tr, na.rm = TRUE) else mean(y_t, na.rm = TRUE)
  sst     <- sum((y_t - mean_tr)^2)
  sse     <- sum((y_t - y_p)^2)
  r2_test <- 1 - (sse / sst)
  
  sd_tr <- if (!is.null(y_tr)) sd(y_tr, na.rm = TRUE) else sd(y_t, na.rm = TRUE)
  nrmse <- rmse / sd_tr
  
  return(c(r = r_val, rmse = rmse, nrmse = nrmse, r2_test = r2_test))
}

find_best_corner <- function(dt_xy, min_val, max_val) {
  df <- data.frame(x = dt_xy$y_pred, y = dt_xy$y_true)
  margin_offset <- 0.04 * (max_val - min_val)
  corners <- list(
    top_left     = list(x = min_val + margin_offset, y = max_val - margin_offset, hjust = 0, vjust = 1),
    top_right    = list(x = max_val - margin_offset, y = max_val - margin_offset, hjust = 1, vjust = 1),
    bottom_left  = list(x = min_val + margin_offset, y = min_val + margin_offset, hjust = 0, vjust = 0),
    bottom_right = list(x = max_val - margin_offset, y = min_val + margin_offset, hjust = 1, vjust = 0)
  )
  calc_min_dist <- function(corner) min(sqrt((df$x - corner$x)^2 + (df$y - corner$y)^2))
  dists <- sapply(corners, calc_min_dist)
  best_name <- names(which.max(dists))
  return(corners[[best_name]])
}

name_obj_by_protein <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(list())
  pnames <- sapply(obj, get_prot_id)
  names(obj) <- pnames
  return(obj[!is.na(names(obj))])
}

load_env <- function(path) {
  e <- new.env()
  if (file.exists(path)) load(path, envir = e)
  return(e)
}

get_list_from_env <- function(e, model_type = "Baseline") {
  if (is.null(e)) return(list())
  objs <- ls(e)
  if (length(objs) == 0) return(list())
  
  m <- if (model_type == "Baseline") {
    grep("rna", objs, value = TRUE)
  } else if (model_type == "Protein") {
    grep("protein", objs, value = TRUE)
  } else if (model_type == "Mastery 50") {
    grep("protein", objs, value = TRUE)
  }
  m <- m[!grepl("plus_rna|combined", m)]
  if (length(m) == 0) return(list())
  return(name_obj_by_protein(e[[m[1]]]))
}

compute_group_medians <- function(obj_list, model_type = "Baseline") {
  if (is.null(obj_list) || length(obj_list) == 0) return(c(r_med = NA, rmse_med = NA, nrmse_med = NA, r2_med = NA))
  metrics_mat <- t(sapply(obj_list, function(item) {
    dt_xy <- safe_extract_pred_y(item, model_type = model_type)
    return(calc_full_metrics_dt(dt_xy))
  }))
  return(c(
    r_med     = median(metrics_mat[, "r"], na.rm = TRUE),
    rmse_med  = median(metrics_mat[, "rmse"], na.rm = TRUE),
    nrmse_med = median(metrics_mat[, "nrmse"], na.rm = TRUE),
    r2_med    = median(metrics_mat[, "r2_test"], na.rm = TRUE)
  ))
}

create_panel <- function(panel_title, dt_xy, group_med_metrics, min_val, max_val, label_size = 3.3, title_size = 14) {
  pair_metrics <- calc_full_metrics_dt(dt_xy)
  
  r_val     <- pair_metrics["r"]
  rmse_val  <- pair_metrics["rmse"]
  nrmse_val <- pair_metrics["nrmse"]
  r2_val    <- pair_metrics["r2_test"]
  
  r_med     <- group_med_metrics["r_med"]
  rmse_med  <- group_med_metrics["rmse_med"]
  nrmse_med <- group_med_metrics["nrmse_med"]
  r2_med    <- group_med_metrics["r2_med"]
  
  ann_parse_str <- sprintf(
    "atop('Group Med:' ~ rho[BP] == '%.3f' ~ '| RMSE =' ~ '%.3f' ~ '| nRMSE =' ~ '%.3f' ~ '|' ~ R^2 == '%.3f', 'This Pair:' ~ rho[BP] == '%.3f' ~ '| RMSE =' ~ '%.3f' ~ '| nRMSE =' ~ '%.3f' ~ '|' ~ R^2 == '%.3f')",
    r_med, rmse_med, nrmse_med, r2_med,
    r_val, rmse_val, nrmse_val, r2_val
  )
  
  corner_pos <- find_best_corner(dt_xy, min_val, max_val)
  break_seq  <- seq(floor(min_val), ceiling(max_val), by = 1)
  
  p <- ggplot(dt_xy, aes(x = y_pred, y = y_true)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "#2C3E50", linetype = "solid", linewidth = 0.9) +
    
    geom_point(aes(color = Group), alpha = 0.85, size = 3.0) +
    scale_color_manual(
      values = c("Control" = "#1F77B4", "Disease" = "#D62728"),
      drop = TRUE
    ) +
    
    annotate("label", x = corner_pos$x, y = corner_pos$y, 
             label = ann_parse_str, parse = TRUE, 
             hjust = corner_pos$hjust, vjust = corner_pos$vjust, 
             family = "sans", size = label_size, label.padding = unit(0.2, "lines"),
             fill = "white", color = "black") +
    scale_x_continuous(breaks = break_seq, limits = c(min_val, max_val)) +
    scale_y_continuous(breaks = break_seq, limits = c(min_val, max_val)) +
    labs(title = panel_title, x = "Predicted Intensity", y = "Observed Intensity", color = "Mouse Status") +
    theme_bw(base_size = 13, base_family = "sans") +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = title_size, face = "bold", color = "black", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text  = element_text(size = 12, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.9),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.1),
      legend.position  = "none"
    )
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 5: LOAD FILES AND GENERATE PROC 3 SUBSET SCATTERPLOTS
# ------------------------------------------------------------------------------
prime_files <- sprintf(file.path(prime_dir, "Models_rf_preselection_full_objects_mergedata_blind_batch_prime_over_%s_24.RData"), dipa_groups)
m50_files   <- sprintf(file.path(prime_dir, "Mastery50/Models_mastery50_mergedata_blind_batch_over_%s_24.RData"), dipa_groups)

base_panels <- list()
panels_12   <- list()
panels_m50  <- list()
panels_prot <- list()

for (i in seq_along(dipa_groups)) {
  dg <- dipa_groups[i]
  dl <- dipa_labels[i]
  
  env_prime <- load_env(prime_files[i])
  env_m50   <- load_env(m50_files[i])
  
  rna_obj  <- get_list_from_env(env_prime, "Baseline")
  m50_obj  <- get_list_from_env(env_m50,   "Mastery 50")
  prot_obj <- get_list_from_env(env_prime, "Protein")
  
  if (length(rna_obj) == 0) next
  
  base_meds <- compute_group_medians(rna_obj,  "Baseline")
  m50_meds  <- compute_group_medians(m50_obj,  "Mastery 50")
  prot_meds <- compute_group_medians(prot_obj, "Protein")
  
  r_vals <- sapply(rna_obj, function(item) {
    dt_xy <- safe_extract_pred_y(item, "Baseline")
    if (is.null(dt_xy) || nrow(dt_xy) < 3) return(NA_real_)
    metrics <- calc_full_metrics_dt(dt_xy)
    return(as.numeric(metrics["r"]))
  })
  valid_r <- r_vals[!is.na(r_vals)]
  if (length(valid_r) == 0) next
  
  target_75 <- as.numeric(quantile(valid_r, 0.75, type = 7))
  prot_75   <- names(valid_r)[which.min(abs(valid_r - target_75))]
  
  dt_base <- safe_extract_pred_y(rna_obj[[prot_75]],  "Baseline")
  dt_m50  <- safe_extract_pred_y(m50_obj[[prot_75]],  "Mastery 50")
  dt_prot <- safe_extract_pred_y(prot_obj[[prot_75]], "Protein")
  
  all_p <- c(dt_base$y_pred, dt_m50$y_pred, dt_prot$y_pred)
  all_t <- c(dt_base$y_true, dt_m50$y_true, dt_prot$y_true)
  
  min_val <- min(c(all_p, all_t), na.rm = TRUE) - 0.3
  max_val <- max(c(all_p, all_t), na.rm = TRUE) + 0.3
  
  p1 <- create_panel(sprintf("%s | Baseline (%s)", dl, prot_75), dt_base, base_meds, min_val, max_val, label_size = 3.3, title_size = 14)
  p2 <- create_panel(sprintf("%s | Mastery 50 (%s)", dl, prot_75), dt_m50, m50_meds, min_val, max_val, label_size = 3.3, title_size = 14)
  p3 <- create_panel(sprintf("%s | Protein (%s)", dl, prot_75), dt_prot, prot_meds, min_val, max_val, label_size = 3.3, title_size = 14)
  
  base_panels[[dl]] <- p1
  panels_m50[[dl]]  <- p2
  panels_prot[[dl]] <- p3
  
  panels_12[[paste0(dl, "_1")]] <- p1
  panels_12[[paste0(dl, "_2")]] <- p2
  panels_12[[paste0(dl, "_3")]] <- p3
}

# --- 1. STRICT 3-PAGE LANDSCAPE 1x4 PDF ---
if (length(base_panels) == 4) {
  pdf_3p <- file.path(out_dir, "Proc3_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf")
  cairo_pdf(pdf_3p, width = 21, height = 6.0)
  
  make_proc3_header <- function(title_expr) {
    draw_point_custom <- function(x, y, color) {
      cowplot::draw_grob(grid::pointsGrob(x = unit(x, "npc"), y = unit(y, "npc"), pch = 16, gp = grid::gpar(col = color, fill = color, fontsize = 18)))
    }
    header_grob <- cowplot::as_grob(
      cowplot::ggdraw() +
        cowplot::draw_label(title_expr, x = 0.470, y = 0.5, hjust = 1, size = 22, fontface = "bold") +
        cowplot::draw_label("|", x = 0.482, y = 0.5, hjust = 0.5, size = 22, fontface = "bold", color = "black") +
        cowplot::draw_label("Mouse status:", x = 0.493, y = 0.5, hjust = 0, size = 20, fontface = "bold", color = "black") +
        draw_point_custom(x = 0.588, y = 0.5, color = "#1F77B4") +
        cowplot::draw_label("Control", x = 0.598, y = 0.5, hjust = 0, size = 18, fontface = "bold", color = "black") +
        draw_point_custom(x = 0.658, y = 0.5, color = "#D62728") +
        cowplot::draw_label("Disease", x = 0.668, y = 0.5, hjust = 0, size = 18, fontface = "bold", color = "black")
    )
    gridExtra::arrangeGrob(
      grid::nullGrob(),
      header_grob,
      nrow = 2,
      heights = grid::unit.c(grid::unit(0.4, "line"), grid::unit(2.4, "line"))
    )
  }
  
  # Page 1: Baseline
  top_expr_pg1 <- bquote(bold("Baseline model: BDL + CCl"[4]))
  g_base       <- gridExtra::arrangeGrob(grobs = base_panels, ncol = 4, nrow = 1, top = make_proc3_header(top_expr_pg1))
  grid::grid.draw(g_base)
  
  # Page 2: Mastery 50
  grid::grid.newpage()
  top_expr_pg2 <- bquote(bold("Mastery 50 model: BDL + CCl"[4]))
  g_m50        <- gridExtra::arrangeGrob(grobs = panels_m50, ncol = 4, nrow = 1, top = make_proc3_header(top_expr_pg2))
  grid::grid.draw(g_m50)
  
  # Page 3: Protein
  grid::grid.newpage()
  top_expr_pg3 <- bquote(bold("Protein model: BDL + CCl"[4]))
  g_prot       <- gridExtra::arrangeGrob(grobs = panels_prot, ncol = 4, nrow = 1, top = make_proc3_header(top_expr_pg3))
  grid::grid.draw(g_prot)
  
  dev.off()
  cat("  -> Saved NEW Proc 3 STRICT 3-PAGE 1x4 PDF (3 Pages Exactly):", pdf_3p, "\n")
}

cat("\nPROCEDURE 3 SUBSET SCATTERPLOTS COMPLETED!\n")
