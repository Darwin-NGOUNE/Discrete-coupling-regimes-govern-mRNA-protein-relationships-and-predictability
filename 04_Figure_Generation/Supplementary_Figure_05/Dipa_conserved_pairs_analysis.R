# ==============================================================================
# SCRIPT: jan_dipa_conserved_pairs_analysis.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/jan_dipa_conserved_pairs_analysis.R
# PURPOSE: Focused Publication-Quality Stacked Barplots (PEARSON ONLY):
#          "Proportion of Conserved Best-Partner Pairs Across DiPa Groups (rho_BP)"
#          EXACT STRUCTURE & COLORS:
#          - Green (#27AE60): High (rho_BP >= 0.8)
#          - Blue  (#2980B9): Moderate (0.5 <= rho_BP < 0.8)
#          - Orange (#F39C12): Weak (0.2 <= rho_BP < 0.5)
#          - Red   (#E74C3C): Lost / Noise (rho_BP < 0.2)
#          - Segment Labels: Dynamically formatted (single line for small slices < 6%) to prevent ANY border clipping
#          - Top Annotations: "Conserved (>= 0.5):\nXX.X% (YY/Total)"
#          - Generates BOTH:
#            1. SUBSET COHORT (12 BDL mice vs 12 CCl4 mice)
#            2. FULL COHORT   (18 BDL mice vs 36 CCl4 mice)
# ==============================================================================

library(data.table)
library(ggplot2)
library(scales)

# ------------------------------------------------------------------------------
# STEP 1: LOAD DATASET & MAP CORE PROTEINS
# ------------------------------------------------------------------------------
cat("1. Loading ComBat batch-corrected dataset...\n")

load_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData"
load(load_path)
Full_DT <- as.data.table(Full_DT)

bdl_all  <- Full_DT[is.na(TreatmentTime)]
ccl4_all <- Full_DT[!is.na(TreatmentTime)]
core_proteins <- intersect(unique(bdl_all$GeneProtein), unique(ccl4_all$GeneProtein))
n_prot <- length(core_proteins)
cat(sprintf("Total shared core proteins: %d\n", n_prot))

# ------------------------------------------------------------------------------
# STEP 2: MAP EXACT DIPA REPLICATION OVERLAP SUBSETS (118, 67, 23, 202, 533)
# ------------------------------------------------------------------------------
prime_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/"
load_gp <- function(fname) {
  e <- new.env(); load(paste0(prime_dir, fname), envir = e)
  vars <- ls(e); v <- vars[grepl("protein", tolower(vars)) & !grepl("combined|plus_rna", tolower(vars))][1]
  sapply(e[[v]], function(x) { m <- if(!is.null(x$prediction.obj)) x$prediction.obj$model else x$model; paste0(m$gene, "_", m$protein) })
}

gp_12 <- load_gp("Models_rf_lasso_full_testing_new_Batch_prime_1_2_sub.RData")
gp_34 <- load_gp("Models_rf_lasso_full_testing_new_Batch_prime_3_4_sub.RData")
gp_56 <- load_gp("Models_rf_lasso_full_testing_new_Batch_prime_5_6_sub.RData")
gp_0  <- load_gp("Models_rf_lasso_full_testing_new_Batch_prime_0_sub.RData")

dipa_map <- data.table(Protein = core_proteins, DiPa_Group = "Other (N = 533)")
dipa_map[Protein %in% gp_12, DiPa_Group := "DiPa 1 & 2 (N = 118)"]
dipa_map[Protein %in% gp_34, DiPa_Group := "DiPa 3 & 4 (N = 67)"]
dipa_map[Protein %in% gp_56, DiPa_Group := "DiPa 5 & 6 (N = 23)"]
dipa_map[Protein %in% gp_0,  DiPa_Group := "DiPa 8 (N = 202)"]

dipa_levels <- c("DiPa 1 & 2 (N = 118)", "DiPa 3 & 4 (N = 67)", "DiPa 5 & 6 (N = 23)", "DiPa 8 (N = 202)", "Other (N = 533)")
dipa_map$DiPa_Group <- factor(dipa_map$DiPa_Group, levels = dipa_levels)

# Helper for tiers
assign_tier <- function(val) {
  fcase(
    val >= 0.8, "High",
    val >= 0.5, "Moderate",
    val >= 0.2, "Low",
    default = "Lost"
  )
}
tier_levels <- c("Lost", "Low", "Moderate", "High")
tier_colors <- c(
  "High"     = "#27AE60", # Green
  "Moderate" = "#2980B9", # Blue
  "Low"      = "#F39C12", # Orange
  "Lost"     = "#E74C3C"  # Red
)

# ------------------------------------------------------------------------------
# STEP 3: PLOTTING FUNCTION WITH ADAPTIVE SINGLE/MULTI-LINE LABELS
# ------------------------------------------------------------------------------
out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/Scatterplots_75th_Percentile"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

build_stacked_plot <- function(dt_master, dir_col, train_name, test_name, cohort_label) {
  
  dt_plot <- copy(dt_master)
  dt_plot[, Tier := factor(assign_tier(get(dir_col)), levels = tier_levels)]
  
  agg_dt <- dt_plot[, .(Count = .N), by = .(DiPa_Group, Tier)]
  tot_dt <- agg_dt[, .(Total = sum(Count)), by = DiPa_Group]
  agg_dt <- merge(agg_dt, tot_dt, by = "DiPa_Group")
  agg_dt[, Prop := Count / Total]
  agg_dt[, Pct := Prop * 100]
  
  # Calculate exact cumulative center inside the stack
  setorder(agg_dt, DiPa_Group, -Tier)
  agg_dt[, CumProp := cumsum(Prop) - (0.5 * Prop), by = DiPa_Group]
  
  # Smart adaptive labels: 2 lines if space >= 6%, 1 line if space between 2.5% and 6%
  agg_dt[, Label := fcase(
    Pct >= 6.0, sprintf("%.1f%%\n(n = %d)", Pct, Count),
    Pct >= 2.5, sprintf("%.1f%% (n = %d)", Pct, Count),
    default = ""
  )]
  agg_dt[, FontSize := fifelse(Pct >= 6.0, 3.8, 3.1)]
  
  # Top annotations for Conserved (>= 0.5)
  top_annot <- dt_plot[, .(
    Cons_Count = sum(get(dir_col) >= 0.5),
    Total = .N
  ), by = DiPa_Group]
  top_annot[, Cons_Pct := Cons_Count / Total * 100]
  top_annot[, TopLabel := sprintf("Conserved (>= 0.5):\n%.1f%% (%d/%d)", Cons_Pct, Cons_Count, Total)]
  
  legend_labels <- c(
    "Lost"     = expression(bold(paste("Lost (", rho[BP] < 0.2, ")"))),
    "Low"      = expression(bold(paste("Low (0.2" <= rho[BP], " < 0.5)"))),
    "Moderate" = expression(bold(paste("Moderate (0.5" <= rho[BP], " < 0.8)"))),
    "High"     = expression(bold(paste("High (", rho[BP] >= 0.8, ")")))
  )
  
  train_expr <- if (train_name %in% c("CCL4", "CCl4")) quote(bold(CCl[4])) else bquote(bold(.(train_name)))
  test_expr  <- if (test_name %in% c("CCL4", "CCl4")) quote(bold(CCl[4])) else bquote(bold(.(test_name)))
  
  p <- ggplot(agg_dt, aes(x = DiPa_Group, y = Prop, fill = Tier)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.65, linewidth = 0.55) +
    
    # Text inside segments with size aesthetic mapped
    geom_text(aes(y = CumProp, label = Label, size = FontSize), 
              family = "sans", fontface = "bold", color = "white", lineheight = 0.85) +
    scale_size_identity() +
    
    # Top annotation above 100%
    geom_text(data = top_annot, aes(x = DiPa_Group, y = 1.025, label = TopLabel, fill = NULL, size = NULL),
              family = "sans", size = 3.6, fontface = "bold", color = "black", lineheight = 0.85, vjust = 0) +
    
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      name = "Proportion of conservation pairs (%)",
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0),
      limits = c(0, 1.11)
    ) +
    coord_cartesian(ylim = c(0, 1.10), clip = "off") +
    scale_fill_manual(
      values = tier_colors,
      labels = legend_labels,
      name = bquote(bold("Best-partner ") * .(train_expr) * bold(" and conservation in ") * .(test_expr) * bold(":")),
      guide = guide_legend(reverse = TRUE, nrow = 1)
    ) +
    theme_bw(base_size = 13, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 16.5, hjust = 0.5, color = "black", margin = margin(b = 10)),
      plot.margin = margin(t = 15, r = 25, b = 15, l = 25),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(face = "bold", size = 12.5, color = "black", margin = margin(t = 6)),
      axis.title.y = element_text(face = "bold", size = 13.5, color = "black"),
      axis.text.y  = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "top",
      legend.title    = element_text(face = "bold", size = 12, color = "black"),
      legend.text     = element_text(face = "bold", size = 10.8, color = "black"),
      legend.spacing.x = unit(0.15, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
      panel.border     = element_rect(color = "#2C3E50", fill = NA, linewidth = 1.0)
    ) +
    labs(
      title = expression(bold("Conservation of protein best-partner pairs (") * bold(rho[BP]) * bold(")"))
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# STEP 4: RUN CALCULATION FOR COHORTS
# ------------------------------------------------------------------------------
run_cohort_analysis <- function(bdl_input, ccl4_input, cohort_label, filename_suffix) {
  cat(sprintf("\n=== RUNNING ANALYSIS FOR: %s ===\n", cohort_label))
  
  bdl_wide  <- dcast(bdl_input[GeneProtein %in% core_proteins], GeneProtein ~ MiceInfo, value.var = "ComBat_Protein_Raw", fun.aggregate = mean)
  ccl4_wide <- dcast(ccl4_input[GeneProtein %in% core_proteins], GeneProtein ~ MiceInfo, value.var = "ComBat_Protein_Raw", fun.aggregate = mean)
  
  mat_bdl  <- as.matrix(bdl_wide[, -1, with = FALSE]); rownames(mat_bdl) <- bdl_wide$GeneProtein
  mat_ccl4 <- as.matrix(ccl4_wide[, -1, with = FALSE]); rownames(mat_ccl4) <- ccl4_wide$GeneProtein
  
  cor_bdl_p  <- cor(t(mat_bdl), use = "pairwise.complete.obs")
  cor_ccl4_p <- cor(t(mat_ccl4), use = "pairwise.complete.obs")
  
  dt_master <- copy(dipa_map)
  dt_master[, `:=`(
    Best_Partner_D1 = character(n_prot),
    rho_BP_Test_D1  = numeric(n_prot),
    Best_Partner_D2 = character(n_prot),
    rho_BP_Test_D2  = numeric(n_prot)
  )]
  
  for (i in seq_len(n_prot)) {
    p <- core_proteins[i]
    # Direction 1: Train BDL -> Best Partner
    r_bdl_p <- cor_bdl_p[p, ]; r_bdl_p[p] <- -Inf
    bp_d1 <- names(which.max(r_bdl_p))
    dt_master[i, `:=`(Best_Partner_D1 = bp_d1, rho_BP_Test_D1 = cor_ccl4_p[p, bp_d1])]
    
    # Direction 2: Train CCl4 -> Best Partner
    r_ccl4_p <- cor_ccl4_p[p, ]; r_ccl4_p[p] <- -Inf
    bp_d2 <- names(which.max(r_ccl4_p))
    dt_master[i, `:=`(Best_Partner_D2 = bp_d2, rho_BP_Test_D2 = cor_bdl_p[p, bp_d2])]
  }
  
  # --- DIRECTION 1: Train BDL -> Test CCl4 ---
  p1 <- build_stacked_plot(dt_master, "rho_BP_Test_D1", "BDL", "CCL4", cohort_label)
  pdf_d1 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_%s.pdf", filename_suffix))
  cairo_pdf(pdf_d1, width = 12.5, height = 7.8)
  print(p1)
  dev.off()
  if (filename_suffix == "Subset") {
    cairo_pdf(file.path(out_dir, "Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson.pdf"), width = 12.5, height = 7.8)
    print(p1)
    dev.off()
  }
  png_d1 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_%s.png", filename_suffix))
  png(png_d1, width = 2500, height = 1560, res = 200)
  print(p1)
  dev.off()
  
  # --- DIRECTION 2: Train CCl4 -> Test BDL ---
  p2 <- build_stacked_plot(dt_master, "rho_BP_Test_D2", "CCL4", "BDL", cohort_label)
  pdf_d2 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_%s.pdf", filename_suffix))
  cairo_pdf(pdf_d2, width = 12.5, height = 7.8)
  print(p2)
  dev.off()
  if (filename_suffix == "Subset") {
    cairo_pdf(file.path(out_dir, "Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson.pdf"), width = 12.5, height = 7.8)
    print(p2)
    dev.off()
  }
  png_d1 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_%s.png", filename_suffix))
  png(png_d1, width = 2300, height = 1560, res = 200)
  print(p1)
  dev.off()
  
  # --- DIRECTION 2: Train CCl4 -> Test BDL ---
  p2 <- build_stacked_plot(dt_master, "rho_BP_Test_D2", "CCL4", "BDL", cohort_label)
  pdf_d2 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_%s.pdf", filename_suffix))
  cairo_pdf(pdf_d2, width = 11.5, height = 7.8)
  print(p2)
  dev.off()
  if (filename_suffix == "Subset") {
    cairo_pdf(file.path(out_dir, "Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson.pdf"), width = 11.5, height = 7.8)
    print(p2)
    dev.off()
  }
  png_d2 <- file.path(out_dir, sprintf("Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_%s.png", filename_suffix))
  png(png_d2, width = 2300, height = 1560, res = 200)
  print(p2)
  dev.off()
  
  cat(sprintf("Successfully generated %s Direction 1 and Direction 2 plots!\n", cohort_label))
}

# ------------------------------------------------------------------------------
# STEP 5: RUN FOR BOTH SUBSET AND FULL COHORTS
# ------------------------------------------------------------------------------

# 1. SUBSET COHORT (12 mice vs 12 mice)
bdl_sub  <- bdl_all[Treatment %in% c("BDL", "control")]
ccl4_sub <- ccl4_all[(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)]
run_cohort_analysis(bdl_sub, ccl4_sub, "| Subset Cohort", "Subset")

# 2. FULL COHORT (18 mice vs 36 mice)
run_cohort_analysis(bdl_all, ccl4_all, "| Full Cohort", "Full")

cat("\nALL CONSERVED PAIRS BARPLOTS (SUBSET & FULL) COMPLETED WITH PERFECT LABEL FORMATTING!\n")
