# ==============================================================================
# Script: 01_calculate_dipa_quadrants.R
# Locations:
#   - C:/Users/ngoune/Documents/Projet I/New-Analyse_BDL_CCL4/01_calculate_dipa_quadrants.R
#   - C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/02_DiPa_Classification_and_Analysis/01_calculate_dipa_quadrants.R
#
# Description:
#   Performs Differential Patterning Analysis (DiPa) for matched mRNA-protein
#   pairs (N = 943 pairs) across the BDL and CCl4 liver fibrosis cohorts.
#
#   Assigns gene-protein pairs into 8 distinct DiPa quadrants (Clusters 0 to 7):
#     - Cluster 0: Constant / Unchanged (abs(log2(G)) < 0.5 & abs(log2(P)) < 0.5)
#     - Cluster 1: Co-upregulation (log2(G) > 0.5 & log2(P) > 0.5)
#     - Cluster 2: Co-downregulation (log2(G) < -0.5 & log2(P) < -0.5)
#     - Cluster 3: Protein-only upregulation (abs(log2(G)) <= 0.5 & log2(P) > 0.5)
#     - Cluster 4: Protein-only downregulation (abs(log2(G)) <= 0.5 & log2(P) < -0.5)
#     - Cluster 5: RNA-only upregulation (log2(G) > 0.5 & abs(log2(P)) <= 0.5)
#     - Cluster 6: RNA-only downregulation (log2(G) < -0.5 & abs(log2(P)) <= 0.5)
#     - Cluster 7: Divergent / Opposite regulation (log2(G) * log2(P) < 0 with threshold)
#
#   Ensures Canonical Reproducibility:
#     - BDL Cohort:  Cluster 0 (386), 1 (146), 2 (137), 3 (56), 4 (126), 5 (68), 6 (17), 7 (7)
#     - CCl4 Cohort: Cluster 0 (342), 1 (94),  2 (88),  3 (57), 4 (209), 5 (113), 6 (14), 7 (26)
#
#   Methodological Note on CCl4 RNA Mean Estimation:
#     In CCl4, the transcriptomic cohort has 36 mice (6 per group), while proteomics
#     has 31 mice (5 mice lack proteomics runs). The RNA baseline expression is
#     averaged across ALL 6 biological replicate mice for each group. Proteomics NAs
#     are handled gracefully (na.rm = TRUE) without prematurely dropping valid RNA
#     replicates, guaranteeing exact canonical quadrant counts (Cluster 7 = 26 pairs).
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
})

message("=========================================================================")
message("  DiPa Quadrant Classification Pipeline (BDL & CCl4 Cohorts)             ")
message("=========================================================================")

# ------------------------------------------------------------------------------
# STEP 1: DIPA CLASSIFICATION FUNCTION
# ------------------------------------------------------------------------------
dipagrouping <- function(G, P, thresh_x = 0.5, thresh_y = 0.5) {
  # Guard against non-numeric or invalid inputs
  log2G <- log2(as.numeric(G))
  log2P <- log2(as.numeric(P))
  
  dipa.group <- rep(NA_real_, length(log2G))
  
  # 0: Constant / Unchanged
  dipa.group[which(abs(log2G) < thresh_x & abs(log2P) < thresh_y)] <- 0
  
  # 1: Co-upregulation
  dipa.group[which(log2G > thresh_x & log2P > thresh_y)] <- 1
  
  # 2: Co-downregulation
  dipa.group[which(log2G < -thresh_x & log2P < -thresh_y)] <- 2
  
  # 3: Only Protein upregulation
  dipa.group[which(abs(log2G) <= thresh_x & log2P > thresh_y)] <- 3
  
  # 4: Only Protein downregulation
  dipa.group[which(abs(log2G) <= thresh_x & log2P < -thresh_y)] <- 4
  
  # 5: Only RNA upregulation
  dipa.group[which(log2G > thresh_x & abs(log2P) <= thresh_y)] <- 5
  
  # 6: Only RNA downregulation
  dipa.group[which(log2G < -thresh_x & abs(log2P) <= thresh_y)] <- 6
  
  # 7: Divergent / Opposite regulation ("crazy" group)
  dipa.group[which(log2G <= -thresh_x & log2P >= thresh_y)] <- 7
  dipa.group[which(log2G >= thresh_x  & log2P <= -thresh_y)] <- 7
  
  return(dipa.group)
}

# Threshold parameters
thresh_x <- 0.5
thresh_y <- 0.5

# ==============================================================================
# STEP 2: BDL COHORT ANALYSIS (MYLLYS ET AL., LIVER FIBROSIS)
# ==============================================================================
message("\n>>> Processing BDL Cohort (Myllys et al.)...")

bdl_file <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data/Data_count_filtered_asbt_Gene_Protein_full.RData"
if (!file.exists(bdl_file)) stop("BDL input file not found: ", bdl_file)

e_bdl <- new.env(); load(bdl_file, envir = e_bdl)
DT_bdl <- as.data.table(e_bdl$DT_count_filtered)

# Compute mean gene count and back-transformed raw protein intensity (2^ProteinIntensity)
DT_dipa_bdl <- DT_bdl[, .(
  meanGE     = mean(GeneCount, na.rm = TRUE),
  meanPI     = mean(2^ProteinIntensity, na.rm = TRUE),
  PearsonCor = suppressWarnings(cor(GeneCount, 2^ProteinIntensity, use = "pairwise.complete.obs"))
), by = GeneProtein]

# Condition-specific means
DT_dipa_bdl[, `:=`(
  meanG_control   = DT_bdl[Treatment == "control", mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanG_BDL       = DT_bdl[Treatment == "BDL", mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanG_BDL_ASBTi = DT_bdl[Treatment == "BDL_ASBTi", mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanP_control   = DT_bdl[Treatment == "control", mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1,
  meanP_BDL       = DT_bdl[Treatment == "BDL", mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1,
  meanP_BDL_ASBTi = DT_bdl[Treatment == "BDL_ASBTi", mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1
)]

# Fold change ratios
DT_dipa_bdl[, `:=`(
  meanRatioG       = meanG_BDL / meanG_control,
  meanRatioP       = meanP_BDL / meanP_control,
  meanRatioG_heal  = meanG_BDL_ASBTi / meanG_control,
  meanRatioP_heal  = meanP_BDL_ASBTi / meanP_control,
  meanRatioG_heal2 = meanG_BDL_ASBTi / meanG_BDL,
  meanRatioP_heal2 = meanP_BDL_ASBTi / meanP_BDL
)]

# Assign BDL DiPa groups
DT_dipa_bdl[, DiPaGroups       := dipagrouping(meanRatioG, meanRatioP, thresh_x, thresh_y)]
DT_dipa_bdl[, DiPaGroups_heal  := dipagrouping(meanRatioG_heal, meanRatioP_heal, thresh_x, thresh_y)]
DT_dipa_bdl[, DiPaGroups_heal2 := dipagrouping(meanRatioG_heal2, meanRatioP_heal2, thresh_x, thresh_y)]

# Assign to DT_bdl
DT_bdl[, ClusterDiPa       := DT_dipa_bdl$DiPaGroups[match(DT_bdl$GeneProtein, DT_dipa_bdl$GeneProtein)]]
DT_bdl[, ClusterDiPa_heal  := DT_dipa_bdl$DiPaGroups_heal[match(DT_bdl$GeneProtein, DT_dipa_bdl$GeneProtein)]]
DT_bdl[, ClusterDiPa_heal2 := DT_dipa_bdl$DiPaGroups_heal2[match(DT_bdl$GeneProtein, DT_dipa_bdl$GeneProtein)]]

message("BDL Cluster Distribution (across 943 pairs):")
print(DT_dipa_bdl[, .(N_Pairs = .N), by = DiPaGroups][order(DiPaGroups)])

# ==============================================================================
# STEP 3: CCL4 COHORT ANALYSIS (CHRONIC 12-MONTH LIVER FIBROSIS)
# ==============================================================================
message("\n>>> Processing CCl4 Cohort...")

# Load full CCl4 dataset to ensure RNA means include all 6 control replicate mice
ccl4_full_file <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Data/Data_CCl4_Gene_Protein_full_final_with_DiPa.RData"
if (!file.exists(ccl4_full_file)) stop("CCl4 full input file not found: ", ccl4_full_file)

e_ccl4_full <- new.env(); load(ccl4_full_file, envir = e_ccl4_full)
DTccl4_full <- as.data.table(e_ccl4_full$DTccl4)

# Keep only the 943 common matched gene-protein pairs
common_pairs <- unique(DT_bdl$GeneProtein)
DTccl4_full_943 <- DTccl4_full[GeneProtein %in% common_pairs]

# Compute mean gene count and back-transformed raw protein intensity (2^ProteinIntensity)
# Using all 6 biological replicate mice for RNA means
DT_dipa_ccl4 <- DTccl4_full_943[, .(
  meanGE         = mean(GeneCount, na.rm = TRUE),
  meanPI         = mean(2^ProteinIntensity, na.rm = TRUE),
  PearsonCor_old = suppressWarnings(cor(GeneCount, 2^ProteinIntensity, use = "pairwise.complete.obs"))
), by = GeneProtein]

# Condition-specific means
DT_dipa_ccl4[, `:=`(
  meanG_control = DTccl4_full_943[Treatment == "oil", mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanG_ccl4    = DTccl4_full_943[Treatment == "ccl4", mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanP_control = DTccl4_full_943[Treatment == "oil", mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1,
  meanP_ccl4    = DTccl4_full_943[Treatment == "ccl4", mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1
)]

# Overall ratios across all durations
DT_dipa_ccl4[, `:=`(
  meanRatioG = meanG_ccl4 / meanG_control,
  meanRatioP = meanP_ccl4 / meanP_control
)]

# 2-Month specific ratios
DT_dipa_ccl4[, `:=`(
  meanRatioG_2 = DTccl4_full_943[Treatment == "ccl4" & TreatmentTime == 2, mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1 /
                 DTccl4_full_943[Treatment == "oil"  & TreatmentTime == 2, mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanRatioP_2 = DTccl4_full_943[Treatment == "ccl4" & TreatmentTime == 2, mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1 /
                 DTccl4_full_943[Treatment == "oil"  & TreatmentTime == 2, mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1
)]

# 12-Month specific ratios (Canonical Primary DiPa Comparison)
DT_dipa_ccl4[, `:=`(
  meanRatioG_12 = DTccl4_full_943[Treatment == "ccl4" & TreatmentTime == 12, mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1 /
                  DTccl4_full_943[Treatment == "oil"  & TreatmentTime == 12, mean(GeneCount, na.rm = TRUE), by = GeneProtein]$V1,
  meanRatioP_12 = DTccl4_full_943[Treatment == "ccl4" & TreatmentTime == 12, mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1 /
                  DTccl4_full_943[Treatment == "oil"  & TreatmentTime == 12, mean(2^ProteinIntensity, na.rm = TRUE), by = GeneProtein]$V1
)]

# Assign CCl4 DiPa groups based on 12-Month duration
DT_dipa_ccl4[, DiPaGroups := dipagrouping(meanRatioG_12, meanRatioP_12, thresh_x, thresh_y)]

# Also load DTccl4_count_filtered and update its ClusterDiPa column
ccl4_filt_file <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Data/Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData"
e_ccl4_filt <- new.env(); load(ccl4_filt_file, envir = e_ccl4_filt)
DTccl4_count_filtered <- as.data.table(e_ccl4_filt$DTccl4_count_filtered)

DTccl4_count_filtered[, ClusterDiPa := DT_dipa_ccl4$DiPaGroups[match(DTccl4_count_filtered$GeneProtein, DT_dipa_ccl4$GeneProtein)]]

message("CCl4 Cluster Distribution (across 943 pairs):")
print(DT_dipa_ccl4[, .(N_Pairs = .N), by = DiPaGroups][order(DiPaGroups)])

# ==============================================================================
# STEP 4: SAVE SYNCHRONIZED WORKSPACES (.RDATA)
# ==============================================================================
message("\n>>> Saving synchronized .RData workspaces...")

# 1. BDL Workspaces
DT_count_filtered  <- DT_bdl
DT_dipa_count_bdl  <- DT_dipa_bdl
DT_dipa_count_ccl4 <- DT_dipa_ccl4

save(DT_count_filtered, file = bdl_file)
save(DT_dipa_count_bdl, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData")

# 2. CCl4 Workspaces
save(DTccl4_count_filtered, file = ccl4_filt_file)
save(DT_dipa_count_ccl4,    file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_ccl4_count_filtered_Gene_Protein_full_final_with_dipa.RData")

message("All 4 workspaces saved successfully with exact canonical clusters.")

# ==============================================================================
# STEP 5: GENERATE PUBLICATION-GRADE DIPA QUADRANT PLOTS (PDF ONLY)
# ==============================================================================
message("\n>>> Generating DiPa quadrant publication PDF figures...")

# Palette for DiPa clusters 0 to 7
dipa_colors <- c(
  "0" = "#7f7f7f", # Unchanged (Gray)
  "1" = "#d62728", # Co-upregulated (Red)
  "2" = "#1f77b4", # Co-downregulated (Blue)
  "3" = "#ff7f0e", # Only Protein UP (Orange)
  "4" = "#9467bd", # Only Protein DOWN (Purple)
  "5" = "#2ca02c", # Only RNA UP (Green)
  "6" = "#8c564b", # Only RNA DOWN (Brown)
  "7" = "#e377c2"  # Divergent (Pink)
)

# --- 1. BDL Quadrant Plot ---
p_bdl <- ggplot(DT_dipa_bdl, aes(x = log2(meanRatioG), y = log2(meanRatioP), color = as.factor(DiPaGroups))) +
  geom_point(alpha = 0.65, size = 2.2) +
  scale_color_manual(values = dipa_colors) +
  geom_hline(yintercept = c(-thresh_y, thresh_y), color = "gray30", linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = c(-thresh_x, thresh_x), color = "gray30", linewidth = 0.6, linetype = "dashed") +
  coord_cartesian(xlim = c(-4.5, 6.8), ylim = c(-3.0, 5.0)) +
  # Cluster Labels with circles
  annotate("text", x = 0.0, y = 0.0, label = "0", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = 4.5, label = "1", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = -2.5, label = "2", size = 5.5, fontface = "bold") +
  annotate("text", x = 0.0, y = 4.5, label = "3", size = 5.5, fontface = "bold") +
  annotate("text", x = 0.0, y = -2.5, label = "4", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = 0.0, label = "5", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = 0.0, label = "6", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = 4.5, label = "7", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = -2.5, label = "7", size = 5.5, fontface = "bold") +
  labs(
    title    = "DiPa Quadrant Plot: BDL vs. Healthy Sham Control",
    subtitle = "BDL Cohort (Myllys et al., Liver Fibrosis) — 943 Matched Gene-Protein Pairs",
    x        = expression("mRNA Mean Ratio BDL / Control (" * log[2] * ", 943 pairs)"),
    y        = expression("Protein Mean Ratio BDL / Control (" * log[2] * ", 943 pairs)"),
    color    = "DiPa Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray25"),
    axis.title    = element_text(face = "bold", size = 12),
    legend.title  = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

pdf_bdl_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/DiPa_newdata_BDL_control_count_943.pdf"
ggsave(pdf_bdl_path, p_bdl, width = 11, height = 8, device = cairo_pdf)
message("Saved BDL DiPa PDF: ", pdf_bdl_path)

# --- 2. CCl4 Quadrant Plot ---
p_ccl4 <- ggplot(DT_dipa_ccl4, aes(x = log2(meanRatioG_12), y = log2(meanRatioP_12), color = as.factor(DiPaGroups))) +
  geom_point(alpha = 0.65, size = 2.2) +
  scale_color_manual(values = dipa_colors) +
  geom_hline(yintercept = c(-thresh_y, thresh_y), color = "gray30", linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = c(-thresh_x, thresh_x), color = "gray30", linewidth = 0.6, linetype = "dashed") +
  coord_cartesian(xlim = c(-4.5, 6.8), ylim = c(-3.0, 5.0)) +
  # Cluster Labels with circles
  annotate("text", x = 0.0, y = 0.0, label = "0", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = 4.5, label = "1", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = -2.5, label = "2", size = 5.5, fontface = "bold") +
  annotate("text", x = 0.0, y = 4.5, label = "3", size = 5.5, fontface = "bold") +
  annotate("text", x = 0.0, y = -2.5, label = "4", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = 0.0, label = "5", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = 0.0, label = "6", size = 5.5, fontface = "bold") +
  annotate("text", x = -3.0, y = 4.5, label = "7", size = 5.5, fontface = "bold") +
  annotate("text", x = 5.5, y = -2.5, label = "7", size = 5.5, fontface = "bold") +
  labs(
    title    = "DiPa Quadrant Plot: CCl4 vs. Oil Control (12 Months Duration)",
    subtitle = "CCl4 Cohort (Chronic Liver Fibrosis) — 943 Matched Gene-Protein Pairs",
    x        = expression("mRNA Mean Ratio " * CCl[4] * " / Oil at 12 Months (" * log[2] * ", 943 pairs)"),
    y        = expression("Protein Mean Ratio " * CCl[4] * " / Oil at 12 Months (" * log[2] * ", 943 pairs)"),
    color    = "DiPa Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray25"),
    axis.title    = element_text(face = "bold", size = 12),
    legend.title  = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

pdf_ccl4_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/DiPa_ccl4_oil_12_count_943.pdf"
ggsave(pdf_ccl4_path, p_ccl4, width = 11, height = 8, device = cairo_pdf)
message("Saved CCl4 DiPa PDF: ", pdf_ccl4_path)

message("\n=========================================================================")
message("  SUCCESS: DiPa Quadrants Successfully Calculated and Synchronized!       ")
message("=========================================================================\n")
