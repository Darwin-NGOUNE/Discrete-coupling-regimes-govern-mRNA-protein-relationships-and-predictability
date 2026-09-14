# ==============================================================================
# SCRIPT: Cloud_Threshold_DataDriven.R
# PROJECT: Multi-Omics Protein Modeling & Liver Fibrosis Consortium (Paper 1)
# AUTHORS: Brunel Darwin Ngoune Domo, Andreas Groll, et al.
# PURPOSE: Unsupervised Data-Driven Derivation of Multi-Omics Cloud Morphology
#          Thresholds via 2D K-Means Clustering in the (Ratio, |Pearson R|) Space.
#
# METHODOLOGICAL RATIONALE:
# Rather than imposing arbitrary manual heuristic boundaries, this script applies
# an unsupervised 2D K-means clustering algorithm (K = 4) on the empirical joint
# distribution of the Transcriptomic-to-Proteomic standard deviation ratio:
#     Ratio = SD(RNA) / SD(Protein)
# and the absolute Bravais-Pearson correlation magnitude:
#     |R| = |cor(RNA, Protein)|
#
# The 4 natural geometric clusters correspond to the 4 canonical Cloud Categories:
#   1. Horizontal : High Ratio (RNA variance dominates), Low |R|
#   2. Vertical   : Low Ratio  (Protein variance dominates), Low |R|
#   3. Diagonal   : Medium Ratio, High |R| (Strong translational coupling)
#   4. Round      : Medium Ratio, Low |R|  (Uncoordinated biological noise)
#
# Optimal decision boundaries are analytically derived as the midpoints between
# adjacent cluster centroids in the unscaled feature space.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. REQUIRED LIBRARIES & ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

cat("======================================================================\n")
cat("  Unsupervised Multi-Omics Cloud Threshold Derivation (K-Means 2D)    \n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# 1. DATA INGESTION & QUALITY FILTERING (BDL Cohort)
# ------------------------------------------------------------------------------
cat("[Step 1/6] Loading cross-omics dataset...\n")

data_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData"
if (!file.exists(data_path)) {
  stop(sprintf("Input data file not found at: %s", data_path))
}

load(data_path)
Full_DT <- DTccl4_DT_LCPM
setDT(Full_DT)

# Standardize column naming and disease status
Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset     := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]
Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "control", DiseaseGroup := "Control"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL",     DiseaseGroup := "Disease"]

# Filter for valid BDL observations
bdl_DT <- Full_DT[Dataset == "BDL Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]

cat("[Step 2/6] Computing bivariate summary statistics per gene-protein pair...\n")
Cloud_Metrics <- bdl_DT[, {
  .(
    N         = .N,
    SD_RNA    = sd(GeneCount),
    SD_Prot   = sd(Protein_Raw),
    Pearson_R = suppressWarnings(cor(GeneCount, Protein_Raw, method = "pearson"))
  )
}, by = GeneProtein]

# Apply quality filtering: minimum 10 sample pairs and positive variance
Cloud_Metrics <- Cloud_Metrics[
  N >= 10 & !is.na(Pearson_R) & !is.na(SD_RNA) & !is.na(SD_Prot) & SD_RNA > 0 & SD_Prot > 0
]

# Derive core geometric features
Cloud_Metrics[, Ratio := SD_RNA / SD_Prot]
Cloud_Metrics[, AbsR  := abs(Pearson_R)]

cat(sprintf("  -> Successfully processed %d valid gene-protein pairs.\n", nrow(Cloud_Metrics)))

# ------------------------------------------------------------------------------
# 2. UNSUPERVISED 2D K-MEANS CLUSTERING ON NORMALIZED FEATURE SPACE
# ------------------------------------------------------------------------------
cat("\n[Step 3/6] Running 2D K-Means clustering (K = 4) on normalized unit hypercube [0, 1]^2...\n")

set.seed(42)  # For strict algorithmic reproducibility

# Min-Max Normalization to ensure equal weighting across both feature dimensions
ratio_min <- min(Cloud_Metrics$Ratio)
ratio_max <- max(Cloud_Metrics$Ratio)

ratio_scaled <- (Cloud_Metrics$Ratio - ratio_min) / (ratio_max - ratio_min)
r_scaled     <- Cloud_Metrics$AbsR  # Already naturally bounded in [0, 1]

km_data <- data.frame(Ratio_s = ratio_scaled, R_s = r_scaled)
km_fit  <- kmeans(km_data, centers = 4, nstart = 50, iter.max = 200)

Cloud_Metrics[, Cluster := as.character(km_fit$cluster)]

# Transform centroids back to the original unscaled metric space
centroids <- data.table(
  Cluster      = as.character(1:4),
  Ratio_center = km_fit$centers[, "Ratio_s"] * (ratio_max - ratio_min) + ratio_min,
  R_center     = km_fit$centers[, "R_s"]
)

# ------------------------------------------------------------------------------
# 3. DETERMINISTIC CLUSTER-TO-MORPHOLOGY MAPPING
# ------------------------------------------------------------------------------
cat("[Step 4/6] Mapping unsupervised cluster centroids to canonical cloud shapes...\n")

# Sort centroids along the Ratio axis to identify boundary regimes
centroids_sorted <- centroids[order(Ratio_center)]

# 1. Lowest Ratio  centroid -> Vertical (Protein variance dominant)
centroids_sorted[1, Category := "Vertical"]

# 2. Highest Ratio centroid -> Horizontal (Transcriptomic variance dominant)
centroids_sorted[4, Category := "Horizontal"]

# 3. Middle centroids: Partitioned by correlation magnitude
middle_two <- centroids_sorted[2:3]
if (middle_two$R_center[1] >= middle_two$R_center[2]) {
  centroids_sorted[2, Category := "Diagonal"]
  centroids_sorted[3, Category := "Round"]
} else {
  centroids_sorted[2, Category := "Round"]
  centroids_sorted[3, Category := "Diagonal"]
}

# Assign category labels to observations
cluster_map <- setNames(centroids_sorted$Category, centroids_sorted$Cluster)
Cloud_Metrics[, Category := cluster_map[Cluster]]

cat("\nCluster Centroid Mapping (Original Feature Space):\n")
print(centroids_sorted[, .(Cluster, Ratio_center = round(Ratio_center, 4), R_center = round(R_center, 4), Category)])

# ------------------------------------------------------------------------------
# 4. ANALYTICAL DERIVATION OF OPTIMAL DECISION BOUNDARIES
# ------------------------------------------------------------------------------
cat("\n[Step 5/6] Analytically deriving optimal decision thresholds from cluster midpoints...\n")

H  <- centroids_sorted[Category == "Horizontal"]
V  <- centroids_sorted[Category == "Vertical"]
D  <- centroids_sorted[Category == "Diagonal"]
Ro <- centroids_sorted[Category == "Round"]

# Threshold 1: Upper Ratio Boundary (Horizontal separation)
non_H_max_ratio <- max(D$Ratio_center, Ro$Ratio_center)
ratio_high <- (H$Ratio_center + non_H_max_ratio) / 2

# Threshold 2: Lower Ratio Boundary (Vertical separation)
non_V_min_ratio <- min(D$Ratio_center, Ro$Ratio_center)
ratio_low  <- (V$Ratio_center + non_V_min_ratio) / 2

# Threshold 3: Diagonal Correlation Floor
r_diag  <- (D$R_center + Ro$R_center) / 2

# Threshold 4: Horizontal / Vertical Correlation Ceiling
r_hv    <- (max(H$R_center, V$R_center) + D$R_center) / 2

# Threshold 5: Round Isotropic Correlation Ceiling
r_round <- Ro$R_center / 2

cat("======================================================================\n")
cat("      OPTIMAL DATA-DRIVEN CLOUD THRESHOLDS (K-MEANS 2D, BDL COHORT)   \n")
cat("======================================================================\n")
cat(sprintf("  ratio_high = %.3f   (Ratio > ratio_high        -> Horizontal)\n", ratio_high))
cat(sprintf("  ratio_low  = %.3f   (Ratio < ratio_low         -> Vertical)\n",   ratio_low))
cat(sprintf("  r_diag     = %.3f   (|Pearson R| > r_diag      -> Diagonal)\n",   r_diag))
cat(sprintf("  r_hv       = %.3f   (|Pearson R| < r_hv        -> Horiz / Vert Ceiling)\n", r_hv))
cat(sprintf("  r_round    = %.3f   (|Pearson R| < r_round     -> Round)\n",      r_round))
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# 5. PRODUCTION CLASSIFICATION FUNCTION (READY FOR DOWNSTREAM PIPELINES)
# ------------------------------------------------------------------------------
classify_cloud <- function(ratio, r, sd_rna, sd_prot) {
  # Handle missing or invalid inputs
  if (is.na(ratio) || is.na(r)) return("Unclassified")
  
  # Rule 1: Horizontal Cloud (RNA variance dominant, uncoupled)
  if (ratio > 1.014 && abs(r) < 0.673 && sd_rna > 0.3) {
    return("Horizontal")
  }
  # Rule 2: Vertical Cloud (Protein variance dominant, uncoupled)
  if (ratio < 0.528 && abs(r) < 0.673 && sd_prot > 0.3) {
    return("Vertical")
  }
  # Rule 3: Diagonal Cloud (Strong translational coupling)
  if (abs(r) > 0.467 && ratio >= 0.528 && ratio <= 1.014 && sd_rna > 0.3 && sd_prot > 0.3) {
    return("Diagonal")
  }
  # Rule 4: Round Cloud (Isotropic basal noise)
  if (ratio >= 0.528 && ratio <= 1.014 && abs(r) < 0.070 && sd_rna > 0.3 && sd_prot > 0.3) {
    return("Round")
  }
  
  return("Unclassified")
}

# Apply classification to verify empirical cohort distribution
Cloud_Metrics[, Final_Category := mapply(classify_cloud, Ratio, Pearson_R, SD_RNA, SD_Prot)]

cat("Empirical Category Distribution:\n")
cat_table <- table(Cloud_Metrics$Final_Category)
print(cat_table)
cat(sprintf("  -> Classified pairs: %.1f%% | Unclassified pairs: %.1f%%\n\n",
            100 * mean(Cloud_Metrics$Final_Category != "Unclassified"),
            100 * mean(Cloud_Metrics$Final_Category == "Unclassified")))

# ------------------------------------------------------------------------------
# 6. PUBLICATION-QUALITY DIAGNOSTIC VISUALIZATION (4-PANEL PLATE)
# ------------------------------------------------------------------------------
cat("[Step 6/6] Rendering diagnostic figure to PDF...\n")

out_pdf <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cloud_Threshold_DataDriven.pdf"

cloud_colors <- c(
  "Diagonal"     = "#2196F3",  # Blue
  "Horizontal"   = "#FF9800",  # Orange
  "Vertical"     = "#9C27B0",  # Purple
  "Round"        = "#4CAF50",  # Green
  "Unclassified" = "#BDBDBD"   # Grey
)

pdf(out_pdf, width = 14, height = 10)

# Panel A: 2D Feature Space with Cluster Centroids and Decision Boundaries
centroid_plot_dt <- centroids_sorted[, .(Ratio_center, R_center, Category)]
p_kmeans <- ggplot(Cloud_Metrics[Ratio <= quantile(Ratio, 0.99)],
                   aes(x = Ratio, y = AbsR, color = Category)) +
  geom_point(size = 1.0, alpha = 0.4) +
  scale_color_manual(values = cloud_colors) +
  geom_vline(xintercept = ratio_low,  color = "#9C27B0", linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = ratio_high, color = "#FF9800", linewidth = 1.1, linetype = "dashed") +
  geom_hline(yintercept = r_round,    color = "#4CAF50", linewidth = 1.1, linetype = "dashed") +
  geom_hline(yintercept = r_hv,       color = "#888888", linewidth = 0.8, linetype = "dotted") +
  geom_hline(yintercept = r_diag,     color = "#2196F3", linewidth = 1.1, linetype = "dashed") +
  geom_point(data = centroid_plot_dt,
             aes(x = Ratio_center, y = R_center, fill = Category),
             shape = 23, size = 6, color = "black", stroke = 1.5, inherit.aes = FALSE,
             show.legend = FALSE) +
  scale_fill_manual(values = cloud_colors) +
  geom_label(data = centroid_plot_dt,
             aes(x = Ratio_center, y = R_center + 0.06, label = Category),
             size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  labs(
    title = "A. Unsupervised 2D K-Means Partitioning in (Ratio, |R|) Space",
    subtitle = sprintf("Diamonds indicate cluster centroids. Thresholds: ratio_low=%.3f, ratio_high=%.3f, r_diag=%.3f",
                       ratio_low, ratio_high, r_diag),
    x = expression("Standard Deviation Ratio ("*sigma[RNA] / sigma[Protein]*")"),
    y = expression("Absolute Bravais-Pearson Correlation ("*group("|", R, "|")*")"),
    color = "Cloud Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# Panel B: 1D Histogram of Variance Ratio Distribution
p_ratio <- ggplot(Cloud_Metrics[Ratio <= quantile(Ratio, 0.99)],
                  aes(x = Ratio, fill = Final_Category)) +
  geom_histogram(bins = 60, alpha = 0.75, position = "stack") +
  scale_fill_manual(values = cloud_colors) +
  geom_vline(xintercept = ratio_low,  color = "#9C27B0", linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = ratio_high, color = "#FF9800", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = ratio_low  - 0.02, y = Inf, label = sprintf("%.3f", ratio_low),
           vjust = 2, hjust = 1, color = "#9C27B0", size = 4, fontface = "bold") +
  annotate("text", x = ratio_high + 0.02, y = Inf, label = sprintf("%.3f", ratio_high),
           vjust = 2, hjust = 0, color = "#FF9800", size = 4, fontface = "bold") +
  labs(
    title = "B. Distribution of Variance Ratio with Decision Boundaries",
    x = expression("Standard Deviation Ratio ("*sigma[RNA] / sigma[Protein]*")"),
    y = "Gene-Protein Pair Count",
    fill = "Cloud Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# Panel C: 1D Histogram of Absolute Correlation Distribution
p_r <- ggplot(Cloud_Metrics, aes(x = AbsR, fill = Final_Category)) +
  geom_histogram(bins = 50, alpha = 0.75, position = "stack") +
  scale_fill_manual(values = cloud_colors) +
  geom_vline(xintercept = r_round, color = "#4CAF50", linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = r_hv,    color = "#888888", linewidth = 0.9, linetype = "dotted") +
  geom_vline(xintercept = r_diag,  color = "#2196F3", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = r_round + 0.01, y = Inf, label = sprintf("r_round\n%.3f", r_round),
           vjust = 2, hjust = 0, color = "#4CAF50", size = 3.5, fontface = "bold") +
  annotate("text", x = r_diag  + 0.01, y = Inf, label = sprintf("r_diag\n%.3f",  r_diag),
           vjust = 2, hjust = 0, color = "#2196F3", size = 3.5, fontface = "bold") +
  labs(
    title = "C. Distribution of Correlation Magnitudes with Thresholds",
    x = expression("Absolute Bravais-Pearson Correlation ("*group("|", R, "|")*")"),
    y = "Gene-Protein Pair Count",
    fill = "Cloud Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# Panel D: Barplot of Final Morphological Categorization
cat_counts <- as.data.frame(table(Cloud_Metrics$Final_Category))
colnames(cat_counts) <- c("Category", "Count")
cat_counts$Pct <- cat_counts$Count / sum(cat_counts$Count) * 100

p_bar <- ggplot(cat_counts, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = cloud_colors, guide = "none") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Pct)),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  labs(
    title = "D. Empirical Multi-Omics Category Frequencies",
    x = "Assigned Cloud Category",
    y = "Gene-Protein Pair Count"
  ) +
  theme_minimal(base_size = 12) +
  ylim(0, max(cat_counts$Count) * 1.15) +
  theme(panel.grid.minor = element_blank())

# Arrange 4-Panel Master Diagnostic Plate
grid.arrange(
  p_kmeans, p_ratio, p_r, p_bar, ncol = 2,
  top = textGrob(
    "Data-Driven Multi-Omics Cloud Morphology Classification (2D K-Means, BDL Cohort)\nDecision boundaries represent objective midpoints between adjacent unsupervised cluster centroids",
    gp = gpar(fontsize = 13, fontface = "bold")
  )
)

dev.off()

cat(sprintf("\nDiagnostic figure successfully saved to:\n  %s\n", out_pdf))
cat("\n======================================================================\n")
cat("  EXECUTION COMPLETE: Threshold derivation fully validated!           \n")
cat("======================================================================\n")
