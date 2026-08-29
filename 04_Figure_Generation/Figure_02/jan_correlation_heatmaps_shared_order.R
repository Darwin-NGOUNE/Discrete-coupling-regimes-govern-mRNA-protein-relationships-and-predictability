# ==============================================================================
# SCRIPT: jan_correlation_heatmaps_shared_order.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/jan_correlation_heatmaps_shared_order.R
# PURPOSE: 1. Compute empirical correlation matrices R^(1) (BDL) and R^(2) (CCl4) [943 x 943]
#          2. Perform joint hierarchical clustering on R^(mean) = (R^(1) + R^(2)) / 2
#          3. Render Heatmap BDL, Heatmap CCl4, and Difference Heatmap Delta R = R^(1) - R^(2)
#             all with IDENTICAL shared protein ordering.
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

cat("1. Loading batch-corrected data and extracting 943 shared core proteins...\n")

load_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData"
load(load_path)
Full_DT <- as.data.table(Full_DT)

# Separate Datasets
bdl_dt  <- Full_DT[is.na(TreatmentTime)]
ccl4_dt <- Full_DT[!is.na(TreatmentTime)]
core_proteins <- intersect(unique(bdl_dt$GeneProtein), unique(ccl4_dt$GeneProtein))
n_prot <- length(core_proteins)
cat(sprintf("Found %d shared core proteins.\n", n_prot))

# Subset and cast to wide matrix (samples x proteins)
bdl_sub  <- bdl_dt[GeneProtein %in% core_proteins]
ccl4_sub <- ccl4_dt[GeneProtein %in% core_proteins]

# Extract wide matrices using ComBat_Protein_Raw
mat_bdl_wide  <- dcast(bdl_sub, Sample_ID ~ GeneProtein, value.var = "ComBat_Protein_Raw")
mat_ccl4_wide <- dcast(ccl4_sub, Sample_ID ~ GeneProtein, value.var = "ComBat_Protein_Raw")

mat_bdl  <- as.matrix(mat_bdl_wide[, -1, with = FALSE])
mat_ccl4 <- as.matrix(mat_ccl4_wide[, -1, with = FALSE])

# Ensure exact identical column ordering
mat_bdl  <- mat_bdl[, core_proteins]
mat_ccl4 <- mat_ccl4[, core_proteins]

cat("2. Computing empirical Pearson correlation matrices R^(1) and R^(2)...\n")
R1 <- cor(mat_bdl,  method = "pearson", use = "pairwise.complete.obs") # BDL (943 x 943)
R2 <- cor(mat_ccl4, method = "pearson", use = "pairwise.complete.obs") # CCl4 (943 x 943)

cat("3. Computing consensus correlation matrix and joint hierarchical clustering...\n")
R_mean <- (R1 + R2) / 2

# Distance matrix based on correlation distance: D = 1 - R_mean
dist_mat <- as.dist(1 - R_mean)
hc <- hclust(dist_mat, method = "ward.D2")
shared_order <- hc$order
ordered_proteins <- core_proteins[shared_order]

# Reorder all matrices identically
R1_ordered <- R1[shared_order, shared_order]
R2_ordered <- R2[shared_order, shared_order]
Delta_R    <- R1_ordered - R2_ordered

cat("Converting matrices to long format for ggplot raster plotting...\n")
# Function to convert matrix to long data.table with numeric coordinate ranks (1 to 943)
mat_to_dt <- function(mat, val_name) {
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("Prot1", "Prot2", val_name))
  dt[, `:=`(
    x_rank = as.integer(factor(Prot2, levels = ordered_proteins)),
    y_rank = as.integer(factor(Prot1, levels = rev(ordered_proteins))) # reverse so top-left is (1,1)
  )]
  return(dt)
}

dt_r1    <- mat_to_dt(R1_ordered, "Correlation")
dt_r2    <- mat_to_dt(R2_ordered, "Correlation")
dt_delta <- mat_to_dt(Delta_R, "Delta")

# ------------------------------------------------------------------------------
# 4. PLOTTING PUBLICATION-GRADE HEATMAPS
# ------------------------------------------------------------------------------
out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/Scatterplots_75th_Percentile"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

theme_heatmap <- function() {
  theme_minimal(base_size = 13, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "#1B4F72", margin = margin(b = 6)),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#566573", margin = margin(b = 10)),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#2C3E50", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11, color = "#1B4F72"),
      legend.text = element_text(size = 9.5),
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.35, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# Heatmap A: BDL
p_bdl <- ggplot(dt_r1, aes(x = x_rank, y = y_rank, fill = Correlation)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1),
    name = expression(bold(rho[BP] ~ "(BDL)"))
  ) +
  labs(
    title = expression(bold("A. Datensatz 1: BDL (") * bold(R^{"(1)"}) * bold(")")),
    subtitle = "Shared Hierarchical Clustering Order"
  ) +
  theme_heatmap()

# Heatmap B: CCl4
p_ccl4 <- ggplot(dt_r2, aes(x = x_rank, y = y_rank, fill = Correlation)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1),
    name = expression(bold(rho[BP] ~ (CCl[4])))
  ) +
  labs(
    title = expression(bold("B. Datensatz 2: ") * bold(CCl[4]) * bold(" (") * bold(R^{"(2)"}) * bold(")")),
    subtitle = "Identical Shared Protein Order"
  ) +
  theme_heatmap()

# Heatmap C: Difference (Delta R)
p_delta <- ggplot(dt_delta, aes(x = x_rank, y = y_rank, fill = Delta)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2980B9", mid = "#FFFFFF", high = "#E74C3C",
    midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
    name = expression(bold(Delta * R[ij] ~ (BDL - CCl[4])))
  ) +
  labs(
    title = expression(bold("C. Differenz-Heatmap: ") * bold(Delta * R == R^{"(1)"} - R^{"(2)"})),
    subtitle = "Blue: CCl4 stronger | White: Equal | Red: BDL stronger"
  ) +
  theme_heatmap() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "#78281F", margin = margin(b = 6))
  )

# ------------------------------------------------------------------------------
# 5. SAVE HIGH-RES TRIPTYCH & INDIVIDUAL PLOTS
# ------------------------------------------------------------------------------
# 1. Triptych Figure (3 Heatmaps Side-by-Side)
g_triptych <- gridExtra::arrangeGrob(
  p_bdl, p_ccl4, p_delta, ncol = 3,
  top = textGrob("Korrelations-Heatmaps mit identischer Proteinreihenfolge (N = 943 Proteine)", 
                 gp = gpar(fontsize = 16, fontface = "bold", col = "#1B4F72"))
)

pdf_triptych <- file.path(out_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_DeltaR.pdf")
cairo_pdf(pdf_triptych, width = 18, height = 7.2)
grid::grid.draw(g_triptych)
dev.off()
cat("Successfully generated:", pdf_triptych, "\n")

png_triptych <- file.path(out_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_DeltaR.png")
png(png_triptych, width = 3600, height = 1500, res = 200)
grid::grid.draw(g_triptych)
dev.off()

cat("\nALL CORRELATION HEATMAPS COMPLETED SUCCESSFULLY!\n")

