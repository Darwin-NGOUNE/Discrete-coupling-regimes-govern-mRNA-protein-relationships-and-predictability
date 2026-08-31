# ==============================================================================
# SCRIPT: correlation_heatmaps_mrna_shared_order.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/correlation_heatmaps_mrna_shared_order.R
# PURPOSE: 1. Compute empirical mRNA correlation matrices R^(1) (BDL) and R^(2) (CCl4) [943 x 943]
#          2. Perform joint hierarchical clustering on R^(mean) = (R^(1) + R^(2)) / 2
#          3. Render Heatmap BDL, Heatmap CCl4, and Difference Heatmap Delta R = R^(1) - R^(2)
#             all with IDENTICAL shared mRNA ordering.
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

cat("1. Loading batch-corrected data and extracting 943 shared core genes/mRNAs...\n")

load_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData"
load(load_path)
Full_DT <- as.data.table(Full_DT)

# Separate Datasets
bdl_dt  <- Full_DT[is.na(TreatmentTime)]
ccl4_dt <- Full_DT[!is.na(TreatmentTime)]
core_genes <- intersect(unique(bdl_dt$GeneProtein), unique(ccl4_dt$GeneProtein))
n_genes <- length(core_genes)
cat(sprintf("Found %d shared core genes/mRNAs.\n", n_genes))

# Subset and cast to wide matrix (samples x genes)
bdl_sub  <- bdl_dt[GeneProtein %in% core_genes]
ccl4_sub <- ccl4_dt[GeneProtein %in% core_genes]

# Extract wide matrices using ComBat_GeneCount
mat_bdl_wide  <- dcast(bdl_sub, Sample_ID ~ GeneProtein, value.var = "ComBat_GeneCount")
mat_ccl4_wide <- dcast(ccl4_sub, Sample_ID ~ GeneProtein, value.var = "ComBat_GeneCount")

mat_bdl  <- as.matrix(mat_bdl_wide[, -1, with = FALSE])
mat_ccl4 <- as.matrix(mat_ccl4_wide[, -1, with = FALSE])

# Ensure exact identical column ordering
mat_bdl  <- mat_bdl[, core_genes]
mat_ccl4 <- mat_ccl4[, core_genes]

cat("2. Computing empirical Pearson correlation matrices R^(1) and R^(2) for mRNA...\n")
R1_mrna <- cor(mat_bdl,  method = "pearson", use = "pairwise.complete.obs") # BDL (943 x 943)
R2_mrna <- cor(mat_ccl4, method = "pearson", use = "pairwise.complete.obs") # CCl4 (943 x 943)

cat("3. Computing consensus correlation matrix and joint hierarchical clustering for mRNA...\n")
R_mean_mrna <- (R1_mrna + R2_mrna) / 2

# Distance matrix based on correlation distance: D = 1 - R_mean
dist_mat_mrna <- as.dist(1 - R_mean_mrna)
hc_mrna <- hclust(dist_mat_mrna, method = "ward.D2")
shared_order_mrna <- hc_mrna$order
ordered_genes <- core_genes[shared_order_mrna]

# Reorder all matrices identically
R1_ordered_mrna <- R1_mrna[shared_order_mrna, shared_order_mrna]
R2_ordered_mrna <- R2_mrna[shared_order_mrna, shared_order_mrna]
Delta_R_mrna    <- R1_ordered_mrna - R2_ordered_mrna

cat("Summary of mRNA Delta R (R^(1) - R^(2)):\n")
cat("Min:", min(Delta_R_mrna), "Median:", median(Delta_R_mrna), "Max:", max(Delta_R_mrna), "\n")

cat("Converting mRNA matrices to long format for ggplot raster plotting...\n")
mat_to_dt <- function(mat, val_name) {
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("Gene1", "Gene2", val_name))
  dt[, `:=`(
    x_rank = as.integer(factor(Gene2, levels = ordered_genes)),
    y_rank = as.integer(factor(Gene1, levels = rev(ordered_genes)))
  )]
  return(dt)
}

dt_r1_mrna    <- mat_to_dt(R1_ordered_mrna, "Correlation")
dt_r2_mrna    <- mat_to_dt(R2_ordered_mrna, "Correlation")
dt_delta_mrna <- mat_to_dt(Delta_R_mrna, "Delta")

# ------------------------------------------------------------------------------
# 4. PLOTTING PUBLICATION-GRADE mRNA HEATMAPS
# ------------------------------------------------------------------------------
out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

theme_heatmap_mrna <- function() {
  theme_minimal(base_size = 14, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black", margin = margin(b = 8)),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
      legend.position = "none", # Legend removed for mRNA as it is combined with Protein
      plot.margin = margin(t = 5, r = 8, b = 5, l = 8)
    )
}

# Heatmap 1: BDL mRNA
p_bdl_mrna <- ggplot(dt_r1_mrna, aes(x = x_rank, y = y_rank, fill = Correlation)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(
    title = "BDL mRNA data",
    subtitle = NULL
  ) +
  theme_heatmap_mrna()

# Heatmap 2: CCl4 mRNA
p_ccl4_mrna <- ggplot(dt_r2_mrna, aes(x = x_rank, y = y_rank, fill = Correlation)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(
    title = expression(bold(CCl[4]~"mRNA data")),
    subtitle = NULL
  ) +
  theme_heatmap_mrna()

# Heatmap 3: Difference (Delta R mRNA)
p_delta_mrna <- ggplot(dt_delta_mrna, aes(x = x_rank, y = y_rank, fill = Delta)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2980B9", mid = "#FFFFFF", high = "#E74C3C",
    midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish
  ) +
  labs(
    title = expression(bold("Difference of BDL on "~CCl[4]~"data")),
    subtitle = NULL
  ) +
  theme_heatmap_mrna()

# ------------------------------------------------------------------------------
# 5. SAVE HIGH-RES TRIPTYCH (NO TOP TITLE, NO BOTTOM LEGEND)
# ------------------------------------------------------------------------------
g_triptych_mrna <- gridExtra::arrangeGrob(
  p_bdl_mrna, p_ccl4_mrna, p_delta_mrna, ncol = 3
)

pdf_triptych_mrna <- file.path(out_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_mRNA_DeltaR.pdf")
cairo_pdf(pdf_triptych_mrna, width = 18, height = 5.8)
grid::grid.draw(g_triptych_mrna)
dev.off()
cat("Successfully generated:", pdf_triptych_mrna, "\n")

png_triptych_mrna <- file.path(out_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_mRNA_DeltaR.png")
png(png_triptych_mrna, width = 3600, height = 1500, res = 200)
grid::grid.draw(g_triptych_mrna)
dev.off()

cat("\nALL mRNA CORRELATION HEATMAPS COMPLETED SUCCESSFULLY!\n")

