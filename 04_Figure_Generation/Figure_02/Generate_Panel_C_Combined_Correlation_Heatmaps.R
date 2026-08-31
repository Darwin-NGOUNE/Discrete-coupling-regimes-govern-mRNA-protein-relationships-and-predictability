# ==============================================================================
# SCRIPT: Generate_Panel_C_Combined_Correlation_Heatmaps.R
# PURPOSE: Build a unified 2x3 matrix for Panel C where RNA (Row 1) and Protein (Row 2)
#          have 100% identical square sizes, perfectly aligned columns, and 3 dedicated
#          colorbars directly under BDL, CCl4, and Difference columns.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

cat("1. Loading batch-corrected data for RNA & Protein...\n")
load_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData"
load(load_path)
Full_DT <- as.data.table(Full_DT)

bdl_dt  <- Full_DT[is.na(TreatmentTime)]
ccl4_dt <- Full_DT[!is.na(TreatmentTime)]
core_features <- intersect(unique(bdl_dt$GeneProtein), unique(ccl4_dt$GeneProtein))
n_feat <- length(core_features)

# ==============================================================================
# 2. COMPUTE CORRELATION MATRICES
# ==============================================================================
# --- RNA ---
mat_bdl_rna_wide  <- dcast(bdl_dt[GeneProtein %in% core_features], Sample_ID ~ GeneProtein, value.var = "ComBat_GeneCount")
mat_ccl4_rna_wide <- dcast(ccl4_dt[GeneProtein %in% core_features], Sample_ID ~ GeneProtein, value.var = "ComBat_GeneCount")
mat_bdl_rna  <- as.matrix(mat_bdl_rna_wide[, -1, with = FALSE])[, core_features]
mat_ccl4_rna <- as.matrix(mat_ccl4_rna_wide[, -1, with = FALSE])[, core_features]

R1_mrna <- cor(mat_bdl_rna,  method = "pearson", use = "pairwise.complete.obs")
R2_mrna <- cor(mat_ccl4_rna, method = "pearson", use = "pairwise.complete.obs")
R_mean_mrna <- (R1_mrna + R2_mrna) / 2
hc_mrna <- hclust(as.dist(1 - R_mean_mrna), method = "ward.D2")
order_mrna <- hc_mrna$order
ordered_genes <- core_features[order_mrna]

R1_mrna_ord <- R1_mrna[order_mrna, order_mrna]
R2_mrna_ord <- R2_mrna[order_mrna, order_mrna]
Delta_mrna  <- R1_mrna_ord - R2_mrna_ord

# --- PROTEIN ---
mat_bdl_prot_wide  <- dcast(bdl_dt[GeneProtein %in% core_features], Sample_ID ~ GeneProtein, value.var = "ComBat_Protein_Raw")
mat_ccl4_prot_wide <- dcast(ccl4_dt[GeneProtein %in% core_features], Sample_ID ~ GeneProtein, value.var = "ComBat_Protein_Raw")
mat_bdl_prot  <- as.matrix(mat_bdl_prot_wide[, -1, with = FALSE])[, core_features]
mat_ccl4_prot <- as.matrix(mat_ccl4_prot_wide[, -1, with = FALSE])[, core_features]

R1_prot <- cor(mat_bdl_prot,  method = "pearson", use = "pairwise.complete.obs")
R2_prot <- cor(mat_ccl4_prot, method = "pearson", use = "pairwise.complete.obs")
R_mean_prot <- (R1_prot + R2_prot) / 2
hc_prot <- hclust(as.dist(1 - R_mean_prot), method = "ward.D2")
order_prot <- hc_prot$order
ordered_prots <- core_features[order_prot]

R1_prot_ord <- R1_prot[order_prot, order_prot]
R2_prot_ord <- R2_prot[order_prot, order_prot]
Delta_prot  <- R1_prot_ord - R2_prot_ord

# ==============================================================================
# 3. HELPER TO CONVERT MATRIX TO GGPLOT HEATMAP
# ==============================================================================
mat_to_dt <- function(mat, val_name, ordered_items) {
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("Item1", "Item2", val_name))
  dt[, `:=`(
    x_rank = as.integer(factor(Item2, levels = ordered_items)),
    y_rank = as.integer(factor(Item1, levels = rev(ordered_items)))
  )]
  return(dt)
}

theme_square <- function() {
  theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black", margin = margin(b = 6)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      plot.margin = margin(t = 4, r = 6, b = 4, l = 6)
    )
}

# --- RNA Plots (Row 1) ---
dt_r1_m <- mat_to_dt(R1_mrna_ord, "Val", ordered_genes)
dt_r2_m <- mat_to_dt(R2_mrna_ord, "Val", ordered_genes)
dt_d_m  <- mat_to_dt(Delta_mrna,  "Val", ordered_genes)

p_rna_1 <- ggplot(dt_r1_m, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#1B4F72", mid = "#FFFFFF", high = "#922B21", midpoint = 0, limits = c(-1, 1), guide = "none") +
  labs(title = "BDL mRNA data") + theme_square()

p_rna_2 <- ggplot(dt_r2_m, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#1B4F72", mid = "#FFFFFF", high = "#922B21", midpoint = 0, limits = c(-1, 1), guide = "none") +
  labs(title = expression(bold(CCl[4]~"mRNA data"))) + theme_square()

p_rna_3 <- ggplot(dt_d_m, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#2980B9", mid = "#FFFFFF", high = "#E74C3C", midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish, guide = "none") +
  labs(title = expression(bold("Difference of BDL on "~CCl[4]~"data"))) + theme_square()

# --- Protein Plots (Row 2) ---
dt_r1_p <- mat_to_dt(R1_prot_ord, "Val", ordered_prots)
dt_r2_p <- mat_to_dt(R2_prot_ord, "Val", ordered_prots)
dt_d_p  <- mat_to_dt(Delta_prot,  "Val", ordered_prots)

p_prot_1 <- ggplot(dt_r1_p, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#1B4F72", mid = "#FFFFFF", high = "#922B21", midpoint = 0, limits = c(-1, 1), guide = "none") +
  labs(title = "BDL protein data") + theme_square()

p_prot_2 <- ggplot(dt_r2_p, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#1B4F72", mid = "#FFFFFF", high = "#922B21", midpoint = 0, limits = c(-1, 1), guide = "none") +
  labs(title = expression(bold(CCl[4]~"protein data"))) + theme_square()

p_prot_3 <- ggplot(dt_d_p, aes(x = x_rank, y = y_rank, fill = Val)) +
  geom_raster() + coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "#2980B9", mid = "#FFFFFF", high = "#E74C3C", midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish, guide = "none") +
  labs(title = expression(bold("Difference of BDL on "~CCl[4]~"data"))) + theme_square()

# --- 3 Dedicated Colorbars (Row 3) ---
# Colorbar 1: BDL rho_BP
dummy_df1 <- data.frame(x = 1:10, y = 1:10, z = seq(-1, 1, length.out = 10))
p_leg1 <- ggplot(dummy_df1, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = expression(bold(rho[BP]))
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18, color = "black", margin = margin(r = 10)),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.45, "cm")
  )
leg1 <- get_legend(p_leg1)

# Colorbar 2: CCl4 rho_BP
dummy_df2 <- data.frame(x = 1:10, y = 1:10, z = seq(-1, 1, length.out = 10))
p_leg2 <- ggplot(dummy_df2, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#1B4F72", mid = "#FFFFFF", high = "#922B21",
    midpoint = 0, limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = expression(bold(rho[BP]))
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18, color = "black", margin = margin(r = 10)),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.45, "cm")
  )
leg2 <- get_legend(p_leg2)

# Colorbar 3: Difference Delta rho_BP
dummy_df3 <- data.frame(x = 1:10, y = 1:10, z = seq(-1.5, 1.5, length.out = 10))
p_leg3 <- ggplot(dummy_df3, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2980B9", mid = "#FFFFFF", high = "#E74C3C",
    midpoint = 0, limits = c(-1.5, 1.5),
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
    name = expression(bold(Delta * rho[BP]))
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 18, color = "black", margin = margin(r = 10)),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.45, "cm")
  )
leg3 <- get_legend(p_leg3)

# ==============================================================================
# 4. ASSEMBLE 2x3 MASTER GRID WITH cowplot (STRICT 1:1 SQUARES + 3 LEGENDS)
# ==============================================================================
# Row 1 (RNA): 3 equal squares
row1 <- plot_grid(p_rna_1, p_rna_2, p_rna_3, nrow = 1, rel_widths = c(1, 1, 1))

# Row 2 (Protein): 3 equal squares
row2 <- plot_grid(p_prot_1, p_prot_2, p_prot_3, nrow = 1, rel_widths = c(1, 1, 1))

# Row 3 (Legends): 3 dedicated colorbars aligned directly under the 3 columns
row3 <- plot_grid(leg1, leg2, leg3, nrow = 1, rel_widths = c(1, 1, 1))

# Combine all 3 rows
panel_c_combined <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1, 0.20))

out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
out_file <- file.path(out_dir, "Panel_C_Combined_Correlation_Heatmaps.pdf")
cairo_pdf(out_file, width = 16, height = 11)
grid::grid.draw(panel_c_combined)
dev.off()

cat(sprintf("SUCCESS! Perfect unified Panel C with 3 legends saved to: %s\n", out_file))
