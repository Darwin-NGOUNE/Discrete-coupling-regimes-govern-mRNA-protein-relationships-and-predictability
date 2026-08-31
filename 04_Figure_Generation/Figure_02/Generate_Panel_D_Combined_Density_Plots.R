# ==============================================================================
# SCRIPT: Generate_Panel_D_Combined_Density_Plots.R
# PURPOSE: Build a unified 1x2 panel for Panel D where BDL (left) and CCl4 (right)
#          density plots have 100% identical dimensions, no individual titles,
#          bold axes, and a single shared legend centered at the bottom.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

cat("Loading data for BDL & CCl4 density curves...\n")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")

BDL_DT  <- DTccl4_DT_LCPM[is.na(TreatmentTime)]
CCL4_DT <- DTccl4_DT_LCPM[!is.na(TreatmentTime)]

extract_top_n_stats <- function(mat, n_val) {
  top_vals <- apply(mat, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE)
    return(sorted_x[n_val])
  })
  return(data.table(Protein = names(top_vals), Max_R = top_vals, Tier = paste0("Top ", n_val)))
}
tiers <- c(1, 2, 5, 10, 15, 20, 50, 70, 100)

# --- BDL ---
dat_wide_bdl <- dcast(BDL_DT, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
mat_bdl <- as.matrix(dat_wide_bdl[, -1, with = FALSE])
rownames(mat_bdl) <- dat_wide_bdl$GeneProtein
mat_bdl_clean <- na.omit(mat_bdl)
abs_cor_bdl <- abs(cor(t(mat_bdl_clean), method = "pearson"))
diag(abs_cor_bdl) <- 0
res_bdl <- rbindlist(lapply(tiers, function(n) extract_top_n_stats(abs_cor_bdl, n)))

# --- CCL4 ---
ccl4_clean_dt <- na.omit(CCL4_DT)
dat_wide_ccl4 <- dcast(ccl4_clean_dt, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
mat_ccl4 <- as.matrix(dat_wide_ccl4[, -1, with = FALSE])
rownames(mat_ccl4) <- dat_wide_ccl4$GeneProtein
mat_ccl4_clean <- na.omit(mat_ccl4)
abs_cor_ccl4 <- abs(cor(t(mat_ccl4_clean), method = "pearson"))
diag(abs_cor_ccl4) <- 0
res_ccl4 <- rbindlist(lapply(tiers, function(n) extract_top_n_stats(abs_cor_ccl4, n)))

# --- GGPLOTS WITHOUT LEGENDS ---
theme_dens <- function() {
  theme_bw(base_size = 18) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 20, face = "bold", color = "black"),
      axis.text = element_text(size = 16, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 1.0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 8, b = 5, l = 8)
    )
}

p_dens_bdl <- ggplot(res_bdl, aes(x = Max_R, fill = factor(Tier, levels = paste0("Top ", tiers)))) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(title = NULL, x = expression(bold("|"*rho[BP]*"|")), y = "Density") +
  theme_dens()

p_dens_ccl4 <- ggplot(res_ccl4, aes(x = Max_R, fill = factor(Tier, levels = paste0("Top ", tiers)))) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(title = NULL, x = expression(bold("|"*rho[BP]*"|")), y = "Density") +
  theme_dens()

# --- DUMMY PLOT FOR SINGLE CENTERED SHARED LEGEND ---
dummy_df <- data.frame(x = 1:length(tiers), y = 1:length(tiers), Tier = factor(paste0("Top ", tiers), levels = paste0("Top ", tiers)))
p_leg_dummy <- ggplot(dummy_df, aes(x, y, fill = Tier)) +
  geom_tile() +
  scale_fill_brewer(palette = "YlOrRd", direction = -1, name = "Correlation rank") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 16, color = "black", margin = margin(r = 10)),
    legend.text = element_text(face = "bold", size = 14, color = "black"),
    legend.key.size = unit(0.55, "cm")
  )
shared_legend <- get_legend(p_leg_dummy)

# Assemble 1x2 row + shared legend beneath
top_plots <- plot_grid(p_dens_bdl, p_dens_ccl4, nrow = 1, rel_widths = c(1, 1))
panel_d_combined <- plot_grid(top_plots, shared_legend, ncol = 1, rel_heights = c(1, 0.15))

out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
out_file <- file.path(out_dir, "Panel_D_Combined_Density_Plots.pdf")
cairo_pdf(out_file, width = 16, height = 8)
grid::grid.draw(panel_d_combined)
dev.off()

cat(sprintf("SUCCESS! Perfect unified Panel D saved to: %s\n", out_file))
