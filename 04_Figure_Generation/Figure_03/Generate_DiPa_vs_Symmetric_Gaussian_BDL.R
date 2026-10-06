# ==============================================================================
# SCRIPT: Generate_DiPa_vs_Symmetric_Gaussian_BDL.R
# OBJECTIVE: Comparative visualization of Observed DiPa vs. Symmetric Bivariate
#            Independence Gaussian Model (Null Model proposed by Prof. Andreas Groll).
#
# PANEL A (Left) : Observed BDL Data (RNA log2FC vs Protein log2FC, n = 943)
# PANEL B (Right): Symmetric Bivariate Gaussian Model with:
#                  mu = (0, 0)'
#                  sd_pooled = sqrt((s_1^2 + s_2^2)/2)
#                  rho_BP = 0 (Independence Copula)
#                  Sample size n = 943
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Fix random seed for exact reproducibility of the Gaussian simulation
set.seed(42)

# ------------------------------------------------------------------------------
# 1. LOAD AND PREPARE OBSERVED BDL DATA (n = 943 common genes)
# ------------------------------------------------------------------------------
data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share"

message("Loading DiPa and expression data...")
load(file.path(data_dir, "New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData"))
load(file.path(data_dir, "Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData"))

# Filter valid DiPa entries for BDL
DT_dipa_valid <- copy(DT_dipa_count_bdl)
setDT(DT_dipa_valid)
DT_dipa_valid[, log2G := log2(meanRatioG)]
DT_dipa_valid[, log2P := log2(meanRatioP)]
DT_dipa_valid <- DT_dipa_valid[is.finite(log2G) & is.finite(log2P) & !is.na(DiPaGroups)]
DT_dipa_valid[, DiPaGroups := as.character(DiPaGroups)]
DT_dipa_valid[DiPaGroups == "0", DiPaGroups := "8"]

# Prepare individual BDL data for variance and correlation
Full_DT <- DTccl4_DT_LCPM
setDT(Full_DT)
Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]
Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "control", DiseaseGroup := "Control"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL",     DiseaseGroup := "Disease"]
bdl_DT <- Full_DT[Dataset == "BDL Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]

Cloud_Metrics_BDL <- bdl_DT[, {
  .(
    N = .N,
    SD_RNA = sd(GeneCount),
    SD_Prot = sd(Protein_Raw),
    Pearson_R = suppressWarnings(cor(GeneCount, Protein_Raw))
  )
}, by = GeneProtein]
Cloud_Metrics_BDL <- Cloud_Metrics_BDL[N >= 10 & !is.na(Pearson_R) & SD_RNA > 0 & SD_Prot > 0]

# Cloud classification function
classify_cloud <- function(ratio, r, sd_rna, sd_prot) {
  if (is.na(ratio) | is.na(r)) return("Unclassified")
  if (ratio > 1.014  & abs(r) < 0.673) return("Horizontal")
  if (ratio < 0.528  & abs(r) < 0.673) return("Vertical")
  if (abs(r) > 0.467 & ratio >= 0.528 & ratio <= 1.014) return("Diagonal")
  if (ratio >= 0.528 & ratio <= 1.014 & abs(r) < 0.070) return("Round")
  return("Unclassified")
}

Cloud_Metrics_BDL[, Ratio := SD_RNA / SD_Prot]
Cloud_Metrics_BDL[, CloudCategory := mapply(classify_cloud, Ratio, Pearson_R, SD_RNA, SD_Prot)]

# Merge to obtain exact 943 genes
Combined <- merge(
  DT_dipa_valid[, .(GeneProtein, log2G, log2P, DiPaGroups, meanRatioG, meanRatioP)],
  Cloud_Metrics_BDL[, .(GeneProtein, CloudCategory, Pearson_R, Ratio, SD_RNA, SD_Prot)],
  by = "GeneProtein"
)

n_genes <- nrow(Combined)
message(sprintf("Loaded %d common genes successfully.", n_genes))

# ------------------------------------------------------------------------------
# 2. STATISTICAL CALCULATIONS: EMPIRICAL & POOLED STANDARD DEVIATION
# ------------------------------------------------------------------------------
mRN_log2_fc  <- Combined$log2G
prot_log2_fc <- Combined$log2P

mean_rna  <- mean(mRN_log2_fc)
mean_prot <- mean(prot_log2_fc)

sd.1 <- sd(mRN_log2_fc)
sd.2 <- sd(prot_log2_fc)

# Pooled standard deviation (Andy's formula with df = n - 1 = 942)
sd.pooled <- sqrt(((n_genes - 1) * (sd.1^2 + sd.2^2)) / (2 * (n_genes - 1)))
obs_corr  <- cor(mRN_log2_fc, prot_log2_fc)

message(sprintf("RNA  Log2FC: Mean = %.4f | SD = %.4f", mean_rna, sd.1))
message(sprintf("Prot Log2FC: Mean = %.4f | SD = %.4f", mean_prot, sd.2))
message(sprintf("Pooled SD  : %.4f", sd.pooled))
message(sprintf("Observed rho_BP: %.4f", obs_corr))

# ------------------------------------------------------------------------------
# 3. SIMULATION: BIVARIATE SYMMETRIC INDEPENDENCE GAUSSIAN COPULA (NULL MODEL)
# ------------------------------------------------------------------------------
# Margins: Normal with mu = 0 and sd = sd.pooled
# Copula: Independence (rho_BP = 0)
x1.vec <- rnorm(n = n_genes, mean = 0, sd = sd.pooled)
x2.vec <- rnorm(n = n_genes, mean = 0, sd = sd.pooled)

sim_DT <- data.table(
  Sim_ID = 1:n_genes,
  log2G = x1.vec,
  log2P = x2.vec
)

# ------------------------------------------------------------------------------
# 4. PLOTTING CONFIGURATION (Strict Symmetry & Publication Aesthetics)
# ------------------------------------------------------------------------------
cloud_colors <- c(
  "Diagonal"     = "#2196F3",  # Blue
  "Horizontal"   = "#FF9800",  # Orange
  "Vertical"     = "#9C27B0",  # Purple
  "Round"        = "#4CAF50",  # Green
  "Unclassified" = "#757575"   # Grey
)

# Common axis limits to guarantee identical geometric scaling
axis_lim <- c(-5.5, 5.5)

# Shared professional theme
theme_pub <- theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(size = 17, face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 13, color = "black", hjust = 0.5, margin = margin(b = 10)),
    axis.title    = element_text(size = 14, face = "bold", color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.ticks    = element_line(color = "black", linewidth = 0.7),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.tag      = element_text(size = 20, face = "bold"),
    legend.position = "bottom",
    legend.title  = element_blank(),
    legend.text   = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(t = -5, b = 5)
  )

# Data frame for clean DiPa quadrant labels with white background pill
quadrant_labels <- data.frame(
  x = c(0, 3.8, -3.8, 0, 0, 3.8, -3.8, -3.8, 3.8),
  y = c(0, 3.8, -3.8, 3.8, -3.8, 0, 0, 3.8, -3.8),
  label = c("8", "1", "2", "3", "4", "5", "6", "7", "7")
)

# ------------------------------------------------------------------------------
# 5. PANEL A: OBSERVED DIPA PLOT
# ------------------------------------------------------------------------------
sub_title_a <- bquote(
  italic(N) == .(n_genes) ~ " | " ~
  italic(s)[RNA] == .(sprintf("%.2f", sd.1)) * "," ~
  italic(s)[Prot] == .(sprintf("%.2f", sd.2)) ~ " | " ~
  bold(rho[BP] == .(sprintf("%.2f", obs_corr)))
)

p_dipa <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.65, size = 1.9) +
  scale_color_manual(values = cloud_colors) +
  # Reference zero-lines (mu = 0)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  # Standard DiPa threshold lines (+/- 0.5)
  geom_hline(yintercept = c(-0.5, 0.5), color = "grey40", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), color = "grey40", linewidth = 0.5, linetype = "dashed") +
  # Sector annotations (1 to 8) with clean white background fill so lines do not clash
  geom_label(
    data = quadrant_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color = "grey30",
    fill = "white",
    linewidth = 0,
    label.padding = unit(0.15, "lines"),
    fontface = "bold",
    size = 5.5
  ) +
  coord_fixed(xlim = axis_lim, ylim = axis_lim) +
  labs(
    tag = "A",
    title = "Observed data (DiPa)",
    subtitle = sub_title_a,
    x = expression(bold(RNA*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)")),
    y = expression(bold(Protein*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"))
  ) +
  theme_pub

# ------------------------------------------------------------------------------
# 6. PANEL B: SYMMETRIC BIVARIATE GAUSSIAN NULL MODEL
# ------------------------------------------------------------------------------
# Theoretical concentric contour circles for bivariate symmetric normal:
radii <- c(1, 2, 3) * sd.pooled
theta <- seq(0, 2*pi, length.out = 200)
circle_df <- rbindlist(lapply(seq_along(radii), function(k) {
  data.table(
    Level = paste0(k, "*sigma"),
    x = radii[k] * cos(theta),
    y = radii[k] * sin(theta)
  )
}))

sub_title_b <- bquote(
  italic(N) == .(n_genes) ~ " | " ~
  bold(mu) == "(0, 0)" ~ " | " ~
  italic(s)[pooled] == .(sprintf("%.2f", sd.pooled)) ~ " | " ~
  bold(rho[BP] == 0)
)

p_gauss <- ggplot(sim_DT, aes(x = log2G, y = log2P)) +
  geom_point(alpha = 0.55, size = 1.9, color = "#37474F") +
  # Theoretical contour circles of the symmetric bivariate normal
  geom_path(data = circle_df, aes(x = x, y = y, group = Level),
            color = "#1976D2", linetype = "dotted", linewidth = 0.7, inherit.aes = FALSE) +
  # Reference zero-lines (mu = 0)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7, linetype = "solid") +
  # Threshold lines (+/- 0.5) matching Panel A
  geom_hline(yintercept = c(-0.5, 0.5), color = "grey40", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), color = "grey40", linewidth = 0.5, linetype = "dashed") +
  coord_fixed(xlim = axis_lim, ylim = axis_lim) +
  labs(
    tag = "B",
    title = "Symmetric bivariate gaussian model",
    subtitle = sub_title_b,
    x = expression(bold(Simulated~RNA~~log[2]~"(fold change)")),
    y = expression(bold(Simulated~Protein~~log[2]~"(fold change)"))
  ) +
  theme_pub

# ------------------------------------------------------------------------------
# 7. ASSEMBLE COMBINED FIGURE WITH PATCHWORK
# ------------------------------------------------------------------------------
layout_combined <- (p_dipa | p_gauss) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# 8. EXPORT FIGURE (PDF ONLY)
# ------------------------------------------------------------------------------
out_pdf_share <- file.path(data_dir, "New_Data/DiPa_vs_Symmetric_Gaussian_BDL.pdf")
out_pdf_fig03 <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Figure_03/DiPa_vs_Symmetric_Gaussian_BDL.pdf"

message("Saving vector PDF outputs...")
ggsave(out_pdf_share, plot = layout_combined, width = 14, height = 7.5, device = cairo_pdf)
ggsave(out_pdf_fig03, plot = layout_combined, width = 14, height = 7.5, device = cairo_pdf)

message("SUCCESS! Only PDF outputs were generated:")
message(paste0("  [1] ", out_pdf_share))
message(paste0("  [2] ", out_pdf_fig03))
