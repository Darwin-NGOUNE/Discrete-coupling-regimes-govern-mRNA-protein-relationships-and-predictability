# ==============================================================================
# SCRIPT: Generate_Figure_5_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Assemble Figure 5 (Procedure 1 Cross-Cohort Validation & 75th Percentile Scatters):
#          - Line 1: Panel A (Train BDL -> Test CCl4 Baseline) & Panel B (Train CCl4 -> Test BDL Baseline)
#          - Line 2: Panel C (Proc1 1x4 Scatterplots 75th Percentile Direction 1)
#          - Line 3: Panel D (Proc1 Full 1x4 Scatterplots 75th Percentile BDL -> CCl4)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Figure 5 Master Plate in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
boxplot_dir  <- file.path(new_data_dir, "Cluster_Modelierung_Hengstler/New_Boxplot")
scatter_dir  <- file.path(boxplot_dir, "Scatterplots_75th_Percentile")
out_pdf      <- file.path(grafiken_dir, "Figure_5.pdf")

# Helper function to convert a PDF page to a raster grob at high resolution
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  img <- image_trim(img)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

cat("Locating and rendering Figure 5 components at 300 DPI...\n")

f_panel_a <- file.path(boxplot_dir, "Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf")
f_panel_b <- file.path(boxplot_dir, "Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf")
f_panel_c <- file.path(scatter_dir, "Proc1_3Pages_1x4_Scatterplot_75thPercentile_Richtung1_Train_BDL_Test_CCl4.pdf")
f_panel_d <- file.path(scatter_dir, "Proc1_Full_3Pages_1x4_Scatterplot_75thPercentile_BDL_CCL4.pdf")

# Copy components to Grafiken_Paper_1
file.copy(f_panel_a, file.path(grafiken_dir, basename(f_panel_a)), overwrite = TRUE)
file.copy(f_panel_b, file.path(grafiken_dir, basename(f_panel_b)), overwrite = TRUE)
file.copy(f_panel_c, file.path(grafiken_dir, basename(f_panel_c)), overwrite = TRUE)
file.copy(f_panel_d, file.path(grafiken_dir, basename(f_panel_d)), overwrite = TRUE)

grob_a <- pdf_to_grob(f_panel_a, page_num = 1)
grob_b <- pdf_to_grob(f_panel_b, page_num = 1)
grob_c <- pdf_to_grob(f_panel_c, page_num = 1)
grob_d <- pdf_to_grob(f_panel_d, page_num = 1)

cat("Building 3-tier master layout canvas (Panels A, B, C, D)...\n")

# Canvas: Width = 18 inches, Height = 18 inches
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TIER 1: Panel A (Top Left) & Panel B (Top Right)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.985, size = 22, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.68, width = 0.47, height = 0.30) +
  
  draw_label("B", x = 0.51, y = 0.985, size = 22, fontface = "bold") +
  draw_grob(grob_b, x = 0.51, y = 0.68, width = 0.47, height = 0.30) +

  # ---------------------------------------------------------------------------
  # TIER 2: Panel C (Middle Full Width)
  # ---------------------------------------------------------------------------
  draw_label("C", x = 0.02, y = 0.66, size = 22, fontface = "bold") +
  draw_grob(grob_c, x = 0.02, y = 0.35, width = 0.96, height = 0.30) +

  # ---------------------------------------------------------------------------
  # TIER 3: Panel D (Bottom Full Width)
  # ---------------------------------------------------------------------------
  draw_label("D", x = 0.02, y = 0.33, size = 22, fontface = "bold") +
  draw_grob(grob_d, x = 0.02, y = 0.02, width = 0.96, height = 0.30)

cat(sprintf("Saving Master Figure 5 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 18, units = "in", dpi = 300)

cat("\nSUCCESS! Master Figure 5 PDF generated and saved to Grafiken_Paper_1!\n")
