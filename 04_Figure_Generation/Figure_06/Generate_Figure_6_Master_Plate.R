# ==============================================================================
# SCRIPT: Generate_Figure_6_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# PURPOSE: Assemble Figure 6 (Procedure 3 4-Models, Proc3 Scatters, Proc1 4-Models):
#          - Tier 1: Panel A (Procedure 3 4-Models Pearson Boxplot, BDL + CCl4)
#          - Tier 2: Panel B (Proc3 1x4 Scatterplots Page 3: Protein model: BDL + CCl4)
#          - Tier 3: Panel C (Proc3 Full 1x4 Scatterplots Page 3: Protein model: BDL + CCl4 (all animals))
#          - Tier 4: Panel D (Procedure 1 4-Models Pearson Boxplot, Train BDL / Test CCl4)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Figure 6 Master Plate in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
boxplot_dir  <- file.path(new_data_dir, "Cluster_Modelierung_Hengstler/New_Boxplot")
scatter_dir  <- file.path(boxplot_dir, "Scatterplots_75th_Percentile")
out_pdf      <- file.path(grafiken_dir, "Figure_6.pdf")

# Helper function to convert a PDF page to a raster grob at high resolution
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

cat("Locating and rendering Figure 6 components at 300 DPI...\n")

f_panel_a <- file.path(boxplot_dir, "Procedure_3_4Models_Pearson_Merged_Batch.pdf")
f_panel_b <- file.path(scatter_dir, "Proc3_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf")
f_panel_c <- file.path(scatter_dir, "Proc3_Full_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf")
f_panel_d <- file.path(boxplot_dir, "Procedure_1_4Models_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf")

# Copy components to Grafiken_Paper_1
file.copy(f_panel_a, file.path(grafiken_dir, basename(f_panel_a)), overwrite = TRUE)
file.copy(f_panel_b, file.path(grafiken_dir, basename(f_panel_b)), overwrite = TRUE)
file.copy(f_panel_c, file.path(grafiken_dir, basename(f_panel_c)), overwrite = TRUE)
file.copy(f_panel_d, file.path(grafiken_dir, basename(f_panel_d)), overwrite = TRUE)

grob_a <- pdf_to_grob(f_panel_a, page_num = 1)
grob_b <- pdf_to_grob(f_panel_b, page_num = 3)
grob_c <- pdf_to_grob(f_panel_c, page_num = 3)
grob_d <- pdf_to_grob(f_panel_d, page_num = 1)

cat("Building 4-tier master layout canvas (Panels A, B, C, D)...\n")

# Canvas: Width = 18 inches, Height = 24 inches
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TIER 1: Panel A (Procedure 3 4-Models Boxplot BDL + CCl4)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.985, size = 28, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.755, width = 0.96, height = 0.225) +
  
  # ---------------------------------------------------------------------------
  # TIER 2: Panel B (Proc3 1x4 Scatterplots Page 3: Protein model: BDL + CCl4)
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.02, y = 0.735, size = 28, fontface = "bold") +
  draw_grob(grob_b, x = 0.02, y = 0.505, width = 0.96, height = 0.225) +

  # ---------------------------------------------------------------------------
  # TIER 3: Panel C (Proc3 Full 1x4 Scatterplots Page 3: Protein model: BDL + CCl4 (all animals))
  # ---------------------------------------------------------------------------
  draw_label("C", x = 0.02, y = 0.485, size = 28, fontface = "bold") +
  draw_grob(grob_c, x = 0.02, y = 0.255, width = 0.96, height = 0.225) +

  # ---------------------------------------------------------------------------
  # TIER 4: Panel D (Procedure 1 4-Models Boxplot Train BDL / Test CCl4)
  # ---------------------------------------------------------------------------
  draw_label("D", x = 0.02, y = 0.235, size = 28, fontface = "bold") +
  draw_grob(grob_d, x = 0.02, y = 0.005, width = 0.96, height = 0.225)

cat(sprintf("Saving Master Figure 6 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 24, units = "in", dpi = 300)

cat("\nSUCCESS! Master Figure 6 PDF generated and saved to Grafiken_Paper_1!\n")
