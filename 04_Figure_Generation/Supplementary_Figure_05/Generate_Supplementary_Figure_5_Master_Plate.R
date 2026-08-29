# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_5_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Supplemental Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Assemble Supplementary Figure 5 (Conserved Pairs Proportions Subset - Vertical Stack):
#          - Panel A (Top): Direction 1 (BDL -> CCl4) Conserved Pairs Subset
#          - Panel B (Bottom): Direction 2 (CCl4 -> BDL) Conserved Pairs Subset
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Supplementary Figure 5 (Vertical Stack) in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
scatter_dir  <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot/Scatterplots_75th_Percentile"
out_pdf      <- file.path(grafiken_dir, "Supplementary_Figure_5.pdf")

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

cat("Rendering PDF components at 300 DPI...\n")

f_d1 <- file.path(scatter_dir, "Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_Subset.pdf")
f_d2 <- file.path(scatter_dir, "Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_Subset.pdf")

file.copy(f_d1, file.path(grafiken_dir, basename(f_d1)), overwrite = TRUE)
file.copy(f_d2, file.path(grafiken_dir, basename(f_d2)), overwrite = TRUE)

grob_a <- pdf_to_grob(f_d1, page_num = 1)
grob_b <- pdf_to_grob(f_d2, page_num = 1)

cat("Building 2-Tier Vertical Stack Canvas (Panel A Top, Panel B Bottom)...\n")

# Canvas: Width = 14 inches, Height = 16 inches (Vertical 2-Tier Stack)
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TIER 1: PANEL A (Direction 1: BDL -> CCl4 - TOP)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.985, size = 24, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.51, width = 0.96, height = 0.47) +

  # ---------------------------------------------------------------------------
  # TIER 2: PANEL B (Direction 2: CCl4 -> BDL - BOTTOM)
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.02, y = 0.485, size = 24, fontface = "bold") +
  draw_grob(grob_b, x = 0.02, y = 0.01, width = 0.96, height = 0.47)

cat(sprintf("Saving Supplementary Figure 5 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 14, height = 16, units = "in", dpi = 300)

cat("\nSUCCESS! Supplementary Figure 5 (Vertical Stack) generated and saved to Grafiken_Paper_1!\n")
