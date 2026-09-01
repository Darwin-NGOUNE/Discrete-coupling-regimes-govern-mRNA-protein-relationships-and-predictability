# ==============================================================================
# SCRIPT: Generate_Figure_4_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Assemble Figure 4 with:
#          - Top: Panel A (DiPa Overlap Analysis)
#          - Bottom: Panel B (Centroid Slope BDL) & Panel C (Centroid Slope CCl4)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Figure 4 (Top: A | Bottom: B & C) in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
out_pdf      <- file.path(grafiken_dir, "Figure_4.pdf")

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

cat("Rendering Figure 4 PDF components at 300 DPI...\n")

# Copy component PDFs to Grafiken_Paper_1
f1 <- file.path(new_data_dir, "DiPa_Overlap_Analysis_Results.pdf")
f2 <- file.path(new_data_dir, "Centroid_Slope_Frequency_BDL.pdf")
f3 <- file.path(new_data_dir, "Centroid_Slope_Frequency_CCL4.pdf")

file.copy(f1, file.path(grafiken_dir, "DiPa_Overlap_Analysis_Results.pdf"), overwrite = TRUE)
file.copy(f2, file.path(grafiken_dir, "Centroid_Slope_Frequency_BDL.pdf"), overwrite = TRUE)
file.copy(f3, file.path(grafiken_dir, "Centroid_Slope_Frequency_CCL4.pdf"), overwrite = TRUE)

grob_overlap    <- pdf_to_grob(file.path(grafiken_dir, "DiPa_Overlap_Analysis_Results.pdf"))
grob_slope_bdl  <- pdf_to_grob(file.path(grafiken_dir, "Centroid_Slope_Frequency_BDL.pdf"))
grob_slope_ccl4 <- pdf_to_grob(file.path(grafiken_dir, "Centroid_Slope_Frequency_CCL4.pdf"))

cat("Building 2-row layout canvas (Top: A | Bottom: B & C)...\n")

# Layout canvas: Width = 16 inches, Height = 14 inches
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TOP ROW: Panel A (DiPa Overlap Analysis)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.985, size = 28, fontface = "bold") +
  draw_grob(grob_overlap,    x = 0.02, y = 0.49, width = 0.96, height = 0.48) +
  
  # ---------------------------------------------------------------------------
  # BOTTOM ROW: Panels B (BDL) & C (CCl4) Centroid Slopes
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.02, y = 0.475, size = 28, fontface = "bold") +
  draw_grob(grob_slope_bdl,  x = 0.03, y = 0.01, width = 0.46, height = 0.45) +
  
  draw_label("C", x = 0.51, y = 0.475, size = 28, fontface = "bold") +
  draw_grob(grob_slope_ccl4, x = 0.52, y = 0.01, width = 0.46, height = 0.45)

cat(sprintf("Saving Master Figure 4 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 16, height = 14, units = "in", dpi = 300)

cat("\nSUCCESS! Figure 4 Master PDF (Top: A | Bottom: B & C) generated!\n")
