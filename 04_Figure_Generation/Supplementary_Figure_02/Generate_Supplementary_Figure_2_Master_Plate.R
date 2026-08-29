# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_2_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Supplemental Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: 1:1 exact reproduction of PowerPoint Slide 2 (Supplementary Figure 2):
#          - Panel A (Left): Protein_Correlation_Density_BDL.pdf (Page 4 - Master Proteins)
#          - Panel B (Right): Protein_Correlation_Density_CCL4.pdf (Page 4 - Master Proteins)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Supplementary Figure 2 Master Plate in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
out_pdf      <- file.path(grafiken_dir, "Supplementary_Figure_2.pdf")

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

cat("Rendering PDF components (Page 4) at 300 DPI...\n")

f_bdl  <- file.path(new_data_dir, "Protein_Correlation_Density_BDL.pdf")
f_ccl4 <- file.path(new_data_dir, "Protein_Correlation_Density_CCL4.pdf")

file.copy(f_bdl,  file.path(grafiken_dir, basename(f_bdl)),  overwrite = TRUE)
file.copy(f_ccl4, file.path(grafiken_dir, basename(f_ccl4)), overwrite = TRUE)

grob_a <- pdf_to_grob(f_bdl,  page_num = 4)
grob_b <- pdf_to_grob(f_ccl4, page_num = 4)

cat("Assembling Exact 16:9 Canvas matching PowerPoint Slide 2...\n")

# Canvas: Width = 18 inches, Height = 9 inches
master_canvas <- ggdraw() +
  # Panel A: BDL Master Proteins Identification (Left)
  draw_label("A", x = 0.02, y = 0.96, size = 22, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.02, width = 0.47, height = 0.92) +
  
  # Panel B: CCl4 Master Proteins Identification (Right)
  draw_label("B", x = 0.51, y = 0.96, size = 22, fontface = "bold") +
  draw_grob(grob_b, x = 0.51, y = 0.02, width = 0.47, height = 0.92)

cat(sprintf("Saving Supplementary Figure 2 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 9, units = "in", dpi = 300)

cat("\nSUCCESS! Supplementary Figure 2 generated and saved to Grafiken_Paper_1!\n")
