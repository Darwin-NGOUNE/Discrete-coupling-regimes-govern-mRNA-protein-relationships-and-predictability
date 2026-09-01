# ==============================================================================
# SCRIPT: Generate_Figure_3_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Exact 1:1 reproduction of PowerPoint Figure 3 containing solely:
#          - Isolated_DiPa_Wolken_BDL.pdf (Left)
#          - Isolated_DiPa_Wolken_CCL4.pdf (Right)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Exact PowerPoint Figure 3 in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
out_pdf      <- file.path(grafiken_dir, "Figure_3.pdf")

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

cat("Rendering DiPa cloud components at 300 DPI...\n")

# Copy clean isolated PDFs to Grafiken_Paper_1
file.copy(file.path(new_data_dir, "Isolated_DiPa_Wolken_BDL.pdf"), file.path(grafiken_dir, "Isolated_DiPa_Wolken_BDL.pdf"), overwrite = TRUE)
file.copy(file.path(new_data_dir, "Isolated_DiPa_Wolken_CCL4.pdf"), file.path(grafiken_dir, "Isolated_DiPa_Wolken_CCL4.pdf"), overwrite = TRUE)

grob_dipa_bdl  <- pdf_to_grob(file.path(grafiken_dir, "Isolated_DiPa_Wolken_BDL.pdf"))
grob_dipa_ccl4 <- pdf_to_grob(file.path(grafiken_dir, "Isolated_DiPa_Wolken_CCL4.pdf"))

# Exact 1:1 PowerPoint slide 3 layout canvas (Widescreen 16:9)
master_canvas <- ggdraw() +
  # Top Headers
  draw_label("BDL", x = 0.25, y = 0.985, size = 24, fontface = "bold") +
  draw_label(expression(bold(CCl[4])), x = 0.75, y = 0.985, size = 24) +
  # Left: Isolated BDL DiPa Wolke + 8 Exemplary Scatterplots
  draw_grob(grob_dipa_bdl,  x = 0.01, y = 0.01, width = 0.485, height = 0.96) +
  # Right: Isolated CCl4 DiPa Wolke + 8 Exemplary Scatterplots
  draw_grob(grob_dipa_ccl4, x = 0.505, y = 0.01, width = 0.485, height = 0.96)

cat(sprintf("Saving Exact Figure 3 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 20, height = 10, units = "in", dpi = 300)

cat("\nSUCCESS! Exact PowerPoint Figure 3 PDF generated!\n")
