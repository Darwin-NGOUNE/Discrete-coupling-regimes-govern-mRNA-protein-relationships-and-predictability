# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_3_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Supplemental Figures)
# PURPOSE: Assemble Supplementary Figure 3 in a 2-Tier Vertical Stack:
#          - Tier 1: Panel A (Procedure 5 Intra BDL 4-Models Boxplot - Top)
#          - Tier 2: Panel B (Procedure 5 Intra CCl4 4-Models Boxplot - Bottom)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Supplementary Figure 3 (2-Tier Vertical Stack: A Top, B Bottom) in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
boxplot_dir  <- file.path(new_data_dir, "Cluster_Modelierung_Hengstler/New_Boxplot")
out_pdf      <- file.path(grafiken_dir, "Supplementary_Figure_3.pdf")

# Helper function to convert a PDF page to a raster grob at high resolution
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

cat("Rendering PDF components at 300 DPI...\n")

f_p5_a <- file.path(boxplot_dir, "Procedure_5_4Models_Pearson_Intra_BDL.pdf")
f_p5_b <- file.path(boxplot_dir, "Procedure_5_4Models_Pearson_Intra_CCl4.pdf")

file.copy(f_p5_a, file.path(grafiken_dir, basename(f_p5_a)), overwrite = TRUE)
file.copy(f_p5_b, file.path(grafiken_dir, basename(f_p5_b)), overwrite = TRUE)

grob_a <- pdf_to_grob(f_p5_a, page_num = 1)
grob_b <- pdf_to_grob(f_p5_b, page_num = 1)

cat("Building 2-Tier Vertical Canvas (A on Top, B at Bottom)...\n")

# Canvas: Width = 18 inches, Height = 18 inches (2-Tier Vertical Stack)
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TIER 1: PANEL A (Procedure 5 Intra BDL 4-Models - TOP)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.985, size = 28, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.505, width = 0.96, height = 0.47) +

  # ---------------------------------------------------------------------------
  # TIER 2: PANEL B (Procedure 5 Intra CCl4 4-Models - BOTTOM)
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.02, y = 0.485, size = 28, fontface = "bold") +
  draw_grob(grob_b, x = 0.02, y = 0.005, width = 0.96, height = 0.47)

cat(sprintf("Saving Supplementary Figure 3 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 18, units = "in", dpi = 300)

cat("\nSUCCESS! Supplementary Figure 3 (2-Tier Vertical Stack) generated!\n")
