# ==============================================================================
# SCRIPT: Generate_Figure_2_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Assemble all 10 individual PDF graphics into a single, publication-ready
#          Master Figure 2 PDF plate matching the PowerPoint structure 100%,
#          with clean labels (A, B, C, D, BDL, CCl4, RNA, Protein) and proper
#          mathematical subscripts (CCl4).
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Paper 1 Figure 2 Master Plate in R ===\n")

# Define directories
grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
out_pdf      <- file.path(grafiken_dir, "Figure_2.pdf")

# Helper function to convert a PDF page to a raster grob at compact high quality
pdf_to_grob <- function(pdf_path, page_num = 1, density = 100) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  img <- image_trim(img)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

cat("Rendering PDF components at compact 100 DPI (~4 MB net)...\n")

# PANEL A: 4 Plots (PCA BDL, PCA CCl4, Volcano BDL, Volcano CCl4)
grob_pca_bdl     <- pdf_to_grob(file.path(grafiken_dir, "Isolated_PCA_BDL.pdf"))
grob_pca_ccl4    <- pdf_to_grob(file.path(grafiken_dir, "Isolated_PCA_CCL4.pdf"))
grob_volc_bdl    <- pdf_to_grob(file.path(grafiken_dir, "Isolated_Volcano_BDL.pdf"))
grob_volc_ccl4   <- pdf_to_grob(file.path(grafiken_dir, "Isolated_Volcano_CCL4.pdf"))

# PANEL B: 4 Heatmaps (Global RNA & Protein Landscapes - Pages 2 & 3)
grob_rna_bdl     <- pdf_to_grob(file.path(grafiken_dir, "Global_RNA_Landscape_Summary.pdf"), page_num = 2)
grob_rna_ccl4    <- pdf_to_grob(file.path(grafiken_dir, "Global_RNA_Landscape_Summary.pdf"), page_num = 3)
grob_prot_bdl    <- pdf_to_grob(file.path(grafiken_dir, "Global_Protein_Landscape_Summary.pdf"), page_num = 2)
grob_prot_ccl4   <- pdf_to_grob(file.path(grafiken_dir, "Global_Protein_Landscape_Summary.pdf"), page_num = 3)

# PANEL C: 2 Correlation Triptychs (mRNA & Protein Delta R)
grob_trip_rna    <- pdf_to_grob(file.path(grafiken_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_mRNA_DeltaR.pdf"))
grob_trip_prot   <- pdf_to_grob(file.path(grafiken_dir, "Correlation_Heatmaps_Triptych_BDL_CCl4_DeltaR.pdf"))

# PANEL D: 2 Correlation Density Plots (BDL & CCl4 - Page 1)
grob_dens_bdl    <- pdf_to_grob(file.path(grafiken_dir, "Protein_Correlation_Density_BDL.pdf"), page_num = 1)
grob_dens_ccl4   <- pdf_to_grob(file.path(grafiken_dir, "Protein_Correlation_Density_CCL4.pdf"), page_num = 1)

cat("Building precise 2x2 multi-panel layout canvas...\n")

# Mathematical expressions for proper bold CCl4 with subscript 4
expr_ccl4_top <- expression(bold(CCl[4]))
expr_ccl4_bot <- expression(bold(CCl[4]))

master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TOP HEADERS: BDL vs CCl4 (with 4 subscript)
  # ---------------------------------------------------------------------------
  draw_label("BDL", x = 0.15, y = 0.975, size = 16, fontface = "bold") +
  draw_label(expr_ccl4_top, x = 0.38, y = 0.975, size = 16) +
  draw_label("BDL", x = 0.65, y = 0.975, size = 16, fontface = "bold") +
  draw_label(expr_ccl4_top, x = 0.88, y = 0.975, size = 16) +

  # ---------------------------------------------------------------------------
  # PANEL A (Top Left): 4 PCA & Volcano plots
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.975, size = 22, fontface = "bold") +
  draw_grob(grob_pca_bdl,   x = 0.04, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_pca_ccl4,  x = 0.27, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_volc_bdl,  x = 0.04, y = 0.51, width = 0.22, height = 0.22) +
  draw_grob(grob_volc_ccl4, x = 0.27, y = 0.51, width = 0.22, height = 0.22) +

  # ---------------------------------------------------------------------------
  # PANEL B (Top Right): 4 RNA & Protein Heatmaps
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.52, y = 0.975, size = 22, fontface = "bold") +
  draw_label("RNA", x = 0.52, y = 0.85, size = 14, fontface = "bold", angle = 90) +
  draw_grob(grob_rna_bdl,   x = 0.54, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_rna_ccl4,  x = 0.77, y = 0.74, width = 0.22, height = 0.22) +
  
  draw_label("Protein", x = 0.52, y = 0.62, size = 14, fontface = "bold", angle = 90) +
  draw_grob(grob_prot_bdl,  x = 0.54, y = 0.51, width = 0.22, height = 0.22) +
  draw_grob(grob_prot_ccl4, x = 0.77, y = 0.51, width = 0.22, height = 0.22) +

  # ---------------------------------------------------------------------------
  # PANEL C (Bottom Left): 2 Correlation Triptychs
  # ---------------------------------------------------------------------------
  draw_label("C", x = 0.02, y = 0.47, size = 22, fontface = "bold") +
  draw_label("RNA", x = 0.02, y = 0.36, size = 13, fontface = "bold", angle = 90) +
  draw_grob(grob_trip_rna,  x = 0.04, y = 0.25, width = 0.45, height = 0.21) +
  
  draw_label("Protein", x = 0.02, y = 0.13, size = 13, fontface = "bold", angle = 90) +
  draw_grob(grob_trip_prot, x = 0.04, y = 0.02, width = 0.45, height = 0.21) +

  # ---------------------------------------------------------------------------
  # PANEL D (Bottom Right): 2 Correlation Density & Centrality plots
  # ---------------------------------------------------------------------------
  draw_label("D", x = 0.52, y = 0.47, size = 22, fontface = "bold") +
  draw_label("BDL", x = 0.65, y = 0.47, size = 16, fontface = "bold") +
  draw_label(expr_ccl4_bot, x = 0.88, y = 0.47, size = 16) +
  draw_grob(grob_dens_bdl,  x = 0.54, y = 0.03, width = 0.22, height = 0.42) +
  draw_grob(grob_dens_ccl4, x = 0.77, y = 0.03, width = 0.22, height = 0.42)

cat(sprintf("Saving Master Figure 2 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 12, units = "in", dpi = 300)

cat("\nSUCCESS! Master Figure 2 PDF assembled in R and saved to Grafiken_Paper_1!\n")
