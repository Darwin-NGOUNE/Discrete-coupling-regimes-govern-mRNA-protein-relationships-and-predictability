# ==============================================================================
# SCRIPT: Generate_Figure_2_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: Assemble all 10 individual PDF graphics into a single, publication-ready
#          Master Figure 2 PDF plate matching the PowerPoint structure 100%,
#          with plot-width underlines under titles, and RNA / Protein lines.
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
source_dir   <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
out_pdf      <- file.path(grafiken_dir, "Figure_2.pdf")

# Helper function to convert a PDF page to a raster grob at high resolution (300 DPI)
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  img <- image_trim(img)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

cat("Rendering PDF components at 300 DPI from updated source directory...\n")

# PANEL A: 4 Plots (PCA BDL, PCA CCl4, Volcano BDL, Volcano CCl4)
grob_pca_bdl     <- pdf_to_grob(file.path(source_dir, "Isolated_PCA_BDL.pdf"))
grob_pca_ccl4    <- pdf_to_grob(file.path(source_dir, "Isolated_PCA_CCL4.pdf"))
grob_volc_bdl    <- pdf_to_grob(file.path(source_dir, "Isolated_Volcano_BDL.pdf"))
grob_volc_ccl4   <- pdf_to_grob(file.path(source_dir, "Isolated_Volcano_CCL4.pdf"))

# PANEL B: 4 Heatmaps (Global RNA & Protein Landscapes - Pages 2 & 3)
grob_rna_bdl     <- pdf_to_grob(file.path(source_dir, "Global_RNA_Landscape_Summary.pdf"), page_num = 2)
grob_rna_ccl4    <- pdf_to_grob(file.path(source_dir, "Global_RNA_Landscape_Summary.pdf"), page_num = 3)
grob_prot_bdl    <- pdf_to_grob(file.path(source_dir, "Global_Protein_Landscape_Summary.pdf"), page_num = 2)
grob_prot_ccl4   <- pdf_to_grob(file.path(source_dir, "Global_Protein_Landscape_Summary.pdf"), page_num = 3)

# PANEL C: Unified Correlation Triptych (RNA + Protein + Single Shared Legend)
grob_panel_c     <- pdf_to_grob(file.path(grafiken_dir, "Panel_C_Combined_Correlation_Heatmaps.pdf"))

# PANEL D: Unified Correlation Density & Centrality plots (BDL + CCl4 + Single Shared Legend)
grob_panel_d     <- pdf_to_grob(file.path(grafiken_dir, "Panel_D_Combined_Density_Plots.pdf"))

cat("Building precise 2x2 multi-panel layout canvas...\n")

# Mathematical expressions for proper bold CCl4 with subscript 4
expr_ccl4_top <- expression(bold(CCl[4]))
expr_ccl4_bot <- expression(bold(CCl[4]))

master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # TOP HEADERS: BDL vs CCl4 (Titles above lines, lines matching graph widths)
  # ---------------------------------------------------------------------------
  # Panel A Top Headers (BDL: x = 0.04 to 0.26; CCl4: x = 0.27 to 0.49)
  draw_label("BDL", x = 0.15, y = 0.982, size = 16, fontface = "bold") +
  draw_line(x = c(0.04, 0.26), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +
  
  draw_label(expr_ccl4_top, x = 0.38, y = 0.982, size = 16) +
  draw_line(x = c(0.27, 0.49), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +
  
  # Panel B Top Headers (BDL: x = 0.54 to 0.76; CCl4: x = 0.77 to 0.99)
  draw_label("BDL", x = 0.65, y = 0.982, size = 16, fontface = "bold") +
  draw_line(x = c(0.54, 0.76), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +
  
  draw_label(expr_ccl4_top, x = 0.88, y = 0.982, size = 16) +
  draw_line(x = c(0.77, 0.99), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +

  # ---------------------------------------------------------------------------
  # PANEL A (Top Left): 4 PCA & Volcano plots
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.015, y = 0.982, size = 22, fontface = "bold") +
  draw_grob(grob_pca_bdl,   x = 0.04, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_pca_ccl4,  x = 0.27, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_volc_bdl,  x = 0.04, y = 0.51, width = 0.22, height = 0.22) +
  draw_grob(grob_volc_ccl4, x = 0.27, y = 0.51, width = 0.22, height = 0.22) +

  # ---------------------------------------------------------------------------
  # PANEL B (Top Right): 4 RNA & Protein Heatmaps
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.51, y = 0.982, size = 22, fontface = "bold") +
  
  # RNA row (y = 0.74 to 0.96)
  draw_label("RNA", x = 0.522, y = 0.85, size = 15, fontface = "bold", angle = 90) +
  draw_line(x = c(0.533, 0.533), y = c(0.74, 0.96), color = "black", linewidth = 1.2) +
  draw_grob(grob_rna_bdl,   x = 0.54, y = 0.74, width = 0.22, height = 0.22) +
  draw_grob(grob_rna_ccl4,  x = 0.77, y = 0.74, width = 0.22, height = 0.22) +
  
  # Protein row (y = 0.51 to 0.73)
  draw_label("Protein", x = 0.522, y = 0.62, size = 15, fontface = "bold", angle = 90) +
  draw_line(x = c(0.533, 0.533), y = c(0.51, 0.73), color = "black", linewidth = 1.2) +
  draw_grob(grob_prot_bdl,  x = 0.54, y = 0.51, width = 0.22, height = 0.22) +
  draw_grob(grob_prot_ccl4, x = 0.77, y = 0.51, width = 0.22, height = 0.22) +

  # ---------------------------------------------------------------------------
  # PANEL C (Bottom Left): Unified 2x3 Correlation Matrix (RNA + Protein)
  # ---------------------------------------------------------------------------
  draw_label("C", x = 0.015, y = 0.47, size = 22, fontface = "bold") +
  
  # RNA row (y = 0.25 to 0.46)
  draw_label("RNA", x = 0.022, y = 0.355, size = 15, fontface = "bold", angle = 90) +
  draw_line(x = c(0.033, 0.033), y = c(0.25, 0.46), color = "black", linewidth = 1.2) +
  
  # Protein row (y = 0.04 to 0.24)
  draw_label("Protein", x = 0.022, y = 0.145, size = 15, fontface = "bold", angle = 90) +
  draw_line(x = c(0.033, 0.033), y = c(0.04, 0.24), color = "black", linewidth = 1.2) +
  
  draw_grob(grob_panel_c, x = 0.04, y = 0.01, width = 0.45, height = 0.46) +

  # ---------------------------------------------------------------------------
  # PANEL D (Bottom Right): Unified Correlation Density & Centrality plots
  # ---------------------------------------------------------------------------
  draw_label("D", x = 0.51, y = 0.47, size = 22, fontface = "bold") +
  
  # BDL plot (x = 0.55 to 0.76)
  draw_label("BDL", x = 0.655, y = 0.455, size = 16, fontface = "bold") +
  draw_line(x = c(0.555, 0.755), y = c(0.440, 0.440), color = "black", linewidth = 1.2) +
  
  # CCl4 plot (x = 0.77 to 0.98)
  draw_label(expr_ccl4_bot, x = 0.875, y = 0.455, size = 16) +
  draw_line(x = c(0.775, 0.975), y = c(0.440, 0.440), color = "black", linewidth = 1.2) +
  
  draw_grob(grob_panel_d, x = 0.54, y = 0.01, width = 0.45, height = 0.43)

cat(sprintf("Saving Master Figure 2 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 18, height = 12, units = "in", dpi = 300)

cat("\nSUCCESS! Master Figure 2 PDF assembled in R and saved to Grafiken_Paper_1!\n")
