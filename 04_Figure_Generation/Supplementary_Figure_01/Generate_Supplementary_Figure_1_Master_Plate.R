# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_1_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Supplemental Figures)
# PURPOSE: Assembles Supplementary Figure 1 with standardized bold styling,
#          full-width graph-aligned underlines, and closely positioned RNA / Protein headers.
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Standardized Supplementary Figure 1 Master Plate ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
out_pdf      <- file.path(grafiken_dir, "Supplementary_Figure_1.pdf")

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

# Vector PDF paths
f_pca_bdl  <- file.path(new_data_dir, "Isolated_PCA_BDL_18Mice_normal_data.pdf")
f_pca_ccl4 <- file.path(new_data_dir, "Isolated_PCA_CCL4_36Mice.pdf")
f_hm_prot  <- file.path(new_data_dir, "Global_Protein_Landscape_Heatmaps_ALL_Mice.pdf")
f_hm_rna   <- file.path(new_data_dir, "Global_RNA_Landscape_Heatmaps_ALL_Mice.pdf")

file.copy(f_pca_bdl,  file.path(grafiken_dir, basename(f_pca_bdl)),  overwrite = TRUE)
file.copy(f_pca_ccl4, file.path(grafiken_dir, basename(f_pca_ccl4)), overwrite = TRUE)
file.copy(f_hm_prot,  file.path(grafiken_dir, basename(f_hm_prot)),  overwrite = TRUE)
file.copy(f_hm_rna,   file.path(grafiken_dir, basename(f_hm_rna)),   overwrite = TRUE)

grob_pca_bdl  <- pdf_to_grob(f_pca_bdl, page_num = 1)
grob_pca_ccl4 <- pdf_to_grob(f_pca_ccl4, page_num = 1)

grob_hm_rna_bdl   <- pdf_to_grob(f_hm_rna,  page_num = 1)
grob_hm_prot_bdl  <- pdf_to_grob(f_hm_prot, page_num = 1)

grob_hm_rna_ccl4  <- pdf_to_grob(f_hm_rna,  page_num = 2)
grob_hm_prot_ccl4 <- pdf_to_grob(f_hm_prot, page_num = 2)

cat("Assembling Clean Canvas matching Publication Standards...\n")

# Canvas: Width = 21 inches, Height = 12 inches
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # BDL ROW (TOP)
  # ---------------------------------------------------------------------------
  # Panel A Label & BDL Header for PCA (x = 0.04 to 0.47)
  draw_label("A", x = 0.015, y = 0.982, size = 24, fontface = "bold") +
  draw_label("BDL", x = 0.255, y = 0.982, size = 20, fontface = "bold") +
  draw_line(x = c(0.04, 0.47), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +
  draw_grob(grob_pca_bdl, x = 0.04, y = 0.51, width = 0.43, height = 0.43) +
  
  # Panel B Label & BDL Header for Heatmaps (x = 0.495 to 0.985)
  draw_label("B", x = 0.485, y = 0.982, size = 24, fontface = "bold") +
  draw_label("BDL", x = 0.74, y = 0.982, size = 20, fontface = "bold") +
  draw_line(x = c(0.495, 0.985), y = c(0.965, 0.965), color = "black", linewidth = 1.2) +
  
  # Sub-headers placed right above the individual heatmaps
  draw_label("RNA", x = 0.62, y = 0.938, size = 16, fontface = "bold") +
  draw_label("Protein", x = 0.87, y = 0.938, size = 16, fontface = "bold") +
  draw_grob(grob_hm_rna_bdl,  x = 0.495, y = 0.505, width = 0.245, height = 0.43) +
  draw_grob(grob_hm_prot_bdl, x = 0.745, y = 0.505, width = 0.245, height = 0.43) +
  
  # ---------------------------------------------------------------------------
  # CCL4 ROW (BOTTOM)
  # ---------------------------------------------------------------------------
  # CCl4 Header for PCA (x = 0.04 to 0.47)
  draw_label(expression(bold(CCl[4])), x = 0.255, y = 0.485, size = 20) +
  draw_line(x = c(0.04, 0.47), y = c(0.468, 0.468), color = "black", linewidth = 1.2) +
  draw_grob(grob_pca_ccl4, x = 0.04, y = 0.02, width = 0.43, height = 0.43) +
  
  # CCl4 Header for Heatmaps (x = 0.495 to 0.985)
  draw_label(expression(bold(CCl[4])), x = 0.74, y = 0.485, size = 20) +
  draw_line(x = c(0.495, 0.985), y = c(0.468, 0.468), color = "black", linewidth = 1.2) +
  
  # Sub-headers placed right above the individual heatmaps
  draw_label("RNA", x = 0.62, y = 0.442, size = 16, fontface = "bold") +
  draw_label("Protein", x = 0.87, y = 0.442, size = 16, fontface = "bold") +
  draw_grob(grob_hm_rna_ccl4,  x = 0.495, y = 0.015, width = 0.245, height = 0.43) +
  draw_grob(grob_hm_prot_ccl4, x = 0.745, y = 0.015, width = 0.245, height = 0.43)

cat(sprintf("Saving Supplementary Figure 1 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 21, height = 12, units = "in", dpi = 300)

cat("\nSUCCESS! Supplementary Figure 1 generated and saved!\n")
