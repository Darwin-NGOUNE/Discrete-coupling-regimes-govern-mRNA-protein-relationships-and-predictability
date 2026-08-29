# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_1_Master_Plate.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Supplemental Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1/
# PURPOSE: 1:1 exact reproduction of PowerPoint Slide 1 (Supplementary Figure 1).
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Exact PowerPoint Supplementary Figure 1 in R ===\n")

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

cat("Assembling Exact 16:9 Canvas matching PowerPoint Slide 1...\n")

# Exact PowerPoint Slide 1 Canvas (16:9 widescreen ratio)
master_canvas <- ggdraw() +
  # ---------------------------------------------------------------------------
  # LEFT COLUMN: PANEL A (PCA PLOTS FOR BDL & CCl4)
  # ---------------------------------------------------------------------------
  draw_label("A", x = 0.02, y = 0.94, size = 24, fontface = "bold") +
  draw_label("BDL", x = 0.23, y = 0.94, size = 18, fontface = "bold") +
  draw_grob(grob_pca_bdl,  x = 0.03, y = 0.50, width = 0.42, height = 0.42) +
  
  draw_label(expression(bold(CCl[4])), x = 0.23, y = 0.46, size = 18) +
  draw_grob(grob_pca_ccl4, x = 0.03, y = 0.03, width = 0.42, height = 0.42) +

  # ---------------------------------------------------------------------------
  # RIGHT COLUMN: PANEL B (HEATMAPS FOR BDL & CCl4)
  # ---------------------------------------------------------------------------
  draw_label("B", x = 0.48, y = 0.94, size = 24, fontface = "bold") +
  
  # Top Row: BDL Heatmaps (RNA left, Protein right)
  draw_grob(grob_hm_rna_bdl,   x = 0.49, y = 0.50, width = 0.24, height = 0.42) +
  draw_grob(grob_hm_prot_bdl,  x = 0.74, y = 0.50, width = 0.24, height = 0.42) +
  
  # Bottom Row: CCl4 Heatmaps (RNA left, Protein right)
  draw_grob(grob_hm_rna_ccl4,  x = 0.49, y = 0.03, width = 0.24, height = 0.42) +
  draw_grob(grob_hm_prot_ccl4, x = 0.74, y = 0.03, width = 0.24, height = 0.42)

cat(sprintf("Saving Exact Supplementary Figure 1 PDF to: %s\n", out_pdf))
ggsave(out_pdf, plot = master_canvas, width = 20, height = 11.25, units = "in", dpi = 300)

cat("\nSUCCESS! Supplementary Figure 1 exactly matching PowerPoint Slide 1 generated!\n")
