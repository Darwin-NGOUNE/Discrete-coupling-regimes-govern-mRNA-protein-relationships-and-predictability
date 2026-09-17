# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_1a_and_1b.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1)
# PURPOSE: Generates two separate, maximized publication figures:
#          - Supplementary_Figure_1a.pdf: Full-cohort PCAs (BDL 18 mice & CCL4 36 mice)
#          - Supplementary_Figure_1b.pdf: Full-cohort 2x2 Heatmaps with generous breathing room
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Standalone Supp Figure 1A & Supp Figure 1B ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
github_dir   <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Supplementary_Figure_01"
new_data_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"

dir.create(grafiken_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(github_dir, recursive = TRUE, showWarnings = FALSE)

# Helper to load vector PDF as high-resolution grob
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  img <- image_trim(img)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

# Vector PDF paths
f_pca_bdl  <- file.path(new_data_dir, "Isolated_PCA_BDL_18Mice_normal_data.pdf")
f_pca_ccl4 <- file.path(new_data_dir, "Isolated_PCA_CCL4_36Mice.pdf")
f_hm_prot  <- file.path(new_data_dir, "Global_Protein_Landscape_Heatmaps_ALL_Mice.pdf")
f_hm_rna   <- file.path(new_data_dir, "Global_RNA_Landscape_Heatmaps_ALL_Mice.pdf")

cat("Rendering vector PDF components at 300 DPI...\n")
grob_pca_bdl      <- pdf_to_grob(f_pca_bdl,  page_num = 1, density = 300)
grob_pca_ccl4     <- pdf_to_grob(f_pca_ccl4, page_num = 1, density = 300)

grob_hm_rna_bdl   <- pdf_to_grob(f_hm_rna,   page_num = 1, density = 300)
grob_hm_prot_bdl  <- pdf_to_grob(f_hm_prot,  page_num = 1, density = 300)
grob_hm_rna_ccl4  <- pdf_to_grob(f_hm_rna,   page_num = 2, density = 300)
grob_hm_prot_ccl4 <- pdf_to_grob(f_hm_prot,  page_num = 2, density = 300)

# ==============================================================================
# 1. SUPPLEMENTARY FIGURE 1A (FULL-COHORT PCAs)
# ==============================================================================
cat("Assembling Supplementary Figure 1A (PCAs)...\n")

# Dimensions: Standard publication A4 / Landscape Canvas (14 x 10 inches)
canvas_1a <- ggdraw() +
  # Top Header (BDL)
  draw_label("BDL", x = 0.50, y = 0.980, size = 22, fontface = "bold") +
  draw_line(x = c(0.04, 0.96), y = c(0.960, 0.960), color = "black", linewidth = 1.3) +
  draw_grob(grob_pca_bdl, x = 0.04, y = 0.505, width = 0.92, height = 0.445) +
  
  # Bottom Header (CCl4)
  draw_label(expression(bold(CCl[4])), x = 0.50, y = 0.485, size = 22) +
  draw_line(x = c(0.04, 0.96), y = c(0.465, 0.465), color = "black", linewidth = 1.3) +
  draw_grob(grob_pca_ccl4, x = 0.04, y = 0.015, width = 0.92, height = 0.445)

out_pdf_1a_share  <- file.path(grafiken_dir, "Supplementary_Figure_1a.pdf")
out_pdf_1a_github <- file.path(github_dir,   "Supplementary_Figure_1a.pdf")

ggsave(out_pdf_1a_share,  plot = canvas_1a, width = 16, height = 12, units = "in", dpi = 300)
ggsave(out_pdf_1a_github, plot = canvas_1a, width = 16, height = 12, units = "in", dpi = 300)
cat("-> Supplementary Figure 1A saved successfully.\n")


# ==============================================================================
# 2. SUPPLEMENTARY FIGURE 1B (FULL-COHORT EXPRESSION HEATMAPS 2x2 WITH BREATHING ROOM)
# ==============================================================================
cat("Assembling Supplementary Figure 1B (2x2 Heatmaps with Generous Whitespace)...\n")

# Dimensions: Large Canvas (16.5 x 12.5 inches)
canvas_1b <- ggdraw() +
  # ---------------------------------------------------------------------------
  # BDL ROW (TOP)
  # ---------------------------------------------------------------------------
  # Main BDL Title
  draw_label("BDL", x = 0.50, y = 0.982, size = 22, fontface = "bold") +
  draw_line(x = c(0.03, 0.97), y = c(0.962, 0.962), color = "black", linewidth = 1.3) +
  
  # RNA Title (Left Column) - Well spaced with clean size
  draw_label("RNA", x = 0.260, y = 0.938, size = 15, fontface = "bold", color = "grey15") +
  
  # Protein Title (Right Column) - Distinct and centered over right heatmap
  draw_label("Protein", x = 0.740, y = 0.938, size = 15, fontface = "bold", color = "grey15") +
  
  # Heatmap grobs with clear top margin below titles (y = 0.510, height = 0.405, top at y = 0.915)
  draw_grob(grob_hm_rna_bdl,  x = 0.030, y = 0.510, width = 0.460, height = 0.405) +
  draw_grob(grob_hm_prot_bdl, x = 0.510, y = 0.510, width = 0.460, height = 0.405) +
  
  # ---------------------------------------------------------------------------
  # CCL4 ROW (BOTTOM)
  # ---------------------------------------------------------------------------
  # Main CCl4 Title
  draw_label(expression(bold(CCl[4])), x = 0.50, y = 0.478, size = 22) +
  draw_line(x = c(0.03, 0.97), y = c(0.458, 0.458), color = "black", linewidth = 1.3) +
  
  # RNA Title (Left Column)
  draw_label("RNA", x = 0.260, y = 0.434, size = 15, fontface = "bold", color = "grey15") +
  
  # Protein Title (Right Column)
  draw_label("Protein", x = 0.740, y = 0.434, size = 15, fontface = "bold", color = "grey15") +
  
  # Heatmap grobs with clear top margin below titles (y = 0.015, height = 0.400, top at y = 0.415)
  draw_grob(grob_hm_rna_ccl4,  x = 0.030, y = 0.015, width = 0.460, height = 0.400) +
  draw_grob(grob_hm_prot_ccl4, x = 0.510, y = 0.015, width = 0.460, height = 0.400)

out_pdf_1b_share  <- file.path(grafiken_dir, "Supplementary_Figure_1b.pdf")
out_pdf_1b_github <- file.path(github_dir,   "Supplementary_Figure_1b.pdf")

ggsave(out_pdf_1b_share,  plot = canvas_1b, width = 18, height = 14, units = "in", dpi = 300)
ggsave(out_pdf_1b_github, plot = canvas_1b, width = 18, height = 14, units = "in", dpi = 300)
cat("-> Supplementary Figure 1B saved successfully.\n")

cat("\n==============================================================================\n")
cat("SUCCESS: Supplementary Figure 1B Updated with Clean Spacing & Margins!\n")
cat("==============================================================================\n")
