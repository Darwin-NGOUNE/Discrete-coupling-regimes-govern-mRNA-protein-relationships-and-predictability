# ==============================================================================
# SCRIPT: Generate_Figure_6ABC_and_6D.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# PURPOSE: Split Figure 6 into two high-readability master figures:
#          1. Figure 6_ABC (Panels A, B, C: Procedure 3 Merged Multi-Omics Models & Scatters)
#          2. Figure 6_D   (Panel D Standalone: Procedure 1 Cross-Cohort Validation Models)
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

cat("=== Generating Figure 6_ABC and Figure 6_D Master Plates in R ===\n")

grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
fig6_dir     <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Figure_06"
boxplot_dir  <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung_Hengstler/New_Boxplot"
scatter_dir  <- file.path(boxplot_dir, "Scatterplots_75th_Percentile")

if (!dir.exists(grafiken_dir)) dir.create(grafiken_dir, recursive = TRUE)

# Helper function to convert a PDF page to a raster grob at high resolution
pdf_to_grob <- function(pdf_path, page_num = 1, density = 300) {
  if (!file.exists(pdf_path)) {
    stop(sprintf("File not found: %s", pdf_path))
  }
  img <- image_read_pdf(pdf_path, pages = page_num, density = density)
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  return(grob)
}

f_panel_a <- file.path(boxplot_dir, "Procedure_3_4Models_Pearson_Merged_Batch.pdf")
f_panel_b <- file.path(scatter_dir, "Proc3_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf")
f_panel_c <- file.path(scatter_dir, "Proc3_Full_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf")
f_panel_d <- file.path(boxplot_dir, "Procedure_1_4Models_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf")

# Verify component files exist
stopifnot(file.exists(f_panel_a), file.exists(f_panel_b), file.exists(f_panel_c), file.exists(f_panel_d))

cat("Rendering component grobs at 300 DPI...\n")
grob_a <- pdf_to_grob(f_panel_a, page_num = 1)
grob_b <- pdf_to_grob(f_panel_b, page_num = 3)
grob_c <- pdf_to_grob(f_panel_c, page_num = 3)
grob_d <- pdf_to_grob(f_panel_d, page_num = 1)

# ==============================================================================
# 1. ASSEMBLE FIGURE 6_ABC (3-TIER MASTER PLATE)
# ==============================================================================
cat("\nBuilding Figure 6_ABC (Panels A, B, C)...\n")

# Canvas: Width = 18 inches, Height = 18.5 inches (Generous vertical space for each tier)
canvas_abc <- ggdraw() +
  # TIER 1: Panel A (Procedure 3 4-Models Pearson Boxplot)
  draw_label("A", x = 0.02, y = 0.985, size = 26, fontface = "bold") +
  draw_grob(grob_a, x = 0.02, y = 0.675, width = 0.96, height = 0.305) +
  
  # TIER 2: Panel B (Proc3 1x4 Scatterplots: Conserved Subset)
  draw_label("B", x = 0.02, y = 0.655, size = 26, fontface = "bold") +
  draw_grob(grob_b, x = 0.02, y = 0.345, width = 0.96, height = 0.305) +

  # TIER 3: Panel C (Proc3 Full 1x4 Scatterplots: All Animals Subset)
  draw_label("C", x = 0.02, y = 0.325, size = 26, fontface = "bold") +
  draw_grob(grob_c, x = 0.02, y = 0.015, width = 0.96, height = 0.305)

out_abc_1 <- file.path(fig6_dir, "Figure_6_ABC.pdf")
out_abc_2 <- file.path(grafiken_dir, "Figure_6_ABC.pdf")

ggsave(out_abc_1, plot = canvas_abc, width = 18, height = 22, units = "in", dpi = 300)
ggsave(out_abc_2, plot = canvas_abc, width = 18, height = 22, units = "in", dpi = 300)
cat(sprintf("  -> Saved Figure 6_ABC to: %s\n", out_abc_1))

# ==============================================================================
# 2. ASSEMBLE FIGURE 6_D (STANDALONE HIGH-VISIBILITY FIGURE)
# ==============================================================================
cat("\nBuilding Figure 6_D (Panel D Standalone)...\n")

# Canvas: Width = 18 inches, Height = 8.5 inches (Full-width, large readable text)
canvas_d <- ggdraw() +
  draw_label("A", x = 0.02, y = 0.96, size = 26, fontface = "bold") +
  draw_grob(grob_d, x = 0.02, y = 0.02, width = 0.96, height = 0.94)

out_d_1 <- file.path(fig6_dir, "Figure_6_D.pdf")
out_d_2 <- file.path(grafiken_dir, "Figure_6_D.pdf")

ggsave(out_d_1, plot = canvas_d, width = 18, height = 8.5, units = "in", dpi = 300)
ggsave(out_d_2, plot = canvas_d, width = 18, height = 8.5, units = "in", dpi = 300)
cat(sprintf("  -> Saved Figure 6_D to: %s\n", out_d_1))

cat("\n==============================================================================\n")
cat("SUCCESS: Both Figure 6_ABC and Figure 6_D generated with optimal readability!\n")
cat("==============================================================================\n")
