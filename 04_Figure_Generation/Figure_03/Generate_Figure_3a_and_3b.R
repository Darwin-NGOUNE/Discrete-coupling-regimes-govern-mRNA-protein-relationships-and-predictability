# ==============================================================================
# SCRIPT: Generate_Figure_3a_and_3b.R
# PROJECT: Liver Fibrosis Protein Modeling Consortium (Paper 1 Figures)
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Figure_03/
# PURPOSE: Generate standalone, full-visibility vector graphics for Figure 3:
#          - Figure_3a.pdf : BDL 2D DiPa Coordinate Cloud + 8 Archetypal Pairs (Panels A & B)
#          - Figure_3b.pdf : CCl4 2D DiPa Coordinate Cloud + 8 Archetypal Pairs (Panels A & B / C & D)
# ==============================================================================

cat("=== Generating Figure 3A and Figure 3B Standalone Vector PDFs ===\n")

source_dir   <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
grafiken_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"
fig03_dir    <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Figure_03"

bdl_pdf_src  <- file.path(source_dir, "Isolated_DiPa_Wolken_BDL.pdf")
ccl4_pdf_src <- file.path(source_dir, "Isolated_DiPa_Wolken_CCL4.pdf")

# Copy to Figure_03 and Grafiken_Paper_1
if (file.exists(bdl_pdf_src)) {
  file.copy(bdl_pdf_src, file.path(fig03_dir, "Figure_3a.pdf"), overwrite = TRUE)
  file.copy(bdl_pdf_src, file.path(grafiken_dir, "Figure_3a.pdf"), overwrite = TRUE)
  cat("SUCCESS: Figure_3a.pdf generated and saved!\n")
} else {
  cat("WARNING: Isolated_DiPa_Wolken_BDL.pdf not found in source_dir\n")
}

if (file.exists(ccl4_pdf_src)) {
  file.copy(ccl4_pdf_src, file.path(fig03_dir, "Figure_3b.pdf"), overwrite = TRUE)
  file.copy(ccl4_pdf_src, file.path(grafiken_dir, "Figure_3b.pdf"), overwrite = TRUE)
  cat("SUCCESS: Figure_3b.pdf generated and saved!\n")
} else {
  cat("WARNING: Isolated_DiPa_Wolken_CCL4.pdf not found in source_dir\n")
}

cat("=== Figure 3 Partition Completed Successfully ===\n")
