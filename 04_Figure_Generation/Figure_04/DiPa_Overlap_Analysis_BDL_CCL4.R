# ==============================================================================
# SCRIPT: DiPa_Overlap_Analysis_BDL_CCL4.R
# OBJECTIVE: Calculate percentage overlap and generate Venn diagram for
#            DiPa groups between BDL and CCL4 datasets.
# ==============================================================================

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading BDL Dipa data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData")
# DT_dipa_count_bdl

print("Loading CCL4 Dipa data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_CCL4_filtered_Gene_Protein_full_final_with_dipa.RData")
# DT_dipa_count_ccl4

# 2. PREPARE DATA
# -------------------------------------------------------------------------
setDT(DT_dipa_count_bdl)
setDT(DT_dipa_count_ccl4)

BDL_short <- DT_dipa_count_bdl[, .(GeneProtein, DiPaGroups)]
BDL_short[DiPaGroups == 0, DiPaGroups := 8]
BDL_short[, DiPaGroups_BDL := factor(DiPaGroups, levels = as.character(1:8))]

CCL4_short <- DT_dipa_count_ccl4[, .(GeneProtein, DiPaGroups)]
CCL4_short[DiPaGroups == 0, DiPaGroups := 8]
CCL4_short[, DiPaGroups_CCL4 := factor(DiPaGroups, levels = as.character(1:8))]

# Merge to find overlapping GeneProtein pairs
Merged_DT <- merge(BDL_short, CCL4_short, by = "GeneProtein")
print(paste("Number of overlapping GeneProtein pairs:", nrow(Merged_DT)))

# 3. PERCENTAGE OVERLAP & OVERLAP INDEX (Di) ANALYSIS
# -------------------------------------------------------------------------
all_groups <- as.character(1:8)
A_background <- 943 # Total proteins/genes in the dataset

# Helper for significance stars
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("ns")
}

Overlap_Stats <- rbindlist(lapply(all_groups, function(g) {
  bdl_genes <- Merged_DT[DiPaGroups_BDL == g, GeneProtein]
  ccl4_genes <- Merged_DT[DiPaGroups_CCL4 == g, GeneProtein]
  common <- intersect(bdl_genes, ccl4_genes)
  
  O <- length(common)
  D <- length(bdl_genes)
  T_val <- length(ccl4_genes)
  A <- A_background
  
  # 1. Calculate Overlap Index (Di)
  di_val <- ifelse(D > 0 & T_val > 0, (O * A) / (T_val * D), 0)
  
  # 2. Fisher's Exact Test
  # Matrix: O, D-O, T-O, A-D-T+O
  contingency_matrix <- matrix(c(O, D - O, T_val - O, A - D - T_val + O), 
                               nrow = 2, byrow = TRUE)
  fisher_res <- fisher.test(contingency_matrix, alternative = "greater")
  p_val <- fisher_res$p.value
  
  list(
    Group       = g,
    N_BDL       = D,
    N_CCL4      = T_val,
    N_Common    = O,
    Pct_Stable  = ifelse(D > 0, round(100 * O / D, 1), 0),
    Overlap_Idx = round(di_val, 2),
    P_Value     = p_val,
    Stars       = get_stars(p_val)
  )
}))

print("Overlap Statistics with Overlap Index (Di):")
print(Overlap_Stats)

# 4. VISUALIZATION FUNCTIONS (Enhanced with Box and Labels)
# -------------------------------------------------------------------------
# Custom function to draw a mini-Venn with a grey box and internal labels
clean_pval <- function(pval, stars) {
  if (pval < 0.001) {
    formatted <- sprintf("%.1e", pval)
    parts <- strsplit(formatted, "e")[[1]]
    base <- parts[1]
    exponent <- as.numeric(parts[2])
    substitute(p == base %*% 10^exponent ~ (st), list(base = as.numeric(base), exponent = exponent, st = stars))
  } else {
    substitute(p == val ~ (st), list(val = round(pval, 3), st = stars))
  }
}

draw_mini_venn_grob <- function(group_name, n_bdl, n_ccl4, n_common, pct_stable, di, pval, stars) {
  formatted_pval_expr <- clean_pval(pval, stars)
  
  grob_out <- gTree(children = gList(
    # 1. Background Universe Box (height = 0.60, center shifted to y = 0.51 so top border is at y = 0.81)
    rectGrob(x = 0.5, y = 0.51, width = 0.94, height = 0.60, 
             gp = gpar(fill = "grey95", col = "grey60", lwd = 2.0)),
    
    # 2. Venn Circles (radius = 0.23, centered at y = 0.53)
    circleGrob(x = 0.35, y = 0.53, r = 0.23, gp = gpar(fill = "yellow", alpha = 0.5, col = "red", lwd = 1.5)),
    circleGrob(x = 0.65, y = 0.53, r = 0.23, gp = gpar(fill = "skyblue", alpha = 0.5, col = "blue", lwd = 1.5)),
    
    # 3. Internal Region Labels (D, O, T)
    textGrob("D", x = 0.22, y = 0.56, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    textGrob("T", x = 0.78, y = 0.56, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    textGrob("O", x = 0.50, y = 0.56, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    
    # 4. Internal Counts (Below Labels)
    textGrob(n_bdl,    x = 0.22, y = 0.48, gp = gpar(fontsize = 18, fontface = "bold")),
    textGrob(n_ccl4,   x = 0.78, y = 0.48, gp = gpar(fontsize = 18, fontface = "bold")),
    textGrob(n_common, x = 0.50, y = 0.48, gp = gpar(fontsize = 18, fontface = "bold")),
    
    # 5. Outer Labels (BDL, CCl4 with subscript 4)
    textGrob("BDL",  x = 0.20, y = 0.75, gp = gpar(fontsize = 20, fontface = "bold", col = "red")),
    textGrob(expression(bold(CCl[4])), x = 0.80, y = 0.75, gp = gpar(fontsize = 20, fontface = "bold", col = "blue")),
    
    # 6. Universe Count (A: 943 cleanly centered at the bottom of the box)
    textGrob("A = 943", x = 0.50, y = 0.26, gp = gpar(fontsize = 18, fontface = "bold", col = "grey30")),
    
    # 7. Title (DiPa group X elevated well above the box with clear whitespace)
    textGrob(paste("DiPa group", group_name), y = 0.93, gp = gpar(fontsize = 22, fontface = "bold")),
    
    # 8. Stats at the bottom (D_i and P-value expressions with clean spacing)
    textGrob(substitute(bold(D[i] == val), list(val = sprintf("%.2f", di))), 
             x = 0.5, y = 0.13, gp = gpar(fontsize = 20, fontface = "bold", col = "darkgreen")),
    textGrob(formatted_pval_expr, 
             x = 0.5, y = 0.05, gp = gpar(fontsize = 17, fontface = "bold"))
  ))
  return(grob_out)
}

# 5. GENERATE PDF
# -------------------------------------------------------------------------
print("Generating PDF results...")
cairo_pdf("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DiPa_Overlap_Analysis_Results.pdf", width = 14, height = 9.2)

# --- Page 3: Grid of Venns (Individual Groups 1-8) ---
venn_list <- lapply(all_groups, function(g) {
  row <- Overlap_Stats[Group == g]
  draw_mini_venn_grob(g, row$N_BDL, row$N_CCL4, row$N_Common, row$Pct_Stable, row$Overlap_Idx, row$P_Value, row$Stars)
})

grid.arrange(
  grobs = venn_list,
  layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)),
  heights = c(1, 1),
  padding = unit(2.2, "line")
)

dev.off()

# Save final stats
# fwrite(Overlap_Stats, "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DiPa_Overlap_Stats.csv")
# print("SUCCESS! File saved: DiPa_Overlap_Analysis_Results.pdf")
