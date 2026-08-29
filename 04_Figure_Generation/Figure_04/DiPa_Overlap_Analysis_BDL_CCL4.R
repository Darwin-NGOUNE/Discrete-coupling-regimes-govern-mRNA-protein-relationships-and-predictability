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
    # 1. Background Universe Box (height increased to 0.78 and center shifted down to y=0.54 for safe margins)
    rectGrob(x = 0.5, y = 0.54, width = 0.95, height = 0.78, 
             gp = gpar(fill = "grey95", col = "grey60", lwd = 2.0)),
    
    # 2. Venn Circles (radius increased to 0.30 and shifted up to y=0.58 to prevent overlap with A at the bottom)
    circleGrob(x = 0.35, y = 0.58, r = 0.30, gp = gpar(fill = "yellow", alpha = 0.5, col = "red", lwd = 1.5)),
    circleGrob(x = 0.65, y = 0.58, r = 0.30, gp = gpar(fill = "skyblue", alpha = 0.5, col = "blue", lwd = 1.5)),
    
    # 3. Internal Region Labels (D, O, T) - tightly grouped vertically
    textGrob("D", x = 0.22, y = 0.61, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    textGrob("T", x = 0.78, y = 0.61, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    textGrob("O", x = 0.50, y = 0.61, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    
    # 4. Internal Counts (Below Labels) - directly centered below labels
    textGrob(n_bdl,    x = 0.22, y = 0.53, gp = gpar(fontsize = 18, fontface = "bold")),
    textGrob(n_ccl4,   x = 0.78, y = 0.53, gp = gpar(fontsize = 18, fontface = "bold")),
    textGrob(n_common, x = 0.50, y = 0.53, gp = gpar(fontsize = 18, fontface = "bold")),
    
    # 5. Outer Labels (BDL, CCl4 with subscript 4)
    textGrob("BDL",  x = 0.20, y = 0.87, gp = gpar(fontsize = 22, fontface = "bold", col = "red")),
    textGrob(expression(CCl[4]), x = 0.80, y = 0.87, gp = gpar(fontsize = 22, fontface = "bold", col = "blue")),
    
    # 6. Universe Count (A and 943 stacked at the bottom of the box, safe from circle overlap and bottom border)
    textGrob("A", x = 0.50, y = 0.23, gp = gpar(fontsize = 20, fontface = "bold", col = "grey20")),
    textGrob("943", x = 0.50, y = 0.18, gp = gpar(fontsize = 18, fontface = "bold")),
    
    # 7. Title (DiPa Group X)
    textGrob(paste("DiPa Group", group_name), y = 0.96, gp = gpar(fontsize = 22, fontface = "bold")),
    
    # 8. Stats at the bottom (D_i and P-value expressions)
    textGrob(substitute(D[i] == val, list(val = sprintf("%.2f", di))), 
             x = 0.5, y = 0.11, gp = gpar(fontsize = 22, fontface = "bold", col="darkgreen")),
    textGrob(formatted_pval_expr, 
             x = 0.5, y = 0.03, gp = gpar(fontsize = 18, fontface = "bold"))
  ))
  return(grob_out)
}

# 5. GENERATE PDF
# -------------------------------------------------------------------------
print("Generating PDF results...")
pdf("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DiPa_Overlap_Analysis_Results.pdf", width = 12, height = 10)

# # --- Page 1: Continuity Matrix (Heatmap) ---
# Matrix_DT <- Merged_DT[, .N, by = .(DiPaGroups_BDL, DiPaGroups_CCL4)]
# p_matrix <- ggplot(Matrix_DT, aes(x = DiPaGroups_BDL, y = DiPaGroups_CCL4, fill = N)) +
#   geom_tile() +
#   geom_text(aes(label = N), color = "black", size = 3.5) +
#   scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
#   theme_minimal() +
#   labs(title = "DiPa Group Assignment Continuity Matrix",
#        subtitle = "Rows show CCL4 group for a given BDL group (and vice versa)",
#        x = "BDL DiPa Category", y = "CCL4 DiPa Category")
# print(p_matrix)

# # --- Page 2: Stability Bar Plot ---
# p_bar <- ggplot(Overlap_Stats, aes(x = factor(Group), y = Pct_Stable, fill = factor(Group))) +
#   geom_bar(stat = "identity", color = "black", width = 0.7) +
#   geom_text(aes(label = paste0(Pct_Stable, "%")), vjust = -0.5, fontface = "bold") +
#   scale_fill_brewer(palette = "Set1") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(
#     title = "Stability of DiPa Category between BDL and CCL4",
#     subtitle = "% of BDL genes remaining in the same category in CCL4",
#     x = "DiPa Group", y = "% Matching Category (relative to BDL)"
#   )
# print(p_bar)

# # --- Page 2: Global Venn (Groups 1-7) ---
# bdl_de  <- Merged_DT[DiPaGroups_BDL %in% as.character(1:7), GeneProtein]
# ccl4_de <- Merged_DT[DiPaGroups_CCL4 %in% as.character(1:7), GeneProtein]
# common_de <- intersect(bdl_de, ccl4_de)
# pct_de  <- round(100 * length(common_de) / length(bdl_de), 1)
# 
# O_de <- length(common_de)
# D_de <- length(bdl_de)
# T_de <- length(ccl4_de)
# A_de <- A_background
# di_de <- (O_de * A_de) / (T_de * D_de)
# 
# # Fisher Global
# con_de <- matrix(c(O_de, D_de - O_de, T_de - O_de, A_de - D_de - T_de + O_de), nrow = 2, byrow = TRUE)
# p_de <- fisher.test(con_de, alternative = "greater")$p.value
# 
# grid.newpage()
# # 1. Background Universe Box
# grid.rect(x = 0.5, y = 0.5, width = 0.7, height = 0.7, 
#           gp = gpar(fill = "grey95", col = "grey60", lwd = 2))
# 
# # 2. Venn Circles
# grid.circle(x = 0.4, y = 0.5, r = 0.22, gp = gpar(fill = "yellow", alpha = 0.5, col = "red"))
# grid.circle(x = 0.6, y = 0.5, r = 0.22, gp = gpar(fill = "skyblue", alpha = 0.5, col = "blue"))
# 
# # 3. Internal Labels (D, O, T)
# grid.text("D", x = 0.3, y = 0.58, gp = gpar(fontsize = 18, fontface = "bold", col = "grey20"))
# grid.text("T", x = 0.7, y = 0.58, gp = gpar(fontsize = 18, fontface = "bold", col = "grey20"))
# grid.text("O", x = 0.5, y = 0.58, gp = gpar(fontsize = 18, fontface = "bold", col = "grey20"))
# 
# # 4. Internal Counts (Below Labels)
# grid.text(D_de, x = 0.3, y = 0.45, gp = gpar(fontsize = 16, fontface = "bold"))
# grid.text(T_de, x = 0.7, y = 0.45, gp = gpar(fontsize = 16, fontface = "bold"))
# grid.text(O_de, x = 0.5, y = 0.45, gp = gpar(fontsize = 16, fontface = "bold"))
# 
# # 5. Outer Labels
# grid.text("BDL (Groups 1-7)", x = 0.3, y = 0.75, gp = gpar(fontsize = 15, fontface = "bold", col="red"))
# grid.text("CCL4 (Groups 1-7)", x = 0.7, y = 0.75, gp = gpar(fontsize = 15, fontface = "bold", col="blue"))
# 
# # 6. Universe Label (A: 943)
# grid.text("A: 943", x = 0.78, y = 0.22, gp = gpar(fontsize = 15, fontface = "bold"))
# 
# # 7. Stats
# grid.text(paste0("Overlap Index (D_i): ", round(di_de, 2)), x = 0.5, y = 0.1, gp = gpar(fontsize = 16, fontface = "bold", col="darkgreen"))
# grid.text(paste0("Fisher P-val: ", format.pval(p_de, digits=3), " (", get_stars(p_de), ")"), x = 0.5, y = 0.04, gp = gpar(fontsize = 12))
# 
# grid.text("Global Overlap of Differentially Regulated Genes", y = 0.94, gp = gpar(fontsize = 18, fontface = "bold"))

# --- Page 3: Grid of Venns (Individual Groups 1-8) ---
venn_list <- lapply(all_groups, function(g) {
  row <- Overlap_Stats[Group == g]
  draw_mini_venn_grob(g, row$N_BDL, row$N_CCL4, row$N_Common, row$Pct_Stable, row$Overlap_Idx, row$P_Value, row$Stars)
})

grid.arrange(
  grobs = venn_list,
  ncol = 4,
  padding = unit(1, "line")
)

# --- Page 5: Stats Table ---
# p_table <- tableGrob(Overlap_Stats)
# grid.newpage()
# grid.draw(p_table)

dev.off()

# Save final stats
# fwrite(Overlap_Stats, "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DiPa_Overlap_Stats.csv")
# print("SUCCESS! File saved: DiPa_Overlap_Analysis_Results.pdf")
