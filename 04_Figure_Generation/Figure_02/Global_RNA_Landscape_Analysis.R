# ==============================================================================
# SCRIPT: Global_RNA_Landscape_Analysis.R
# OBJECTIVE: Characerize the global population response for RNA (Step 1)
#            Focus: Fold-change distributions and proportions of regulation.
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)

# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")

# Separate Datasets
# BDL has NA for TreatmentTime in this combined dataset
# CCL4 has numeric TreatmentTime (0, 2, 6, 12 months)
BDL_DT <- DTccl4_DT_LCPM[is.na(TreatmentTime)]
CCL4_DT <- DTccl4_DT_LCPM[!is.na(TreatmentTime)]

# 2. CALCULATE FOLD-CHANGES (GENERALIZED)
# -------------------------------------------------------------------------

calc_fc_bdl_rna <- function(dt) {
  # Looking for exactly "BDL" and "control" in Treatment column
  dt_sub <- dt[Treatment %in% c("BDL", "control")]
  dt_sub[, Condition := ifelse(Treatment == "BDL", "Treat", "Control")]
  
  # Use GeneCount for RNA
  stats <- dt_sub[, .(Med = median(GeneCount, na.rm = TRUE)), by = .(GeneProtein, Condition)]
  stats_wide <- dcast(stats, GeneProtein ~ Condition, value.var = "Med")
  
  if ("Treat" %in% colnames(stats_wide) & "Control" %in% colnames(stats_wide)) {
    stats_wide[, Log2FC := Treat - Control]
  } else {
    stats_wide[, Log2FC := NA_real_]
  }
  stats_wide[, Dataset := "BDL"]
  return(stats_wide[, .(GeneProtein, Log2FC, Dataset)])
}

calc_fc_ccl4_rna <- function(dt) {
  # Requested comparison: Treatment (ccl4 at time 12) vs Control (oil at time 0)
  dt_sub <- dt[(Treatment == "ccl4" & TreatmentTime == 12) | 
               (Treatment == "oil" & TreatmentTime == 0)]
  dt_sub[, Condition := ifelse(Treatment == "ccl4", "Treat", "Control")]
  
  # Use GeneCount for RNA
  stats <- dt_sub[, .(Med = median(GeneCount, na.rm = TRUE)), by = .(GeneProtein, Condition)]
  stats_wide <- dcast(stats, GeneProtein ~ Condition, value.var = "Med")
  
  if ("Treat" %in% colnames(stats_wide) & "Control" %in% colnames(stats_wide)) {
    stats_wide[, Log2FC := Treat - Control]
  } else {
    stats_wide[, Log2FC := NA_real_]
  }
  
  stats_wide[, Dataset := "CCL4"]
  return(stats_wide[, .(GeneProtein, Log2FC, Dataset)])
}

bdl_fc_rna  <- calc_fc_bdl_rna(BDL_DT)
ccl4_fc_rna <- calc_fc_ccl4_rna(CCL4_DT)
combined_fc_rna <- rbind(bdl_fc_rna, ccl4_fc_rna)

# Categorize (Biological Order: Blue at bottom, Red at top)
threshold <- 0.5
combined_fc_rna[, Status := "Stable"]
combined_fc_rna[Log2FC > threshold, Status := "Up-regulated"]
combined_fc_rna[Log2FC < -threshold, Status := "Down-regulated"]
combined_fc_rna[, Status := factor(Status, levels = c("Up-regulated", "Stable", "Down-regulated"))]

# -------------------------------------------------------------------------
# 3. VISUALIZATION
# -------------------------------------------------------------------------
pdf_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Global_RNA_Landscape_Summary.pdf"
cairo_pdf(pdf_path, width = 10, height = 8)

# Map Dataset to plotmath expression names for facet rendering
combined_fc_rna[, Dataset_Plot := ifelse(Dataset == "BDL", "BDL", "CCl[4]")]

# --- PANEL A: Density Plot (Distributions) ---
p_dist <- ggplot(combined_fc_rna, aes(x = Log2FC, fill = Dataset_Plot)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("BDL" = "skyblue", "CCl[4]" = "pink"), labels = c("BDL" = "BDL", "CCl[4]" = expression(CCl[4]))) +
  labs(title = "Global Distribution of RNA Changes (Median-based Log2FC)",
       subtitle = "Calculated as: Median(Disease RNA) - Median(Control RNA) | Vertical lines at |Log2FC| = 0.5",
       x = "Log2 Fold-Change (Median-based Delta)", 
       y = "Density",
       fill = "Dataset") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  facet_wrap(~Dataset_Plot, ncol = 2, labeller = label_parsed) # SIDE-BY-SIDE with parsed CCl4 label

# --- PANEL B: Proportional Bar Chart ---
prop_dt <- combined_fc_rna[, .N, by = .(Dataset, Status)]
prop_dt[, Percentage := round(100 * N / sum(N), 1), by = Dataset]
prop_dt[, Dataset_Plot := ifelse(Dataset == "BDL", "BDL", "CCl[4]")]

p_bar <- ggplot(prop_dt, aes(x = Dataset_Plot, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) + # Default: Level 1 (Down) is bottom, Level 3 (Up) is top
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), fontface="bold") +
  scale_fill_manual(values = c("Down-regulated" = "#2196F3", "Stable" = "grey80", "Up-regulated" = "#F44336")) +
  scale_x_discrete(labels = c("BDL" = "BDL", "CCl[4]" = expression(CCl[4]))) +
  labs(title = "Proportion of Regulated RNA (Median-based)",
       subtitle = "Status based on |Log2 Fold-Change| > 0.5",
       x = "Dataset", y = "Percentage of Total RNA Population (%)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(reverse = FALSE)) # Reverse legend to match Top-Down visual (Up at top)

# Arrange Page 1
grid.arrange(p_dist, p_bar, ncol = 1, heights = c(1.2, 1))

# --- PAGE 2: UNSUPERVISED HEATMAPS (ALL RNA) ---
# -------------------------------------------------------------------------
library(pheatmap)
library(gtable)
library(grid)

# --- Define custom pheatmap drawing overrides ---

# Custom draw_colnames: draw vertically (rot = 90) inside the annotation bar at y = 0.5 npc
draw_colnames_inside_annotation <- function (coln, gaps, vjust_col, hjust_col, angle_col, ...) 
{
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = grid::textGrob(coln, x = x, y = grid::unit(0.5, "npc"), 
                         vjust = 0.5, hjust = 0.5, rot = 90, 
                         gp = grid::gpar(fontface = "bold", ...))
    return(res)
}
assignInNamespace("draw_colnames", draw_colnames_inside_annotation, ns = "pheatmap")

# Custom draw_annotation_legend: support plotmath expressions (subscripts) in legend labels
draw_annotation_legend_parsed <- function(annotation, annotation_colors, border_color, ...) {
    y = unit(1, "npc")
    text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))
    res = gList()
    for (i in names(annotation)) {
        res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, 
            gp = gpar(fontface = "bold", ...))
        y = y - 1.5 * text_height
        if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
            n = length(annotation_colors[[i]])
            yy = y - (1:n - 1) * 2 * text_height
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), 
                y = yy, hjust = 0, vjust = 1, height = 2 * text_height, 
                width = 2 * text_height, gp = gpar(col = border_color, 
                  fill = annotation_colors[[i]]))
            
            labels_raw <- names(annotation_colors[[i]])
            labels_parsed <- lapply(labels_raw, function(x) {
                if (grepl("\\[|~", x)) parse(text = x)[[1]] else x
            })
            
            for(j in 1:length(labels_parsed)) {
                res[[paste(i, "t", j)]] = textGrob(labels_parsed[[j]], 
                    x = text_height * 2.4, y = yy[j] - text_height, 
                    hjust = 0, vjust = 0.5, gp = gpar(...))
            }
            y = y - n * 2 * text_height
        }
        else {
            yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 
                8 * text_height
            h = 8 * text_height * 0.25
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), 
                y = yy, hjust = 0, vjust = 1, height = h, width = 2 * 
                  text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
            res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), 
                y = y, hjust = 0, vjust = 1, height = 8 * text_height, 
                width = 2 * text_height, gp = gpar(col = border_color, 
                  fill = NA))
            txt = rev(range(grid.pretty(range(annotation[[i]], 
                na.rm = TRUE))))
            yy = y - c(1, 7) * text_height
            res[[paste(i, "t")]] = textGrob(txt, x = text_height * 
                2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
            y = y - 8 * text_height
        }
        y = y - 1.5 * text_height
    }
    res = gTree(children = res)
    return(res)
}
assignInNamespace("draw_annotation_legend", draw_annotation_legend_parsed, ns = "pheatmap")


# --- BDL UNSUPERVISED HEATMAP (n=943) ---
# -------------------------------------------------------------------------

# 1. Prepare BDL Matrix
bdl_heat_dt_rna <- BDL_DT[Treatment %in% c("BDL", "control"), .(GeneProtein, GeneCount, MiceInfo, Condition = Treatment)]
mat_wide_bdl_rna <- dcast(bdl_heat_dt_rna, GeneProtein ~ MiceInfo, value.var = "GeneCount")
mat_matrix_bdl_rna <- as.matrix(mat_wide_bdl_rna[, -1, with = FALSE])
rownames(mat_matrix_bdl_rna) <- mat_wide_bdl_rna$GeneProtein

# RNA is 100% complete
mat_matrix_bdl_rna <- na.omit(mat_matrix_bdl_rna)

# Map BDL MiceInfo to sequential numbers (Control: 1-6, Disease: 7-12)
bdl_mice_mapping <- c(
  "ShamvehicleM1" = "No. 1",
  "ShamvehicleM3" = "No. 2",
  "ShamvehicleM4" = "No. 3",
  "ShamvehicleM6" = "No. 4",
  "ShamvehicleM7" = "No. 5",
  "ShamvehicleM8" = "No. 6",
  "BDLvehicleM5"  = "No. 7",
  "BDLvehicleM7"  = "No. 8",
  "BDLvehicleM10" = "No. 9",
  "BDLvehicleM11" = "No. 10",
  "BDLvehicleM12" = "No. 11",
  "BDLvehicleM15" = "No. 12"
)
colnames(mat_matrix_bdl_rna) <- bdl_mice_mapping[colnames(mat_matrix_bdl_rna)]

# Annotation
anno_bdl_rna <- data.frame(Condition = unique(bdl_heat_dt_rna[, .(MiceInfo, Condition)])$Condition)
rownames(anno_bdl_rna) <- unique(bdl_heat_dt_rna[, .(MiceInfo, Condition)])$MiceInfo
rownames(anno_bdl_rna) <- bdl_mice_mapping[rownames(anno_bdl_rna)]
anno_bdl_rna$Condition <- factor(anno_bdl_rna$Condition, levels = c("BDL", "control"), labels = c("BDL", "Control"))

ann_colors <- list(Condition = c(BDL = "#F44336", Control = "#2196F3"))

# Generate BDL heatmap silently
res_bdl_rna <- pheatmap(mat_matrix_bdl_rna, 
                        cluster_rows = TRUE, cluster_cols = TRUE, 
                        show_rownames = FALSE, show_colnames = TRUE,
                        scale = "row", annotation_col = anno_bdl_rna,
                        annotation_colors = ann_colors,
                        main = "Individual Mice",
                        fontsize = 12, fontsize_col = 12,
                        color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100),
                        silent = TRUE)

gt_bdl_rna <- res_bdl_rna$gtable

# Shift column names (mouse numbers) to Row 3 (annotation bar row)
gt_bdl_rna$layout[gt_bdl_rna$layout$name == "col_names", c("t", "b")] <- 3
# Ensure col_names are drawn on top of the annotation colors (raise z-index)
gt_bdl_rna$layout[gt_bdl_rna$layout$name == "col_names", "z"] <- max(gt_bdl_rna$layout$z) + 1

# Make Row 3 (annotation row) taller to fit vertical text labels
gt_bdl_rna$heights[3] <- grid::unit(1.2, "cm")

# Modify col_annotation grob to fill the entire height of Row 3
anno_idx_bdl_rna <- which(gt_bdl_rna$layout$name == "col_annotation")
if (length(anno_idx_bdl_rna) > 0) {
    gt_bdl_rna$grobs[[anno_idx_bdl_rna]]$height <- grid::unit(1, "npc")
    gt_bdl_rna$grobs[[anno_idx_bdl_rna]]$y <- grid::unit(0.5, "npc")
    gt_bdl_rna$grobs[[anno_idx_bdl_rna]]$vjust <- 0.5
}

# Collapse Row 5 (original col_names row at bottom) height to 0
gt_bdl_rna$heights[5] <- grid::unit(0, "lines")

# Add a dedicated margin column on the far left for the Y-axis title
gt_bdl_rna <- gtable_add_cols(gt_bdl_rna, widths = grid::unit(1.2, "cm"), pos = 0)

# Subtract the new column's width (1.2 cm) from the matrix column (now Col 4)
# to keep the total width of the heatmap at exactly 1 npc (100% of page width)
gt_bdl_rna$widths[[4]] <- gt_bdl_rna$widths[[4]] - grid::unit(1.2, "cm")

y_label_bdl_rna <- grid::textGrob(paste0("Individual genes N = ", nrow(mat_matrix_bdl_rna)), 
                                  rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold"))
# The new column is Col 1. The matrix row in the updated gtable is Row 4.
gt_bdl_rna <- gtable_add_grob(gt_bdl_rna, y_label_bdl_rna, t = 4, l = 1, b = 4, r = 1, name = "y_axis_title")

# Draw the modified BDL heatmap
grid::grid.newpage()
grid::grid.draw(gt_bdl_rna)


# --- CCL4 UNSUPERVISED HEATMAP (n=943) ---
# -------------------------------------------------------------------------

# 1. Prepare CCL4 Matrix
ccl4_heat_dt_rna <- CCL4_DT[(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0), 
                             .(GeneProtein, GeneCount, MiceInfo, Condition = Treatment)]

# No exclusions for RNA analysis as data is already 100% complete
mat_wide_ccl4_rna <- dcast(ccl4_heat_dt_rna, GeneProtein ~ MiceInfo, value.var = "GeneCount")
mat_matrix_ccl4_rna <- as.matrix(mat_wide_ccl4_rna[, -1, with = FALSE])
rownames(mat_matrix_ccl4_rna) <- mat_wide_ccl4_rna$GeneProtein

# Complete cases
mat_matrix_ccl4_rna <- na.omit(mat_matrix_ccl4_rna)

# Map CCL4 MiceInfo to sequential numbers (Control: 1-6, Disease: 7-12)
ccl4_mice_mapping <- c(
  "month0_oil_rep1"    = "No. 1",
  "month0_oil_rep2"    = "No. 2",
  "month0_oil_rep3"    = "No. 3",
  "month0_oil_rep4"    = "No. 4",
  "month0_oil_rep5"    = "No. 5",
  "month0_oil_rep6"    = "No. 6",
  "month12_ccl4_rep1"  = "No. 7",
  "month12_ccl4_rep2"  = "No. 8",
  "month12_ccl4_rep3"  = "No. 9",
  "month12_ccl4_rep4"  = "No. 10",
  "month12_ccl4_rep5"  = "No. 11",
  "month12_ccl4_rep6"  = "No. 12"
)
colnames(mat_matrix_ccl4_rna) <- ccl4_mice_mapping[colnames(mat_matrix_ccl4_rna)]

# Annotation with plotmath-compatible strings for legend parsing
anno_ccl4_rna <- data.frame(Condition = unique(ccl4_heat_dt_rna[, .(MiceInfo, Condition)])$Condition)
rownames(anno_ccl4_rna) <- unique(ccl4_heat_dt_rna[, .(MiceInfo, Condition)])$MiceInfo
rownames(anno_ccl4_rna) <- ccl4_mice_mapping[rownames(anno_ccl4_rna)]
anno_ccl4_rna$Condition <- factor(anno_ccl4_rna$Condition, 
                                  levels = c("ccl4", "oil"), 
                                  labels = c("CCl[4]~(Month~12)", "Oil~(Month~0)"))

ann_colors_ccl4 <- list(Condition = c("CCl[4]~(Month~12)" = "#F44336", "Oil~(Month~0)" = "#2196F3"))

# Generate CCL4 heatmap silently
res_ccl4_rna <- pheatmap(mat_matrix_ccl4_rna, 
                         cluster_rows = TRUE, cluster_cols = TRUE, 
                         show_rownames = FALSE, show_colnames = TRUE,
                         scale = "row", annotation_col = anno_ccl4_rna,
                         annotation_colors = ann_colors_ccl4,
                         main = "Individual Mice",
                         fontsize = 12, fontsize_col = 12,
                         color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100),
                         silent = TRUE)

gt_ccl4_rna <- res_ccl4_rna$gtable

# Shift column names (mouse numbers) to Row 3 (annotation bar row)
gt_ccl4_rna$layout[gt_ccl4_rna$layout$name == "col_names", c("t", "b")] <- 3
# Ensure col_names are drawn on top of the annotation colors (raise z-index)
gt_ccl4_rna$layout[gt_ccl4_rna$layout$name == "col_names", "z"] <- max(gt_ccl4_rna$layout$z) + 1

# Make Row 3 (annotation row) taller to fit vertical text labels
gt_ccl4_rna$heights[3] <- grid::unit(1.2, "cm")

# Modify col_annotation grob to fill the entire height of Row 3
anno_idx_ccl4_rna <- which(gt_ccl4_rna$layout$name == "col_annotation")
if (length(anno_idx_ccl4_rna) > 0) {
    gt_ccl4_rna$grobs[[anno_idx_ccl4_rna]]$height <- grid::unit(1, "npc")
    gt_ccl4_rna$grobs[[anno_idx_ccl4_rna]]$y <- grid::unit(0.5, "npc")
    gt_ccl4_rna$grobs[[anno_idx_ccl4_rna]]$vjust <- 0.5
}

# Collapse Row 5 (original col_names row at bottom) height to 0
gt_ccl4_rna$heights[5] <- grid::unit(0, "lines")

# Add a dedicated margin column on the far left for the Y-axis title
gt_ccl4_rna <- gtable_add_cols(gt_ccl4_rna, widths = grid::unit(1.2, "cm"), pos = 0)

# Subtract the new column's width (1.2 cm) from the matrix column (now Col 4)
# to keep the total width of the heatmap at exactly 1 npc (100% of page width)
gt_ccl4_rna$widths[[4]] <- gt_ccl4_rna$widths[[4]] - grid::unit(1.2, "cm")

y_label_ccl4_rna <- grid::textGrob(paste0("Individual genes N = ", nrow(mat_matrix_ccl4_rna)), 
                                    rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold"))
# The new column is Col 1. The matrix row in the updated gtable is Row 4.
gt_ccl4_rna <- gtable_add_grob(gt_ccl4_rna, y_label_ccl4_rna, t = 4, l = 1, b = 4, r = 1, name = "y_axis_title")

# Draw the modified CCL4 heatmap
grid::grid.newpage()
grid::grid.draw(gt_ccl4_rna)


dev.off()

print(paste("SUCCESS! Global RNA Landscape Report saved to:", pdf_path))

