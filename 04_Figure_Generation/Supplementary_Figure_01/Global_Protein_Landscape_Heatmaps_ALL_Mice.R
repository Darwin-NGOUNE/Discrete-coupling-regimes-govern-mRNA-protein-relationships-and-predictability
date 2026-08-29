# ==============================================================================
# SCRIPT: Global_Protein_Landscape_Heatmaps_ALL_Mice.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/
# OBJECTIVE: Generate a clean 2-Page PDF report containing ONLY the Unsupervised Heatmaps
#            Page 1: BDL Cohort (All 18 Mice including BDL_ASBTi Zwischenmäuse)
#            Page 2: CCL4 Cohort (All Mice including Oil M2, Oil M12, CCl4 M2, CCl4 M6)
# ==============================================================================

library(data.table)
library(ggplot2)
library(pheatmap)
library(gtable)
library(grid)

# Output directory & PDF path
out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"
pdf_path <- file.path(out_dir, "Global_Protein_Landscape_Heatmaps_ALL_Mice.pdf")

# -------------------------------------------------------------------------
# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
Full_DT <- as.data.table(DTccl4_DT_LCPM)

BDL_DT  <- Full_DT[is.na(TreatmentTime)]
CCL4_DT <- Full_DT[!is.na(TreatmentTime)]

# -------------------------------------------------------------------------
# 2. DEFINE CUSTOM PHEATMAP OVERRIDES (Exact formatting from landscape script)
# -------------------------------------------------------------------------

# Custom draw_colnames: draw vertically inside the annotation bar
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

# Open PDF device (2 Pages)
cairo_pdf(pdf_path, width = 11, height = 8.5)

# ==============================================================================
# PAGE 1: BDL UNSUPERVISED HEATMAP (ALL 18 MICE WITH ZWISCHENMÄUSE)
# ==============================================================================
print("Generating Page 1: BDL Heatmap (18 Mice)...")

bdl_heatmap_dt <- BDL_DT[!is.na(ProteinIntensity), .(GeneProtein, ProteinIntensity, MiceInfo, Treatment)]

# Map 3 Treatment groups
bdl_heatmap_dt[, Condition := fcase(
  Treatment == "control",   "Control",
  Treatment == "BDL_ASBTi", "BDL_ASBTi",
  Treatment == "BDL",       "BDL"
)]

mat_wide_bdl <- dcast(bdl_heatmap_dt, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_matrix_bdl <- as.matrix(mat_wide_bdl[, -1, with = FALSE])
rownames(mat_matrix_bdl) <- mat_wide_bdl$GeneProtein
mat_matrix_bdl <- na.omit(mat_matrix_bdl)

# Map 18 BDL mice to sequential numbers No. 1 to No. 18
bdl_mice_mapping <- c(
  "ShamvehicleM1" = "No. 1",
  "ShamvehicleM3" = "No. 2",
  "ShamvehicleM4" = "No. 3",
  "ShamvehicleM6" = "No. 4",
  "ShamvehicleM7" = "No. 5",
  "ShamvehicleM8" = "No. 6",
  "BDLASBTiM8"    = "No. 7",
  "BDLASBTiM9"    = "No. 8",
  "BDLASBTiM10"   = "No. 9",
  "BDLASBTiM11"   = "No. 10",
  "BDLASBTiM12"   = "No. 11",
  "BDLASBTiM14"   = "No. 12",
  "BDLvehicleM5"  = "No. 13",
  "BDLvehicleM7"  = "No. 14",
  "BDLvehicleM10" = "No. 15",
  "BDLvehicleM11" = "No. 16",
  "BDLvehicleM12" = "No. 17",
  "BDLvehicleM15" = "No. 18"
)
colnames(mat_matrix_bdl) <- bdl_mice_mapping[colnames(mat_matrix_bdl)]

# Annotations for BDL
meta_bdl <- unique(bdl_heatmap_dt[, .(MiceInfo, Condition)])
anno_bdl <- data.frame(Condition = meta_bdl$Condition)
rownames(anno_bdl) <- bdl_mice_mapping[meta_bdl$MiceInfo]
anno_bdl$Condition <- factor(anno_bdl$Condition, levels = c("Control", "BDL_ASBTi", "BDL"))

ann_colors_bdl <- list(
  Condition = c(
    "Control"   = "#2196F3", # Blue
    "BDL_ASBTi" = "#4CAF50", # Green (Zwischenmäuse)
    "BDL"       = "#F44336"  # Red
  )
)

res_bdl <- pheatmap(mat_matrix_bdl, 
                    cluster_rows = TRUE, cluster_cols = TRUE, 
                    show_rownames = FALSE, show_colnames = TRUE,
                    scale = "row", annotation_col = anno_bdl,
                    annotation_colors = ann_colors_bdl,
                    main = "Individual Mice",
                    fontsize = 12, fontsize_col = 11,
                    color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100),
                    silent = TRUE)

gt_bdl <- res_bdl$gtable

# Shift column names (mouse numbers) to Row 3 (annotation bar row)
gt_bdl$layout[gt_bdl$layout$name == "col_names", c("t", "b")] <- 3
gt_bdl$layout[gt_bdl$layout$name == "col_names", "z"] <- max(gt_bdl$layout$z) + 1

# Make Row 3 taller to fit vertical text labels
gt_bdl$heights[3] <- grid::unit(1.2, "cm")

# Modify col_annotation grob to fill the entire height of Row 3
anno_idx_bdl <- which(gt_bdl$layout$name == "col_annotation")
if (length(anno_idx_bdl) > 0) {
    gt_bdl$grobs[[anno_idx_bdl]]$height <- grid::unit(1, "npc")
    gt_bdl$grobs[[anno_idx_bdl]]$y <- grid::unit(0.5, "npc")
    gt_bdl$grobs[[anno_idx_bdl]]$vjust <- 0.5
}

# Collapse Row 5 (original col_names row at bottom)
gt_bdl$heights[5] <- grid::unit(0, "lines")

# Add a dedicated margin column on the far left for the Y-axis title
gt_bdl <- gtable_add_cols(gt_bdl, widths = grid::unit(1.2, "cm"), pos = 0)
gt_bdl$widths[[4]] <- gt_bdl$widths[[4]] - grid::unit(1.2, "cm")

y_label_bdl <- grid::textGrob(paste0("Individual proteins N = ", nrow(mat_matrix_bdl)), 
                              rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold"))
gt_bdl <- gtable_add_grob(gt_bdl, y_label_bdl, t = 4, l = 1, b = 4, r = 1, name = "y_axis_title")

# Draw Page 1
grid::grid.newpage()
grid::grid.draw(gt_bdl)


# ==============================================================================
# PAGE 2: CCL4 UNSUPERVISED HEATMAP (ALL MICE WITH ZWISCHENMÄUSE)
# ==============================================================================
print("Generating Page 2: CCL4 Heatmap (All Mice with Zwischenmäuse)...")

ccl4_heatmap_dt <- CCL4_DT[!is.na(ProteinIntensity) & MiceInfo != "month12_ccl4_rep6", 
                           .(GeneProtein, ProteinIntensity, MiceInfo, Treatment, TreatmentTime)]

# Map all 6 Treatment Groups for CCl4
ccl4_heatmap_dt[, Condition := fcase(
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 0,  "Oil~(Month~0)",
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 2,  "Oil~(Month~2)",
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 12, "Oil~(Month~12)",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 2,  "CCl[4]~(Month~2)",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 6,  "CCl[4]~(Month~6)",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 12, "CCl[4]~(Month~12)"
)]

mat_wide_ccl4 <- dcast(ccl4_heatmap_dt, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_matrix_ccl4 <- as.matrix(mat_wide_ccl4[, -1, with = FALSE])
rownames(mat_matrix_ccl4) <- mat_wide_ccl4$GeneProtein
mat_matrix_ccl4 <- na.omit(mat_matrix_ccl4)

# Sequential mouse numbering for CCL4
all_ccl4_mice <- unique(ccl4_heatmap_dt$MiceInfo)
meta_ccl4 <- unique(ccl4_heatmap_dt[, .(MiceInfo, Condition, Treatment, TreatmentTime)])
setorder(meta_ccl4, Treatment, TreatmentTime, MiceInfo)
meta_ccl4[, MouseNum := paste0("No. ", .I)]

ccl4_mice_mapping <- setNames(meta_ccl4$MouseNum, meta_ccl4$MiceInfo)
colnames(mat_matrix_ccl4) <- ccl4_mice_mapping[colnames(mat_matrix_ccl4)]

anno_ccl4 <- data.frame(Condition = meta_ccl4$Condition)
rownames(anno_ccl4) <- meta_ccl4$MouseNum
anno_ccl4$Condition <- factor(anno_ccl4$Condition, levels = c(
  "Oil~(Month~0)", "Oil~(Month~2)", "Oil~(Month~12)",
  "CCl[4]~(Month~2)", "CCl[4]~(Month~6)", "CCl[4]~(Month~12)"
))

# 6 Highly distinct colors
ann_colors_ccl4 <- list(
  Condition = c(
    "Oil~(Month~0)"    = "#1E88E5", # Pure Blue (Control)
    "Oil~(Month~2)"    = "#00ACC1", # Cyan / Teal
    "Oil~(Month~12)"   = "#8E24AA", # Purple
    "CCl[4]~(Month~2)"  = "#FF9800", # Orange
    "CCl[4]~(Month~6)"  = "#43A047", # Green
    "CCl[4]~(Month~12)" = "#E53935"  # Pure Red (Disease)
  )
)

res_ccl4 <- pheatmap(mat_matrix_ccl4, 
                     cluster_rows = TRUE, cluster_cols = TRUE, 
                     show_rownames = FALSE, show_colnames = TRUE,
                     scale = "row", annotation_col = anno_ccl4,
                     annotation_colors = ann_colors_ccl4,
                     main = "Individual Mice",
                     fontsize = 12, fontsize_col = 10,
                     color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100),
                     silent = TRUE)

gt_ccl4 <- res_ccl4$gtable

# Shift column names (mouse numbers) to Row 3 (annotation bar row)
gt_ccl4$layout[gt_ccl4$layout$name == "col_names", c("t", "b")] <- 3
gt_ccl4$layout[gt_ccl4$layout$name == "col_names", "z"] <- max(gt_ccl4$layout$z) + 1

# Make Row 3 taller to fit vertical text labels
gt_ccl4$heights[3] <- grid::unit(1.2, "cm")

# Modify col_annotation grob to fill the entire height of Row 3
anno_idx_ccl4 <- which(gt_ccl4$layout$name == "col_annotation")
if (length(anno_idx_ccl4) > 0) {
    gt_ccl4$grobs[[anno_idx_ccl4]]$height <- grid::unit(1, "npc")
    gt_ccl4$grobs[[anno_idx_ccl4]]$y <- grid::unit(0.5, "npc")
    gt_ccl4$grobs[[anno_idx_ccl4]]$vjust <- 0.5
}

# Collapse Row 5 (original col_names row at bottom)
gt_ccl4$heights[5] <- grid::unit(0, "lines")

# Add a dedicated margin column on the far left for the Y-axis title
gt_ccl4 <- gtable_add_cols(gt_ccl4, widths = grid::unit(1.2, "cm"), pos = 0)
gt_ccl4$widths[[4]] <- gt_ccl4$widths[[4]] - grid::unit(1.2, "cm")

y_label_ccl4 <- grid::textGrob(paste0("Individual proteins N = ", nrow(mat_matrix_ccl4)), 
                                rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold"))
gt_ccl4 <- gtable_add_grob(gt_ccl4, y_label_ccl4, t = 4, l = 1, b = 4, r = 1, name = "y_axis_title")

# Draw Page 2
grid::grid.newpage()
grid::grid.draw(gt_ccl4)

dev.off()

print(paste("SUCCESS! 2-Page Heatmaps Report saved to:", pdf_path))
