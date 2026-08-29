# ==============================================================================
# SCRIPT: Generate_A4_MASTER_BDL_Schwerpunkt_ALL_18Mice.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/
# OBJECTIVE: Strictly preserve 100% of the original layout, styling, and aesthetics
#            from Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R
#            The ONLY difference is the integration of the 6 Zwischenmäuse (18 BDL mice total).
#            PDF-ONLY outputs (no PNG).
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(ggforce)
library(ggrepel)

# Text sizes strictly identical to original script
T_TITLE = 20
T_SUB   = 14
T_AXIS  = 16
T_TICKS = 14

out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"

# -------------------------------------------------------------------------
# 1. LOAD DATA GLOBALLY (ORIGINAL RAW DATASET WITH 18 MICE)
# -------------------------------------------------------------------------
print("Loading raw DT for PCA, Volcano, DiPa, and Clouds (18 BDL Mice)...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
Full_DT <- as.data.table(DTccl4_DT_LCPM)
Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]

# Map all 3 groups for BDL
Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "control",   DiseaseGroup := "Control"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL_ASBTi", DiseaseGroup := "BDL_ASBTi"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL",       DiseaseGroup := "Disease"]

bdl_DT <- Full_DT[Dataset == "BDL Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]
bdl_DT$DiseaseGroup <- factor(bdl_DT$DiseaseGroup, levels = c("Control", "BDL_ASBTi", "Disease"))

cat(sprintf("Loaded %d rows across %d mice for BDL cohort.\n", nrow(bdl_DT), uniqueN(bdl_DT$MiceInfo)))

# Compute DiPa Ratios directly on BDL Control vs Disease
dipa_stats <- bdl_DT[, .(
  meanG_C = mean(GeneCount[DiseaseGroup == "Control"], na.rm = TRUE),
  meanG_D = mean(GeneCount[DiseaseGroup == "Disease"], na.rm = TRUE),
  meanP_C = mean(Protein_Raw[DiseaseGroup == "Control"], na.rm = TRUE),
  meanP_D = mean(Protein_Raw[DiseaseGroup == "Disease"], na.rm = TRUE),
  DiPaGroups = as.character(ClusterDiPa[1])
), by = GeneProtein]

dipa_stats[DiPaGroups == "0", DiPaGroups := "8"]
dipa_stats[, log2G := log2(meanG_D / meanG_C)]
dipa_stats[, log2P := log2(meanP_D / meanP_C)]
dipa_stats <- dipa_stats[is.finite(log2G) & is.finite(log2P)]

# 3-Group Color Palette
group_colors <- c(
  "Control"   = "#2196F3", # Blue (Sham Vehicle)
  "BDL_ASBTi" = "#4CAF50", # Green (BDL ASBTi Zwischenmäuse)
  "Disease"   = "#F44336"  # Red (BDL Vehicle)
)
group_labels <- c(
  "Control"   = "Sham Vehicle",
  "BDL_ASBTi" = "BDL ASBTi",
  "Disease"   = "BDL Vehicle"
)

# -------------------------------------------------------------------------
# 2. STEP 1: PCA PLOTS (18 MICE, IDENTICAL STRUCTURE)
# -------------------------------------------------------------------------
generate_master_pca <- function(dt, value_col, title, ctrl_label, asbti_label, dis_label) {
  wide_data <- dcast(dt, GeneProtein ~ MiceInfo, value.var = value_col, fun.aggregate = mean)
  mat <- as.matrix(wide_data[, -1, with = FALSE])
  rownames(mat) <- wide_data$GeneProtein
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 3) return(ggplot() + theme_void() + annotate("text", x=0, y=0, label="Omitted\nSparsity"))
  mat <- mat[apply(mat, 1, var, na.rm=TRUE) > 0, ]
  pca_res <- prcomp(t(mat), scale. = TRUE)
  var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  pca_df <- data.table(MiceInfo = rownames(pca_res$x), PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2])
  meta <- unique(dt[, .(MiceInfo, DiseaseGroup)])
  pca_df <- merge(pca_df, meta, by = "MiceInfo")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = DiseaseGroup)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(
      values = group_colors, 
      labels = c("Control" = ctrl_label, "BDL_ASBTi" = asbti_label, "Disease" = dis_label), 
      name = ""
    ) +
    labs(title = title, x = paste0("PC1: ", var_explained[1], "% variance"), y = paste0("PC2: ", var_explained[2], "% variance")) +
    theme_bw(base_size = T_AXIS) +
    theme(
      plot.title = element_text(size = T_TITLE, hjust = 0.5), axis.title = element_text(size = T_AXIS), axis.text = element_text(size = T_TICKS, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "bottom", legend.title = element_text(size = T_AXIS), legend.text = element_text(size = T_TICKS), legend.key.size = unit(0.4, "cm"),
      legend.margin=margin(t=-10)
    )
}

p1_pca_rna  <- generate_master_pca(bdl_DT, "GeneCount", "RNA", "Sham Vehicle", "BDL ASBTi", "BDL Vehicle") + 
  labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 20), legend.position = "none")
p1_pca_prot <- generate_master_pca(bdl_DT, "Protein_Raw", "Protein", "Sham Vehicle", "BDL ASBTi", "BDL Vehicle") + 
  theme(legend.position = "none")

# -------------------------------------------------------------------------
# 3. STEP 2: VOLCANO PLOTS (EXACT SAME STRUCTURE AND MATH)
# -------------------------------------------------------------------------
calculate_volcano <- function(dt, value_col) {
  dt_sub <- dt[!is.na(get(value_col)) & DiseaseGroup %in% c("Control", "Disease")]
  stats <- dt_sub[, {
    v_c <- get(value_col)[DiseaseGroup == "Control"]; v_c <- v_c[!is.na(v_c)]
    v_d <- get(value_col)[DiseaseGroup == "Disease"]; v_d <- v_d[!is.na(v_d)]
    if(length(v_c)>=2 && length(v_d)>=2) {
      logFC <- mean(v_d) - mean(v_c)
      pval  <- tryCatch(t.test(v_d, v_c)$p.value, error = function(e) NA_real_)
    } else { logFC <- NA_real_; pval <- NA_real_ }
    .(LogFC_Val = logFC, Pval_Val = pval)
  }, by = GeneProtein]
  stats <- stats[!is.na(Pval_Val)]
  stats[, padj_Val := p.adjust(Pval_Val, method = "BH")]
  stats[, Status := "Not Sig"]
  stats[padj_Val < 0.05 & LogFC_Val > 1,  Status := "Up"]
  stats[padj_Val < 0.05 & LogFC_Val < -1, Status := "Down"]
  return(stats)
}

plot_master_volcano <- function(stats_dt, title_str, text_size = 4.5) {
  n_up   <- sum(stats_dt$Status == "Up")
  n_down <- sum(stats_dt$Status == "Down")
  
  x_max_val <- max(abs(stats_dt$LogFC_Val), na.rm = TRUE)
  y_max_val <- max(-log10(stats_dt$padj_Val), na.rm = TRUE)
  
  x_text_pos <- x_max_val * 0.95
  y_text_pos <- y_max_val * 1.03
  
  ggplot(stats_dt, aes(x = LogFC_Val, y = -log10(padj_Val), color = Status)) +
    geom_point(aes(size = Status, alpha = Status)) +
    scale_size_manual(values = c("Not Sig" = 0.4, "Down" = 1.2, "Up" = 1.2)) +
    scale_alpha_manual(values = c("Not Sig" = 0.4, "Down" = 0.8, "Up" = 0.8)) +
    scale_color_manual(values = c("Down" = "red", "Not Sig" = "black", "Up" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.6) +
    annotate("text", x = -x_text_pos, y = y_text_pos, 
             label = paste0("Down: ", format(n_down, big.mark=",")), 
             color = "black", fontface = "bold", size = text_size, hjust = 0) +
    annotate("text", x = x_text_pos, y = y_text_pos, 
             label = paste0("Up: ", format(n_up, big.mark=",")), 
             color = "black", fontface = "bold", size = text_size, hjust = 1) +
    labs(
      title = title_str, 
      subtitle = NULL, 
      x = expression(Log[2]~"(fold change)"), 
      y = expression(-Log[10]~"(P value)")
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = T_AXIS) +
    theme(
      plot.title = element_text(size = T_TITLE, hjust = 0.5, margin = margin(b = 15)), 
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
      axis.title = element_text(size = T_AXIS), axis.text = element_text(size = T_TICKS, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none"
    )
}

p2_volc_rna  <- plot_master_volcano(calculate_volcano(bdl_DT, "GeneCount"), "RNA") + 
  labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 20), legend.position = "none")
p2_volc_prot <- plot_master_volcano(calculate_volcano(bdl_DT, "Protein_Raw"), "Protein") + 
  theme(legend.position = "none")

# -------------------------------------------------------------------------
# 4. STEP 3: DIPA & CLOUDS — WITH 3 SCHWERPUNKTE ACROSS 18 MICE
# -------------------------------------------------------------------------
Cloud_Metrics_BDL <- bdl_DT[, {
  rna_c  <- GeneCount[DiseaseGroup == "Control"]
  prot_c <- Protein_Raw[DiseaseGroup == "Control"]
  
  rna_z  <- GeneCount[DiseaseGroup == "BDL_ASBTi"]
  prot_z <- Protein_Raw[DiseaseGroup == "BDL_ASBTi"]
  
  rna_d  <- GeneCount[DiseaseGroup == "Disease"]
  prot_d <- Protein_Raw[DiseaseGroup == "Disease"]
  
  cor_c   <- if (length(rna_c) >= 3) cor(rna_c, prot_c) else NA_real_
  cor_z   <- if (length(rna_z) >= 3) cor(rna_z, prot_z) else NA_real_
  cor_d   <- if (length(rna_d) >= 3) cor(rna_d, prot_d) else NA_real_
  cor_all <- suppressWarnings(cor(GeneCount, Protein_Raw))
  
  .(
    N = .N,
    SD_RNA = sd(GeneCount),
    SD_Prot = sd(Protein_Raw),
    Pearson_R = cor_all,
    Cor_Control = cor_c,
    Cor_ASBTi   = cor_z,
    Cor_Disease = cor_d,
    Schwerpunkt_RNA_C = mean(rna_c),
    Schwerpunkt_Prot_C = mean(prot_c),
    Schwerpunkt_RNA_Z = mean(rna_z),
    Schwerpunkt_Prot_Z = mean(prot_z),
    Schwerpunkt_RNA_D = mean(rna_d),
    Schwerpunkt_Prot_D = mean(prot_d),
    Ratio = sd(GeneCount) / sd(Protein_Raw)
  )
}, by = GeneProtein]

Cloud_Metrics_BDL <- Cloud_Metrics_BDL[
  N >= 10 & !is.na(Pearson_R) & SD_RNA > 0 & SD_Prot > 0
]

classify_cloud <- function(ratio, r, sd_rna, sd_prot) {
  if (is.na(ratio) | is.na(r)) return("Unclassified")
  if (ratio > 1.014  & abs(r) < 0.673) return("Horizontal")
  if (ratio < 0.528  & abs(r) < 0.673) return("Vertical")
  if (abs(r) > 0.467 & ratio >= 0.528 & ratio <= 1.014) return("Diagonal")
  if (ratio >= 0.528 & ratio <= 1.014 & abs(r) < 0.070) return("Round")
  return("Unclassified")
}

Cloud_Metrics_BDL[, CloudCategory := mapply(classify_cloud, Ratio, Pearson_R, SD_RNA, SD_Prot)]

Combined <- merge(dipa_stats, Cloud_Metrics_BDL, by = "GeneProtein")

cloud_colors <- c("Diagonal"="#2196F3", "Horizontal"="#FF9800", "Vertical"="#9C27B0",
                  "Round"="#4CAF50", "Unclassified"="#BDBDBD")

extreme_pts <- rbindlist(lapply(1:8, function(g) {
  sub <- Combined[DiPaGroups == as.character(g)]
  if(nrow(sub) == 0) return(NULL)
  idx <- switch(as.character(g),
                "8"=which.min(sqrt(sub$log2G^2 + sub$log2P^2)),
                "1"=which.max(sub$log2G + sub$log2P),
                "2"=which.min(sub$log2G + sub$log2P),
                "3"=which.max(sub$log2P),
                "4"=which.min(sub$log2P),
                "5"=which.max(sub$log2G),
                "6"=which.min(sub$log2G),
                "7"=which.max(abs(sub$log2G - sub$log2P)))
  sub[idx]
}), fill=TRUE)

extreme_pts[, PlotMarker := paste0("Pair ", DiPaGroups)]
extreme_pts[, LegendText := paste0("Pair ", DiPaGroups, ": ", GeneProtein)]

# --- Grand plot DiPa (EXACT ORIGINAL STRUCTURE & GREY SECTOR LABELS) ---
p3_dipa <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = cloud_colors) +
  geom_hline(yintercept = c(0.5, -0.5), color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(0.5, -0.5), color = "black", linewidth = 0.3) +
  xlim(-6, 6) + ylim(-5, 5) +
  labs(
    title = "",
    tag = "C",
    x = expression(RNA*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"),
    y = expression(Protein*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"),
    color = "", caption = "DiPa 1&2 → Diagonal | 5&6 → Horizontal | 3&4 → Vertical | 8 → Round | 7 → Discordant"
  ) +
  theme_bw(base_size = T_AXIS) +
  theme(
    plot.title = element_text(size = T_TITLE, hjust=0.5),
    axis.title = element_text(size = T_AXIS),
    axis.text = element_text(size = T_TICKS, color="black"),
    axis.ticks = element_line(color="black", linewidth=0.8),
    panel.border = element_rect(color="black", fill=NA, linewidth=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = T_TICKS),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(t=-10)
  ) +
  annotate("text", x =  0,    y =  0,    label = "8", size = 7, color = "grey50") +
  annotate("text", x =  4.0,  y =  3.5,  label = "1", size = 7, color = "grey50") +
  annotate("text", x = -4.0,  y = -3.5,  label = "2", size = 7, color = "grey50") +
  annotate("text", x =  0,    y =  3.5,  label = "3", size = 7, color = "grey50") +
  annotate("text", x =  0,    y = -3.5,  label = "4", size = 7, color = "grey50") +
  annotate("text", x =  4.0,  y =  0,    label = "5", size = 7, color = "grey50") +
  annotate("text", x = -4.0,  y =  0,    label = "6", size = 7, color = "grey50") +
  annotate("text", x = -4.0,  y =  3.5,  label = "7", size = 7, color = "grey50") +
  annotate("text", x =  4.0,  y = -3.5,  label = "7", size = 7, color = "grey50") +
  geom_point(data = extreme_pts, aes(x = log2G, y = log2P),
             shape = 18, size = 3, color = "black", inherit.aes=FALSE) +
  geom_text_repel(
    data = extreme_pts,
    aes(label = PlotMarker),
    size = 5.0,
    color="black",
    bg.color="white",
    bg.r=0.15,
    box.padding=0.1,
    max.overlaps=30,
    segment.color="grey40"
  )

# --- 8 Mini-clouds avec 3 Schwerpunkte (EXACT ORIGINAL BOX COLORING & LAYOUT) ---
clouds_list <- list()
clouds_list_iso <- list()

for(g in 1:8) {
  
  ep <- extreme_pts[DiPaGroups == as.character(g)]
  if(nrow(ep)==0) { 
    clouds_list[[as.character(g)]] <- ggplot()+theme_void()
    next 
  }
  
  gp   <- ep$GeneProtein[1]
  ccat <- ep$CloudCategory[1]
  p_dat <- bdl_DT[GeneProtein == gp]
  met   <- Cloud_Metrics_BDL[GeneProtein == gp]
  
  col <- cloud_colors[ccat]; if(is.na(col)) col <- "#BDBDBD"
  
  all_x <- c(p_dat$GeneCount, met$Schwerpunkt_RNA_C, met$Schwerpunkt_RNA_Z, met$Schwerpunkt_RNA_D)
  all_y <- c(p_dat$Protein_Raw, met$Schwerpunkt_Prot_C, met$Schwerpunkt_Prot_Z, met$Schwerpunkt_Prot_D)
  
  cx <- mean(range(all_x, na.rm=TRUE))
  cy <- mean(range(all_y, na.rm=TRUE))
  ps <- max(diff(range(all_x, na.rm=TRUE)), diff(range(all_y, na.rm=TRUE))) * 1.3
  
  df_sp <- data.frame(
    Group = c("Control", "BDL_ASBTi", "Disease"),
    RNA  = c(met$Schwerpunkt_RNA_C[1], met$Schwerpunkt_RNA_Z[1], met$Schwerpunkt_RNA_D[1]),
    Prot = c(met$Schwerpunkt_Prot_C[1], met$Schwerpunkt_Prot_Z[1], met$Schwerpunkt_Prot_D[1])
  )
  
  fit_all <- lm(Protein_Raw ~ GeneCount, data = p_dat)
  
  # Plot Master Grid
  cp <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=DiseaseGroup)) +
    geom_point(size=1.2, alpha=0.8) +
    scale_color_manual(values=group_colors) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    geom_abline(intercept = coef(fit_all)[1], slope = coef(fit_all)[2], color = "black", linewidth = 0.8) +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 2.5, stroke = 1.2, 
               color = c("Control" = "#0D47A1", "BDL_ASBTi" = "#1B5E20", "Disease" = "#B71C1C")[df_sp$Group], 
               inherit.aes = FALSE) +
    theme_bw(base_size=T_AXIS) +
    theme(
      legend.position="none",
      plot.title=element_text(color=col, size=T_TITLE - 4, hjust=0.5, face="bold"),
      plot.subtitle=element_text(size=8.5, color="black", hjust=0.5, lineheight=0.9),
      axis.title=element_text(size=T_TICKS - 2),
      axis.text=element_text(size=8, color="black"),
      axis.ticks=element_line(color="black", linewidth=0.8),
      panel.border=element_rect(color=col, fill=NA, linewidth=1.2),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    ) +
    labs(
      title = paste0("Pair ", g, ":\n", gp),
      subtitle = sprintf(
        "Cor=%.2f\nC=%.2f | ASBTi=%.2f | D=%.2f",
        met$Pearson_R[1], met$Cor_Control[1], met$Cor_ASBTi[1], met$Cor_Disease[1]
      ),
      x="RNA", y="Protein"
    )
  
  clouds_list[[as.character(g)]] <- cp
  
  # Plot Isolated
  cp_iso <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=DiseaseGroup)) +
    geom_point(size=2.2, alpha=0.8) +
    scale_color_manual(values=group_colors) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    geom_abline(intercept = coef(fit_all)[1], slope = coef(fit_all)[2], color = "black", linewidth = 1.0) +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 4.0, stroke = 1.5, 
               color = c("Control" = "#0D47A1", "BDL_ASBTi" = "#1B5E20", "Disease" = "#B71C1C")[df_sp$Group], 
               inherit.aes = FALSE) +
    theme_bw(base_size=14) +
    theme(
      legend.position="none",
      plot.title=element_text(color=col, size=18, hjust=0.5, face="bold"),
      plot.subtitle=element_text(size=13, color="black", hjust=0.5, lineheight=0.9),
      axis.title=element_text(size=14, face="bold"),
      axis.text=element_text(size=12, color="black"),
      axis.ticks=element_line(color="black", linewidth=0.8),
      panel.border=element_rect(color=col, fill=NA, linewidth=1.5),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    ) +
    labs(
      title = paste0("Pair ", g, ": ", gp),
      subtitle = substitute(
        paste(rho[BP] == r_all, " | ", rho[C] == r_ctrl, " | ", rho[Z] == r_asbt, " | ", rho[D] == r_dise),
        list(
          r_all  = sprintf("%.2f", met$Pearson_R[1]),
          r_ctrl = sprintf("%.2f", met$Cor_Control[1]),
          r_asbt = sprintf("%.2f", met$Cor_ASBTi[1]),
          r_dise = sprintf("%.2f", met$Cor_Disease[1])
        )
      ),
      x="RNA", y="Protein"
    )
  clouds_list_iso[[as.character(g)]] <- cp_iso
}

# -------------------------------------------------------------------------
# 5. ASSEMBLE ENTIRE MASTER PLATE (EXACT ORIGINAL DESIGN & HEADERS)
# -------------------------------------------------------------------------
design <- "
  AAABBBCCCDDD
  EEEEEEFFFGGG
  EEEEEEHHHIII
  EEEEEEJJJKKK
  EEEEEEMMMNNN
"

layout_master <- p1_pca_rna + p1_pca_prot + p2_volc_rna + p2_volc_prot +
  p3_dipa + 
  clouds_list[["1"]] + clouds_list[["2"]] + clouds_list[["3"]] + 
  clouds_list[["4"]] + clouds_list[["5"]] + clouds_list[["6"]] + 
  clouds_list[["7"]] + clouds_list[["8"]] +
  plot_layout(design = design) +
  plot_annotation(
    title = "BDL DATASET (18 Mice with BDL ASBTi)",
    subtitle = paste0(
      "Using 943 overlapping genes between BDL and CCl4 datasets\n",
      "** Points:  Blue = Control (Sham) | Green = BDL ASBTi | Red = Disease (BDL) **\n",
      "** CROSSES (+) = SCHWERPUNKTE (Group Centroids) in Dark Blue / Green / Red **\n",
      paste(extreme_pts$LegendText[1:4], collapse = " | "), "\n", 
      paste(extreme_pts$LegendText[5:8], collapse = " | ")
    ),
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), 
                  plot.subtitle = element_text(size = 14, color = "black", hjust = 0.5, margin = margin(b = 10)))
  )

master_pdf <- file.path(out_dir, "A4_MASTER_BDL_Combined_Schwerpunkt_18Mice_normal_data.pdf")
ggsave(master_pdf, plot = layout_master, width = 11.7, height = 14.0, device = cairo_pdf)
print("SUCCESS! File A4_MASTER_BDL_Combined_Schwerpunkt_18Mice_normal_data.pdf was saved!")

# -------------------------------------------------------------------------
# 6. SAVE ISOLATED PLOTS (3 SEPARATE PDF FILES)
# -------------------------------------------------------------------------
print("Saving isolated plots for 18 BDL mice...")

# Isolated PCA
layout_pca <- ((p1_pca_rna + labs(title = "RNA", tag = NULL)) | (p1_pca_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16)
  )
ggsave(file.path(out_dir, "Isolated_PCA_BDL_18Mice_normal_data.pdf"), plot = layout_pca, width = 10, height = 5, device = cairo_pdf)

# Isolated Volcano
layout_volcano <- ((p2_volc_rna + labs(title = "RNA", tag = NULL)) | (p2_volc_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black")
  )
ggsave(file.path(out_dir, "Isolated_Volcano_BDL_18Mice_normal_data.pdf"), plot = layout_volcano, width = 10, height = 5, device = cairo_pdf)

# Isolated DiPa + Clouds
p3_dipa_iso <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.5, size = 1.8) + scale_color_manual(values = cloud_colors) +
  geom_hline(yintercept = c(0.5, -0.5), color = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(0.5, -0.5), color = "black", linewidth = 0.4) +
  xlim(-6, 6) + ylim(-5, 5) +
  labs(
    title = "DiPa Fold Changes",
    x = expression(RNA*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"),
    y = expression(Protein*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"),
    color = "Cloud Type:"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 28, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "bottom", legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18), legend.key.size = unit(0.8, "cm")
  ) +
  annotate("text", x =  0,    y =  0,    label = "8", size = 10, color = "grey50") +
  annotate("text", x =  4.0,  y =  3.5,  label = "1", size = 10, color = "grey50") +
  annotate("text", x = -4.0,  y = -3.5,  label = "2", size = 10, color = "grey50") +
  annotate("text", x =  0,    y =  3.5,  label = "3", size = 10, color = "grey50") +
  annotate("text", x =  0,    y = -3.5,  label = "4", size = 10, color = "grey50") +
  annotate("text", x =  4.0,  y =  0,    label = "5", size = 10, color = "grey50") +
  annotate("text", x = -4.0,  y =  0,    label = "6", size = 10, color = "grey50") +
  annotate("text", x = -4.0,  y =  3.5,  label = "7", size = 10, color = "grey50") +
  annotate("text", x =  4.0,  y = -3.5,  label = "7", size = 10, color = "grey50") +
  geom_point(data = extreme_pts, aes(x = log2G, y = log2P), shape = 18, size = 5, color = "black", inherit.aes = FALSE) +
  geom_text_repel(
    data = extreme_pts,
    aes(label = PlotMarker),
    size = 7.0,
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.2,
    max.overlaps = 30,
    segment.color = "grey40"
  )

design_isolated_dipa <- "
  EEEEEEFFFGGG
  EEEEEEHHHIII
  EEEEEEJJJKKK
  EEEEEEMMMNNN
"
layout_dipa_isolated <- (p3_dipa_iso + labs(tag="A") + theme(plot.tag = element_text(face = "bold", size = 32))) + 
  (clouds_list_iso[["1"]] + labs(tag="B") + theme(plot.tag = element_text(face = "bold", size = 32))) + 
  clouds_list_iso[["2"]] + clouds_list_iso[["3"]] + 
  clouds_list_iso[["4"]] + clouds_list_iso[["5"]] + clouds_list_iso[["6"]] + 
  clouds_list_iso[["7"]] + clouds_list_iso[["8"]] +
  plot_layout(design = design_isolated_dipa) 

ggsave(file.path(out_dir, "Isolated_DiPa_Wolken_BDL_18Mice_normal_data.pdf"), plot = layout_dipa_isolated, width = 19, height = 17.5, device = cairo_pdf)

print("SUCCESS! Isolated plots for 18 BDL mice were saved in PDF only.")
