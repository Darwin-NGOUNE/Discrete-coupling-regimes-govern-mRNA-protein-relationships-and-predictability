# ==============================================================================
# SCRIPT: Generate_A4_MASTER_CCL4_Schwerpunkt_ALL_36Mice.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/
# OBJECTIVE: Generate a unified 13-Panel A4 Master Plate combining
#            PCA (Step 1), Volcano (Step 2), and DiPa/Clouds (Step 3) for CCL4
#            INTEGRATING ALL 36 MICE WITH 6 HIGHLY DISTINCT CONTRASTED COLORS.
#            (Blue = Control, Red = Disease, and 4 completely distinct intermediate colors).
#            PDF-ONLY outputs (no PNG).
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(ggforce)
library(ggrepel)

# Text sizes strictly matching the standard master format
T_TITLE = 20
T_SUB   = 14
T_AXIS  = 16
T_TICKS = 14

out_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data"

# -------------------------------------------------------------------------
# 1. LOAD DATA GLOBALLY (ORIGINAL RAW DATASET FOR ALL 36 CCL4 MICE)
# -------------------------------------------------------------------------
print("Loading Global Raw Data for CCL4 (All 36 Mice with Zwischenmäuse)...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
Full_DT <- as.data.table(DTccl4_DT_LCPM)

Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]

# Map all 6 Treatment Groups for CCl4
ccl4_DT <- Full_DT[Dataset == "CCl4 Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw)]

ccl4_DT[, TreatmentGroup := fcase(
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 0,  "Oil M0",
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 2,  "Oil M2",
  Treatment == "oil"  & as.numeric(as.character(TreatmentTime)) == 12, "Oil M12",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 2,  "CCl4 M2",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 6,  "CCl4 M6",
  Treatment == "ccl4" & as.numeric(as.character(TreatmentTime)) == 12, "CCl4 M12"
)]

ccl4_DT <- ccl4_DT[!is.na(TreatmentGroup)]
ccl4_DT$TreatmentGroup <- factor(ccl4_DT$TreatmentGroup, 
                                 levels = c("Oil M0", "Oil M2", "Oil M12", "CCl4 M2", "CCl4 M6", "CCl4 M12"))

cat(sprintf("Loaded %d rows across %d mice for CCl4 cohort (6 conditions x 6 reps = 36 mice).\n", 
            nrow(ccl4_DT), uniqueN(ccl4_DT$MiceInfo)))

# 6 HIGHLY DISTINCT COLORS (Blue=Control, Red=Disease, Cyan, Purple, Orange, Green)
group_colors <- c(
  "Oil M0"   = "#1E88E5", # Pure Bright Blue (Baseline Control)
  "Oil M2"   = "#00ACC1", # Cyan / Teal (Oil Month 2 Zwischenmäuse)
  "Oil M12"  = "#8E24AA", # Vivid Purple (Oil Month 12 Zwischenmäuse)
  "CCl4 M2"  = "#FF9800", # Vivid Orange (CCl4 Month 2 Zwischenmäuse)
  "CCl4 M6"  = "#43A047", # Vivid Green (CCl4 Month 6 Zwischenmäuse)
  "CCl4 M12" = "#E53935"  # Pure Vivid Red (CCl4 Month 12 Endpoint Disease)
)

group_labels <- c(
  "Oil M0"   = "Oil M0 (Control)",
  "Oil M2"   = "Oil M2",
  "Oil M12"  = "Oil M12",
  "CCl4 M2"  = "CCl4 M2",
  "CCl4 M6"  = "CCl4 M6",
  "CCl4 M12" = "CCl4 M12 (Disease)"
)

# Compute DiPa Ratios (Month 12 CCl4 vs Month 0 Oil)
dipa_stats_ccl4 <- ccl4_DT[, .(
  meanG_C = mean(GeneCount[TreatmentGroup == "Oil M0"], na.rm = TRUE),
  meanG_D = mean(GeneCount[TreatmentGroup == "CCl4 M12"], na.rm = TRUE),
  meanP_C = mean(Protein_Raw[TreatmentGroup == "Oil M0"], na.rm = TRUE),
  meanP_D = mean(Protein_Raw[TreatmentGroup == "CCl4 M12"], na.rm = TRUE),
  DiPaGroups = as.character(ClusterDiPa[1])
), by = GeneProtein]

dipa_stats_ccl4[DiPaGroups == "0", DiPaGroups := "8"]
dipa_stats_ccl4[, log2G := log2(meanG_D / meanG_C)]
dipa_stats_ccl4[, log2P := log2(meanP_D / meanP_C)]
DT_dipa_valid <- dipa_stats_ccl4[is.finite(log2G) & is.finite(log2P) & !is.na(DiPaGroups)]

# -------------------------------------------------------------------------
# 2. STEP 1: PCA PLOTS (ALL 36 MICE)
# -------------------------------------------------------------------------
generate_master_pca <- function(dt, value_col, title) {
  wide_data <- dcast(dt, GeneProtein ~ MiceInfo, value.var = value_col, fun.aggregate = mean)
  mat <- as.matrix(wide_data[, -1, with = FALSE])
  rownames(mat) <- wide_data$GeneProtein
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 3) return(ggplot() + theme_void() + annotate("text", x=0, y=0, label="Omitted\nSparsity"))
  mat <- mat[apply(mat, 1, var, na.rm=TRUE) > 0, ]
  pca_res <- prcomp(t(mat), scale. = TRUE)
  var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  pca_df <- data.table(MiceInfo = rownames(pca_res$x), PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2])
  meta <- unique(dt[, .(MiceInfo, TreatmentGroup)])
  pca_df <- merge(pca_df, meta, by = "MiceInfo")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = TreatmentGroup)) +
    geom_point(size = 2.4, alpha = 0.9) +
    scale_color_manual(values = group_colors, labels = group_labels, name = "") +
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

p1_pca_rna  <- generate_master_pca(ccl4_DT, "GeneCount", "RNA") + 
  labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 20), legend.position = "none")
p1_pca_prot <- generate_master_pca(ccl4_DT, "Protein_Raw", "Protein") + 
  theme(legend.position = "none")

# -------------------------------------------------------------------------
# 3. STEP 2: VOLCANO PLOTS (CCl4 Month 12 vs Oil Month 0)
# -------------------------------------------------------------------------
calculate_volcano <- function(dt, value_col) {
  dt_sub <- dt[!is.na(get(value_col)) & TreatmentGroup %in% c("Oil M0", "CCl4 M12")]
  stats <- dt_sub[, {
    v_c <- get(value_col)[TreatmentGroup == "Oil M0"];   v_c <- v_c[!is.na(v_c)]
    v_d <- get(value_col)[TreatmentGroup == "CCl4 M12"]; v_d <- v_d[!is.na(v_d)]
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
    scale_color_manual(values = c("Down" = "#E53935", "Not Sig" = "black", "Up" = "#E53935")) +
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

p2_volc_rna  <- plot_master_volcano(calculate_volcano(ccl4_DT, "GeneCount"), "RNA") + 
  labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 20), legend.position = "none")
p2_volc_prot <- plot_master_volcano(calculate_volcano(ccl4_DT, "Protein_Raw"), "Protein") + 
  theme(legend.position = "none")

# -------------------------------------------------------------------------
# 4. STEP 3: DIPA & CLOUDS (ACROSS ALL 36 MICE WITH 6 SCHWERPUNKTE)
# -------------------------------------------------------------------------
Cloud_Metrics_CCL4 <- ccl4_DT[, {
  rna_o0  <- GeneCount[TreatmentGroup == "Oil M0"]
  prot_o0 <- Protein_Raw[TreatmentGroup == "Oil M0"]
  
  rna_o2  <- GeneCount[TreatmentGroup == "Oil M2"]
  prot_o2 <- Protein_Raw[TreatmentGroup == "Oil M2"]
  
  rna_o12  <- GeneCount[TreatmentGroup == "Oil M12"]
  prot_o12 <- Protein_Raw[TreatmentGroup == "Oil M12"]
  
  rna_c2  <- GeneCount[TreatmentGroup == "CCl4 M2"]
  prot_c2 <- Protein_Raw[TreatmentGroup == "CCl4 M2"]
  
  rna_c6  <- GeneCount[TreatmentGroup == "CCl4 M6"]
  prot_c6 <- Protein_Raw[TreatmentGroup == "CCl4 M6"]
  
  rna_c12  <- GeneCount[TreatmentGroup == "CCl4 M12"]
  prot_c12 <- Protein_Raw[TreatmentGroup == "CCl4 M12"]
  
  cor_all <- suppressWarnings(cor(GeneCount, Protein_Raw, use = "complete.obs"))
  cor_ctrl <- if (length(rna_o0) >= 3) cor(rna_o0, prot_o0) else NA_real_
  cor_dise <- if (length(rna_c12) >= 3) cor(rna_c12, prot_c12) else NA_real_
  
  .(
    N = .N,
    SD_RNA  = sd(GeneCount, na.rm = TRUE),
    SD_Prot = sd(Protein_Raw, na.rm = TRUE),
    Pearson_R = cor_all,
    Cor_Control = cor_ctrl,
    Cor_Disease = cor_dise,
    
    # 6 Schwerpunkte
    SP_RNA_Oil0   = mean(rna_o0, na.rm = TRUE),
    SP_Prot_Oil0  = mean(prot_o0, na.rm = TRUE),
    SP_RNA_Oil2   = mean(rna_o2, na.rm = TRUE),
    SP_Prot_Oil2  = mean(prot_o2, na.rm = TRUE),
    SP_RNA_Oil12  = mean(rna_o12, na.rm = TRUE),
    SP_Prot_Oil12 = mean(prot_o12, na.rm = TRUE),
    SP_RNA_CCl4_2 = mean(rna_c2, na.rm = TRUE),
    SP_Prot_CCl4_2= mean(prot_c2, na.rm = TRUE),
    SP_RNA_CCl4_6 = mean(rna_c6, na.rm = TRUE),
    SP_Prot_CCl4_6= mean(prot_c6, na.rm = TRUE),
    SP_RNA_CCl4_12= mean(rna_c12, na.rm = TRUE),
    SP_Prot_CCl4_12=mean(prot_c12, na.rm = TRUE)
  )
}, by = GeneProtein]

Cloud_Metrics_CCL4 <- Cloud_Metrics_CCL4[
  N >= 20 & !is.na(Pearson_R) & !is.na(SD_RNA) & !is.na(SD_Prot) & SD_RNA > 0 & SD_Prot > 0
]

Cloud_Metrics_CCL4[, Ratio := SD_RNA / SD_Prot]

classify_cloud <- function(ratio, r, sd_rna, sd_prot) {
  if (is.na(ratio) | is.na(r)) return("Unclassified")
  if (ratio > 1.014  & abs(r) < 0.673) return("Horizontal")
  if (ratio < 0.528  & abs(r) < 0.673) return("Vertical")
  if (abs(r) > 0.467 & ratio >= 0.528 & ratio <= 1.014) return("Diagonal")
  if (ratio >= 0.528 & ratio <= 1.014 & abs(r) < 0.070) return("Round")
  return("Unclassified")
}

Cloud_Metrics_CCL4[, CloudCategory := mapply(classify_cloud, Ratio, Pearson_R, SD_RNA, SD_Prot)]

Combined <- merge(
  DT_dipa_valid[, .(GeneProtein, log2G, log2P, DiPaGroups)],
  Cloud_Metrics_CCL4,
  by = "GeneProtein"
)

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

p3_dipa <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.5, size = 1) + scale_color_manual(values = cloud_colors) +
  geom_hline(yintercept = 0.5, color = "black", linewidth = 0.3) + geom_hline(yintercept = -0.5, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, color = "black", linewidth = 0.3) + geom_vline(xintercept = -0.5, color = "black", linewidth = 0.3) +
  xlim(-6, 6) + ylim(-5, 5) +
  labs(
    title = "", 
    tag = "C", 
    x = expression(RNA*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), 
    y = expression(Protein*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), 
    color = "", 
    caption = "DiPa 1&2 → Diagonal | 5&6 → Horizontal | 3&4 → Vertical | 8 → Round | 7 → Discordant"
  ) + 
  theme_bw(base_size = T_AXIS) +
  theme(
    plot.title = element_text(size = T_TITLE, hjust=0.5), axis.title = element_text(size = T_AXIS), axis.text = element_text(size = T_TICKS, color="black"),
    axis.ticks = element_line(color="black", linewidth=0.8), panel.border = element_rect(color="black", fill=NA, linewidth=1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 20),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = T_TICKS), legend.key.size = unit(0.3, "cm"), legend.margin=margin(t=-10)
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
  geom_point(data = extreme_pts, aes(x = log2G, y = log2P), shape = 18, size = 3, color = "black", inherit.aes=FALSE) +
  geom_text_repel(data = extreme_pts, aes(label = PlotMarker), size = 5.0, color="black", bg.color="white", bg.r=0.15, box.padding=0.1, max.overlaps=30, segment.color="grey40")

# --- 8 Mini-clouds CCL4 (36 Mice with 6 Schwerpunkte) ---
clouds_list <- list()
clouds_list_iso <- list()

for(g in 1:8) {
  
  ep <- extreme_pts[DiPaGroups == as.character(g)]
  if(nrow(ep)==0) { 
    clouds_list[[as.character(g)]] <- ggplot()+theme_void()
    clouds_list_iso[[as.character(g)]] <- ggplot()+theme_void()
    next 
  }
  
  gp <- ep$GeneProtein[1]
  ccat <- ep$CloudCategory[1]
  p_dat <- ccl4_DT[GeneProtein == gp]
  met <- Cloud_Metrics_CCL4[GeneProtein == gp]
  
  col <- cloud_colors[ccat]; if(is.na(col)) col <- "#BDBDBD"
  
  all_x <- p_dat$GeneCount
  all_y <- p_dat$Protein_Raw
  
  cx <- mean(range(all_x, na.rm=TRUE))
  cy <- mean(range(all_y, na.rm=TRUE))
  ps <- max(diff(range(all_x, na.rm=TRUE)), diff(range(all_y, na.rm=TRUE))) * 1.3
  
  # 6 Schwerpunkte with darker high-contrast tones
  df_sp <- data.frame(
    Group = c("Oil M0", "Oil M2", "Oil M12", "CCl4 M2", "CCl4 M6", "CCl4 M12"),
    RNA = c(met$SP_RNA_Oil0[1], met$SP_RNA_Oil2[1], met$SP_RNA_Oil12[1],
            met$SP_RNA_CCl4_2[1], met$SP_RNA_CCl4_6[1], met$SP_RNA_CCl4_12[1]),
    Prot = c(met$SP_Prot_Oil0[1], met$SP_Prot_Oil2[1], met$SP_Prot_Oil12[1],
             met$SP_Prot_CCl4_2[1], met$SP_Prot_CCl4_6[1], met$SP_Prot_CCl4_12[1])
  )
  
  sp_colors <- c(
    "Oil M0"   = "#0D47A1", # Dark Blue
    "Oil M2"   = "#006064", # Dark Teal
    "Oil M12"  = "#4A148C", # Dark Purple
    "CCl4 M2"  = "#E65100", # Dark Orange
    "CCl4 M6"  = "#1B5E20", # Dark Green
    "CCl4 M12" = "#B71C1C"  # Dark Red
  )
  
  fit_all <- lm(Protein_Raw ~ GeneCount, data = p_dat)
  
  # Plot Master Grid
  cp <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=TreatmentGroup)) +
    geom_point(size=1.3, alpha=0.85) +
    scale_color_manual(values=group_colors) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    geom_abline(intercept = coef(fit_all)[1], slope = coef(fit_all)[2], color = "black", linewidth = 0.8) +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 2.5, stroke = 1.3, 
               color = sp_colors[df_sp$Group], 
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
      panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    ) +
    labs(
      title = paste0("Pair ", g, ":\n", gp),
      subtitle = sprintf(
        "Cor=%.2f\nCtrl(M0)=%.2f | Dis(M12)=%.2f",
        met$Pearson_R[1], met$Cor_Control[1], met$Cor_Disease[1]
      ),
      x="RNA", y="Protein"
    )
  
  clouds_list[[as.character(g)]] <- cp
  
  # Plot Isolated
  cp_iso <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=TreatmentGroup)) +
    geom_point(size=2.4, alpha=0.85) +
    scale_color_manual(values=group_colors) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    geom_abline(intercept = coef(fit_all)[1], slope = coef(fit_all)[2], color = "black", linewidth = 1.0) +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 4.2, stroke = 1.6, 
               color = sp_colors[df_sp$Group], 
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
        paste(rho[BP] == r_all, " | ", rho[C] == r_ctrl, " | ", rho[D] == r_dise),
        list(
          r_all  = sprintf("%.2f", met$Pearson_R[1]),
          r_ctrl = sprintf("%.2f", met$Cor_Control[1]),
          r_dise = sprintf("%.2f", met$Cor_Disease[1])
        )
      ),
      x="RNA", y="Protein"
    )
  
  clouds_list_iso[[as.character(g)]] <- cp_iso
}

# -------------------------------------------------------------------------
# 5. ASSEMBLE ENTIRE MASTER PLATE
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
    title = "CCL4 DATASET (Full Cohort: 36 Mice with Zwischenmäuse)",
    subtitle = paste0(
      "Using 943 overlapping genes between BDL and CCl4 datasets\n",
      "** Points: Oil M0 (Blue), Oil M2 (Cyan), Oil M12 (Purple) | CCl4 M2 (Orange), CCl4 M6 (Green), CCl4 M12 (Red) **\n",
      "** CROSSES (+) = SCHWERPUNKTE (Group Centroids) for the 6 Experimental Conditions **\n",
      paste(extreme_pts$LegendText[1:4], collapse = " | "), "\n", 
      paste(extreme_pts$LegendText[5:8], collapse = " | ")
    ),
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5), 
      plot.subtitle = element_text(size = 13, color = "black", hjust = 0.5, margin = margin(b = 10))
    )
  )

master_pdf <- file.path(out_dir, "A4_MASTER_CCL4_Combined_Schwerpunkt_36Mice.pdf")
ggsave(master_pdf, plot = layout_master, width = 11.7, height = 14.0, device = cairo_pdf)
print("SUCCESS! File A4_MASTER_CCL4_Combined_Schwerpunkt_36Mice.pdf was saved!")

# -------------------------------------------------------------------------
# 6. SAVE ISOLATED PLOTS (3 SEPARATE PDF FILES)
# -------------------------------------------------------------------------
print("Saving isolated plots for CCL4 (All 36 Mice)...")

# Isolated PCA
layout_pca <- ((p1_pca_rna + labs(title = "RNA", tag = NULL)) | (p1_pca_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16)
  )
ggsave(file.path(out_dir, "Isolated_PCA_CCL4_36Mice.pdf"), plot = layout_pca, width = 10, height = 5, device = cairo_pdf)

# Isolated Volcano
layout_volcano <- ((p2_volc_rna + labs(title = "RNA", tag = NULL)) | (p2_volc_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black")
  )
ggsave(file.path(out_dir, "Isolated_Volcano_CCL4_36Mice.pdf"), plot = layout_volcano, width = 10, height = 5, device = cairo_pdf)

# Isolated DiPa + Clouds
p3_dipa_iso <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.5, size = 1.8) + scale_color_manual(values = cloud_colors) +
  geom_hline(yintercept = 0.5, color = "black", linewidth = 0.4) + geom_hline(yintercept = -0.5, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0.5, color = "black", linewidth = 0.4) + geom_vline(xintercept = -0.5, color = "black", linewidth = 0.4) +
  xlim(-6, 6) + ylim(-5, 5) +
  labs(
    title = "DiPa Fold Changes",
    x = expression(RNA*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"),
    y = expression(Protein*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"),
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

ggsave(file.path(out_dir, "Isolated_DiPa_Wolken_CCL4_36Mice.pdf"), plot = layout_dipa_isolated, width = 19, height = 17.5, device = cairo_pdf)

print("SUCCESS! Isolated plots for CCL4 (All 36 Mice) were saved in PDF only.")
