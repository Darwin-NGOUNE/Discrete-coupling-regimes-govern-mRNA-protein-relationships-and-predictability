# ==============================================================================
# SCRIPT: Generate_A4_MASTER_CCL4.R
# OBJECTIVE: Generate a strictly scaled, unified 13-Panel A4 Master Plate combining
#            PCA (Step 1), Volcano (Step 2), and DiPa/Clouds (Step 3).
# AESTHETICS: Derived strictly from the finalized BDL Master Script.
# ==============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(ggforce)
library(ggrepel)

# Set base text size suitable for 13 miniaturized plots on one A4 (Increased for readability)
T_TITLE = 20
T_SUB   = 14
T_AXIS  = 16
T_TICKS = 14

# -------------------------------------------------------------------------
# 1. LOAD DATA GLOBALLY
# -------------------------------------------------------------------------
print("Loading Global Data for PCA and Volcano...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
Full_DT <- DTccl4_DT_LCPM
setDT(Full_DT)

# Format Global DT exactly as we did for BDL
Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]
Full_DT[, DiseaseGroup := NA_character_]

# Disease Group Logic for CCl4
Full_DT[Dataset == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) == 0, DiseaseGroup := "Control"]
#Full_DT[Dataset == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) > 0 & as.character(Treatment) == "ccl4", DiseaseGroup := "Disease"]

Full_DT[Dataset == "CCl4 Dataset" &  as.numeric(as.character(TreatmentTime)) == 12 & as.character(Treatment) == "ccl4", DiseaseGroup := "Disease"]


# Isolate clean CCl4 DT
ccl4_DT <- Full_DT[Dataset == "CCl4 Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]

# The main structure used for PCA and Volcano
DT_CCL4_BatchCorrected <- ccl4_DT[DiseaseGroup %in% c("Control", "Disease")]

print("Loading DiPa Base Data for CCl4...")
env_dipa <- new.env()
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_CCL4_filtered_Gene_Protein_full_final_with_dipa.RData", envir = env_dipa)
# Dynamically extract the single data.frame/data.table loaded
dipa_var <- ls(env_dipa)[sapply(ls(env_dipa), function(x) inherits(env_dipa[[x]], "data.frame"))][1]
DT_dipa_valid <- copy(env_dipa[[dipa_var]])
setDT(DT_dipa_valid)

DT_dipa_valid[, log2G := log2(meanRatioG_12)]
DT_dipa_valid[, log2P := log2(meanRatioP_12)]
DT_dipa_valid <- DT_dipa_valid[is.finite(log2G) & is.finite(log2P) & !is.na(DiPaGroups)]
DT_dipa_valid[, DiPaGroups := as.character(DiPaGroups)]
DT_dipa_valid[DiPaGroups == "0", DiPaGroups := "8"]

# -------------------------------------------------------------------------
# 2. STEP 1: PCA PLOTS
# -------------------------------------------------------------------------
generate_master_pca <- function(dt, value_col, title, ctrl_label, dis_label) {
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
    scale_color_manual(values = c("Control" = "#2196F3", "Disease" = "#F44336"), labels = c("Control" = ctrl_label, "Disease" = dis_label), name = "") +
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
p1_pca_rna  <- generate_master_pca(DT_CCL4_BatchCorrected, "GeneCount", "RNA", "Control (Oil/Baseline)", "Disease (CCl4)") + labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 20)) + theme(legend.position = "none")
p1_pca_prot <- generate_master_pca(DT_CCL4_BatchCorrected, "Protein_Raw", "Protein", "Control (Oil/Baseline)", "Disease (CCl4)") + theme(legend.position = "none")

# -------------------------------------------------------------------------
# 3. STEP 2: VOLCANO PLOTS
# -------------------------------------------------------------------------
calculate_combat_volcano <- function(dt, value_col) {
  dt <- dt[!is.na(get(value_col))]
  stats <- dt[, {
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
# MODIFIED: Place Down/Up labels far-left and far-right above the grid, and capitalize axis labels
plot_master_volcano <- function(stats_dt, title_str, text_size = 4.5) {
  n_up   <- sum(stats_dt$Status == "Up")
  n_down <- sum(stats_dt$Status == "Down")
  
  # Calculate dynamic text placement coordinates above the grid
  x_max_val <- max(abs(stats_dt$LogFC_Val), na.rm = TRUE)
  y_max_val <- max(-log10(stats_dt$padj_Val), na.rm = TRUE)
  
  # Far left and far right (0.95 of max width) and just above the top border (1.03 of max height)
  x_text_pos <- x_max_val * 0.95
  y_text_pos <- y_max_val * 1.03
  
  ggplot(stats_dt, aes(x = LogFC_Val, y = -log10(padj_Val), color = Status)) +
    geom_point(aes(size = Status, alpha = Status)) +
    scale_size_manual(values = c("Not Sig" = 0.4, "Down" = 1.2, "Up" = 1.2)) +
    scale_alpha_manual(values = c("Not Sig" = 0.4, "Down" = 0.8, "Up" = 0.8)) +
    scale_color_manual(values = c("Down" = "red", "Not Sig" = "black", "Up" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.6) +
    
    # Far left (hjust = 0 at -x_text_pos) and Far right (hjust = 1 at x_text_pos) above the grid
    annotate("text", x = -x_text_pos, y = y_text_pos, 
             label = paste0("Down: ", format(n_down, big.mark=",")), 
             color = "black", fontface = "bold", size = text_size, hjust = 0) +
    annotate("text", x = x_text_pos, y = y_text_pos, 
             label = paste0("Up: ", format(n_up, big.mark=",")), 
             color = "black", fontface = "bold", size = text_size, hjust = 1) +
             
    # MODIFIED: Capitalized Log2 and -Log10 (P value)
    labs(
      title = title_str, 
      subtitle = NULL, 
      x = expression(Log[2]~"(fold change)"), 
      y = expression(-Log[10]~"(P value)")
    ) +
    coord_cartesian(clip = "off") + # Allows annotations outside the grid
    theme_bw(base_size = T_AXIS) +
    theme(
      plot.title = element_text(size = T_TITLE, hjust = 0.5, margin = margin(b = 15)), 
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10), # Extra top margin for labels
      axis.title = element_text(size = T_AXIS), axis.text = element_text(size = T_TICKS, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none"
    )
}
p2_volc_rna  <- plot_master_volcano(calculate_combat_volcano(DT_CCL4_BatchCorrected, "GeneCount"), "RNA") + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 20)) + theme(legend.position = "none")
p2_volc_prot <- plot_master_volcano(calculate_combat_volcano(DT_CCL4_BatchCorrected, "Protein_Raw"), "Protein") + theme(legend.position = "none")

# -------------------------------------------------------------------------
# 4. STEP 3: DIPA & CLOUDS
# -------------------------------------------------------------------------

# ---- 1. Calcul des métriques Cloud ----
Cloud_Metrics_CCL4 <- ccl4_DT[, {
  
  # Modèle global
  mod <- tryCatch(lm(Protein_Raw ~ GeneCount), error = function(e) NULL)
  
  # Modèle inverse
  mod_inv <- tryCatch(lm(GeneCount ~ Protein_Raw), error = function(e) NULL)
  
  # Sous-ensembles
  rna_c <- GeneCount[DiseaseGroup == "Control"]
  prot_c <- Protein_Raw[DiseaseGroup == "Control"]
  
  rna_d  <- GeneCount[DiseaseGroup == "Disease"]
  prot_d <- Protein_Raw[DiseaseGroup == "Disease"]
  
  # --- Corrélations par groupe ---
  cor_c <- if (length(rna_c) >= 3) cor(rna_c, prot_c) else NA_real_
  cor_d <- if (length(rna_d) >= 3) cor(rna_d, prot_d) else NA_real_
  
  # --- Schwerpunkt Control ---
  mean_rna_c  <- mean(rna_c)
  mean_prot_c <- mean(prot_c)
  
  # --- Schwerpunkt Disease ---
  mean_rna_d  <- mean(rna_d)
  mean_prot_d <- mean(prot_d)
  
  .(
    N = .N,
    SD_RNA = sd(GeneCount),
    SD_Prot = sd(Protein_Raw),
    Pearson_R = suppressWarnings(cor(GeneCount, Protein_Raw)),
    Cor_Control = cor_c,
    Cor_Disease = cor_d,
    Schwerpunkt_RNA_C = mean_rna_c,
    Schwerpunkt_Prot_C = mean_prot_c,
    Schwerpunkt_RNA_D = mean_rna_d,
    Schwerpunkt_Prot_D = mean_prot_d,
    Schwerpunkt_Prot_D = mean_prot_d
  )
}, by = GeneProtein]

Cloud_Metrics_CCL4 <- Cloud_Metrics_CCL4[
  N >= 10 & !is.na(Pearson_R) & !is.na(SD_RNA) & !is.na(SD_Prot) & SD_RNA > 0 & SD_Prot > 0
]

Cloud_Metrics_CCL4[, Ratio := SD_RNA / SD_Prot]

# ---- 2. Classification des clouds ----
classify_cloud <- function(ratio, r, sd_rna, sd_prot) {
  if (is.na(ratio) | is.na(r)) return("Unclassified")
  if (ratio > 1.014  & abs(r) < 0.673) return("Horizontal")
  if (ratio < 0.528  & abs(r) < 0.673) return("Vertical")
  if (abs(r) > 0.467 & ratio >= 0.528 & ratio <= 1.014) return("Diagonal")
  if (ratio >= 0.528 & ratio <= 1.014 & abs(r) < 0.070) return("Round")
  return("Unclassified")
}

Cloud_Metrics_CCL4[, CloudCategory := mapply(classify_cloud, Ratio, Pearson_R, SD_RNA, SD_Prot)]

# ---- 3. Fusion avec DiPa ----
Combined <- merge(
  DT_dipa_valid[, .(GeneProtein, log2G, log2P, DiPaGroups, meanRatioG, meanRatioP)],
  Cloud_Metrics_CCL4[, .(GeneProtein, CloudCategory, Pearson_R, Cor_Control, Cor_Disease,
                         Schwerpunkt_RNA_C, Schwerpunkt_Prot_C, Schwerpunkt_RNA_D, Schwerpunkt_Prot_D,
                         Ratio, SD_RNA, SD_Prot)],
  by = "GeneProtein"
)

# ---- 4. Points extrêmes ----
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
  labs(title = "", tag = "C", x = expression(RNA*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), y = expression(Protein*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), color = "", caption = "DiPa 1&2 → Diagonal | 5&6 → Horizontal | 3&4 → Vertical | 8 → Round | 7 → Discordant") + theme_bw(base_size = T_AXIS) +
  theme(plot.title = element_text(size = T_TITLE, hjust=0.5), axis.title = element_text(size = T_AXIS), axis.text = element_text(size = T_TICKS, color="black"),
        axis.ticks = element_line(color="black", linewidth=0.8), panel.border = element_rect(color="black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 20),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = T_TICKS), legend.key.size = unit(0.3, "cm"), legend.margin=margin(t=-10)) +
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


# ---- 5. Mini-clouds ----
clouds_list <- list()

for(g in 1:8) {
  
  ep <- extreme_pts[DiPaGroups == as.character(g)]
  if(nrow(ep)==0) { clouds_list[[as.character(g)]] <- ggplot()+theme_void(); next }
  
  gp <- ep$GeneProtein[1]
  ccat <- ep$CloudCategory[1]
  p_dat <- ccl4_DT[GeneProtein == gp]
  met <- Cloud_Metrics_CCL4[GeneProtein == gp]
  
  col <- cloud_colors[ccat]; if(is.na(col)) col <- "#BDBDBD"
  
  # --- Inclure aussi les Schwerpunkt dans le calcul du cadre ---
  all_x <- c(p_dat$GeneCount, met$Schwerpunkt_RNA_C, met$Schwerpunkt_RNA_D)
  all_y <- c(p_dat$Protein_Raw, met$Schwerpunkt_Prot_C, met$Schwerpunkt_Prot_D)
  
  cx <- mean(range(all_x))
  cy <- mean(range(all_y))
  ps <- max(diff(range(all_x)), diff(range(all_y))) * 1.5
  
  # --- Create a local Schwerpunkt data frame ---
  df_sp <- data.frame(
    Group = c("Control", "Disease"),
    RNA = c(met$Schwerpunkt_RNA_C[1], met$Schwerpunkt_RNA_D[1]),
    Prot = c(met$Schwerpunkt_Prot_C[1], met$Schwerpunkt_Prot_D[1])
  )
  
  x1 <- df_sp$RNA[1]; y1 <- df_sp$Prot[1]
  x2 <- df_sp$RNA[2]; y2 <- df_sp$Prot[2]
  dx <- x2 - x1; dy <- y2 - y1
  
  cp <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=DiseaseGroup)) +
    geom_point(size=1.2, alpha=0.8) +
    scale_color_manual(values=c("Control"="#2196F3", "Disease"="#F44336")) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    theme_bw(base_size=T_AXIS) +
    theme(
      legend.position="none",
      plot.title=element_text(color=col, size=T_TITLE - 4, hjust=0.5, face="bold"),
      plot.subtitle=element_text(size=T_SUB - 2, color="black", hjust=0.5),
      axis.title=element_text(size=T_TICKS - 2),
      axis.text=element_text(size=8, color="black"),
      axis.ticks=element_line(color="black", linewidth=0.8),
      panel.border=element_rect(color=col, fill=NA, linewidth=1.2),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    )
  
  # --- Ligne Noire Longue traversant les Schwerpunkte ---
  if (abs(dx) < 1e-6) {
    cp <- cp + geom_vline(xintercept = x1, color = "black", linewidth = 0.8)
  } else {
    slope <- dy / dx
    intercept <- y1 - (slope * x1)
    cp <- cp + geom_abline(slope = slope, intercept = intercept, color = "black", linewidth = 0.8)
  }
  
  # --- Points Schwerpunkt par-dessus (Refined Visualization) ---
  cp <- cp +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 2.5, stroke = 1.0, 
               color = c("Control" = "#0D47A1", "Disease" = "#B71C1C")[df_sp$Group], 
               inherit.aes = FALSE) +
    labs(
      title = paste0("Pair ", g, ":\n", gp),
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
  
  clouds_list[[as.character(g)]] <- cp
}
# -------------------------------------------------------------------------
# 5. ASSEMBLE ENTIRE MASTER PLATE
# -------------------------------------------------------------------------
# Top Row: PCA & Volcano (2x2 structure)
top_row <- (p1_pca_rna | p1_pca_prot | p2_volc_rna | p2_volc_prot) + plot_layout(ncol = 4)

# Bottom half: DiPa and Clouds (1 big DiPa + 8 small clouds)
# Using nested design with a gigantic central DiPa:
design <- "
  AAABBBCCCDDD
  EEEEEEFFFGGG
  EEEEEEHHHIII
  EEEEEEJJJKKK
  EEEEEEMMMNNN
"
# A=PCA_RNA, B=PCA_Prot, C=Volc_RNA, D=Volc_Prot
# E=DiPa (spans 4 rows, massive width & height)
# F=1, G=2, H=3, I=4, J=5, K=6, M=7, N=8

layout_master <- p1_pca_rna + p1_pca_prot + p2_volc_rna + p2_volc_prot +
  p3_dipa + 
  clouds_list[["1"]] + clouds_list[["2"]] + clouds_list[["3"]] + 
  clouds_list[["4"]] + clouds_list[["5"]] + clouds_list[["6"]] + 
  clouds_list[["7"]] + clouds_list[["8"]] +
  plot_layout(design = design) +
  plot_annotation(
    title = "CCL4 DATASET",
    subtitle = paste0(
      "Using 943 overlapping genes between BDL and CCl4 datasets\n",
      "** Points:  Blue = Control (Oil) | Red = Disease (CCl4) **\n",
      "** CROSSES (+) = SCHWERPUNKTE (Group Centroids) in Dark Blue/Red **\n",
      paste(extreme_pts$LegendText[1:4], collapse = " | "), "\n", 
      paste(extreme_pts$LegendText[5:8], collapse = " | ")
    ),
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5), 
                  plot.subtitle = element_text(size = 14, color = "black", hjust = 0.5, margin = margin(b = 10)))
  )

ggsave("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/A4_MASTER_CCL4_Combined_Schwerpunkt_normal_data.pdf", plot = layout_master, width = 11.7, height = 14.0, device = cairo_pdf)
print("SUCCESS! File A4_MASTER_CCL4_Combined_Schwerpunkt_normal_data.pdf was saved!")

# -------------------------------------------------------------------------
# 6. SAVE ISOLATED PLOTS (3 SEPARATE FILES)
# -------------------------------------------------------------------------
print("Saving isolated plots for CCL4...")

# MODIFIED: Removed global titles, set subplot titles to "RNA" and "Protein", and increased text sizes for isolated PCA
layout_pca <- ((p1_pca_rna + labs(title = "RNA", tag = NULL)) | (p1_pca_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Large title
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16)
  )
ggsave("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Isolated_PCA_CCL4.pdf", plot = layout_pca, width = 10, height = 5, device = cairo_pdf)

# MODIFIED: Removed global titles, set subplot titles to "RNA" and "Protein", and increased text sizes for isolated Volcano
layout_volcano <- ((p2_volc_rna + labs(title = "RNA", tag = NULL)) | (p2_volc_prot + labs(title = "Protein"))) & 
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5), # Large title
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 18, color = "black")
  )
ggsave("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Isolated_Volcano_CCL4.pdf", plot = layout_volcano, width = 10, height = 5, device = cairo_pdf)



# Plot 3: DiPa + Wolken only (Rebuilt with larger fonts, points, and labels for the isolated PDF)
p3_dipa_iso <- ggplot(Combined, aes(x = log2G, y = log2P, color = CloudCategory)) +
  geom_point(alpha = 0.5, size = 1.8) + scale_color_manual(values = cloud_colors) +
  geom_hline(yintercept = 0.5, color = "black", linewidth = 0.4) + geom_hline(yintercept = -0.5, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0.5, color = "black", linewidth = 0.4) + geom_vline(xintercept = -0.5, color = "black", linewidth = 0.4) +
  xlim(-6, 6) + ylim(-5, 5) +
  labs(title = "", tag = NULL, x = expression(RNA*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), y = expression(Protein*":"~CCl[4]~(Month~12)~vs.~Oil~(Month~0)~~log[2]~"(fold change)"), color = "", caption = "DiPa 1&2 → Diagonal | 5&6 → Horizontal | 3&4 → Vertical | 8 → Round | 7 → Discordant") + theme_bw(base_size = 20) +
  theme(plot.title = element_text(size = 24, hjust=0.5), axis.title = element_text(size = 22, face="bold"), axis.text = element_text(size = 18, color="black"),
        axis.ticks = element_line(color="black", linewidth=0.8), panel.border = element_rect(color="black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 28),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 18), legend.key.size = unit(0.5, "cm"), legend.margin=margin(t=-10)) +
  annotate("text", x =  0,    y =  0,    label = "8", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x =  4.0,  y =  3.5,  label = "1", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x = -4.0,  y = -3.5,  label = "2", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x =  0,    y =  3.5,  label = "3", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x =  0,    y = -3.5,  label = "4", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x =  4.0,  y =  0,    label = "5", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x = -4.0,  y =  0,    label = "6", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x = -4.0,  y =  3.5,  label = "7", size = 10, color = "grey50", fontface = "bold") +
  annotate("text", x =  4.0,  y = -3.5,  label = "7", size = 10, color = "grey50", fontface = "bold") +
  geom_point(data = extreme_pts, aes(x = log2G, y = log2P), shape = 18, size = 4.5, color = "black", inherit.aes=FALSE) +
  geom_text_repel(data = extreme_pts, aes(label = PlotMarker), size = 6.5, color="black", bg.color="white", bg.r=0.15, box.padding=0.1, max.overlaps=30, segment.color="grey40")

clouds_list_iso <- list()
for(g in 1:8) {
  ep <- extreme_pts[DiPaGroups == as.character(g)]
  if(nrow(ep)==0) { clouds_list_iso[[as.character(g)]] <- ggplot()+theme_void(); next }
  
  gp <- ep$GeneProtein[1]
  ccat <- ep$CloudCategory[1]
  p_dat <- ccl4_DT[GeneProtein == gp]
  met <- Cloud_Metrics_CCL4[GeneProtein == gp]
  
  col <- cloud_colors[ccat]; if(is.na(col)) col <- "#BDBDBD"
  
  all_x <- c(p_dat$GeneCount, met$Schwerpunkt_RNA_C, met$Schwerpunkt_RNA_D)
  all_y <- c(p_dat$Protein_Raw, met$Schwerpunkt_Prot_C, met$Schwerpunkt_Prot_D)
  
  cx <- mean(range(all_x))
  cy <- mean(range(all_y))
  ps <- max(diff(range(all_x)), diff(range(all_y))) * 1.5
  
  df_sp <- data.frame(
    Group = c("Control", "Disease"),
    RNA = c(met$Schwerpunkt_RNA_C[1], met$Schwerpunkt_RNA_D[1]),
    Prot = c(met$Schwerpunkt_Prot_C[1], met$Schwerpunkt_Prot_D[1])
  )
  
  x1 <- df_sp$RNA[1]; y1 <- df_sp$Prot[1]
  x2 <- df_sp$RNA[2]; y2 <- df_sp$Prot[2]
  dx <- x2 - x1; dy <- y2 - y1
  
  cp <- ggplot(p_dat, aes(x=GeneCount, y=Protein_Raw, color=DiseaseGroup)) +
    geom_point(size=2.2, alpha=0.8) +
    scale_color_manual(values=c("Control"="#2196F3", "Disease"="#F44336")) +
    coord_cartesian(xlim=c(cx-ps/2, cx+ps/2), ylim=c(cy-ps/2, cy+ps/2)) +
    theme_bw(base_size=14) +
    theme(
      legend.position="none",
      plot.title=element_text(color=col, size=18, hjust=0.5, face="bold"),
      plot.subtitle=element_text(size=14, color="black", hjust=0.5),
      axis.title=element_text(size=14, face="bold"),
      axis.text=element_text(size=12, color="black"),
      axis.ticks=element_line(color="black", linewidth=0.8),
      panel.border=element_rect(color=col, fill=NA, linewidth=1.5),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank()
    )
  
  if (abs(dx) < 1e-6) {
    cp <- cp + geom_vline(xintercept = x1, color = "black", linewidth = 1.0)
  } else {
    slope <- dy / dx
    intercept <- y1 - (slope * x1)
    cp <- cp + geom_abline(slope = slope, intercept = intercept, color = "black", linewidth = 1.0)
  }
  
  cp <- cp +
    geom_point(data = df_sp, aes(x = RNA, y = Prot), 
               shape = 3, size = 4.0, stroke = 1.5, 
               color = c("Control" = "#0D47A1", "Disease" = "#B71C1C")[df_sp$Group], 
               inherit.aes = FALSE) +
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
  
  clouds_list_iso[[as.character(g)]] <- cp
}

design_isolated_dipa <- "
  EEEEEEFFFGGG
  EEEEEEHHHIII
  EEEEEEJJJKKK
  EEEEEEMMMNNN
"
layout_dipa_isolated <- (p3_dipa_iso + labs(tag="C") + theme(plot.tag = element_text(face = "bold", size = 32))) + 
  (clouds_list_iso[["1"]] + labs(tag="D") + theme(plot.tag = element_text(face = "bold", size = 32))) + clouds_list_iso[["2"]] + clouds_list_iso[["3"]] + 
  clouds_list_iso[["4"]] + clouds_list_iso[["5"]] + clouds_list_iso[["6"]] + 
  clouds_list_iso[["7"]] + clouds_list_iso[["8"]] +
  plot_layout(design = design_isolated_dipa)

ggsave("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Isolated_DiPa_Wolken_CCL4.pdf", plot = layout_dipa_isolated, width = 19, height = 17.5, device = cairo_pdf)

print("SUCCESS! Isolated plots for CCL4 were saved.")
