# ==============================================================================
# SCRIPT: Schwerpunkt_Slope_Distribution_BDL.R
# OBJECTIVE: Calculate transition slopes between control/disease centroids
#            and visualize their distribution for Pair 1 and Pair 2 (BDL).
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

# -------------------------------------------------------------------------
# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading BDL raw data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
Full_DT <- DTccl4_DT_LCPM
setDT(Full_DT)

# Load DiPa categorization (BDL)
print("Loading BDL DiPa data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData")
dipa_dt <- DT_dipa_count_bdl
setDT(dipa_dt)

# -------------------------------------------------------------------------
# 2. CALCULATE CENTROID SLOPES (SP_SLOPE)
# -------------------------------------------------------------------------
# Preprocess BDL data
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]
Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "control", DiseaseGroup := "Control"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL",     DiseaseGroup := "Disease"]

bdl_subset <- Full_DT[Dataset == "BDL Dataset" & !is.na(GeneCount) & !is.na(ProteinIntensity) & !is.na(DiseaseGroup)]

#length(unique(bdl_subset$MiceInfo)) #12

# Calculate Group Centroids per Gene
sp_metrics <- bdl_subset[, {
  rna_c  <- mean(GeneCount[DiseaseGroup == "Control"], na.rm=TRUE)
  prot_c <- mean(ProteinIntensity[DiseaseGroup == "Control"], na.rm=TRUE)
  
  rna_d  <- mean(GeneCount[DiseaseGroup == "Disease"], na.rm=TRUE)
  prot_d <- mean(ProteinIntensity[DiseaseGroup == "Disease"], na.rm=TRUE)
  
  dx <- rna_d - rna_c
  dy <- prot_d - prot_c
  
  .(
    RNA_C = rna_c, Prot_C = prot_c,
    RNA_D = rna_d, Prot_D = prot_d,
    dX = dx, dY = dy,
    SP_Slope = dy / dx
  )
}, by = GeneProtein]

# Merge with DiPa groups
final_dt <- merge(sp_metrics, dipa_dt[, .(GeneProtein, DiPaGroups)], by = "GeneProtein")
final_dt[, DiPaGroups := as.character(DiPaGroups)]

# -------------------------------------------------------------------------
# 3. VISUALIZATION (MATCHING HARMONIZED STYLE)
# -------------------------------------------------------------------------

# plot_sp_slope_dist <- function(data, group_id, title_suffix, color_fill) {
#   df <- data[DiPaGroups == as.character(group_id)]
#   if(nrow(df) == 0) return(ggplot() + theme_void())
#   
#   pct <- round(100 * sum(df$SP_Slope >= 0.5 & df$SP_Slope <= 1.5, na.rm=TRUE) / nrow(df), 1)
#   
#   ggplot(df, aes(x = SP_Slope)) +
#     geom_histogram(fill = color_fill, alpha = 0.6, binwidth = 0.1, color = "black", linewidth = 0.2) +
#     geom_vline(xintercept = c(0.5, 1.5), linetype = "dashed", color = "black", linewidth = 1) +
#     annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = Inf, alpha = 0.1, fill = "green") +
#     xlim(-1, 3) +
#     labs(title = paste("Centroid Slopes:", title_suffix),
#          subtitle = paste0(pct, "% fall in [0.5, 1.5] range\nN = ", nrow(df)),
#          x = "Schwerpunkt Slope (m)", y = "Frequency") +
#     theme_minimal() +
#     theme(plot.title = element_text(face="bold", size=14),
#           plot.subtitle = element_text(size=10, color="grey30"))
# }
# 
# p1 <- plot_sp_slope_dist(final_dt, "1", "DiPa Group 1", "blue")
# p2 <- plot_sp_slope_dist(final_dt, "2", "DiPa Group 2", "red")

# -------------------------------------------------------------------------
# 4. COMBINED DENSITY PLOT (GROUPS 1 & 2)
# -------------------------------------------------------------------------
df_comb <- final_dt[DiPaGroups %in% c("1", "2")]
df_comb[, DiPaGroups := factor(DiPaGroups, levels = c("1", "2"), labels = c("1", "2"))]

n1 <- nrow(df_comb[DiPaGroups == "1"])
n2 <- nrow(df_comb[DiPaGroups == "2"])
pct1 <- 100 * sum(df_comb$SP_Slope[df_comb$DiPaGroups == "1"] >= 0.5 & df_comb$SP_Slope[df_comb$DiPaGroups == "1"] <= 1.5, na.rm=TRUE) / n1
pct2 <- 100 * sum(df_comb$SP_Slope[df_comb$DiPaGroups == "2"] >= 0.5 & df_comb$SP_Slope[df_comb$DiPaGroups == "2"] <= 1.5, na.rm=TRUE) / n2
max_y <- max(density(df_comb$SP_Slope[df_comb$DiPaGroups == "1"], na.rm=TRUE)$y,
             density(df_comb$SP_Slope[df_comb$DiPaGroups == "2"], na.rm=TRUE)$y)

p_comb <- ggplot(df_comb, aes(x = SP_Slope, fill = DiPaGroups)) +
  geom_density(alpha = 0.40) +
  geom_vline(xintercept = c(0.5, 1.5), linetype = "dashed", color = "black", linewidth = 1.0) +
  scale_fill_manual(
    values = c("1" = "blue", "2" = "red"),
    labels = c("1" = "DiPa group 1", "2" = "DiPa group 2")
  ) +
  xlim(-1, 3) +
  labs(
    title = "BDL",
    subtitle = NULL,
    x = "Centroid slope",
    y = "Density",
    fill = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_blank(),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 18, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 1.0),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, face = "bold", color = "black"),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1.7), breaks = c(0, 0.5, 1.0, 1.5)) +
  annotate("text", x = 2.2, y = 1.7 * 0.76, label = paste0("bold(N[1] == ", n1, ")"), parse = TRUE, size = 8, color = "blue", hjust = 0) +
  annotate("text", x = 2.2, y = 1.7 * 0.67, label = paste0("bold(N[2] == ", n2, ")"), parse = TRUE, size = 8, color = "red", hjust = 0)

# -------------------------------------------------------------------------
# 5. SAVE RESULTS
# -------------------------------------------------------------------------
pdf_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Centroid_Slope_Frequency_BDL.pdf"
cairo_pdf(pdf_path, width = 10, height = 10)
print(p_comb)
dev.off()

# stats <- final_dt[DiPaGroups %in% c("1", "2"), .(
#   N = .N,
#   Med_Slope = median(SP_Slope, na.rm=TRUE),
#   Perc_in_Range = round(100 * sum(SP_Slope >= 0.5 & SP_Slope <= 1.5, na.rm=TRUE) / .N, 1)
# ), by = DiPaGroups]
# 
# csv_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Schwerpunkt_Slope_Stats_BDL.csv"
# fwrite(stats, csv_path)
# 
# print(paste("SUCCESS! PDF saved to:", pdf_path))
