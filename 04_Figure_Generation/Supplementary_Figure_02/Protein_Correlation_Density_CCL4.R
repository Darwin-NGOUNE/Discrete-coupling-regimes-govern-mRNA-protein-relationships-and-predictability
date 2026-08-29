# ==============================================================================
# SCRIPT: Protein_Correlation_Density_CCL4.R
# OBJECTIVE: Quantify the 'Correlatability' of the protein landscape.
#            Calculate the distribution of Top N absolute correlations (|r|)
#            among 861 proteins in the CCL4 dataset.
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)

# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading CCL4 data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData")
CCL4_DT <- DTccl4_DT_LCPM[!is.na(TreatmentTime)] 
#CCL4_DT <- Full_DT[!is.na(TreatmentTime)] 

CCL4_DT <- na.omit(CCL4_DT) 
unique(CCL4_DT$MiceInfo) #31 mice

# Filter for Month 12 CCl4 vs Month 0 Oil and exclude outlier rep6
#sub_dt <- CCL4_DT[(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)]
#sub_dt <- sub_dt[MiceInfo != "month12_ccl4_rep6"]

# Create a mapping for ProteinID (Uniprot Accessions)
#protein_map <- unique(sub_dt[, .(GeneProtein, ProteinID)])
protein_map <- unique(CCL4_DT[, .(GeneProtein, ProteinID)])
setkey(protein_map, GeneProtein)

# 2. PREPARE MATRIX (Complete cases only)
# -------------------------------------------------------------------------
#dat_wide <- dcast(sub_dt, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
dat_wide <- dcast(CCL4_DT, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
mat <- as.matrix(dat_wide[, -1, with = FALSE])
rownames(mat) <- dat_wide$GeneProtein

# Complete cases only (n=801)
mat_clean <- na.omit(mat)
print(paste("Proteins analyzed:", nrow(mat_clean))) 

# 3. COMPUTE CORRELATION MATRIX (Pearson)
# -------------------------------------------------------------------------
print("Calculating correlation matrix...")
cor_mat <- cor(t(mat_clean), method = "pearson")
abs_cor_mat <- abs(cor_mat)

# 4. EXTRACT TOP N CORRELS PER PROTEIN
# -------------------------------------------------------------------------
diag(abs_cor_mat) <- 0 # Ignore self

extract_top_n_stats <- function(mat, n_val) {
  top_vals <- apply(mat, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE)
    return(sorted_x[n_val])
  })
  return(data.table(Protein = names(top_vals), Max_R = top_vals, Tier = paste0("Top ", n_val)))
}

tiers <- c(1, 2, 5, 10, 15, 20, 50, 70, 100)
results_list <- lapply(tiers, function(n) extract_top_n_stats(abs_cor_mat, n))
combined_results <- rbindlist(results_list)

# 5. VISUALIZATION
# -------------------------------------------------------------------------
pdf_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Protein_Correlation_Density_CCL4.pdf"
cairo_pdf(pdf_path, width = 12, height = 8)

p_dens <- ggplot(combined_results, aes(x = Max_R, fill = factor(Tier, levels = paste0("Top ", tiers)))) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(title = expression(CCl[4]),
       x = expression(Absolute~Correlation~Coefficient~(abs(rho[BP]))), 
       y = "Density",
       fill = "Correlation Rank") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

p_facet <- ggplot(combined_results, aes(x = Max_R, fill = Tier)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.1) +
  facet_wrap(~factor(Tier, levels = paste0("Top ", tiers)), scales = "free_y") +
  labs(title = expression(CCl[4]),
       x = expression(Absolute~Correlation~Coefficient~(abs(rho[BP]))),
       y = "Frequency") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")

# Page 1: Overlaid Densities
print(p_dens)

# Page 2: Faceted Histograms
print(p_facet)

# -------------------------------------------------------------------------
# 6. GLOBAL CENTRALITY (MASTER DRIVERS)
# -------------------------------------------------------------------------
centrality_stats <- data.table(
  GeneProtein = rownames(abs_cor_mat),
  Min_Abs_R   = apply(abs_cor_mat, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat, 1, mean),
  Thresh_80pct = apply(abs_cor_mat, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat, 1, function(x) quantile(x, probs = 0.10))
)

# Merge with IDs for plotting
centrality_stats <- merge(centrality_stats, protein_map, by = "GeneProtein")
setnames(centrality_stats, "ProteinID", "Protein")

# Page 3: Top 10 Master Proteins (saved) and Top 30 Master Proteins (plotted)
top_10_min <- centrality_stats[order(-Min_Abs_R)][1:10]
top_30_min <- centrality_stats[order(-Min_Abs_R)][1:30]

p_top_min <- ggplot(top_30_min, aes(x = reorder(Protein, Min_Abs_R), y = Min_Abs_R)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  labs(title = expression(CCl[4]),
       x = "Protein", y = expression(Minimum~abs(rho[BP]))) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"))

# Page 3: Top 50 Master Proteins (Mean |rho_BP|)
top_50_mean <- centrality_stats[order(-Mean_Abs_R)][1:50]

p_top_mean <- ggplot(top_50_mean, aes(x = reorder(Protein, Mean_Abs_R), y = Mean_Abs_R)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = expression(CCl[4]),
       x = "Protein", y = expression(Mean~abs(rho[BP]))) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 8.5, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

# Print Page 3 (Top 50 Hubs)
print(p_top_mean)

# Page 4: Top 50 Mastery (80% Coverage)
top_50_80pct <- centrality_stats[order(-Thresh_80pct)][1:50]

p_top_80 <- ggplot(top_50_80pct, aes(x = reorder(Protein, Thresh_80pct), y = Thresh_80pct)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  labs(title = expression(CCl[4]),
       x = "Protein", y = expression(Threshold~abs(rho[BP])~(80*'% Coverage'))) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 8.5, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

# Page 5: Top 50 Mastery (90% Coverage)
top_50_90pct <- centrality_stats[order(-Thresh_90pct)][1:50]

p_top_90 <- ggplot(top_50_90pct, aes(x = reorder(Protein, Thresh_90pct), y = Thresh_90pct)) +
  geom_bar(stat = "identity", fill = "brown") +
  coord_flip() +
  labs(title = expression(CCl[4]),
       x = "Protein", y = expression(Threshold~abs(rho[BP])~(90*'% Coverage'))) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 8.5, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

# Print Page 4 (Top 50 Mastery 80% Coverage)
print(p_top_80)

# Print Page 5 (Top 50 Mastery 90% Coverage)
print(p_top_90)

dev.off()

#rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_CCL4_Raw.RData"
#save(top_10_mean, top_10_80pct, top_10_90pct, file = rdata_path)


top_10_mean  <- centrality_stats[order(-Mean_Abs_R)][1:10]
top_10_80pct <- centrality_stats[order(-Thresh_80pct)][1:10]
top_10_90pct <- centrality_stats[order(-Thresh_90pct)][1:10]

rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_CCL4_Batch.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = rdata_path)


# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_CCL4_Batch_Full.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))
# 
# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_CCL4_Batch_Subset.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))
# 
# 
# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_CCL4_Full.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))
# 
# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_CCL4_Subset.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))



print(paste("SUCCESS! CCL4 Correlation Density & Centrality report saved to:", pdf_path))
