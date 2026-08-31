# ==============================================================================
# SCRIPT: Protein_Correlation_Density_BDL.R
# OBJECTIVE: Quantify the 'Correlatability' of the protein landscape.
#            Calculate the distribution of Top N absolute correlations (|r|)
#            among 943 proteins in the BDL dataset.
# ==============================================================================

library(data.table)
library(ggplot2)
library(gridExtra)

# 1. LOAD DATA
# -------------------------------------------------------------------------
print("Loading BDL data...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData")
BDL_DT <- DTccl4_DT_LCPM[is.na(TreatmentTime)] # BDL samples have NA TreatmentTime
#BDL_DT <- Full_DT[is.na(TreatmentTime)] # BDL samples have NA TreatmentTime

#sub_dt <- BDL_DT[Treatment %in% c("BDL", "control")]

# Create a mapping for ProteinID (Uniprot Accessions)
#protein_map <- unique(sub_dt[, .(GeneProtein, ProteinID)])
protein_map <- unique(BDL_DT[, .(GeneProtein, ProteinID)])
setkey(protein_map, GeneProtein)

# 2. PREPARE MATRIX (943 Proteins x 12 Samples)
# -------------------------------------------------------------------------
#dat_wide <- dcast(sub_dt, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
dat_wide <- dcast(BDL_DT, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity") 
mat <- as.matrix(dat_wide[, -1, with = FALSE])
rownames(mat) <- dat_wide$GeneProtein

# Ensure no missing values for correlation
mat_clean <- na.omit(mat)
print(paste("Proteins analyzed:", nrow(mat_clean)))

# 3. COMPUTE CORRELATION MATRIX (Pearson)
# -------------------------------------------------------------------------
print("Calculating 943x943 correlation matrix...")
cor_mat <- cor(t(mat_clean), method = "pearson")
# Take absolute values to handle anti-correlations
abs_cor_mat <- abs(cor_mat)

# 4. EXTRACT TOP N CORRELS PER PROTEIN
# -------------------------------------------------------------------------
# Diagonal is 1.0 (self), so we must ignore it
diag(abs_cor_mat) <- 0 

extract_top_n_stats <- function(mat, n_val) {
  # Sort each row and pick the n-th largest value
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
pdf_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Protein_Correlation_Density_BDL.pdf"
cairo_pdf(pdf_path, width = 12, height = 8)

# Figure 1: Overlaid Densities (Clean for Figure 2 Master Plate)
p_dens <- ggplot(combined_results, aes(x = Max_R, fill = factor(Tier, levels = paste0("Top ", tiers)))) +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(title = NULL,
       x = expression(bold("|"*rho[BP]*"|")), 
       y = "Density",
       fill = "Correlation Rank") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 22, face = "bold", color = "black"),
        axis.text = element_text(size = 18, face = "bold", color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 1.0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Figure 2: Faceted View for Clarity
p_facet <- ggplot(combined_results, aes(x = Max_R, fill = Tier)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.1) +
  facet_wrap(~factor(Tier, levels = paste0("Top ", tiers)), scales = "free_y") +
  labs(title = "BDL",
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
# Identify proteins that correlate well with EVERYTHING
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
top_10_min <- centrality_stats[order(-Min_Abs_R)][1:10] #10
top_30_min <- centrality_stats[order(-Min_Abs_R)][1:30]

p_top_min <- ggplot(top_30_min, aes(x = reorder(Protein, Min_Abs_R), y = Min_Abs_R)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  labs(title = "BDL",
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
  labs(title = "BDL",
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
  labs(title = "BDL",
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
  labs(title = "BDL",
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

#rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_BDL_Raw.RData"
#save(top_10_mean, top_10_80pct, top_10_90pct, file = rdata_path)



top_10_mean  <- centrality_stats[order(-Mean_Abs_R)][1:10]
top_10_80pct <- centrality_stats[order(-Thresh_80pct)][1:10]
top_10_90pct <- centrality_stats[order(-Thresh_90pct)][1:10]

rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_BDL_Batch.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = rdata_path)



# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_BDL_Batch_Full.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))

# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_BDL_Batch_Subset.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))


# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_BDL_Full.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))

# rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top30_Mastery_Hubs_BDL_Subset.RData"
# save(top_30_mean, top_30_80pct, top_30_90pct, file = rdata_path)
# print(paste("Saved Top 30 RData to:", rdata_path))

print(paste("SUCCESS! BDL Correlation Density & Centrality report saved to:", pdf_path))



# ==============================================================================
# OPTION 3: MERGED DATASET MASTERY HUB PROTEINS CALCULATION (Darwin's addition)
# ==============================================================================
# This section calculates the Top 10 Mastery Hubs directly on the combined
# BDL + CCl4 merged dataset (after ComBat batch correction) so they can be
# used for Procedure 3 (merged data modeling).

print("Loading Merged Batch-Corrected dataset...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData")
Full_DT <- as.data.table(Full_DT)

# Use ComBat corrected values for protein intensities
Full_DT$ProteinIntensity <- Full_DT$ComBat_Protein_Raw

# Identify the 943 core overlap proteins in the merged dataset
prot_bdl <- unique(Full_DT[Treatment %in% c("BDL", "control", "BDL_ASBTi"), GeneProtein])
prot_ccl4 <- unique(Full_DT[Treatment %in% c("ccl4", "oil"), GeneProtein])
core_merged_proteins <- intersect(prot_bdl, prot_ccl4)

Full_DT_core <- Full_DT[GeneProtein %in% core_merged_proteins]

# Create a mapping for ProteinID (Uniprot Accessions)
protein_map_merged <- unique(Full_DT_core[, .(GeneProtein, ProteinID)])
setkey(protein_map_merged, GeneProtein)

# Prepare Matrix (943 Proteins x all mice in merged dataset)
dat_wide_merged <- dcast(Full_DT_core, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_merged <- as.matrix(dat_wide_merged[, -1, with = FALSE])
rownames(mat_merged) <- dat_wide_merged$GeneProtein

# Omit missing values for correlation (excluding mice with no protein data)
non_empty_cols_merged <- colnames(mat_merged)[colSums(!is.na(mat_merged)) > 0]
mat_merged_valid <- mat_merged[, non_empty_cols_merged, drop = FALSE]
mat_merged_clean <- na.omit(mat_merged_valid)
print(paste("Merged Proteins analyzed:", nrow(mat_merged_clean)))

# Compute Correlation Matrix (Pearson)
print("Calculating 943x943 merged correlation matrix...")
cor_mat_merged <- cor(t(mat_merged_clean), method = "pearson")
abs_cor_mat_merged <- abs(cor_mat_merged)
diag(abs_cor_mat_merged) <- 0

# Compute Centrality Statistics
print("Calculating merged centrality statistics...")
centrality_stats_merged <- data.table(
  GeneProtein = rownames(abs_cor_mat_merged),
  Min_Abs_R   = apply(abs_cor_mat_merged, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_merged, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_merged, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_merged, 1, function(x) quantile(x, probs = 0.10))
)

centrality_stats_merged <- merge(centrality_stats_merged, protein_map_merged, by = "GeneProtein")
setnames(centrality_stats_merged, "ProteinID", "Protein")

# Extract Top 10 Mastery Hubs
top_10_mean  <- centrality_stats_merged[order(-Mean_Abs_R)][1:10]
top_10_80pct <- centrality_stats_merged[order(-Thresh_80pct)][1:10]
top_10_90pct <- centrality_stats_merged[order(-Thresh_90pct)][1:10]

# Save Merged Mastery RData file
merged_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_Merged_Batch.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = merged_rdata_path)

print(paste("SUCCESS! Merged Mastery Hubs saved to:", merged_rdata_path))
print("Selected Top 10 Merged Mastery Hubs (80% Coverage threshold):")
print(top_10_80pct[, .(Protein, Thresh_80pct)])


# ==============================================================================
# BATCH-CORRECTED BDL & CCL4 MASTERY HUB PROTEINS CALCULATION (Darwin's addition)
# ==============================================================================
# This section calculates the Top 10 Mastery Hubs for BDL and CCL4 separately
# using the batch-corrected dataset (Full_DT with ComBat_Protein_Raw).
# This corrects the original bug where BDL_Batch and CCL4_Batch files were
# calculated on raw, uncorrected data.

# A. BDL BATCH-CORRECTED HUBS
print("\nCalculating BDL batch-corrected centrality statistics...")
BDL_DT_batch <- Full_DT_core[is.na(TreatmentTime)] # BDL has NA TreatmentTime

dat_wide_bdl_batch <- dcast(BDL_DT_batch, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_bdl_batch <- as.matrix(dat_wide_bdl_batch[, -1, with = FALSE])
rownames(mat_bdl_batch) <- dat_wide_bdl_batch$GeneProtein

non_empty_cols_bdl_batch <- colnames(mat_bdl_batch)[colSums(!is.na(mat_bdl_batch)) > 0]
mat_bdl_batch_valid <- mat_bdl_batch[, non_empty_cols_bdl_batch, drop = FALSE]
mat_bdl_batch_clean <- na.omit(mat_bdl_batch_valid)
print(paste("BDL Batch Proteins analyzed:", nrow(mat_bdl_batch_clean)))

cor_mat_bdl <- cor(t(mat_bdl_batch_clean), method = "pearson")
abs_cor_mat_bdl <- abs(cor_mat_bdl)
diag(abs_cor_mat_bdl) <- 0

centrality_stats_bdl <- data.table(
  GeneProtein = rownames(abs_cor_mat_bdl),
  Min_Abs_R   = apply(abs_cor_mat_bdl, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_bdl, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_bdl, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_bdl, 1, function(x) quantile(x, probs = 0.10))
)
centrality_stats_bdl <- merge(centrality_stats_bdl, protein_map_merged, by = "GeneProtein")
setnames(centrality_stats_bdl, "ProteinID", "Protein")

top_10_mean_bdl  <- centrality_stats_bdl[order(-Mean_Abs_R)][1:10]
top_10_80pct_bdl <- centrality_stats_bdl[order(-Thresh_80pct)][1:10]
top_10_90pct_bdl <- centrality_stats_bdl[order(-Thresh_90pct)][1:10]

# Assign to target names and save
top_10_mean  <- top_10_mean_bdl
top_10_80pct <- top_10_80pct_bdl
top_10_90pct <- top_10_90pct_bdl

bdl_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_BDL_Batch.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = bdl_rdata_path)
print(paste("SUCCESS! BDL Batch-Corrected Mastery Hubs saved to:", bdl_rdata_path))


# B. CCL4 BATCH-CORRECTED HUBS
print("\nCalculating CCL4 batch-corrected centrality statistics...")
CCL4_DT_batch <- Full_DT_core[!is.na(TreatmentTime)] # CCL4 has non-NA TreatmentTime

dat_wide_ccl4_batch <- dcast(CCL4_DT_batch, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_ccl4_batch <- as.matrix(dat_wide_ccl4_batch[, -1, with = FALSE])
rownames(mat_ccl4_batch) <- dat_wide_ccl4_batch$GeneProtein

non_empty_cols_ccl4_batch <- colnames(mat_ccl4_batch)[colSums(!is.na(mat_ccl4_batch)) > 0]
mat_ccl4_batch_valid <- mat_ccl4_batch[, non_empty_cols_ccl4_batch, drop = FALSE]
mat_ccl4_batch_clean <- na.omit(mat_ccl4_batch_valid)
print(paste("CCL4 Batch Proteins analyzed:", nrow(mat_ccl4_batch_clean)))

cor_mat_ccl4 <- cor(t(mat_ccl4_batch_clean), method = "pearson")
abs_cor_mat_ccl4 <- abs(cor_mat_ccl4)
diag(abs_cor_mat_ccl4) <- 0

centrality_stats_ccl4 <- data.table(
  GeneProtein = rownames(abs_cor_mat_ccl4),
  Min_Abs_R   = apply(abs_cor_mat_ccl4, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_ccl4, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_ccl4, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_ccl4, 1, function(x) quantile(x, probs = 0.10))
)
centrality_stats_ccl4 <- merge(centrality_stats_ccl4, protein_map_merged, by = "GeneProtein")
setnames(centrality_stats_ccl4, "ProteinID", "Protein")

top_10_mean_ccl4  <- centrality_stats_ccl4[order(-Mean_Abs_R)][1:10]
top_10_80pct_ccl4 <- centrality_stats_ccl4[order(-Thresh_80pct)][1:10]
top_10_90pct_ccl4 <- centrality_stats_ccl4[order(-Thresh_90pct)][1:10]

# Assign to target names and save
top_10_mean  <- top_10_mean_ccl4
top_10_80pct <- top_10_80pct_ccl4
top_10_90pct <- top_10_90pct_ccl4

ccl4_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_CCL4_Batch.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = ccl4_rdata_path)
print(paste("SUCCESS! CCL4 Batch-Corrected Mastery Hubs saved to:", ccl4_rdata_path))


# ==============================================================================
# RAW (UNCORRECTED) BDL, CCL4 & MERGED MASTERY HUB PROTEINS CALCULATION
# ==============================================================================
# This section calculates the Top 10 Mastery Hubs for BDL, CCL4, and Merged
# separately using the raw uncorrected dataset (DTccl4_DT_LCPM).

print("\nLoading Raw Uncorrected dataset...")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")
DT_raw <- as.data.table(DTccl4_DT_LCPM)

# Identify the 943 core overlap proteins in the raw dataset
prot_bdl_raw <- unique(DT_raw[Treatment %in% c("BDL", "control", "BDL_ASBTi"), GeneProtein])
prot_ccl4_raw <- unique(DT_raw[Treatment %in% c("ccl4", "oil"), GeneProtein])
core_raw_proteins <- intersect(prot_bdl_raw, prot_ccl4_raw)

DT_raw_core <- DT_raw[GeneProtein %in% core_raw_proteins]

# Create a mapping for ProteinID (Uniprot Accessions)
protein_map_raw <- unique(DT_raw_core[, .(GeneProtein, ProteinID)])
setkey(protein_map_raw, GeneProtein)

# A. BDL RAW HUBS
print("Calculating BDL raw centrality statistics...")
BDL_DT_raw <- DT_raw_core[is.na(TreatmentTime)]

dat_wide_bdl_raw <- dcast(BDL_DT_raw, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_bdl_raw <- as.matrix(dat_wide_bdl_raw[, -1, with = FALSE])
rownames(mat_bdl_raw) <- dat_wide_bdl_raw$GeneProtein

non_empty_cols_bdl_raw <- colnames(mat_bdl_raw)[colSums(!is.na(mat_bdl_raw)) > 0]
mat_bdl_raw_valid <- mat_bdl_raw[, non_empty_cols_bdl_raw, drop = FALSE]
mat_bdl_raw_clean <- na.omit(mat_bdl_raw_valid)
print(paste("BDL Raw Proteins analyzed:", nrow(mat_bdl_raw_clean)))

cor_mat_bdl_raw <- cor(t(mat_bdl_raw_clean), method = "pearson")
abs_cor_mat_bdl_raw <- abs(cor_mat_bdl_raw)
diag(abs_cor_mat_bdl_raw) <- 0

centrality_stats_bdl_raw <- data.table(
  GeneProtein = rownames(abs_cor_mat_bdl_raw),
  Min_Abs_R   = apply(abs_cor_mat_bdl_raw, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_bdl_raw, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_bdl_raw, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_bdl_raw, 1, function(x) quantile(x, probs = 0.10))
)
centrality_stats_bdl_raw <- merge(centrality_stats_bdl_raw, protein_map_raw, by = "GeneProtein")
setnames(centrality_stats_bdl_raw, "ProteinID", "Protein")

top_10_mean_bdl_raw  <- centrality_stats_bdl_raw[order(-Mean_Abs_R)][1:10]
top_10_80pct_bdl_raw <- centrality_stats_bdl_raw[order(-Thresh_80pct)][1:10]
top_10_90pct_bdl_raw <- centrality_stats_bdl_raw[order(-Thresh_90pct)][1:10]

top_10_mean  <- top_10_mean_bdl_raw
top_10_80pct <- top_10_80pct_bdl_raw
top_10_90pct <- top_10_90pct_bdl_raw

bdl_raw_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_BDL_Raw.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = bdl_raw_rdata_path)
print(paste("SUCCESS! BDL Raw Mastery Hubs saved to:", bdl_raw_rdata_path))


# B. CCL4 RAW HUBS
print("Calculating CCL4 raw centrality statistics...")
CCL4_DT_raw <- DT_raw_core[!is.na(TreatmentTime)]

dat_wide_ccl4_raw <- dcast(CCL4_DT_raw, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_ccl4_raw <- as.matrix(dat_wide_ccl4_raw[, -1, with = FALSE])
rownames(mat_ccl4_raw) <- dat_wide_ccl4_raw$GeneProtein

non_empty_cols_ccl4_raw <- colnames(mat_ccl4_raw)[colSums(!is.na(mat_ccl4_raw)) > 0]
mat_ccl4_raw_valid <- mat_ccl4_raw[, non_empty_cols_ccl4_raw, drop = FALSE]
mat_ccl4_raw_clean <- na.omit(mat_ccl4_raw_valid)
print(paste("CCL4 Raw Proteins analyzed:", nrow(mat_ccl4_raw_clean)))

cor_mat_ccl4_raw <- cor(t(mat_ccl4_raw_clean), method = "pearson")
abs_cor_mat_ccl4_raw <- abs(cor_mat_ccl4_raw)
diag(abs_cor_mat_ccl4_raw) <- 0

centrality_stats_ccl4_raw <- data.table(
  GeneProtein = rownames(abs_cor_mat_ccl4_raw),
  Min_Abs_R   = apply(abs_cor_mat_ccl4_raw, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_ccl4_raw, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_ccl4_raw, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_ccl4_raw, 1, function(x) quantile(x, probs = 0.10))
)
centrality_stats_ccl4_raw <- merge(centrality_stats_ccl4_raw, protein_map_raw, by = "GeneProtein")
setnames(centrality_stats_ccl4_raw, "ProteinID", "Protein")

top_10_mean_ccl4_raw  <- centrality_stats_ccl4_raw[order(-Mean_Abs_R)][1:10]
top_10_80pct_ccl4_raw <- centrality_stats_ccl4_raw[order(-Thresh_80pct)][1:10]
top_10_90pct_ccl4_raw <- centrality_stats_ccl4_raw[order(-Thresh_90pct)][1:10]

top_10_mean  <- top_10_mean_ccl4_raw
top_10_80pct <- top_10_80pct_ccl4_raw
top_10_90pct <- top_10_90pct_ccl4_raw

ccl4_raw_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_CCL4_Raw.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = ccl4_raw_rdata_path)
print(paste("SUCCESS! CCL4 Raw Mastery Hubs saved to:", ccl4_raw_rdata_path))


# C. MERGED RAW HUBS
print("Calculating Merged raw centrality statistics...")
Merged_DT_raw <- DT_raw_core

dat_wide_merged_raw <- dcast(Merged_DT_raw, GeneProtein ~ MiceInfo, value.var = "ProteinIntensity")
mat_merged_raw <- as.matrix(dat_wide_merged_raw[, -1, with = FALSE])
rownames(mat_merged_raw) <- dat_wide_merged_raw$GeneProtein

non_empty_cols_merged_raw <- colnames(mat_merged_raw)[colSums(!is.na(mat_merged_raw)) > 0]
mat_merged_raw_valid <- mat_merged_raw[, non_empty_cols_merged_raw, drop = FALSE]
mat_merged_raw_clean <- na.omit(mat_merged_raw_valid)
print(paste("Merged Raw Proteins analyzed:", nrow(mat_merged_raw_clean)))

cor_mat_merged_raw <- cor(t(mat_merged_raw_clean), method = "pearson")
abs_cor_mat_merged_raw <- abs(cor_mat_merged_raw)
diag(abs_cor_mat_merged_raw) <- 0

centrality_stats_merged_raw <- data.table(
  GeneProtein = rownames(abs_cor_mat_merged_raw),
  Min_Abs_R   = apply(abs_cor_mat_merged_raw, 1, min),
  Mean_Abs_R  = apply(abs_cor_mat_merged_raw, 1, mean),
  Thresh_80pct = apply(abs_cor_mat_merged_raw, 1, function(x) quantile(x, probs = 0.20)),
  Thresh_90pct = apply(abs_cor_mat_merged_raw, 1, function(x) quantile(x, probs = 0.10))
)
centrality_stats_merged_raw <- merge(centrality_stats_merged_raw, protein_map_raw, by = "GeneProtein")
setnames(centrality_stats_merged_raw, "ProteinID", "Protein")

top_10_mean_merged_raw  <- centrality_stats_merged_raw[order(-Mean_Abs_R)][1:10]
top_10_80pct_merged_raw <- centrality_stats_merged_raw[order(-Thresh_80pct)][1:10]
top_10_90pct_merged_raw <- centrality_stats_merged_raw[order(-Thresh_90pct)][1:10]

top_10_mean  <- top_10_mean_merged_raw
top_10_80pct <- top_10_80pct_merged_raw
top_10_90pct <- top_10_90pct_merged_raw

merged_raw_rdata_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top10_Mastery_Hubs_Protein_Merged_Raw.RData"
save(top_10_mean, top_10_80pct, top_10_90pct, file = merged_raw_rdata_path)
print(paste("SUCCESS! Merged Raw Mastery Hubs saved to:", merged_raw_rdata_path))
