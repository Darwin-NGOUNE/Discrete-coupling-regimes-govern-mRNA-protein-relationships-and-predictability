library(data.table)
library(writexl)

cat("=== Building Jan's Requested Supplementary Excel Table ===\n")

# 1. Load Data
e <- new.env()
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData", envir = e)
dt <- as.data.table(e$Full_DT)

# 2. Define ordered 54 mice metadata
# BDL cohort (18 mice)
bdl_sham <- paste0("ShamvehicleM", c(1, 3, 4, 6, 7, 8))
bdl_veh  <- paste0("BDLvehicleM", c(5, 7, 10, 11, 12, 15))
bdl_as0  <- paste0("BDLASBTiM", c(8, 9, 10, 11, 12, 14))

# CCl4 cohort (36 mice)
ccl4_m0_oil  <- paste0("month0_oil_rep", 1:6)
ccl4_m2_oil  <- paste0("month2_oil_rep", 1:6)
ccl4_m2_ccl4 <- paste0("month2_ccl4_rep", 1:6)
ccl4_m6_ccl4 <- paste0("month6_ccl4_rep", 1:6)
ccl4_m12_oil <- paste0("month12_oil_rep", 1:6)
ccl4_m12_ccl4 <- paste0("month12_ccl4_rep", 1:6)

ordered_sample_ids <- c(
  bdl_sham, bdl_veh, bdl_as0,
  ccl4_m0_oil, ccl4_m2_oil, ccl4_m2_ccl4,
  ccl4_m6_ccl4, ccl4_m12_oil, ccl4_m12_ccl4
)

treatment_labels <- c(
  rep("Sham control", 6),
  rep("BDL", 6),
  rep("BDL + ASBTi (AS0369)", 6),
  rep("0 months, oil control", 6),
  rep("2 months, oil control", 6),
  rep("2 months, CCl4", 6),
  rep("6 months, CCl4", 6),
  rep("12 months, oil control", 6),
  rep("12 months, CCl4", 6)
)

intervention_groups <- c(
  rep("BDL", 18),
  rep("CCl4", 36)
)

mice_metadata_df <- data.table(
  `Mouse number` = 1:54,
  `Intervention group` = intervention_groups,
  `Treatment` = treatment_labels,
  `Sample_ID` = ordered_sample_ids
)

cat("Mice count:", nrow(mice_metadata_df), "\n")

# 3. Reshape GeneProtein columns (Wide format)
# We have 943 unique GeneProtein pairs
unique_pairs <- unique(dt[, .(GeneProtein, GeneSyn, ProteinID)])
setorder(unique_pairs, GeneSyn, ProteinID)

cat("Number of unique Gene-Protein pairs:", nrow(unique_pairs), "\n")

# Create wide matrices for each of the 4 variables
d_prot_raw <- dcast(dt, Sample_ID ~ GeneSyn, value.var = "Protein_Raw")
d_prot_cb  <- dcast(dt, Sample_ID ~ GeneSyn, value.var = "ComBat_Protein_Raw")
d_rna_raw  <- dcast(dt, Sample_ID ~ GeneSyn, value.var = "GeneCount")
d_rna_cb   <- dcast(dt, Sample_ID ~ GeneSyn, value.var = "ComBat_GeneCount")

# Reorder rows according to ordered_sample_ids
setkey(d_prot_raw, Sample_ID)
setkey(d_prot_cb, Sample_ID)
setkey(d_rna_raw, Sample_ID)
setkey(d_rna_cb, Sample_ID)

d_prot_raw <- d_prot_raw[ordered_sample_ids]
d_prot_cb  <- d_prot_cb[ordered_sample_ids]
d_rna_raw  <- d_rna_raw[ordered_sample_ids]
d_rna_cb   <- d_rna_cb[ordered_sample_ids]

# Interleave the 4 columns per GeneSyn
data_cols_list <- list()

for (g in unique_pairs$GeneSyn) {
  c_p_raw <- paste0("Protein ", g, ", without batch correction")
  c_p_cb  <- paste0("Protein ", g, ", with batch correction")
  c_r_raw <- paste0("mRNA ", g, ", without batch correction")
  c_r_cb  <- paste0("mRNA ", g, ", with batch correction")
  
  data_cols_list[[c_p_raw]] <- d_prot_raw[[g]]
  data_cols_list[[c_p_cb]]  <- d_prot_cb[[g]]
  data_cols_list[[c_r_raw]] <- d_rna_raw[[g]]
  data_cols_list[[c_r_cb]]  <- d_rna_cb[[g]]
}

data_matrix_dt <- as.data.table(data_cols_list)

# Combine metadata and data columns
final_wide_table <- cbind(
  mice_metadata_df[, .(`Mouse number`, `Intervention group`, `Treatment`)],
  data_matrix_dt
)

cat("Final Table Dimensions: ", nrow(final_wide_table), "rows x", ncol(final_wide_table), "columns\n")

# 4. Save Excel files safely
out_file1 <- "C:/Users/ngoune/Documents/Paper I/Supplementary_Table_1_RNA_Protein_Pairs.xlsx"
out_file2 <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/01_Data_Harmonization_and_Batch_Correction/Supplementary_Table_1_RNA_Protein_Pairs.xlsx"

safe_save <- function(df, path) {
  tryCatch({
    write_xlsx(df, path = path)
    cat("Successfully written:", path, "\n")
  }, error = function(e) {
    alt_path <- gsub("\\.xlsx$", "_unrounded.xlsx", path)
    cat("Warning: Original file is open/locked in Excel. Saving as:", alt_path, "\n")
    write_xlsx(df, path = alt_path)
    cat("Successfully written:", alt_path, "\n")
  })
}

safe_save(final_wide_table, out_file1)
safe_save(final_wide_table, out_file2)

cat("\nSUCCESS: Excel tables with exact unrounded values generated!\n")
