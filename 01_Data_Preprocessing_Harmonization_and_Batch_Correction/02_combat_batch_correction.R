# ==============================================================================
# SCRIPT: Batch_Effect_Correction_OilControl.R
#
# IMPLEMENTATION OF BATCH CORRECTION WITH A DIFFERENT CONTROL DEFINITION
#
# HYPOTHESIS TEST: Prof. Hengstler's idea. 
# "What if we use the mice that ACTUALLY received the Oil injection (Months 2 & 12) 
# as the CCl4 Control group for Batch Correction, rather than the Naive Month 0 mice?"
#
# This script applies the same ComBat correction but reassigns the "DiseaseGroup" 
# metadata exactly according to this hypothesis.
# ==============================================================================

library(data.table)

# -------------------------------------------------------------------------
# 1. LOAD DATA & DEFINE CONTROLS
# -------------------------------------------------------------------------
file_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData"
load(file_path)
Full_DT <- DTccl4_DT_LCPM
setDT(Full_DT)

Full_DT[, Protein_Raw := ProteinIntensity]
unique_genes <- unique(Full_DT$GeneProtein)

# Define Batch Column Robustly
if ("Dataset" %in% names(Full_DT)) {
    Full_DT[, Batch := Dataset]
} else if ("experiment" %in% names(Full_DT)) {
    Full_DT[, Batch := experiment]
    Full_DT[, Dataset := experiment]
} else {
    Full_DT[, Batch := "Unknown"]
    Full_DT[Treatment %in% c("control", "BDL","BDL_ASBTi"), Batch := "BDL Dataset"]
    Full_DT[Treatment %in% c("oil", "ccl4") | TreatmentTime == 0, Batch := "CCl4 Dataset"]
    Full_DT[, Dataset := Batch]
}

# Define Disease Groups (Crucial for all 3 methods)
Full_DT[, DiseaseGroup := "Unknown"]

# ----------- PRECISE GROUP ASSIGNMENT (THE HYPOTHESIS SWITCH) ------------

# 1. BDL Dataset
Full_DT[Batch == "BDL Dataset" & as.character(Treatment) == "control", DiseaseGroup := "Control"]
Full_DT[Batch == "BDL Dataset" & as.character(Treatment) == "BDL", DiseaseGroup := "BDL"]
Full_DT[Batch == "BDL Dataset" & as.character(Treatment) == "BDL_ASBTi", DiseaseGroup := "BDL_ASBTi"]

# 2. CCl4 Dataset
# -> Control is Month 0 Naive
Full_DT[Batch == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) == 0 & as.character(Treatment) == "oil", DiseaseGroup := "Control"]

# -> Disease is ccl4 at times > 0
Full_DT[Batch == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) > 0 & as.character(Treatment) == "ccl4", DiseaseGroup := "Ccl4_Month_2_6_12"]

# -> Remedy is oil at times > 0 (Month 2, 12, etc.)
Full_DT[Batch == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) > 0 & as.character(Treatment) == "oil", DiseaseGroup := "Oil_Month_2_12"]

print("Data Loaded and NEW Disease Groups defined (Control, Disease, Remedy).")


# =========================================================================
# STRATEGY 3: CONTROL-CENTERING
# =========================================================================
print("Executing STRATEGY 3: Control-Centering...")
Control_Baselines <- Full_DT[DiseaseGroup == "Control", .(
    Baseline_RNA = mean(GeneCount, na.rm = TRUE),
    Baseline_Prot = mean(Protein_Raw, na.rm = TRUE)
), by = .(GeneProtein, Batch)]

Full_DT <- merge(Full_DT, Control_Baselines, by = c("GeneProtein", "Batch"), all.x = TRUE)
Full_DT[, Norm_GeneCount := GeneCount - Baseline_RNA]
Full_DT[, Norm_Protein_Raw := Protein_Raw - Baseline_Prot]


# =========================================================================
# PREPARATION FOR COMBAT: WIDE MATRIX CREATION
# =========================================================================
print("Preparing Wide Matrices for ComBat...")

if (!"Sample_ID" %in% names(Full_DT)) {
    Full_DT[, Sample_ID := MiceInfo]
}

# Create Wide Matrices
Wide_Matrix_RNA <- dcast(Full_DT[!is.na(GeneCount)], GeneProtein ~ Sample_ID, value.var = "GeneCount", fun.aggregate = mean)
Matrix_RNA <- as.matrix(Wide_Matrix_RNA[, -1, with = FALSE])
rownames(Matrix_RNA) <- Wide_Matrix_RNA$GeneProtein

Wide_Matrix_Prot <- dcast(Full_DT[!is.na(Protein_Raw)], GeneProtein ~ Sample_ID, value.var = "Protein_Raw", fun.aggregate = mean)
Matrix_Prot <- as.matrix(Wide_Matrix_Prot[, -1, with = FALSE])
rownames(Matrix_Prot) <- Wide_Matrix_Prot$GeneProtein

# Create Metadata
colData <- unique(Full_DT[, .(Sample_ID, DiseaseGroup, Batch)])
colData_RNA <- colData[match(colnames(Matrix_RNA), colData$Sample_ID)]
colData_Prot <- colData[match(colnames(Matrix_Prot), colData$Sample_ID)]
rownames(colData_RNA) <- colData_RNA$Sample_ID
rownames(colData_Prot) <- colData_Prot$Sample_ID


# =========================================================================
# STRATEGY 2: COMBAT (Empirical Bayes via 'sva')
# =========================================================================
print("Executing STRATEGY 2: ComBat Setup (Oil as Control)")

if (require(sva, quietly = TRUE)) {
    tryCatch({
        # 1. Filter Metadata to keep Control, Disease, and Remedy
        colData_RNA_cb <- colData_RNA[DiseaseGroup %in% c("Control", "BDL", "BDL_ASBTi", "Ccl4_Month_2_6_12", "Oil_Month_2_12"), ]
        colData_Prot_cb <- colData_Prot[DiseaseGroup %in% c("Control", "BDL", "BDL_ASBTi", "Ccl4_Month_2_6_12", "Oil_Month_2_12"), ]

        colData_RNA_cb$DiseaseGroup <- droplevels(as.factor(colData_RNA_cb$DiseaseGroup))
        colData_Prot_cb$DiseaseGroup <- droplevels(as.factor(colData_Prot_cb$DiseaseGroup))

        # 2. Filter Matrices
        Matrix_RNA_cb <- Matrix_RNA[, colData_RNA_cb$Sample_ID, drop = FALSE]
        Matrix_Prot_cb <- Matrix_Prot[, colData_Prot_cb$Sample_ID, drop = FALSE]

        # 3. Process RNA
        mod_rna <- model.matrix(~DiseaseGroup, data = colData_RNA_cb)
        Matrix_RNA_ComBat <- ComBat(dat = Matrix_RNA_cb, batch = colData_RNA_cb$Batch, mod = mod_rna, par.prior = TRUE, prior.plots = FALSE)

        # 4. Process Protein
        mod_prot <- model.matrix(~DiseaseGroup, data = colData_Prot_cb)
        Matrix_Prot_ComBat <- ComBat(dat = Matrix_Prot_cb, batch = colData_Prot_cb$Batch, mod = mod_prot, par.prior = TRUE, prior.plots = FALSE)

        # 5. Melt and merge back
        dt_rna_cb <- as.data.table(Matrix_RNA_ComBat, keep.rownames = "GeneProtein")
        dt_rna_cb <- melt(dt_rna_cb, id.vars = "GeneProtein", variable.name = "Sample_ID", value.name = "ComBat_GeneCount")

        dt_prot_cb <- as.data.table(Matrix_Prot_ComBat, keep.rownames = "GeneProtein")
        dt_prot_cb <- melt(dt_prot_cb, id.vars = "GeneProtein", variable.name = "Sample_ID", value.name = "ComBat_Protein_Raw")

        Full_DT <- merge(Full_DT, dt_rna_cb, by = c("GeneProtein", "Sample_ID"), all.x = TRUE)
        Full_DT <- merge(Full_DT, dt_prot_cb, by = c("GeneProtein", "Sample_ID"), all.x = TRUE)

        print("ComBat execution complete.")
    }, error = function(e) { print(paste("ComBat Error:", e$message)) })
} else {
    print("Skipping ComBat: 'sva' not installed.")
}

# -------------------------------------------------------------------------
# SAVE OUTPUTS (Under new names to prevent overwriting original)
# -------------------------------------------------------------------------
out_data_path <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData"
save(Full_DT, file = out_data_path)



# -------------------------------------------------------------------------
# VERIFICATION: mRNA-PROTEIN CORRELATION BEFORE VS AFTER COMBAT
# -------------------------------------------------------------------------
cat("\n=== Verifying mRNA-Protein Correlations (Before vs After ComBat) ===\n")

cor_check <- Full_DT[!is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(ComBat_GeneCount) & !is.na(ComBat_Protein_Raw), .(
  r_before = cor(GeneCount, Protein_Raw, use = "pairwise.complete.obs"),
  r_after  = cor(ComBat_GeneCount, ComBat_Protein_Raw, use = "pairwise.complete.obs")
), by = GeneProtein]

cor_check[, delta_r := r_after - r_before]

cat(sprintf("Average correlation BEFORE ComBat : %.4f\n", mean(cor_check$r_before, na.rm = TRUE)))
cat(sprintf("Average correlation AFTER ComBat  : %.4f\n", mean(cor_check$r_after, na.rm = TRUE)))
cat(sprintf("Mean change in correlation (delta): %.4f\n", mean(cor_check$delta_r, na.rm = TRUE)))
cat(sprintf("Proportion of proteins with dropped correlation: %.2f%%\n", 
            mean(cor_check$delta_r < 0, na.rm = TRUE) * 100))






