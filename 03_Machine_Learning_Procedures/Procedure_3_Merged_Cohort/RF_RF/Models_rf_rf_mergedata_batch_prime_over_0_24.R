#### RF LASSO MODEL WITH MERGED BATCH DATA - PRIME (OVERLAP ONLY) ##################
# DiPa Cluster 0
# This script runs the Prime modeling (943 candidate predictors) for DiPa Cluster 0
# on ComBat batch corrected data.
################################################################################

library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

setwd("/work/smbrngo1")

# 1. Load the full merged 3-groups batch-corrected dataset
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups_24.RData")
load("DTccl4_DT_LCPM_BatchCorrected_3_Groups_24.RData")
Full_DT <- as.data.table(Full_DT_24)

# Set ComBat corrected counts to active columns
Full_DT$GeneCount <- Full_DT$ComBat_GeneCount
Full_DT$ProteinIntensity <- Full_DT$ComBat_Protein_Raw

# Source the modeling functions
source("Function_protein_modelling_full_new_version_Meg.R")
source("Function_design_filtered.R")
source("Function_model_measures.R")


# OPTION C FILTERING: Extract target proteins for DiPa Cluster 0
cluster_DT <- Full_DT[ClusterDiPa == 0]
all_proteins <- unique(cluster_DT$GeneProtein)

treatments_ccl4 <- c("ccl4", "oil")
treatments_bdl <- c("control", "BDL", "BDL_ASBTi")

overlap_status <- sapply(all_proteins, function(p) {
  sub_p <- cluster_DT[GeneProtein == p]
  has_ccl4 <- any(!is.na(sub_p[Treatment %in% treatments_ccl4, ProteinIntensity]))
  has_bdl <- any(!is.na(sub_p[Treatment %in% treatments_bdl, ProteinIntensity]))
  return(has_ccl4 && has_bdl)
})
pairs.list.final <- all_proteins[overlap_status]
cat("Number of target proteins in overlap (Option C) :", length(pairs.list.final), "\n")


rf.rf.pre.full.mergedata_blind_protein_batch_prime_over_0_24 <- pblapply(FUN = protein.regression.complete,
                                                                            pairs.list.final,
                                                                            dataset = Full_DT,
                                                                            weighted = FALSE,
                                                                            grouping = "all",
                                                                            design.m = "proteins",
                                                                            scaled = FALSE,
                                                                            RNA.log.scale = FALSE,
                                                                            na.process = "max.n",
                                                                            include.treatment = FALSE,
                                                                            treatment.info = "full",
                                                                            duration.scale = "factor",
                                                                            PE.RNA.penalty = FALSE,
                                                                            treatment.penalty = TRUE,
                                                                            baseline.interactions = FALSE,
                                                                            prediction = TRUE,
                                                                            model.method = "rf.rf.pre",
                                                                            prediction.method = "loo",
                                                                            lasso.fam = "gaussian",
                                                                            type.measure = "deviance",
                                                                            alpha = 1,
                                                                            intercept = TRUE,
                                                                            nfolds = 10,
                                                                            top.n = 30,
                                                                            enforce.treatment = FALSE,
                                                                            tune.ranger = FALSE,
                                                                            manual.weights = FALSE,
                                                                            weights.vector = NULL,
                                                                            seed.set = 123)



save(rf.rf.pre.full.mergedata_blind_protein_batch_prime_over_0_24,
     file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP/RF_RF/Models_rf_rf_preselection_full_objects_mergedata_blind_batch_prime_over_0_24.RData")



