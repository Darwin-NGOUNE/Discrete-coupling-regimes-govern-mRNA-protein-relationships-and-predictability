#### RF LASSO MODEL WITH MERGED BATCH DATA - PRIME (OVERLAP ONLY) ##################
# DiPa Cluster 3_4
# This script runs the Prime modeling (943 candidate predictors) for DiPa Cluster 3_4
# on ComBat batch corrected data.
################################################################################

library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

setwd("/work/smbrngo1")

# 1. Load the full merged 3-groups batch-corrected dataset
# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups_24.RData")
load("DTccl4_DT_LCPM_BatchCorrected_3_Groups_24.RData")
Full_DT <- as.data.table(Full_DT_24)

# Set ComBat corrected counts to active columns
Full_DT$GeneCount <- Full_DT$ComBat_GeneCount
Full_DT$ProteinIntensity <- Full_DT$ComBat_Protein_Raw

# Source the modeling functions
source("Function_protein_modelling_full_new_version_Meg.R")
source("Function_design_filtered.R")
source("Function_model_measures.R")

# OPTION C FILTERING: Extract target proteins for DiPa Cluster 3_4
cluster_DT <- Full_DT[ClusterDiPa %in% c(3, 4)]
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

# RUN MODELS (RNA, Protein, Combined)
lasso.rf.pre.full.mergedata_blind_rna_batch_prime_over_3_4_24 <- pblapply(FUN = protein.regression.complete,
                                                           pairs.list.final,
                                                           dataset = Full_DT,
                                                           weighted = FALSE,
                                                           grouping = "all",
                                                           design.m = "genes",
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
                                                           model.method = "lasso.rf.pre",
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

lasso.rf.pre.full.mergedata_blind_protein_batch_prime_over_3_4_24 <- pblapply(FUN = protein.regression.complete,
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
                                                           model.method = "lasso.rf.pre",
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

lasso.rf.pre.full.mergedata_blind_combined_batch_prime_over_3_4_24 <- pblapply(FUN = protein.regression.complete,
                                                           pairs.list.final,
                                                           dataset = Full_DT,
                                                           weighted = FALSE,
                                                           grouping = "all",
                                                           design.m = "all",
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
                                                           model.method = "lasso.rf.pre",
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

# SAVE RESULTS
save(lasso.rf.pre.full.mergedata_blind_rna_batch_prime_over_3_4_24,
     lasso.rf.pre.full.mergedata_blind_protein_batch_prime_over_3_4_24,
     lasso.rf.pre.full.mergedata_blind_combined_batch_prime_over_3_4_24,
     file = "Models_rf_preselection_full_objects_mergedata_blind_batch_prime_over_3_4_24.RData")

# ==============================================================================
# 4. EVALUATION AND SUMMARY STATISTICS
# ==============================================================================
# # load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP/Models_rf_preselection_full_objects_mergedata_blind_batch_prime_over_3_4_24.RData")
# 
# res.rna <- gof.all.complete(model = lasso.rf.pre.full.mergedata_blind_rna_batch_prime_over_3_4_24, prediction.method = "lasso")
# res.prot <- gof.all.complete(model = lasso.rf.pre.full.mergedata_blind_protein_batch_prime_over_3_4_24, prediction.method = "lasso")
# res.comb <- gof.all.complete(model = lasso.rf.pre.full.mergedata_blind_combined_batch_prime_over_3_4_24, prediction.method = "lasso")
# 
# cat("\n=== RNA CORRELATIONS ===\n")
# cor_rna <- sapply(lasso.rf.pre.full.mergedata_blind_rna_batch_prime_over_3_4_24, function(x) as.numeric(x$prediction.obj$model$correlation.pearson))
# print(summary(cor_rna))
# cat("Count < 0.5:", sum(cor_rna < 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.5:", sum(cor_rna >= 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.8:", sum(cor_rna >= 0.8, na.rm=TRUE), "\n")
# 
# cat("\n=== PROTEIN CORRELATIONS ===\n")
# cor_prot <- sapply(lasso.rf.pre.full.mergedata_blind_protein_batch_prime_over_3_4_24, function(x) as.numeric(x$prediction.obj$model$correlation.pearson))
# print(summary(cor_prot))
# cat("Count < 0.5:", sum(cor_prot < 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.5:", sum(cor_prot >= 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.8:", sum(cor_prot >= 0.8, na.rm=TRUE), "\n")
# 
# cat("\n=== COMBINED CORRELATIONS ===\n")
# cor_comb <- sapply(lasso.rf.pre.full.mergedata_blind_combined_batch_prime_over_3_4_24, function(x) as.numeric(x$prediction.obj$model$correlation.pearson))
# print(summary(cor_comb))
# cat("Count < 0.5:", sum(cor_comb < 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.5:", sum(cor_comb >= 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.8:", sum(cor_comb >= 0.8, na.rm=TRUE), "\n")

