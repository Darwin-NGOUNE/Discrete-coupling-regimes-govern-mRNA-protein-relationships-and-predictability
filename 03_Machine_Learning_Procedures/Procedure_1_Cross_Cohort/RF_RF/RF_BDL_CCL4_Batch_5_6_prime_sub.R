library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)


setwd("/work/smbrngo1")

# 1. Load the 3-groups batch-corrected dataset (Full_DT)
# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_analysis/DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData")
# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_analysis/DTccl4_filtered_final_Gene_Protein_full_Batch.RData" )

load("DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData")
load("DTccl4_filtered_final_Gene_Protein_full_Batch.RData")

BDL_DT <- as.data.table(BDL_DT)
CCL4_DTccl4 <- as.data.table(CCL4_DTccl4)

# Separate into BDL_DT and CCL4_DTccl4 (Subset Cohort)
BDL_DT <- BDL_DT[Treatment %in% c("BDL", "control")]
CCL4_DTccl4 <- CCL4_DTccl4[(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)]

# Set ComBat corrected counts to active columns
BDL_DT$GeneCount <- BDL_DT$ComBat_GeneCount
BDL_DT$ProteinIntensity <- BDL_DT$ComBat_Protein_Raw

CCL4_DTccl4$GeneCount <- CCL4_DTccl4$ComBat_GeneCount
CCL4_DTccl4$ProteinIntensity <- CCL4_DTccl4$ComBat_Protein_Raw

# Find the 943 core overlap proteins for the predictor space
prot_full_bdl <- unique(BDL_DT$GeneProtein)
prot_full_ccl4 <- unique(CCL4_DTccl4$GeneProtein)
core_full_proteins <- intersect(prot_full_bdl, prot_full_ccl4)

BDL_Full <- BDL_DT[GeneProtein %in% core_full_proteins]
CCL4_Full <- CCL4_DTccl4[GeneProtein %in% core_full_proteins]

# 2. Extract the 118 target proteins directly from the split datasets using ClusterDiPa
prot_bdl_5_6 <- unique(BDL_Full[ClusterDiPa %in% c(5, 6), GeneProtein])
prot_ccl4_5_6 <- unique(CCL4_Full[ClusterDiPa %in% c(5, 6), GeneProtein])
core_5_6_proteins <- intersect(prot_bdl_5_6, prot_ccl4_5_6)

cat("Number of target proteins (DiPa Groups 5 & 6):", length(core_5_6_proteins), "\n")
cat("Number of candidate predictors in dataset:", length(core_full_proteins), "\n")

# Target proteins to run:
pairs.list <- core_5_6_proteins

# 3. Set working directory to parent for sourced relative paths
# old_wd <- getwd()
# setwd("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing")

################################################################################

# Source modeling functions
source("Function_protein_modelling_testing_full_new.R")
source("Function_design.R")
source("Function_design_newdata_filtered.R")
source("Function_model_measures.R")

# 4. Run Lasso + RF preselection modeling
# cat("\n--> Running RNA-only Lasso (preselected)...")
# rf_rf_full_rna_testing_new_Batch_prime_5_6_sub <- pblapply(FUN = protein.regression.function,
#                                                               pairs.list,
#                                                               dataset   = BDL_Full,
#                                                               newdata   = CCL4_Full,
#                                                               weighted = FALSE,
#                                                               grouping = "all",
#                                                               scaled = FALSE,
#                                                               RNA.log.scale = FALSE,
#                                                               design.m = "genes", 
#                                                               na.process = "max.n",
#                                                               include.treatment = FALSE,
#                                                               treatment.info = "full",
#                                                               duration.scale = "factor",
#                                                               treatment.penalty = FALSE, 
#                                                               PE.RNA.penalty = FALSE,
#                                                               baseline.interactions = FALSE,
#                                                               collinearity.reduce = FALSE,
#                                                               rho.cutoff = 0.75,
#                                                               manual.weights = FALSE,
#                                                               weights.vector = NULL,
#                                                               top.n = 30,
#                                                               enforce.treatment = FALSE,
#                                                               tune.ranger = FALSE,
#                                                               model.method = "rf.rf.pre",
#                                                               alpha = 1,
#                                                               lasso.fam = "gaussian",
#                                                               nfolds = 10,
#                                                               type.measure = "deviance",
#                                                               intercept = TRUE,
#                                                               seed.set = 123)

cat("\n--> Running Protein-only Lasso (preselected)...")
rf_rf_full_proteins_testing_new_Batch_prime_5_6_sub  <- pblapply(FUN = protein.regression.function,
                                                                    pairs.list,
                                                                    dataset   = BDL_Full,
                                                                    newdata   = CCL4_Full,
                                                                    weighted = FALSE,
                                                                    grouping = "all",
                                                                    scaled = FALSE,
                                                                    RNA.log.scale = FALSE,
                                                                    design.m = "proteins", 
                                                                    na.process = "max.n",
                                                                    include.treatment = FALSE,
                                                                    treatment.info = "full",
                                                                    duration.scale = "factor",
                                                                    treatment.penalty = FALSE, 
                                                                    PE.RNA.penalty = FALSE,
                                                                    baseline.interactions = FALSE,
                                                                    collinearity.reduce = FALSE,
                                                                    rho.cutoff = 0.75,
                                                                    manual.weights = FALSE,
                                                                    weights.vector = NULL,
                                                                    top.n = 30,
                                                                    enforce.treatment = FALSE,
                                                                    tune.ranger = TRUE,
                                                                    model.method = "rf.rf.pre",
                                                                    alpha = 1,
                                                                    lasso.fam = "gaussian",
                                                                    nfolds = 10,
                                                                    type.measure = "deviance",
                                                                    intercept = TRUE,
                                                                    seed.set = 123)

# cat("\n--> Running Combined Lasso (preselected)...")
# rf_rf_full_combined_testing_new_Batch_prime_5_6_sub  <- pblapply(FUN = protein.regression.function,
#                                                                     pairs.list,
#                                                                     dataset   = BDL_Full,
#                                                                     newdata   = CCL4_Full,
#                                                                     weighted = FALSE,
#                                                                     grouping = "all",
#                                                                     scaled = FALSE,
#                                                                     RNA.log.scale = FALSE,
#                                                                     design.m = "all", 
#                                                                     na.process = "max.n",
#                                                                     include.treatment = FALSE,
#                                                                     treatment.info = "full",
#                                                                     duration.scale = "factor",
#                                                                     treatment.penalty = FALSE, 
#                                                                     PE.RNA.penalty = FALSE,
#                                                                     baseline.interactions = FALSE,
#                                                                     collinearity.reduce = FALSE,
#                                                                     rho.cutoff = 0.75,
#                                                                     manual.weights = FALSE,
#                                                                     weights.vector = NULL,
#                                                                     top.n = 30,
#                                                                     enforce.treatment = FALSE,
#                                                                     tune.ranger = FALSE,
#                                                                     model.method = "rf.rf.pre",
#                                                                     alpha = 1,
#                                                                     lasso.fam = "gaussian",
#                                                                     nfolds = 10,
#                                                                     type.measure = "deviance",
#                                                                     intercept = TRUE,
#                                                                     seed.set = 123)

# Restore working directory
# setwd(old_wd)
# 
# # Save the models
# output_file <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/RF_RF/Models_rf_rf_full_testing_new_Batch_prime_5_6_sub.RData"
# # save(rf_rf_full_rna_testing_new_Batch_prime_5_6_sub ,
# #      rf_rf_full_proteins_testing_new_Batch_prime_5_6_sub ,
# #      rf_rf_full_combined_testing_new_Batch_prime_5_6_sub ,
# #      file = output_file)
# # cat("\nModels saved successfully to:", output_file, "\n")
# save(rf_rf_full_proteins_testing_new_Batch_prime_5_6_sub ,
#      file = output_file)
# cat("\nModels saved successfully to:", output_file, "\n")

save(rf_rf_full_proteins_testing_new_Batch_prime_5_6_sub,
     file = "Models_rf_rf_rf_full_testing_new_Batch_prime_5_6_sub.RData"
)


################################################################################
##### Analysing the result of the modeling RF PRESELECTION
################################################################################

# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/Models_rf_lasso_full_testing_new_Batch_prime_5_6_sub.RData")
# 
# 
# ######################rna#######################################################
# 
# cor_vals <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub, function(x) {
#   as.numeric(x$model$correlation.pearson)
# })
# 
# 
# sum(cor_vals < 0.5) #18
# sum(cor_vals >= 0.5) #5
# sum(cor_vals >= 0.8) #1
# 
# 
# 
# 
# cor_vals1 <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub, function(x) {
#   as.numeric(x$baseline.model$Pearson)
# })
# 
# sum(cor_vals1 < 0.5) #21
# sum(cor_vals1 >= 0.5) #2
# sum(cor_vals1 >= 0.8) #0
# 
# ################proteins########################################################
# 
# cor_vals <- sapply(rf_lasso_full_proteins_testing_new_Batch_prime_5_6_sub, function(x) {
#   as.numeric(x$model$correlation.pearson)
# })
# 
# sum(cor_vals < 0.5) #19
# sum(cor_vals >= 0.5) #4
# sum(cor_vals >= 0.8) #1
# 
# ##################combined######################################################
# 
# cor_vals <- sapply(rf_lasso_full_combined_testing_new_Batch_prime_5_6_sub, function(x) {
#   as.numeric(x$model$correlation.pearson)
# })
# 
# sum(cor_vals < 0.5) #18
# sum(cor_vals >= 0.5) #5
# sum(cor_vals >= 0.8) #1
# 
# ################################################################################
# 
# # 5. Analyze and print the result summary
# rmse_bas <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$baseline.model$RMSE)
# mae_bas <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$baseline.model$MAE)
# cor_bas  <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$baseline.model$Pearson)
# 
# rmse_rna <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$model$rmse)
# mae_rna <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$model$mae)
# cor_rna  <- sapply(rf_lasso_full_rna_testing_new_Batch_prime_5_6_sub , function(x) x$model$correlation.pearson)
# 
# rmse_prot <- sapply(rf_lasso_full_proteins_testing_new_Batch_prime_5_6_sub , function(x) x$model$rmse)
# mae_prot <- sapply(rf_lasso_full_proteins_testing_new_Batch_prime_5_6_sub , function(x) x$model$mae)
# cor_prot  <- sapply(rf_lasso_full_proteins_testing_new_Batch_prime_5_6_sub , function(x) x$model$correlation.pearson)
# 
# rmse_comb <- sapply(rf_lasso_full_combined_testing_new_Batch_prime_5_6_sub , function(x) x$model$rmse)
# mae_comb <- sapply(rf_lasso_full_combined_testing_new_Batch_prime_5_6_sub , function(x) x$model$mae)
# cor_comb  <- sapply(rf_lasso_full_combined_testing_new_Batch_prime_5_6_sub , function(x) x$model$correlation.pearson)
# 
# results_table <- data.frame(
#   Design = c("Baseline", "RNA", "Protein", "Combined"),
#   Mean_RMSE = c(mean(rmse_bas, na.rm = TRUE), mean(rmse_rna, na.rm = TRUE), mean(rmse_prot, na.rm = TRUE), mean(rmse_comb, na.rm = TRUE)),
#   Mean_MAE = c(mean(mae_bas, na.rm = TRUE), mean(mae_rna, na.rm = TRUE), mean(mae_prot, na.rm = TRUE), mean(mae_comb, na.rm = TRUE)),
#   Mean_Pearson = c(mean(cor_bas, na.rm = TRUE), mean(cor_rna, na.rm = TRUE), mean(cor_prot, na.rm = TRUE), mean(cor_comb, na.rm = TRUE))
# )
# 
# cat("\n=== SUMMARY TABLE (PROCEDURE 1-PRIME - SUBSET COHORT) ===\n")
# print(results_table)
# 
# Design Mean_RMSE  Mean_MAE Mean_Pearson
# 1 Baseline 0.4608389 0.3611135  -0.10914522
# 2      RNA 0.7406871 0.5791999   0.03622782
# 3  Protein 0.5070503 0.4066799   0.04327138
# 4 Combined 0.5193121 0.4060925   0.11455512


