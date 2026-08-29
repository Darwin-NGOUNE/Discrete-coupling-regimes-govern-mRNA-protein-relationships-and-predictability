library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_analysis/DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_analysis/DTccl4_filtered_final_Gene_Protein_full_Batch.RData" )

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

# 2. Extract the target proteins directly from the split datasets using ClusterDiPa
prot_bdl_cluster <- unique(BDL_Full[ClusterDiPa %in% c(3, 4), GeneProtein])
prot_ccl4_cluster <- unique(CCL4_Full[ClusterDiPa %in% c(3, 4), GeneProtein])
core_cluster_proteins <- intersect(prot_bdl_cluster, prot_ccl4_cluster)

cat("Number of target proteins (DiPa Cluster 3_4):", length(core_cluster_proteins), "\n")
cat("Number of candidate predictors in dataset:", length(core_full_proteins), "\n")

# Target proteins to run:
pairs.list <- core_cluster_proteins

# 3. Load Mastery Hub Proteins (CCL4)
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top50_Mastery_Hubs_Protein_CCL4_Batch.RData")
# Extract exactly the top 10 from 80% coverage
top_50_data <- top_50_80pct[order(top_50_80pct$Thresh_80pct, decreasing = TRUE), ][1:50, ]
top50_proteins <- str_split(top_50_data$GeneProtein, "_", simplify = TRUE)[, 2]
my_top_50_list <- paste0(top50_proteins, "_P")

cat("Selected top 10 mastery hub covariates:\n")
print(my_top_50_list)

# 4. Set working directory to parent for sourced relative paths
old_wd <- getwd()
setwd("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing")

# Source modeling functions
source("Function_protein_modelling_testing_full.R")
source("Function_design.R")
source("Function_design_newdata_filtered.R")
source("Function_model_measures.R")

# 5. Run Lasso model with predefined mastery proteins
cat("\n--> Running Mastery 10 Protein Lasso...")
rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub <- pblapply(FUN = protein.regression.function,
                                                           pairs.list,
                                                           dataset   = CCL4_Full,
                                                           newdata   = BDL_Full,
                                                           weighted = FALSE,
                                                           grouping = "all",
                                                           scaled = FALSE,
                                                           RNA.log.scale = FALSE,
                                                           design.m = "proteins", 
                                                           manual.covar = my_top_50_list,  # <-- MASTERY 10
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
                                                           top.n = 50,
                                                           enforce.treatment = FALSE,
                                                           tune.ranger = FALSE,
                                                           model.method = "lasso.predefined",
                                                           alpha = 1,
                                                           lasso.fam = "gaussian",
                                                           nfolds = 10,
                                                           type.measure = "deviance",
                                                           intercept = TRUE,
                                                           seed.set = 123)

# 5b. Run Lasso model with all RNA + predefined mastery proteins (RNA design + 10 mastery)
cat("\n--> Running Mastery 10 Protein + All RNA Lasso...")
rf_lasso_mastery50_plus_rna_CCL4_BDL_3_4_prime_sub <- pblapply(FUN = protein.regression.function,
                                                           pairs.list,
                                                           dataset   = CCL4_Full,
                                                           newdata   = BDL_Full,
                                                           weighted = FALSE,
                                                           grouping = "all",
                                                           scaled = FALSE,
                                                           RNA.log.scale = FALSE,
                                                           design.m = "all", 
                                                           manual.covar = my_top_50_list,  # <-- MASTERY 10
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
                                                           tune.ranger = FALSE,
                                                           model.method = "lasso.rf.pre",
                                                           alpha = 1,
                                                           lasso.fam = "gaussian",
                                                           nfolds = 10,
                                                           type.measure = "deviance",
                                                           intercept = TRUE,
                                                           seed.set = 123)

# Restore working directory
setwd(old_wd)

# Save the models
output_file <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/Batch_corrected_prime_analysis/Mastery50/Models_mastery50_CCL4_BDL_3_4_prime_sub.RData"
save(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub,
     rf_lasso_mastery50_plus_rna_CCL4_BDL_3_4_prime_sub,
     file = output_file)
cat("\nModels saved successfully to:", output_file, "\n")

# 6. Analyze and print the result summary
rmse_bas <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$RMSE)
mae_bas <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$MAE)
cor_bas  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$Pearson)

rmse_mod <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$rmse)
mae_mod <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$mae)
cor_mod  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$correlation.pearson)


# Clean analysis and printing of results
rmse_bas <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$RMSE)
mae_bas  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$MAE)
cor_bas  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$baseline.model$Pearson)

rmse_mod <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$rmse)
mae_mod  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$mae)
cor_mod  <- sapply(rf_lasso_mastery50_proteins_CCL4_BDL_3_4_prime_sub, function(x) x$model$correlation.pearson)

results_table <- data.frame(
  Design = c("Baseline", "Mastery50_Proteins"),
  Mean_RMSE = c(mean(rmse_bas, na.rm = TRUE), mean(rmse_mod, na.rm = TRUE)),
  Mean_MAE = c(mean(mae_bas, na.rm = TRUE), mean(mae_mod, na.rm = TRUE)),
  Mean_Pearson = c(mean(cor_bas, na.rm = TRUE), mean(cor_mod, na.rm = TRUE))
)

cat("\n=== SUMMARY TABLE (CCL4_BDL_Mastery50_3_4_prime_sub.R) ===\n")
print(results_table)

cor_vals <- cor_mod
cat("\nCorrelations summary:\n")
cat("sum(cor_vals < 0.5):", sum(cor_vals < 0.5, na.rm = TRUE), "\n")
cat("sum(cor_vals >= 0.5):", sum(cor_vals >= 0.5, na.rm = TRUE), "\n")
cat("sum(cor_vals >= 0.8):", sum(cor_vals >= 0.8, na.rm = TRUE), "\n")
