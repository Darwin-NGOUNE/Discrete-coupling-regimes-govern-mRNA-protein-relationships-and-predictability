#### RF LASSO MODEL WITH MERGED BATCH DATA - MASTERY 10 (OVERLAP ONLY) ##############
# DiPa Cluster 0
# This script runs the Mastery 10 modeling (predefined top 10 Mastery Hubs) for DiPa Cluster 0
# on ComBat batch corrected data.
################################################################################

library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

setwd("/work/smbrngo1")

# 1. Load the full merged 3-groups batch-corrected dataset
# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/Batch_corrected_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData")
load("DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData")
Full_DT <- as.data.table(Full_DT)

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

# Load Mastery Hub Proteins (using CCL4 batch hubs as reference)
# load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top50_Mastery_Hubs_Protein_Merged_Batch.RData")
load("Top50_Mastery_Hubs_Protein_Merged_Batch.RData")
# Extract exactly the top 10 from 80% coverage
top_50_data <- top_50_80pct[order(top_50_80pct$Thresh_80pct, decreasing = TRUE), ][1:50, ]
top50_proteins <- str_split(top_50_data$GeneProtein, "_", simplify = TRUE)[, 2]
my_top_50_list <- paste0(top50_proteins, "_P")

cat("Selected top 10 mastery hub covariates:\n")
print(my_top_50_list)

# RUN MODELS (Protein only, using predefined mastery hubs covariates)
lasso.rf.pre.full.mergedata_blind_mastery50_protein_batch_over_0 <- pblapply(FUN = protein.regression.complete,
                                                           pairs.list.final,
                                                           dataset = Full_DT,
                                                           weighted = FALSE,
                                                           grouping = "all",          # <-- grouping all to build full design matrix
                                                           manual.covar = my_top_50_list, # <-- Mastery 10 variables
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
                                                           model.method = "lasso.predefined", # <-- lasso.predefined mode
                                                           prediction.method = "loo",
                                                           lasso.fam = "gaussian",
                                                           type.measure = "deviance",
                                                           alpha = 1,
                                                           intercept = TRUE,
                                                           nfolds = 10,
                                                           top.n = 50,
                                                           enforce.treatment = FALSE,
                                                           tune.ranger = FALSE,
                                                           manual.weights = FALSE,
                                                           weights.vector = NULL,
                                                           seed.set = 123)

# RUN MODELS (RNA + Protein, using predefined mastery hubs covariates + all RNA)
lasso.rf.pre.full.mergedata_blind_mastery50_plus_rna_batch_over_0 <- pblapply(FUN = protein.regression.complete,
                                                           pairs.list.final,
                                                           dataset = Full_DT,
                                                           weighted = FALSE,
                                                           grouping = "all",          # <-- grouping all to build full design matrix
                                                           manual.covar = my_top_50_list, # <-- Mastery 10 variables
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
                                                           model.method = "lasso.rf.pre", # <-- lasso.predefined mode
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
save(lasso.rf.pre.full.mergedata_blind_mastery50_protein_batch_over_0,
     lasso.rf.pre.full.mergedata_blind_mastery50_plus_rna_batch_over_0,
     file = "Models_mastery50_mergedata_blind_batch_over_0.RData")

# ==============================================================================
# 4. EVALUATION AND SUMMARY STATISTICS
# ==============================================================================
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_Merged_Data/BATCH_PRIME_OVERLAP/Mastery50/Models_mastery50_mergedata_blind_batch_over_0.RData")

res.m10 <- gof.all.complete(model = lasso.rf.pre.full.mergedata_blind_mastery50_protein_batch_over_0, prediction.method = "lasso")

cat("\n=== MASTERY 10 PROTEIN CORRELATIONS ===\n")
cor_m10 <- sapply(lasso.rf.pre.full.mergedata_blind_mastery50_protein_batch_over_0, function(x) as.numeric(x$prediction.obj$model$correlation.pearson))
print(summary(cor_m10))
cat("Count < 0.5:", sum(cor_m10 < 0.5, na.rm=TRUE), "\n")
cat("Count >= 0.5:", sum(cor_m10 >= 0.5, na.rm=TRUE), "\n")
cat("Count >= 0.8:", sum(cor_m10 >= 0.8, na.rm=TRUE), "\n")

