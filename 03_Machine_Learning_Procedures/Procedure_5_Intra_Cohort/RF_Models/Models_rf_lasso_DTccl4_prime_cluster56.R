library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

setwd("/work/smbrngo1")

# 1. LOAD DATA
# Load full uncorrected dataset (containing all 943 candidate covariates)
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/Funktionen_CD_Testing/raw_corrected_prime_analysis/DTccl4_filtered_final_Gene_Protein_full.RData")
load("DTccl4_filtered_final_Gene_Protein_full.RData")
CCL4_Full <- as.data.table(DTccl4_filtered_final)

# Filter for the subset cohort samples if applicable
# No subset filtering needed

# Extract target proteins (pairs) directly from ClusterDiPa column
pairs.list.5.6 <- unique(CCL4_Full[ClusterDiPa %in% c(5, 6), GeneProtein])

source("Function_protein_modelling_full_new_version_CD.R")
source("Function_design.R")
source("Function_model_measures.R")

# 2. RUN MODELS (Global pool of 943 covariates)

# A. RNA-only model
lasso.rf.pre.full.CD.new_rna_DTccl4_prime_cluster_5_6 <- pblapply(FUN = protein.regression.complete,
                                                          pairs.list.5.6,
                                                          dataset = CCL4_Full, # <-- Full dataset for Prime
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
                                                          model.method = "lasso.rf.pre", # <-- RF preselection
                                                          prediction.method = "loo",
                                                          lasso.fam = "gaussian",
                                                          type.measure = "deviance",
                                                          alpha = 1,
                                                          intercept = TRUE,
                                                          nfolds = NULL,
                                                          top.n = 30,
                                                          enforce.treatment = FALSE,
                                                          tune.ranger = FALSE,
                                                          manual.weights = FALSE,
                                                          weights.vector = NULL,
                                                          seed.set = 23)

# B. Protein-only model
lasso.rf.pre.full.CD.new_protein_DTccl4_prime_cluster_5_6 <- pblapply(FUN = protein.regression.complete,
                                                              pairs.list.5.6,
                                                              dataset = CCL4_Full,
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
                                                              nfolds = NULL,
                                                              top.n = 30,
                                                              enforce.treatment = FALSE,
                                                              tune.ranger = FALSE,
                                                              manual.weights = FALSE,
                                                              weights.vector = NULL,
                                                              seed.set = 23)

# C. Combined model
lasso.rf.pre.full.CD.new_combined_DTccl4_prime_cluster_5_6 <- pblapply(FUN = protein.regression.complete,
                                                               pairs.list.5.6,
                                                               dataset = CCL4_Full,
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
                                                               nfolds = NULL,
                                                               top.n = 30,
                                                               enforce.treatment = FALSE,
                                                               tune.ranger = FALSE,
                                                               manual.weights = FALSE,
                                                               weights.vector = NULL,
                                                               seed.set = 123)

# 3. SAVE RESULTS
save(lasso.rf.pre.full.CD.new_rna_DTccl4_prime_cluster_5_6,
     lasso.rf.pre.full.CD.new_protein_DTccl4_prime_cluster_5_6,
     lasso.rf.pre.full.CD.new_combined_DTccl4_prime_cluster_5_6,
     file = "Models_rf_preselection_full_objects_CD_new_DTccl4_prime_cluster_5_6.RData")
