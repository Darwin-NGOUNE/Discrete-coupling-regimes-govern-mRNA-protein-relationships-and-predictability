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
CCL4_Full <- CCL4_Full[(Treatment == 'ccl4' & TreatmentTime == 12) | (Treatment == 'oil' & TreatmentTime == 0)]

# Extract target proteins (pairs) directly from ClusterDiPa column
pairs.list.5.6.sub <- unique(CCL4_Full[ClusterDiPa %in% c(5, 6), GeneProtein])

source("Function_protein_modelling_full_new_version_CD.R")
source("Function_design.R")
source("Function_model_measures.R")

# 2. RUN MODELS (Global pool of 943 covariates)


# B. Protein-only model
rf.rf.pre.full.CD.new_protein_DTccl4_prime_cluster_5_6_sub <- pblapply(FUN = protein.regression.complete,
                                                                          pairs.list.5.6.sub,
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
                                                                          model.method = "rf.rf.pre",
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




# 3. SAVE RESULTS
save(rf.rf.pre.full.CD.new_protein_DTccl4_prime_cluster_5_6_sub,
     file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung/Intra_RF_RF/Models_rf_rf_preselection_full_objects_CD_new_DTccl4_prime_cluster_5_6_sub.RData")


#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Cluster_Modelierung/Intra_new/Models_rf_preselection_full_objects_CD_new_DTccl4_prime_cluster_5_6_sub.RData")


