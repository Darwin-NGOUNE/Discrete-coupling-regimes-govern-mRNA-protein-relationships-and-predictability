library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library(pbapply)

setwd("/work/smbrngo1")

# 1. LOAD DATA
# Load full uncorrected dataset (so that the columns for the 10 Mastery proteins are present)
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

# 2. LOAD MASTERY HUB PROTEINS (CCL4 Raw)
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Top50_Mastery_Hubs_Protein_CCL4_Raw.RData")
load("Top50_Mastery_Hubs_Protein_CCL4_Raw.RData")
# Extract exactly the top 10 from 80% coverage
top_50_data <- top_50_80pct[order(top_50_80pct$Thresh_80pct, decreasing = TRUE), ][1:50, ]
top50_proteins <- str_split(top_50_data$GeneProtein, "_", simplify = TRUE)[, 2]
my_top_50_list <- paste0(top50_proteins, "_P")

cat("Selected top 10 mastery hub covariates:\n")
print(my_top_50_list)

# 3. RUN MODEL (Mastery 10 Proteins on Full Dataset)
lasso.rf.pre.full.CD.new_mastery50_DTccl4_cluster_5_6 <- pblapply(FUN = protein.regression.complete,
                          pairs.list.5.6,
                          dataset = CCL4_Full, # <-- Full dataset to access all 10 hubs
                          weighted = FALSE,
                          grouping = "all",
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
                          nfolds = NULL,
                          top.n = 50,
                          enforce.treatment = FALSE,
                          tune.ranger = FALSE,
                          manual.weights = FALSE,
                          weights.vector = NULL,
                          seed.set = 123)

# 3b. RUN MODEL (Mastery 10 Proteins + RNA on Full Dataset)
lasso.rf.pre.full.CD.new_mastery50_plus_rna_DTccl4_cluster_5_6 <- pblapply(FUN = protein.regression.complete,
                          pairs.list.5.6,
                          dataset = CCL4_Full, # <-- Full dataset to access all 10 hubs
                          weighted = FALSE,
                          grouping = "all",
                          manual.covar = my_top_50_list,
                          design.m = "all",                 # <-- both RNA and proteins
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
                          model.method = "lasso.rf.pre",    # <-- lasso.rf.pre mode
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

# 4. SAVE RESULTS
save(lasso.rf.pre.full.CD.new_mastery50_DTccl4_cluster_5_6,
     lasso.rf.pre.full.CD.new_mastery50_plus_rna_DTccl4_cluster_5_6,
     file = "Models_mastery50_DTccl4_cluster_5_6.RData")

# 5. EVALUATION AND SUMMARY STATISTICS
# res.m10 <- gof.all.complete(model = lasso.rf.pre.full.CD.new_mastery50_DTccl4_cluster_5_6, prediction.method = "lasso")
# cat("\n=== MASTERY 10 PROTEIN CORRELATIONS ===\n")
# cor_m10 <- sapply(lasso.rf.pre.full.CD.new_mastery50_DTccl4_cluster_5_6, function(x) as.numeric(x$prediction.obj$model$correlation.pearson))
# print(summary(cor_m10))
# cat("Count < 0.5:", sum(cor_m10 < 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.5:", sum(cor_m10 >= 0.5, na.rm=TRUE), "\n")
# cat("Count >= 0.8:", sum(cor_m10 >= 0.8, na.rm=TRUE), "\n")
