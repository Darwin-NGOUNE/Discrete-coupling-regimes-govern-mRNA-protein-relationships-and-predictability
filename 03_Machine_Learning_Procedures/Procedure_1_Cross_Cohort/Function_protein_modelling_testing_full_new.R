# load required packages
library("mgcv")
library("glmnet")
library("stringr")
library(data.table)
library("pbapply")
library(LBI)
library(corrr)
library(caret)
library("ranger")
library("tuneRanger")
library("mlr")

# implement external functions for the design matrix and adaptive lasso

source("Function_design_newdata_filtered.R")
source("Function_design.R")
source("Function_adaptive_lasso.R")



# function to manually compute the log-likelihood of a LASSO model
log_likelihood <- function(y, y_pred) {
  n <- length(y)
  residuals <- y - y_pred
  sigma2 <- sum(residuals^2) / n
  logL <- -n / 2 * log(2 * pi * sigma2) - 1 / (2 * sigma2) * sum(residuals^2)
  return(logL)
}


# FUNCTION FOR MULTICOLLINEARITY REDUCTION
collinearityReduce <- function(mRNA, Rho_cutoff) {
  gene <- mRNA[, apply(mRNA != 0, 2, all)]
  row.names(gene) <- gene[, 1]
  gene <- gene[, -c(1, 2)]

  data_cor <- cor(gene, method = "pearson")
  data_highcor1 <- findCorrelation(data_cor, cutoff = Rho_cutoff)
  gene.temp <- gene[, -data_highcor1]

  return(gene.temp)
}

#-------------------------------------------------------------------------------

protein.regression.baseline <- function(dataset,
                                        geneprotein,
                                        include.treatment = FALSE,
                                        intercept = TRUE,
                                        newdata,
                                        seed.set = NULL) {
  # 1. Préparer les données d'entraînement
  sub_data <- dataset[GeneProtein == geneprotein, ]
  if (include.treatment) {
    sub_data <- sub_data[!is.na(ProteinIntensity) & !is.na(GeneCount) & !is.na(Treatment), ]
  } else {
    sub_data <- sub_data[!is.na(ProteinIntensity) & !is.na(GeneCount), ]
  }

  if (nrow(sub_data) < 3) {
    warning(paste("Not enough complete observations in the baseline dataset for", geneprotein, "- Skipping..."))
    return(NULL)
  }

  gene <- str_replace(unique(sub_data$GeneSyn), "-", ".")
  protein <- unique(sub_data$ProteinID)
  y <- sub_data$ProteinIntensity

  sub_data$TreatmentTotal <- factor(sub_data$TreatmentTotal,
    levels = c(
      "month0_oil", "month2_oil", "month12_oil",
      "month2_ccl4", "month6_ccl4", "month12_ccl4"
    )
  )

  # 2. Formule du modèle
  if (include.treatment) {
    form.obj <- ProteinIntensity ~ Treatment + GeneCount
  } else {
    form.obj <- ProteinIntensity ~ GeneCount
  }

  # 3. Entraîner le modèle
  gam.obj <- gam(form.obj,
    data = sub_data,
    family = "gaussian", drop.intercept = !(intercept)
  )

  # R^2 et deviance expliquée
  rsq <- summary(gam.obj)$r.sq
  dev.expl <- summary(gam.obj)$dev.expl

  # AIC, BIC, logLik
  gam.aic <- gam.obj$aic
  gam.bic <- BIC(gam.obj)
  ll_gam <- logLik(gam.obj)
  k <- dim(model.matrix(gam.obj))[2] - 1

  ##############################################################################
  # Subset to the target protein first (BUG FIXED)
  newdata_sub <- newdata[newdata$GeneProtein == geneprotein, ]
  
  # Strict filtration des NAs pour les variables de test (comme pour l'entraînement)
  if (include.treatment) {
    newdata_sub <- newdata_sub[!is.na(ProteinIntensity) & !is.na(GeneCount) & !is.na(Treatment), ]
  } else {
    newdata_sub <- newdata_sub[!is.na(ProteinIntensity) & !is.na(GeneCount), ]
  }

  if (nrow(newdata_sub) == 0) {
    warning(paste("No complete observations in the test dataset for", geneprotein, "- Skipping..."))
    return(NULL)
  }
  ##############################################################################

  # 4. Prédictions sur newdata
  newdat <- get_all_vars(form.obj, newdata_sub)
  prediction <- predict(gam.obj, newdat)

  # y = valeurs observées dans newdata
  y_test <- newdata_sub$ProteinIntensity

  # 5. Calcul des métriques
  rmse <- sqrt(mean((y_test - prediction)^2, na.rm = TRUE))
  mae <- mean(abs(y_test - prediction), na.rm = TRUE)
  pearson_cor <- cor(y_test, prediction, method = "pearson", use = "complete.obs")
  spearman_cor <- cor(y_test, prediction, method = "spearman", use = "complete.obs")

  sd_y_train   <- sd(y, na.rm = TRUE)
  mean_y_train <- mean(y, na.rm = TRUE)
  nrmse_trainSD <- if (!is.na(sd_y_train) && sd_y_train > 0) rmse / sd_y_train else NA
  ss_res <- sum((y_test - prediction)^2, na.rm = TRUE)
  ss_tot <- sum((y_test - mean_y_train)^2, na.rm = TRUE)
  r2_test <- if (!is.na(ss_tot) && ss_tot > 0) 1 - (ss_res / ss_tot) else NA

  # 6. Retourner les résultats
  gam_list <- list(
    regression.object = gam.obj,
    gene = gene,
    protein = protein,
    logLikelihood = ll_gam,
    k = k,
    adj.R.squared = rsq,
    y = y,
    y_test = y_test,
    prediction = prediction,
    deviance.explained = dev.expl,
    AIC = gam.aic,
    BIC = gam.bic,
    RMSE = rmse,
    MAE = mae,
    nrmse_trainSD = nrmse_trainSD,
    r2_test = r2_test,
    Pearson = pearson_cor,
    Spearman = spearman_cor
  )

  return(gam_list)
}



# Lasso ------------------------------------------------------------------------

lasso.mod <- function(Design,
                      y,
                      y.t,
                      w,
                      n,
                      gene = NULL,
                      protein = NULL,
                      newdata, # <-- obligatoire : jeu de données externe
                      lasso.fam = "gaussian",
                      alpha = 1,
                      nfolds = 10,
                      type.measure = "deviance",
                      intercept = TRUE,
                      seed.set = NULL,
                      na.pro = "drop.cols") {
  # 1. Validation croisée pour choisir lambda optimal --------------------------
  set.seed(seed.set)
  lasso.obj <- cv.glmnet(
    y = y,
    x = as.matrix(Design),
    alpha = alpha,
    family = lasso.fam,
    nfolds = nfolds,
    penalty.factor = w,
    type.measure = type.measure,
    relax = TRUE,
    gamma = 0,
    intercept = intercept
  )

  set.seed(seed.set)

  # 3. Ajustement final sur toutes les données d'entraînement ---------
  final.lasso <- do.call(glmnet,
    args = list(
      y = y,
      x = as.matrix(Design),
      alpha = alpha,
      family = as.character(lasso.fam),
      penalty.factor = w,
      relax = TRUE,
      gamma = 0,
      intercept = intercept,
      lambda = lasso.obj$lambda
    )
  )

  # RELAXED LASSO
  # coef.glmnet_relaxed <- coef(final.lasso$relaxed)
  # #opt.glmnet_relaxed <- match(lasso.obj$relaxed$lambda.min, lasso.obj$lambda)
  # opt.glmnet_relaxed <- which(abs(lasso.obj$relaxed$lambda.min - final.lasso$lambda) ==
  #                               min(abs(lasso.obj$relaxed$lambda.min - final.lasso$lambda)))


  lambda.opt_relaxed <- lasso.obj$relaxed$lambda.min

  # ## coefficients at optimal lambda
  # coef.glmnet.opt_relaxed <- coef.glmnet_relaxed[, opt.glmnet_relaxed]
  #
  # # deviance ratio (for elastic net this is the R^2)
  # final.lasso$relaxed$dev.ratio[opt.glmnet_relaxed]

  ###################################################### only for RF PRESELECTION
  # ✅ **En résumé** : ce code garantit que `newdata` possède les mêmes variables (colonnes) que `Design`.

  common_vars <- colnames(Design)
  newdata <- newdata[, common_vars, drop = FALSE]

  common_vars_new <- colnames(newdata)
  selected_vars <- as.list(common_vars)
  selected_vars_new <- as.list(common_vars_new)

  # common_vars <- intersect(colnames(Design), colnames(newdata))
  # newdata <- newdata[, common_vars, drop = FALSE]
  # selected_vars <- as.list(common_vars)

  ####################################################

  # prediction based on model
  # prediction.vector_relaxed <- predict(final.lasso$relaxed, as.matrix(newdata), s = lambda.opt_relaxed)

  ##############################################################################
  # THIS WORKS
  # directly via the CV object
  coef.relaxed <- coef(lasso.obj, s = "lambda.min")[which(coef(lasso.obj, s = "lambda.min") != 0)]
  names(coef.relaxed) <- rownames(coef(lasso.obj, s = "lambda.min"))[which(coef(lasso.obj, s = "lambda.min") != 0)]

  prediction.vector_relaxed <- predict(lasso.obj, newx = as.matrix(newdata), s = "lambda.min")
  ##############################################################################

  # UNRELAXED LASSO
  coef.glmnet <- coef(final.lasso)

  # opt.glmnet <- match(lasso.obj$lambda.min, lasso.obj$lambda)
  # workaround because of some internal glmnet discrepancies between lambda object and lambda.min
  opt.glmnet <- which(abs(lasso.obj$lambda.min - final.lasso$lambda) == min(abs(lasso.obj$lambda.min - final.lasso$lambda)))

  lambda.opt <- lasso.obj$lambda.min

  ## coefficients at optimal lambda
  coef.glmnet.opt <- coef.glmnet[, opt.glmnet]

  # deviance ratio (for elastic net this is the R^2)
  final.lasso$dev.ratio[opt.glmnet]

  # 5. Prédictions UNIQUEMENT sur le nouveau jeu ----------------------
  prediction.vector <- predict(final.lasso, as.matrix(newdata), s = lambda.opt)

  ######
  # 6. Retourner les résultats ----------------------------------------
  variable.selection <- coef.glmnet.opt[which(coef.glmnet.opt != 0)]
  # variable.selection_relaxed <- coef.glmnet.opt_relaxed[which(coef.glmnet.opt_relaxed != 0)]
  variable.selection_relaxed <- coef.relaxed[which(coef.relaxed != 0)]


  # for gaussian family adj. R2 is:
  r2 <- lasso.obj$glmnet.fit$dev.ratio[opt.glmnet]
  # r2_relaxed <- lasso.obj$glmnet.fit$relaxed$dev.ratio[opt.glmnet_relaxed]


  # (1 - lasso.obj$cvm/var(y))[opt.glmnet]

  ## compute adj. R^2, AIC and BIC
  k <- final.lasso$df[opt.glmnet] + 1
  n <- final.lasso$nobs

  # TODO adj. R2 ? (double check) ????
  adj.r2 <- 1 - (1 - r2) * (n - 1) / (n - k + 1)

  # k_relaxed <- final.lasso$relaxed$df[opt.glmnet_relaxed] + 1
  n_relaxed <- final.lasso$relaxed$nobs

  # TODO adj. R2 ? (double check) ????
  adj.r2_relaxed <- 1 - (1 - r2) * (n - 1) / (n - k + 1)

  # NOT CORRECT:
  # tLL <- final.lasso$nulldev - deviance(final.lasso)[opt.glmnet]
  # lasso.aic <- -tLL+2*k#+2*k*(k+1)/(n-k-1) #correction term (AICc)
  # lasso.bic <-log(n)*k - tLL

  correlation.pearson <- cor(prediction.vector, y.t, method = "pearson", use = "complete.obs")
  correlation.spearman <- cor(prediction.vector, y.t, method = "spearman", use = "complete.obs")

  correlation.pearson_relaxed <- cor(prediction.vector_relaxed, y.t, method = "pearson", use = "complete.obs")
  correlation.spearman_relaxed <- cor(prediction.vector_relaxed, y.t, method = "spearman", use = "complete.obs")

  rmse <- sqrt(mean((prediction.vector - y.t)^2, na.rm = TRUE))
  mae <- mean(abs(prediction.vector - y.t), na.rm = TRUE)

  rmse_relaxed <- sqrt(mean((prediction.vector_relaxed - y.t)^2, na.rm = TRUE))
  mae_relaxed <- mean(abs(prediction.vector_relaxed - y.t), na.rm = TRUE)

  sd_y_train   <- sd(y, na.rm = TRUE)
  mean_y_train <- mean(y, na.rm = TRUE)

  nrmse_trainSD <- if (!is.na(sd_y_train) && sd_y_train > 0) rmse / sd_y_train else NA
  ss_res <- sum((y.t - prediction.vector)^2, na.rm = TRUE)
  ss_tot <- sum((y.t - mean_y_train)^2, na.rm = TRUE)
  r2_test <- if (!is.na(ss_tot) && ss_tot > 0) 1 - (ss_res / ss_tot) else NA

  nrmse_trainSD_relaxed <- if (!is.na(sd_y_train) && sd_y_train > 0) rmse_relaxed / sd_y_train else NA
  ss_res_relaxed <- sum((y.t - prediction.vector_relaxed)^2, na.rm = TRUE)
  r2_test_relaxed <- if (!is.na(ss_tot) && ss_tot > 0) 1 - (ss_res_relaxed / ss_tot) else NA

  # likelihood for the LASSO model
  ll_lasso <- log_likelihood(y.t, prediction.vector)
  ll_lasso_relaxed <- log_likelihood(y.t, prediction.vector_relaxed)

  k_ll <- length(variable.selection) + 1
  k_ll_relaxed <- length(variable.selection_relaxed) + 1

  lasso.aic <- 2 * k_ll - 2 * ll_lasso
  lasso.bic <- log(n) * k_ll - 2 * ll_lasso

  lasso.aic_relaxed <- 2 * k_ll_relaxed - 2 * ll_lasso_relaxed
  lasso.bic_relaxed <- log(n) * k_ll_relaxed - 2 * ll_lasso_relaxed


  lasso_list <- list(
    final.lasso = final.lasso,
    cv.lasso.obj = lasso.obj,
    lambda.opt = lasso.obj$lambda.min,
    gene = gene,
    protein = protein,
    y = y,
    y.t = y.t,
    Design = Design,
    newdata = newdata,
    selected_vars = selected_vars,
    selected_vars_new = selected_vars_new,


    # unrelaxed
    variable.selection = variable.selection,
    all.covar = coef.glmnet.opt,
    k = k,
    predictions.insample = prediction.vector,
    correlation.pearson = correlation.pearson,
    correlation.spearman = correlation.spearman,
    logLikelihood = ll_lasso,
    deviance.explained = final.lasso$dev.ratio[opt.glmnet],
    adj.R.squared = adj.r2,
    AIC = lasso.aic,
    BIC = lasso.bic,
    rmse = rmse,
    mae = mae,
    nrmse_trainSD = nrmse_trainSD,
    r2_test = r2_test,

    # relaxed
    variable.selection_relaxed = variable.selection_relaxed,
    all.covar_relaxed = coef.relaxed,
    predictions.insample_relaxed = prediction.vector_relaxed,
    correlation.pearson_relaxed = correlation.pearson_relaxed,
    correlation.spearman_relaxed = correlation.spearman_relaxed,
    logLikelihood_relaxed = ll_lasso_relaxed,
    adj.R.squared_relaxed = adj.r2_relaxed,
    AIC_relaxed = lasso.aic_relaxed,
    BIC_relaxed = lasso.bic_relaxed,
    rmse_relaxed = rmse_relaxed,
    mae_relaxed = mae_relaxed,
    nrmse_trainSD_relaxed = nrmse_trainSD_relaxed,
    r2_test_relaxed = r2_test_relaxed
  )

  return(lasso_list)
}

# RF----------------------------------------------------------------------------


# Random Forest model avec prédiction uniquement sur newdata
rf.mod <- function(Design,
                   y,
                   y.t,
                   w,
                   n,
                   gene,
                   protein,
                   newdata, # <-- obligatoire : jeu de données externe
                   tune.ranger = FALSE,
                   n.cv = 10,
                   seed.set = 23) {
  # 1. Préparer les données d'entraînement ----------------------------
  DTrf <- data.frame(Design, y = y)
  rf.task <- makeRegrTask(data = DTrf, target = "y")

  # 2. Entraîner le modèle Random Forest -------------------------------
  if (tune.ranger) {
    # Avec tuning automatique
    set.seed(seed.set)
    rf.tune <- tuneRanger(rf.task,
      measure = list(mse), num.trees = 1000,
      num.threads = 2, iters = 70, show.info = TRUE
    )

    fit.ranger <- ranger(
      x = Design, y = y, data = DTrf,
      importance = "impurity_corrected",
      probability = FALSE, classification = FALSE,
      num.trees = 500, oob.error = TRUE,
      mtry = rf.tune$recommended.pars$mtry,
      min.node.size = rf.tune$recommended.pars$min.node.size,
      sample.fraction = rf.tune$recommended.pars$sample.fraction
    )
  } else {
    # Sans tuning (paramètres par défaut)
    set.seed(seed.set)
    fit.ranger <- ranger(
      x = Design, y = y, data = DTrf,
      importance = "impurity_corrected",
      probability = FALSE, classification = FALSE,
      num.trees = 500, oob.error = TRUE,
      mtry = ncol(Design) / 3,
      min.node.size = 5
    )
  }

  # 3. Prédictions UNIQUEMENT sur le nouveau jeu -----------------------
  prediction <- predict(fit.ranger, newdata)
  prediction.vector <- prediction$predictions

  # 4. Importance des variables ----------------------------------------
  variable.importance <- importance(fit.ranger)

  # 5. Calcul des métriques --------------------------------------------
  r2 <- fit.ranger$r.squared
  n <- length(y)
  p <- ncol(Design)
  adj.r2 <- 1 - ((1 - r2) * (n - 1) / (n - p - 1))

  # Corrélations (si newdata contient aussi y_test)
  correlation.pearson <- cor(prediction.vector, y.t, method = "pearson", use = "complete.obs")
  correlation.spearman <- cor(prediction.vector, y.t, method = "spearman", use = "complete.obs")

  rmse <- sqrt(mean((prediction.vector - y.t)^2, na.rm = TRUE))
  mae  <- mean(abs(prediction.vector - y.t), na.rm = TRUE)

  sd_y_train   <- sd(y, na.rm = TRUE)
  mean_y_train <- mean(y, na.rm = TRUE)
  nrmse_trainSD <- if (!is.na(sd_y_train) && sd_y_train > 0) rmse / sd_y_train else NA
  ss_res <- sum((y.t - prediction.vector)^2, na.rm = TRUE)
  ss_tot <- sum((y.t - mean_y_train)^2, na.rm = TRUE)
  r2_test <- if (!is.na(ss_tot) && ss_tot > 0) 1 - (ss_res / ss_tot) else NA

  # likelihood for the LASSO model
  ll_rf <- log_likelihood(y.t, prediction.vector)

  # Varianz der Residuen schätzen
  sigma2 <- mean((y.t - prediction.vector)^2, na.rm = TRUE)

  # Anzahl der Parameter (p + 1 für Intercept)
  k <- p + 1

  # AIC und BIC berechnen
  aic <- n * log(sigma2) + 2 * k
  bic <- n * log(sigma2) + log(n) * k

  # 6. Retourner les résultats -----------------------------------------
  rf_list <- list(
    rf.model = fit.ranger,
    gene = gene,
    protein = protein,
    variable.importance = variable.importance,
    k = k,
    predictions.insample = prediction.vector,
    correlation.pearson = correlation.pearson,
    correlation.spearman = correlation.spearman,
    rmse = rmse,
    mae = mae,
    nrmse_trainSD = nrmse_trainSD,
    r2_test = r2_test,
    y = y,
    y.t = y.t,
    logLikelihood = ll_rf,
    adj.R.squared = adj.r2,
    AIC = aic,
    BIC = bic
  )

  return(rf_list)
}


# FUNCTION TO PRESELECT COVARIATES FROM RF FEATURE IMPORTANCE ------------------

rf.preselect <- function(Design, # matrice des covariables (X)
                         y, # vecteur réponse
                         w, # poids éventuels (non utilisés ici)
                         n, # nombre d’observations
                         gene = gene, # nom du gène (info)
                         protein, # nom de la protéine (info)
                         seed.set = 23, # graine aléatoire
                         tune.ranger = FALSE, # tuning automatique (non utilisé ici)
                         top.n = 30) { # nombre de variables à garder

  # 1. Créer un data.frame avec X et y -------------------------------
  DTrf <- data.frame(Design, y = y)

  # 2. Définir une tâche de régression pour mlr (optionnel) ----------
  rf.task <- makeRegrTask(data = DTrf, target = "y")

  # 3. Entraîner un modèle Random Forest -----------------------------
  set.seed(seed.set)
  fit.ranger <- ranger(
    x = Design, y = y, data = DTrf,
    importance = "impurity_corrected", # méthode d’importance
    probability = FALSE, classification = FALSE,
    num.trees = 2500, # nombre d’arbres
    oob.error = TRUE, # erreur OOB
    mtry = floor(ncol(Design) / 3), # nb de variables testées à chaque split
    min.node.size = 5
  ) # taille min des feuilles

  # 4. Extraire l’importance des variables ---------------------------
  # importance(fit.ranger) renvoie un score pour chaque covariable
  # On trie par ordre décroissant et on garde les top.n
  top.important.features <- sort(importance(fit.ranger), decreasing = TRUE)[1:top.n]

  # 5. Retourner les variables les plus importantes ------------------
  return(top.important.features)
}


################################################################################


protein.regression.function <- function(
  dataset,
  geneprotein,
  newdata, # <-- obligatoire : jeu de données externe
  weighted = FALSE,
  weight.method = NULL,
  weight.matrix = NULL,
  weight.matrix.p = NULL,
  grouping = "all",
  scaled = TRUE,
  RNA.log.scale = TRUE,
  manual.covar = NULL,
  design.m = "genes",
  na.process = "max.n",
  na.pro = "drop.cols",
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
  model.method = "lasso", # choix: baseline.only, lasso, lasso.rf.pre, rf
  lasso.fam = "gaussian",
  type.measure = "deviance",
  alpha = 1,
  top.n = 30,
  enforce.treatment = FALSE,
  nfolds = 10,
  tune.ranger = FALSE,
  n.cv = 10,
  seed.set = 23,
  sample.seed = NULL,
  intercept = TRUE
) {
  # 1. Créer l’objet Design pour train -------------------------------
  Design.obj <- design.function.new(
    dataset = dataset,
    geneprotein = geneprotein,
    weighted = weighted,
    weight.method = weight.method,
    weight.matrix = weight.matrix,
    weight.matrix.p = weight.matrix.p,
    grouping = grouping,
    manual.covar = manual.covar,
    design.m = design.m,
    scaled = scaled,
    RNA.log.scale = RNA.log.scale,
    na.process = na.process,
    include.treatment = include.treatment,
    treatment.info = treatment.info,
    duration.scale = duration.scale,
    PE.RNA.penalty = PE.RNA.penalty,
    treatment.penalty = treatment.penalty,
    baseline.interactions = baseline.interactions
  )

  # 2. Créer l’objet Design pour test -------------------------------
  Design.test.obj <- design.function(
    dataset = newdata,
    geneprotein = geneprotein,
    weighted = weighted,
    weight.method = weight.method,
    weight.matrix = weight.matrix,
    weight.matrix.p = weight.matrix.p,
    grouping = grouping,
    manual.covar = manual.covar,
    design.m = design.m,
    scaled = scaled,
    RNA.log.scale = RNA.log.scale,
    na.process = na.process,
    include.treatment = include.treatment,
    treatment.info = treatment.info,
    duration.scale = duration.scale,
    PE.RNA.penalty = PE.RNA.penalty,
    treatment.penalty = treatment.penalty,
    baseline.interactions = baseline.interactions
  )

  # Entraînement
  gene.train <- Design.obj$gene
  protein.train <- Design.obj$protein
  y.train <- Design.obj$y
  Design.train <- Design.obj$Design
  n.train <- Design.obj$n
  protein.vector.train <- Design.obj$protein.vector
  mice.train <- Design.obj$mice
  experiments.train <- Design.obj$experiments

  # Test
  gene.test <- Design.test.obj$gene
  protein.test <- Design.test.obj$protein
  y.test <- Design.test.obj$y
  Design.test <- Design.test.obj$Design
  n.test <- Design.test.obj$n
  protein.vector.test <- Design.test.obj$protein.vector
  mice.test <- Design.test.obj$mice
  experiments.test <- Design.test.obj$experiments


  # Poids
  if (manual.weights == FALSE) {
    w.train <- Design.obj$w
    w.test <- Design.test.obj$w
  } else {
    w.train <- weights.vector
    w.test <- weights.vector
  }

  # 4. Réduction de colinéarité --------------------------------------
  if (collinearity.reduce) {
    Design_temp <- collinearityReduce(Design.train, rho.cutoff)
    w_temp <- w.train[match(colnames(Design_temp), names(w.train))]
    i.g <- which(colnames(Design.train) == paste0(gene.train, "_G"))
    Design.train <- cbind(Design.train[, i.g, drop = FALSE], Design_temp)
    w.train <- c(w.train[i.g], w_temp)
    colnames(Design.train)[1] <- paste0(gene.train, "_G")
    names(w.train)[1] <- paste0(gene.train, "_G")
  }

  # 5. Harmonisation GLOBALE des matrices Train et Test ----------------
  # On s'assure que TOUS les modèles (Lasso, RF-pre, RF) reçoivent exactement
  # les mêmes variables (colonnes) pour éviter le crash "undefined columns selected".
  ############################################################################
  common_vars <- intersect(colnames(Design.train), colnames(Design.test))

  Design.train <- Design.train[, common_vars, drop = FALSE]
  Design.test <- Design.test[, common_vars, drop = FALSE]

  # Subset the weights accordingly
  w.train <- w.train[common_vars]
  w.test <- w.test[common_vars]
  ############################################################################
  rf_design_cols <- ncol(Design.train)

  # 6. Baseline model ------------------------------------------------
  baseline.model <- protein.regression.baseline(
    dataset = dataset,
    geneprotein = geneprotein,
    include.treatment = include.treatment,
    intercept = intercept,
    newdata = newdata
  )

  # 7. Choix du modèle principal -------------------------------------
  if (model.method == "baseline.only") {
    model <- baseline.model
  } else if (model.method == "lasso") {
    model <- lasso.mod(
      Design = Design.train,
      y = y.train,
      y.t = y.test,
      w = w.train,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      newdata = Design.test, # prédiction externe
      alpha = alpha,
      lasso.fam = lasso.fam,
      nfolds = nfolds,
      type.measure = type.measure,
      intercept = intercept,
      seed.set = seed.set,
      na.pro = na.pro
    )
  } else if (model.method == "lasso.predefined") {
    # =========================================================================
    # NOUVEAU MODE : LASSO SUR UNE LISTE PREDÉFINIE
    # =========================================================================

    if (is.null(manual.covar)) {
      stop("Pour utiliser le modele 'lasso.predefined', vous devez fournir une liste de variables dans l'argument 'manual.covar'.")
    }

    # 1. On s'assure que le gène propre de la protéine est forcé dans les variables sélectionnées
    target_gene_name <- paste0(gene.train, "_G")
    if (design.m == "all") {
      rna_covars <- colnames(Design.train)[grep("_G$", colnames(Design.train))]
      variables_to_keep <- c(manual.covar, rna_covars)
    } else {
      variables_to_keep <- c(manual.covar, target_gene_name)
    }

    # enforce treatment = on inclut aussi la variable de traitement si demandé (où penalité w == 0)
    if (enforce.treatment == TRUE & include.treatment == TRUE) {
      unpenalized.variables <- names(w.train)[which(w.train == 0)]
      variables_to_keep <- c(variables_to_keep, unpenalized.variables)
    }

    # 2. On trouve l'intersection stricte entre train, test et notre liste
    common_all <- intersect(colnames(Design.train), colnames(Design.test))
    cols_to_keep <- intersect(variables_to_keep, common_all)
    
    Design.predef.train <- Design.train[, cols_to_keep, drop = FALSE]
    Design.predef.test <- Design.test[, cols_to_keep, drop = FALSE]
    
    # glmnet crash prevention: it requires >= 2 columns
    if (ncol(Design.predef.train) < 2) {
      model <- "Aucune Mastery Protein trouvee dans ce cluster."
    } else {
      # On s'assure de garder les pénalités correspondantes pour ces mêmes colonnes
      w.predef <- w.train[cols_to_keep]

      # 3. On passe ces données au modèle Lasso classique
      model <- lasso.mod(
        Design = Design.predef.train,
        y = y.train,
        y.t = y.test,
        w = w.predef,
        n = n.train,
        gene = gene.train,
        protein = protein.train,
        newdata = Design.predef.test,
        alpha = alpha,
        lasso.fam = lasso.fam,
        nfolds = nfolds,
        type.measure = type.measure,
        intercept = intercept,
        seed.set = seed.set,
        na.pro = na.pro
      )
    }
  } else if (model.method == "lasso.rf.pre") {
    if (design.m == "all" & !is.null(manual.covar)) {
      rna_covars <- colnames(Design.train)[grep("_G$", colnames(Design.train))]
      variables_to_keep <- c(manual.covar, rna_covars)
      common_all <- intersect(colnames(Design.train), colnames(Design.test))
      cols_to_keep <- intersect(variables_to_keep, common_all)
      Design.rf.input <- Design.train[, cols_to_keep, drop = FALSE]
      w.rf.input <- w.train[cols_to_keep]
      rf_design_cols <- ncol(Design.rf.input)
    } else {
      Design.rf.input <- Design.train
      w.rf.input <- w.train
    }

    rf.pre.model <- rf.preselect(
      Design = Design.rf.input,
      y = y.train,
      w = w.rf.input,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      top.n = top.n,
      seed.set = seed.set
    )

    if (enforce.treatment == TRUE & include.treatment == TRUE) {
      unpenalized.variables <- names(w.train)[which(w.train == 0)]
      Design.rf <- Design.train[, colnames(Design.train) %in% names(rf.pre.model) | colnames(Design.train) %in% unpenalized.variables]
      w.rf <- w.train[names(w.train) %in% names(rf.pre.model) | names(w.train) %in% unpenalized.variables]
    } else {
      Design.rf <- Design.train[, colnames(Design.train) %in% names(rf.pre.model) | colnames(Design.train) %in% paste0(gene.train, "_G")]
      w.rf <- w.train[names(w.train) %in% names(rf.pre.model) | names(w.train) %in% paste0(gene.train, "_G")]
    }

    model <- lasso.mod(
      Design = Design.rf,
      y = y.train,
      y.t = y.test,
      w = w.rf,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      newdata = Design.test, # prédiction externe
      alpha = alpha,
      lasso.fam = lasso.fam,
      nfolds = nfolds,
      type.measure = type.measure,
      intercept = intercept,
      seed.set = seed.set,
      na.pro = na.pro
    )
  } else if (model.method == "rf") {
    model <- rf.mod(
      Design = Design.train,
      y = y.train,
      y.t = y.test,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      newdata = Design.test, # prédiction externe
      tune.ranger = tune.ranger,
      seed.set = seed.set
    )
  } else if (model.method == "rf.rf.pre") {
    # =========================================================================
    # TWO-STEP RANDOM FOREST PIPELINE (RF + RF):
    #   Step 1: Random Forest Feature Pre-selection (top.n most important covariates)
    #   Step 2: Random Forest Regression Model on the pre-selected feature space
    # =========================================================================

    # 1. Prepare candidate predictor matrices for feature pre-selection
    if (design.m == "all" & !is.null(manual.covar)) {
      rna_covars <- colnames(Design.train)[grep("_G$", colnames(Design.train))]
      variables_to_keep <- c(manual.covar, rna_covars)
      common_all <- intersect(colnames(Design.train), colnames(Design.test))
      cols_to_keep <- intersect(variables_to_keep, common_all)
      Design.rf.input <- Design.train[, cols_to_keep, drop = FALSE]
      w.rf.input <- w.train[cols_to_keep]
      rf_design_cols <- ncol(Design.rf.input)
    } else {
      Design.rf.input <- Design.train
      w.rf.input <- w.train
    }

    # 2. Execute Step 1: Pre-select top.n features using Random Forest importance
    rf.pre.model <- rf.preselect(
      Design = Design.rf.input,
      y = y.train,
      w = w.rf.input,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      top.n = top.n,
      seed.set = seed.set
    )

    # 3. Subset Training and Testing matrices to the pre-selected top.n features (plus target gene)
    if (enforce.treatment == TRUE & include.treatment == TRUE) {
      unpenalized.variables <- names(w.train)[which(w.train == 0)]
      selected_cols <- colnames(Design.train)[colnames(Design.train) %in% names(rf.pre.model) | colnames(Design.train) %in% unpenalized.variables]
    } else {
      selected_cols <- colnames(Design.train)[colnames(Design.train) %in% names(rf.pre.model) | colnames(Design.train) %in% paste0(gene.train, "_G")]
    }

    Design.rf.train <- Design.train[, selected_cols, drop = FALSE]
    Design.rf.test  <- Design.test[, selected_cols, drop = FALSE]

    # 4. Execute Step 2: Fit final Random Forest regression model on the pre-selected subset
    model <- rf.mod(
      Design = Design.rf.train,
      y = y.train,
      y.t = y.test,
      n = n.train,
      gene = gene.train,
      protein = protein.train,
      newdata = Design.rf.test, # External cross-cohort test evaluation
      tune.ranger = tune.ranger,
      seed.set = seed.set
    )
  }

  # 7. Retourner les résultats ----------------------------------------
  # return(list(model = model,
  #             baseline.model = baseline.model,
  #             experiments = experiments))
  return(list(
    model = model,
    baseline.model = baseline.model,
    experiments.train = Design.obj$experiments,
    experiments.test = Design.test.obj$experiments,
    rf_mtry = if (model.method %in% c("lasso.rf.pre", "rf.rf.pre", "rf", "lasso.predefined")) floor(rf_design_cols / 3) else NULL
  ))
}

################################################################################
protein.regression.complete <- function(
  dataset,
  geneprotein,
  newdata, # <-- obligatoire
  weighted = FALSE,
  weight.method = NULL,
  weight.matrix = NULL,
  weight.matrix.p = NULL,
  grouping = "all",
  scaled = TRUE,
  RNA.log.scale = TRUE,
  manual.covar = NULL,
  design.m = "genes",
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
  model.method = "lasso", # choix: baseline.only, lasso, lasso.rf.pre, rf
  lasso.fam = "gaussian",
  type.measure = "deviance",
  alpha = 1,
  top.n = 30,
  enforce.treatment = FALSE,
  nfolds = 10,
  tune.ranger = FALSE,
  n.cv = 10,
  seed.set = 23,
  intercept = TRUE
) {
  # 1. Entraînement du modèle principal -------------------------------
  model.obj <- protein.regression.function(
    dataset = dataset,
    geneprotein = geneprotein,
    newdata = newdata, # toujours fourni
    weighted = weighted,
    weight.method = weight.method,
    weight.matrix = weight.matrix,
    weight.matrix.p = weight.matrix.p,
    grouping = grouping,
    manual.covar = manual.covar,
    design.m = design.m,
    scaled = scaled,
    RNA.log.scale = RNA.log.scale,
    na.process = na.process,
    include.treatment = include.treatment,
    treatment.info = treatment.info,
    duration.scale = duration.scale,
    PE.RNA.penalty = PE.RNA.penalty,
    treatment.penalty = treatment.penalty,
    baseline.interactions = baseline.interactions,
    collinearity.reduce = collinearity.reduce,
    rho.cutoff = rho.cutoff,
    manual.weights = manual.weights,
    weights.vector = weights.vector,
    model.method = model.method,
    lasso.fam = lasso.fam,
    alpha = alpha,
    type.measure = type.measure,
    top.n = top.n,
    enforce.treatment = enforce.treatment,
    nfolds = nfolds,
    tune.ranger = tune.ranger,
    n.cv = n.cv,
    seed.set = seed.set,
    intercept = intercept
  )

  # 2. Prédiction externe ---------------------------------------------
  prediction.obj <- protein.regression.function(
    dataset = dataset,
    geneprotein = geneprotein,
    newdata = newdata, # toujours fourni
    weighted = weighted,
    weight.method = weight.method,
    weight.matrix = weight.matrix,
    weight.matrix.p = weight.matrix.p,
    grouping = grouping,
    manual.covar = manual.covar,
    design.m = design.m,
    scaled = scaled,
    RNA.log.scale = RNA.log.scale,
    na.process = na.process,
    include.treatment = include.treatment,
    treatment.info = treatment.info,
    duration.scale = duration.scale,
    PE.RNA.penalty = PE.RNA.penalty,
    treatment.penalty = treatment.penalty,
    baseline.interactions = baseline.interactions,
    collinearity.reduce = collinearity.reduce,
    rho.cutoff = rho.cutoff,
    manual.weights = manual.weights,
    weights.vector = weights.vector,
    model.method = model.method,
    lasso.fam = lasso.fam,
    alpha = alpha,
    type.measure = type.measure,
    top.n = top.n,
    enforce.treatment = enforce.treatment,
    nfolds = nfolds,
    tune.ranger = tune.ranger,
    n.cv = n.cv,
    seed.set = seed.set,
    intercept = intercept
  )

  # 3. Retourner les résultats ----------------------------------------
  return(list(
    model.obj = model.obj,
    prediction.obj = prediction.obj
  ))
}
