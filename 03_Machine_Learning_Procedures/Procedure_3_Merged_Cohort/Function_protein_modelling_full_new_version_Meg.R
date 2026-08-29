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
#source("Function_design.R")
source("Function_design_filtered.R")
#source("Function_design_newdata.R")
source("Function_adaptive_lasso.R")


# separate functions for either lasso, svr, rf model based on design (either full or loo, or cv)
# all 3 functions should return the same list 


# function to manually compute the log-likelihood of a LASSO model
log_likelihood <- function(y, y_pred) {
  n <- length(y)
  residuals <- y - y_pred
  sigma2 <- sum(residuals^2) / n
  logL <- -n/2 * log(2 * pi * sigma2) - 1/(2 * sigma2) * sum(residuals^2)
  return(logL)
}


# FUNCTION FOR MULTICOLLINEARITY REDUCTION
collinearityReduce <- function (mRNA, Rho_cutoff) {
  gene <- mRNA[,apply(mRNA!=0, 2, all)]
  row.names(gene) <- gene[,1]
  gene <- gene[,-c(1,2)]

  data_cor <- cor(gene, method = "pearson")
  data_highcor1 <- findCorrelation(data_cor, cutoff=Rho_cutoff)
  gene.temp <- gene[,-data_highcor1]
  
  return(gene.temp)
  
}



# BASELINE MODEL ---------------------------------------------------------------

protein.regression.baseline <- function(dataset,
                                      geneprotein,
                                      #form.obj.manual = NULL, #manual formula for model
                                      include.treatment = FALSE,
                                      intercept = TRUE,
                                      prediction.method = "loo", # loo, cv, full, random
                                      seed.set = NULL
                                      #loo = TRUE,
                                      #sample.seed = 23,
                                      # prediction.sample = "random", #prediction measures based on which mice (random or timepoint)
                                      # prediction.group = "month12_ccl4", 
                                      # test.size.strata = 1,
                                      # test.size.random = 6
                                      ){
  
  # sub-data set including the investigated gene
  #sub_data <- dataset[GeneProtein == geneprotein, ][complete.cases(dataset[GeneProtein == geneprotein, ]),]
  
  ############################################################################## #pour le baseline modele de merged data, car new data a aucun na
  sub_data <- dataset[GeneProtein == geneprotein, ]
  # Retirer uniquement les NA pour oil et ccl4
  sub_data <- sub_data[!(Treatment %in% c("oil", "ccl4") & is.na(ProteinIntensity)), ]
  ##############################################################################
  
  # extract gene, protein names and target vector
  gene <- str_replace(unique(sub_data$GeneSyn), "-",  ".")
  protein <- unique(sub_data$ProteinID)
  
  y <- sub_data$ProteinIntensity
  
  sub_data$TreatmentTotal <- factor(sub_data$TreatmentTotal, levels = c("month0_oil", 
                                                                        "month2_oil",  
                                                                        "month12_oil", 
                                                                        "month2_ccl4",
                                                                        "month6_ccl4",
                                                                        "month12_ccl4"))
  
  
  # argument to include the treatment variable
  if(include.treatment == TRUE){
    #form.obj <- ProteinIntensity ~ TreatmentTotal + GeneCount
    form.obj <- ProteinIntensity ~ Treatment + GeneCount
  }else{
    form.obj <- ProteinIntensity ~ GeneCount
  }
  

  
  if(prediction.method == "full"){
    # GAM/GLM
    gam.obj <- gam(form.obj, data = sub_data,
                   family = "gaussian", drop.intercept = !(intercept))
    
    # R^2
    rsq <- summary(gam.obj)$r.sq
    
    # explained deviance 
    #dev.expl <- 1-gam.obj$deviance / gam.obj$null.deviance
    dev.expl <- summary(gam.obj)$dev.expl
    
    # AIC and BIC
    gam.aic <- gam.obj$aic
    gam.bic <- BIC(gam.obj)
    
    ll_gam <- logLik(gam.obj) 
    
    k <- dim(model.matrix(gam.obj))[2] - 1
    
    
    newdat <- get_all_vars(form.obj, sub_data)
    
    prediction <- predict(gam.obj, newdat)
    
    #tLL <- gam.obj$null.deviance - gam.obj$deviance
    
    gam_list <- list(regression.object = gam.obj, 
                     gene = gene,
                     protein = protein,
                     logLikelihood = ll_gam, 
                     k = k,
                     adj.R.squared = rsq,
                     y = y,
                     prediction = prediction,
                     deviance.explained = dev.expl,
                     AIC = gam.aic,
                     BIC = gam.bic
    )
    
    return(gam_list)
    
  }else if(prediction.method == "loo"){ ## Train-Test-Split of the data with resulting prediction error -------
    
      loo.pred <- sapply(1:nrow(sub_data), function(x){
        
        testDT <- sub_data[x,]
        trainDT <- sub_data[-x,]
        
        n.test <- nrow(testDT)
        
        # PROBLEM: BEI OHNEHIN SCHON KLEINEN STICHPROBEN WIRD DER 
        # TRAIN-DATESATZ NOCH KLEINER, ERGO VERLIERT ER WEITER AN AUSSAGEKRAFT
        
        ## compute GAM/GLM from train dataset
        gam.obj <- gam(form.obj, data = trainDT,
                       family = "gaussian", drop.intercept = !(intercept))
        
        # test data for prediction
        #newdat <- data.frame(TreatmentTime = newDT$TreatmentTime, Treatment = newDT$Treatment)
        newdat <- get_all_vars(form.obj, testDT)
        
        prediction <- predict(gam.obj, newdat)
        
        return(prediction)
      })
      
      # compute the prediction error (RMSE)
      rmse <- sqrt(mean((loo.pred-y)^2))
      
      mae <- 1/length(y) * sum(abs(loo.pred-y))
      
      sd.y <- sd(y)
      
      spearman.correlation <- cor(loo.pred, y, method = "spearman", use = "complete.obs")
      pearson.correlation <- cor(loo.pred, y, method = "pearson", use = "complete.obs")
      
      #print(geneprotein)
      
      return(list(rmse = rmse,
                  mae = mae,
                  rmse.scaled = rmse/sd.y,
                  mae.scaled = mae/sd.y,
                  correlation.spearman = spearman.correlation,
                  correlation.pearson = pearson.correlation,
                  gene = gene,
                  protein = protein,
                  y = y,
                  prediction = loo.pred))
      
      
    }else if(prediction.method == "cv"){ #randomly drawn train-test-splits
      
     
      # function to split a vector of indices into n equally large groups
      chunk.fun <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
      
      # sample indices of the sub-data
      set.seed(seed.set)
      i.sample <- sample(1:nrow(sub_data), nrow(sub_data), replace = FALSE)
      
      split.list <- chunk.fun(i.sample, n = n.cv) 
      
      
      cv.pred <- lapply(split.list, function(x){
        
        test.data <- sub_data[x,]
        train.data <- sub_data[-x,]
        
        y.test <- y[x]
        #colnames(y.test) <- colnames(protein.vector)
        y.train <- y[-x]
        #colnames(y.train) <- colnames(y)
        
        gam.obj <- gam(form.obj, data = train.data,
                       family = "gaussian", drop.intercept = !(intercept))
        
        # prediction based on model
        newdat <- get_all_vars(form.obj, test.data)
        
        prediction <- predict(gam.obj, newdat)
        
        return(prediction)
        
      })
      
      # full prediction vector (reordered to original order)
      prediction.vector <- unlist(cv.pred)[order(i.sample)]
      
      # compute the prediction error (RMSE)
      rmse <- sqrt(mean((prediction.vector-y)^2))
      
      mae <- 1/length(y) * sum(abs(prediction.vector-y))
      
      sd.y <- sd(y)
      
      correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
      correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
      
      #print(geneprotein)
      
      return_list <- list(rmse = rmse,
                         mae = mae,
                         rmse.scaled = rmse/sd.y,
                         mae.scaled = mae/sd.y,
                         correlation.pearson = correlation.pearson, 
                         correlation.spearman = correlation.spearman, 
                         y = y,
                         gene = gene,
                         protein = protein,
                         prediction = prediction.vector)
      
      #names(lasso_list) <- protein 
      
      return(return_list)
    
    }else if(prediction.method == "random"){
      
      #TODO
      
    }
  
}




# LASSO MODEL ------------------------------------------------------------------

lasso.mod <- function(Design,
                      y,
                      w,
                      n,
                      gene = gene,
                      protein,
                      prediction.method = "loo", # loo, cv, full, random
                      lasso.fam = "gaussian",
                      alpha = 1, 
                      nfolds = 10,
                      #relaxed = TRUE,
                      #gamma = 0,
                      type.measure = "deviance",
                      n.cv = 10,
                      intercept = TRUE,
                      seed.set = NULL){
  
  #print(seed.set)
  
  
  if(prediction.method == "full"){
    
    ## Calculate lambda path (first get lambda_max):
    mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
    sx <- scale(Design, scale = apply(Design, 2, mysd))
    sy <- as.vector(scale(y, scale=mysd(y)))
    
    # own lambda sequence (currently not used)
    lambda_max <- max(abs(colSums(sx*sy)))/length(y)
    epsilon <- .0001
    K <- 100
    lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                                length.out = K)), digits = 10)
    
    # lasso CV
    set.seed(seed.set)
    lasso.obj <- cv.glmnet(y = y, 
                           x = as.matrix(Design), 
                           alpha = alpha, 
                           family = lasso.fam,
                           #lambda = lambdapath,
                           nfolds = nfolds, 
                           penalty.factor = w, 
                           type.measure = type.measure,
                           relax= TRUE,
                           gamma = 0,
                           intercept = intercept)
    
    
    ## fit lasso again on whole data set
    set.seed(seed.set)
    
    # final.lasso <- glmnet(y = y, 
    #                       x = as.matrix(Design), 
    #                       alpha = alpha, 
    #                       family = match.arg(lasso.fam), 
    #                       penalty.factor = w, 
    #                       relax= TRUE,
    #                       gamma = 0,
    #                       intercept = intercept,
    #                       lambda = lasso.obj$lambda)
    
    final.lasso <- do.call(glmnet,
                           args = list(y = y, 
                                       x = as.matrix(Design), 
                                       alpha = alpha, 
                                       family = as.character(lasso.fam), 
                                       penalty.factor = w, 
                                       relax= TRUE,
                                       gamma = 0,
                                       intercept = intercept,
                                       lambda = lasso.obj$lambda))
    
    
    # RELAXED LASSO 
    # (the commented part leads to errors regarding the lambda sequence for the relaxed lasso)
    # coef.glmnet_relaxed <- coef(final.lasso$relaxed)
    # #opt.glmnet_relaxed <- match(lasso.obj$relaxed$lambda.min, lasso.obj$lambda)
    # opt.glmnet_relaxed <- which(abs(lasso.obj$relaxed$lambda.min - final.lasso$lambda) == 
    #                               min(abs(lasso.obj$relaxed$lambda.min - final.lasso$lambda)))
    # 
    # 
    lambda.opt_relaxed <-  lasso.obj$relaxed$lambda.min
    # 
    # ## coefficients at optimal lambda
    # coef.glmnet.opt_relaxed <- coef.glmnet_relaxed[, opt.glmnet_relaxed]
    # 
    # # deviance ratio (for elastic net this is the R^2)
    # final.lasso$relaxed$dev.ratio[opt.glmnet_relaxed]
    # 
    # # prediction based on model
    # prediction.vector_relaxed <- predict(final.lasso$relaxed, as.matrix(Design), s = lambda.opt_relaxed)
    
    # THIS WORKS
    # directly via the CV object
    coef.relaxed <- coef(lasso.obj, s = "lambda.min")[which(coef(lasso.obj, s = "lambda.min") != 0)]
    names(coef.relaxed) <- rownames(coef(lasso.obj, s = "lambda.min"))[which(coef(lasso.obj, s = "lambda.min") != 0)]
    
    prediction.vector_relaxed <- predict(lasso.obj, newx = as.matrix(Design), s = "lambda.min")
    
     
    
    # UNRELAXED LASSO
    coef.glmnet <- coef(final.lasso)
    
    #opt.glmnet <- match(lasso.obj$lambda.min, lasso.obj$lambda)
    # workaround because of some internal glmnet discrepancies between lambda object and lambda.min
    opt.glmnet <- which(abs(lasso.obj$lambda.min - final.lasso$lambda) == min(abs(lasso.obj$lambda.min - final.lasso$lambda)))
    
    lambda.opt <-  lasso.obj$lambda.min
    
    ## coefficients at optimal lambda
    coef.glmnet.opt <- coef.glmnet[, opt.glmnet]
    
    # deviance ratio (for elastic net this is the R^2)
    final.lasso$dev.ratio[opt.glmnet]
    
    # prediction based on model
    prediction.vector <- predict(final.lasso, as.matrix(Design), s = lambda.opt)
    
    
  }else if(prediction.method == "loo"){
    
    loo.pred <- sapply(1:nrow(Design), function(x){
      
      test.data <- Design[x, , drop = FALSE]
      train.data <- Design[-x, , drop = FALSE] # ici j ai ajouter drop egale false pour le 3_4_24
      
      y.test <- y[x]
      names(y.test) <- protein
      y.train <- y[-x]
      #names(y.train) <- protein
      
      # Create lambda path (currently not used)
      mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
      sx <- scale(train.data, scale = apply(train.data, 2, mysd))
      sy <- as.vector(scale(y.train, scale=mysd(y.train)))
      
      lambda_max <- max(abs(colSums(sx*sy)), na.rm = TRUE)/length(y.train)
      epsilon <- .0001
      K <- 100
      lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                                  length.out = K)), digits = 10)
      
      # lasso CV object
      set.seed(seed.set)
      lasso.obj <- cv.glmnet(y = y.train, 
                             x = as.matrix(train.data), 
                             alpha = alpha, 
                             family = lasso.fam,
                             penalty.factor = w, 
                             #lambda = lambdapath,
                             nfolds = nfolds, 
                             type.measure = type.measure,
                             relax= TRUE,
                             gamma = 0,
                             intercept = intercept)
      
      # UNRELAXED LASSO
      ## fit lasso again on whole data set
      set.seed(seed.set)
      # final.lasso <- glmnet(y = y.train, 
      #                       x = as.matrix(train.data), 
      #                       alpha = alpha, 
      #                       family = lasso.fam, 
      #                       penalty.factor = w, 
      #                       relax= TRUE,
      #                       gamma = 0,
      #                       #thresh = 1e-12, maxit = 1e5, lambda.min.ratio = NULL, standardize = FALSE,
      #                       #dfmax = ncol(train.data), 
      #                       #pmax = ncol(train.data),
      #                       intercept = intercept,
      #                       lambda = lasso.obj$lambda)
      
      
      final.lasso <- do.call(glmnet,
                             args = list(y = y.train, 
                                         x = as.matrix(train.data), 
                                         alpha = alpha, 
                                         family = as.character(lasso.fam), 
                                         penalty.factor = w, 
                                         relax= TRUE,
                                         gamma = 0,
                                         intercept = intercept,
                                         lambda = lasso.obj$lambda))
      
      coef.glmnet <- coef(final.lasso)
      #opt.glmnet <- match(lasso.obj$lambda.min, final.lasso$lambda)
      opt.glmnet <- which(abs(lasso.obj$lambda.min - final.lasso$lambda) == min(abs(lasso.obj$lambda.min - final.lasso$lambda)))
      
      lambda.opt <- lasso.obj$lambda.min
      
      ## coefficients at optimal lambda
      coef.glmnet.opt <- coef.glmnet[, opt.glmnet]
      
      # deviance ratio (for elastic net this is the R^2)
      final.lasso$dev.ratio[opt.glmnet]
      
      # prediction based on model
      prediction <- predict(final.lasso, as.matrix(test.data), s = lambda.opt)
      
      
      # RELAXED LASSO
      
      # does not work:
      # coef.glmnet_relaxed <- coef(final.lasso$relaxed)
      # opt.glmnet_relaxed <- match(lasso.obj$relaxed$lambda.min, final.lasso$lambda)
      # 
      # lambda.opt_relaxed <- lasso.obj$relaxed$lambda.min
      # 
      # ## coefficients at optimal lambda
      # coef.glmnet.opt_relaxed <- coef.glmnet_relaxed[, opt.glmnet_relaxed]
      # 
      # # deviance ratio (for elastic net this is the R^2)
      # final.lasso_relaxed$relaxed$dev.ratio[opt.glmnet_relaxed]
      # 
      # # prediction based on model
      # prediction_relaxed <- predict(final.lasso_relaxed$relaxed, as.matrix(test.data), s = lambda.opt_relaxed)
      
      # directly via the CV object
      coef.relaxed <- coef(lasso.obj, s = "lambda.min")[which(coef(lasso.obj, s = "lambda.min") != 0)]
      names(coef.relaxed) <- rownames(coef(lasso.obj, s = "lambda.min"))[which(coef(lasso.obj, s = "lambda.min") != 0)]
      
      prediction_relaxed <- predict(lasso.obj, newx = as.matrix(test.data), s = "lambda.min")
      
      return(list(prediction = prediction,
                  prediction_relaxed = prediction_relaxed))
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(loo.pred[1,])
    prediction.vector_relaxed <- unlist(loo.pred[2,])
    
    
  }else if(prediction.method == "cv"){
    #TODO MUSS NOCH ANGEPASST WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # function to split a vector of indices into n equally large groups
    chunk.fun <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
    
    # sample indices of the Design data
    set.seed(seed.set)
    i.sample <- sample(1:nrow(Design), nrow(Design), replace = FALSE)
    
    split.list <- chunk.fun(i.sample, n = n.cv) 
    
    
    cv.pred <- lapply(split.list, function(x){
      
      test.data <- Design[x, , drop = FALSE] #
      train.data <- Design[-x, , drop = FALSE] # ici j ai ajouter drop egale false pour le 3_4_24
      
      y.test <- y[x]
      #colnames(y.test) <- colnames(protein.vector)
      y.train <- y[-x]
      #colnames(y.train) <- colnames(y)
      
      
      # Create lambda path
      mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
      sx <- scale(train.data, scale = apply(train.data, 2, mysd))
      sy <- as.vector(scale(y.train, scale=mysd(y.train)))
      
      lambda_max <- max(abs(colSums(sx*sy)))/length(y.train)
      epsilon <- .0001
      K <- 100
      lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                                  length.out = K)), digits = 10)
      
      
      set.seed(seed.set)
      lasso.obj <- cv.glmnet(y = y.train, 
                             x = as.matrix(train.data), 
                             alpha = alpha, 
                             family = lasso.fam,
                             #lambda = lambdapath,
                             nfolds = nfolds,
                             relax= TRUE,
                             gamma = 0,
                             type.measure = type.measure)
      
      
      ## fit lasso again on whole data set
      set.seed(seed.set)
      # final.lasso <- glmnet(y = y.train, 
      #                       x = as.matrix(train.data), 
      #                       alpha = alpha, 
      #                       family = lasso.fam, 
      #                       relax= TRUE,
      #                       gamma = 0,
      #                       lambda = lasso.obj$lambda)
      
      final.lasso <- do.call(glmnet,
                             args = list(y = y.train, 
                                         x = as.matrix(train.data), 
                                         alpha = alpha, 
                                         family = as.character(lasso.fam), 
                                         penalty.factor = w, 
                                         relax= TRUE,
                                         gamma = 0,
                                         intercept = intercept,
                                         lambda = lasso.obj$lambda))
      
      #opt.glmnet <- match(lasso.obj$lambda.min, lasso.obj$lambda)
      opt.glmnet <- which(abs(lasso.obj$lambda.min - final.lasso$lambda) == min(abs(lasso.obj$lambda.min - final.lasso$lambda)))
      
      lambda.opt <-  lasso.obj$lambda.min
      
      # deviance ratio (for elastic net this is the R^2)
      final.lasso$dev.ratio[opt.glmnet]
      
      # prediction based on model
      prediction <- predict(final.lasso, as.matrix(test.data), s = lambda.opt)
      
      return(prediction)
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(cv.pred)[order(i.sample)]
    
    
  }else if(prediction.method == "random"){
    
    #TODO
    
  }
  
  
  # what to return --------------------------------------------------------------
  
  
  if(prediction.method == "full"){
    
    variable.selection <- coef.glmnet.opt[which(coef.glmnet.opt != 0)]
    variable.selection_relaxed <- coef.relaxed[which(coef.relaxed != 0)]
    
    # for gaussian family adj. R2 is:
    r2 <- lasso.obj$glmnet.fit$dev.ratio[opt.glmnet]
    #r2_relaxed <- lasso.obj$glmnet.fit$relaxed$dev.ratio[opt.glmnet_relaxed]
    
    
    #(1 - lasso.obj$cvm/var(y))[opt.glmnet]
    
    ## compute adj. R^2, AIC and BIC
    k <- final.lasso$df[opt.glmnet] + 1
    n <- final.lasso$nobs
    
    #TODO adj. R2 ? (double check) ????
    adj.r2 <- 1 - (1-r2) * (n-1)/(n-k+1)
    
    #k_relaxed <- final.lasso$relaxed$df[opt.glmnet_relaxed] + 1
    n_relaxed <- final.lasso$relaxed$nobs
    
    #TODO adj. R2 ? (double check) ????
    adj.r2_relaxed <- 1 - (1-r2) * (n-1)/(n-k+1)
    
    # NOT CORRECT: 
    #tLL <- final.lasso$nulldev - deviance(final.lasso)[opt.glmnet]
    #lasso.aic <- -tLL+2*k#+2*k*(k+1)/(n-k-1) #correction term (AICc)
    #lasso.bic <-log(n)*k - tLL
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    correlation.pearson_relaxed <- cor(prediction.vector_relaxed, y, method = "pearson", use = "complete.obs")
    correlation.spearman_relaxed <- cor(prediction.vector_relaxed, y, method = "spearman", use = "complete.obs")
    
    # likelihood for the LASSO model
    ll_lasso <- log_likelihood(y, prediction.vector)
    ll_lasso_relaxed <- log_likelihood(y, prediction.vector_relaxed)
    
    k_ll <- length(variable.selection) + 1
    k_ll_relaxed <- length(variable.selection_relaxed) + 1
    
    lasso.aic <- 2*k_ll - 2 * ll_lasso
    lasso.bic <- log(n)*k_ll - 2*ll_lasso 
    
    lasso.aic_relaxed <- 2*k_ll_relaxed - 2 * ll_lasso_relaxed
    lasso.bic_relaxed <- log(n)*k_ll_relaxed - 2*ll_lasso_relaxed 
    
    
    lasso_list <- list(final.lasso = final.lasso,
                       cv.lasso.obj = lasso.obj, 
                       lambda.opt = lasso.obj$lambda.min,
                       gene = gene,
                       protein = protein,
                       y = y,
                       #Design = Design,
                       
                       #unrelaxed
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
                       
                       #relaxed
                       variable.selection_relaxed = variable.selection_relaxed, 
                       all.covar_relaxed = coef.relaxed,
                       #k_relaxed = k_relaxed,
                       predictions.insample_relaxed = prediction.vector_relaxed,
                       correlation.pearson_relaxed = correlation.pearson_relaxed, 
                       correlation.spearman_relaxed = correlation.spearman_relaxed, 
                       logLikelihood_relaxed = ll_lasso_relaxed,
                       #deviance.explained = final.lasso$relaxed$dev.ratio[opt.glmnet_relaxed],
                       adj.R.squared_relaxed = adj.r2_relaxed,
                       AIC_relaxed = lasso.aic_relaxed,
                       BIC_relaxed = lasso.bic_relaxed
                       
                       )
    
    #names(lasso_list) <- protein 
    
    return(lasso_list)
    
    
  }else if(prediction.method == "loo"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    rmse_relaxed <- sqrt(mean((prediction.vector_relaxed-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    mae_relaxed <- 1/length(y) * sum(abs(prediction.vector_relaxed-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    correlation.pearson_relaxed <- cor(prediction.vector_relaxed, y, method = "pearson", use = "complete.obs")
    correlation.spearman_relaxed <- cor(prediction.vector_relaxed, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    lasso_list <- list(y = y,
                       gene = gene,
                       protein = protein,
                       
                       # unrelaxed
                       rmse = rmse,
                       mae = mae,
                       rmse.scaled = rmse/sd.y,
                       mae.scaled = mae/sd.y,
                       correlation.pearson = correlation.pearson, 
                       correlation.spearman = correlation.spearman,
                       prediction = prediction.vector,
                       
                       # relaxed
                       rmse_relaxed = rmse_relaxed,
                       mae_relaxed = mae_relaxed,
                       rmse.scaled_relaxed = rmse_relaxed/sd.y,
                       mae.scaled_relaxed = mae_relaxed/sd.y,
                       correlation.pearson_relaxed = correlation.pearson_relaxed, 
                       correlation.spearman_relaxed = correlation.spearman_relaxed, 
                       prediction_relaxed = prediction.vector_relaxed
                       )
    
    #names(lasso_list) <- protein 
    
    return(lasso_list)
    
    
  }else if(prediction.method == "cv"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    rmse_relaxed <- sqrt(mean((prediction.vector_relaxed-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    mae_relaxed <- 1/length(y) * sum(abs(prediction.vector_relaxed-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    correlation.pearson_relaxed <- cor(prediction.vector_relaxed, y, method = "pearson", use = "complete.obs")
    correlation.spearman_relaxed <- cor(prediction.vector_relaxed, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    lasso_list <- list(y = y,
                       gene = gene,
                       protein = protein,
                       
                       # unrelaxed
                       rmse = rmse,
                       mae = mae,
                       rmse.scaled = rmse/sd.y,
                       mae.scaled = mae/sd.y,
                       correlation.pearson = correlation.pearson, 
                       correlation.spearman = correlation.spearman,
                       prediction = prediction.vector,
                       
                       # relaxed
                       rmse_relaxed = rmse_relaxed,
                       mae_relaxed = mae_relaxed,
                       rmse.scaled_relaxed = rmse_relaxed/sd.y,
                       mae.scaled_relaxed = mae_relaxed/sd.y,
                       correlation.pearson_relaxed = correlation.pearson_relaxed, 
                       correlation.spearman_relaxed = correlation.spearman_relaxed, 
                       prediction_relaxed = prediction.vector_relaxed
    )
    
    #names(lasso_list) <- protein 
    
    return(lasso_list)
    
    
  }else if(prediction.method == "random"){
    #TODO
  }
  
  
}


# Support Vector regression (comparable to Prabahar et al.)
svr.mod <- function(Design,
                    y,
                    w,
                    n,
                    gene,
                    protein,
                    prediction.method = "loo", # loo, cv, full, random
                    nfolds = 5,
                    n.cv = 10,
                    seed.set = NULL){
  
  
  if(prediction.method == "full"){
    
    df.train <- data.frame(Design, y = y)
    ctrl <- trainControl(method = "CV", number = nfolds, savePred=TRUE) 
    
    model.obj <- caret::train(y ~., data = df.train, method = "svmRadial", trControl = ctrl)
    
    #model.obj <- trainOpt(data_out = mRNA.data, protein.vector = protein.vector)
    
    # prediction based on model
    prediction.vector <- predict(model.obj, Design)
    
    
  }else if(prediction.method == "loo"){
    
    loo.pred <- sapply(1:nrow(Design), function(x){
      
      test.data <- Design[x,]
      train.data <- Design[-x,]
      
      y.test <- y[x]
      #names(y.test) <- protein
      y.train <- y[-x]
      #names(y.train) <- protein
      
      
      df.train <- data.frame(train.data, y = y.train)
      ctrl <- trainControl(method = "CV", number = nfolds, savePred=TRUE) 
      
      model.cv <- caret::train(y ~., data = df.train, method = "svmRadial", trControl = ctrl)
      #model.cv <- trainOpt(data_out = train.data, protein.vector = y.train)
      
      # prediction based on model
      prediction <- predict(model.cv, test.data)
      
      return(prediction)
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(loo.pred)
    
    
  }else if(prediction.method == "cv"){
    
    # function to split a vector of indices into n equally large groups
    chunk.fun <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
    
    # sample indices of the Design data
    set.seed(seed.set)
    i.sample <- sample(1:nrow(Design), nrow(Design), replace = FALSE)
    
    split.list <- chunk.fun(i.sample, n = n.cv) 
    
    
    cv.pred <- lapply(split.list, function(x){
      
      test.data <- Design[x,]
      train.data <- Design[-x,]
      
      y.test <- y[x]
      #colnames(y.test) <- colnames(protein.vector)
      y.train <- y[-x]
      #colnames(y.train) <- colnames(y)
      
      
      df.train <- data.frame(train.data, y = y.train)
      ctrl <- trainControl(method = "CV", number = nfolds, savePred=TRUE) 
      
      set.seed(seed.set)
      model.cv <- caret::train(y ~., data = df.train, method = "svmRadial", trControl = ctrl)
      
      # prediction based on model
      prediction <- predict(model.cv, test.data)
      
      return(prediction)
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(cv.pred)[order(i.sample)]
    
    
  }else if(prediction.method == "random"){
    
    #TODO
    
  }
  
  
  # what to return --------------------------------------------------------------
  
  if(prediction.method == "full"){
    
    #variable.selection <- coef.glmnet.opt[which(coef.glmnet.opt != 0)]
    
    # for gaussian family adj. R2 is:
    i.tuned <- which(model.obj$results$C == model.obj$bestTune$C & model.obj$results$sigma == model.obj$bestTune$sigma)
    
    model.results <- model.obj$results[i.tuned,]
    
    r2 <- model.results$Rsquared
    
    rmse <- sqrt(mean((prediction.vector-y)^2))
    #rmse <- model.results$rmse
    
    mae <- model.results$MAE
    
    # Berechnung der Residuen
    residuals <- y - prediction.vector
    rss <- sum(residuals^2)  # Residual Sum of Squares
    
    # Schätzung der effektiven Anzahl der Parameter
    sigma2 <- var(residuals)  # Varianz der Residuen
    edf <- sum((residuals / sigma2)^2)  # Effektive Freiheitsgrade (grobe Schätzung)
    
    # AIC und BIC berechnen
    n <- length(y)  # Anzahl der Beobachtungen
    svr.aic <- n * log(rss/n) + 2 * edf
    svr.bic <- n * log(rss/n) + log(n) * edf
    
    
    #R2(prediction.vector, y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    # likelihood for the LASSO model
    ll_svr <- log_likelihood(y, prediction.vector)
    
    
    svr_list <- list(gene = gene,
                     protein = protein,
                     #Design = Design,
                     predictions.insample = prediction.vector,
                     correlation.pearson = correlation.pearson, 
                     correlation.spearman = correlation.spearman, 
                     y = y,
                     logLikelihood = ll_svr,
                     #deviance.explained = final.lasso$dev.ratio[opt.glmnet],
                     adj.R.squared = r2,
                     AIC = svr.aic,
                     BIC = svr.bic)
    
    #names(lasso_list) <- protein 
    
    return(svr_list)
    
    
  }else if(prediction.method == "loo"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    svr_list <- list(rmse = rmse,
                     mae = mae,
                     rmse.scaled = rmse/sd.y,
                     mae.scaled = mae/sd.y,
                     correlation.pearson = correlation.pearson, 
                     correlation.spearman = correlation.spearman, 
                     y = y,
                     gene = gene,
                     protein = protein,
                     prediction = prediction.vector)
    
    #names(lasso_list) <- protein 
    
    return(svr_list)
    
    
  }else if(prediction.method == "cv"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    svr_list <- list(rmse = rmse,
                     mae = mae,
                     rmse.scaled = rmse/sd.y,
                     mae.scaled = mae/sd.y,
                     correlation.pearson = correlation.pearson, 
                     correlation.spearman = correlation.spearman, 
                     y = y,
                     gene = gene,
                     protein = protein,
                     prediction = prediction.vector)
    
    #names(lasso_list) <- protein 
    
    return(svr_list)
    
    
  }else if(prediction.method == "random"){
    #TODO
  }
  
  
}


# Random forest model
rf.mod <- function(Design,
                   y,
                   w,
                   n,
                   gene,
                   protein,
                   prediction.method = "loo", # loo, cv, full, random
                   tune.ranger = FALSE,
                   #num.trees,
                   #mtry,
                   #min.node.size,
                   #sample.fraction,
                   n.cv = 10,
                   seed.set = 23){
  
  
  
  
  if(prediction.method == "full"){
    
    DTrf <- data.frame(Design, y = y)
  
    rf.task <- makeRegrTask(data = DTrf, target = "y")
    
    if(tune.ranger == TRUE){
      
      set.seed(seed.set)
      rf.tune = tuneRanger(rf.task, measure = list(mse), num.trees = 1000, 
                           num.threads = 2, iters = 70, save.file.path = NULL,
                           show.info = TRUE)
      
      # ranger with tuned params
      fit.ranger <- ranger(x = Design, y = y, data = DTrf, 
                           importance = "impurity_corrected", 
                           probability = FALSE, classification = FALSE,
                           num.trees = 500, oob.error = TRUE,
                           mtry = rf.tune$recommended.pars$mtry, 
                           min.node.size = rf.tune$recommended.pars$min.node.size, 
                           sample.fraction = rf.tune$recommended.pars$sample.fraction)
      
    }else{
      
      # ranger without tuned params
      set.seed(seed.set)
      fit.ranger <- ranger(x = Design, y = y, data = DTrf, 
                           importance = "impurity_corrected", 
                           probability = FALSE, classification = FALSE,
                           num.trees = 500, oob.error = TRUE,
                           mtry = ncol(Design)/3, 
                           min.node.size = 5#, 
                           #sample.fraction = rf.tune$recommended.pars$sample.fraction
                           )
      
    }
    
    
    #print(fit.ranger)
    #summary(fit.ranger)
    
    prediction <- predict(fit.ranger, Design)
    
    prediction.vector <- prediction$predictions
    
    
  }else if(prediction.method == "loo"){
    
    loo.pred <- sapply(1:nrow(Design), function(x){
      
      test.data <- Design[x,]
      train.data <- Design[-x,]
      
      y.test <- y[x]
      names(y.test) <- protein
      y.train <- y[-x]
      #names(y.train) <- protein
      
      DTrf <- data.frame(train.data, y = y.train)
      
      rf.task <- makeRegrTask(data = DTrf, target = "y")
      
      if(tune.ranger == TRUE){
        
        set.seed(seed.set)
        rf.tune = tuneRanger(rf.task, measure = list(mse), num.trees = 1000, 
                             num.threads = 2, iters = 70, save.file.path = NULL,
                             show.info = TRUE)
        
        # ranger with tuned params
        fit.ranger <- ranger(x = train.data, y = y.train, data = DTrf, 
                             importance = "impurity_corrected", 
                             probability = FALSE, classification = FALSE,
                             num.trees = 500, oob.error = TRUE,
                             mtry = rf.tune$recommended.pars$mtry, 
                             min.node.size = rf.tune$recommended.pars$min.node.size, 
                             sample.fraction = rf.tune$recommended.pars$sample.fraction)
        
      }else{
        
        # ranger without tuned params
        fit.ranger <- ranger(x = train.data, y = y.train, data = DTrf, 
                             importance = "impurity_corrected", 
                             probability = FALSE, classification = FALSE,
                             num.trees = 500, oob.error = TRUE,
                             mtry = ncol(Design)/3, 
                             min.node.size = 5#, 
                             #sample.fraction = rf.tune$recommended.pars$sample.fraction
        )
        
      }
      
      
      #print(fit.ranger)
      #summary(fit.ranger)
      
      prediction <- predict(fit.ranger, test.data)
      
      prediction.value <- prediction$predictions
      
      return(prediction.value)
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(loo.pred)
    
    
  }else if(prediction.method == "cv"){
    
    # function to split a vector of indices into n equally large groups
    chunk.fun <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
    
    # sample indices of the Design data
    set.seed(seed.set)
    i.sample <- sample(1:nrow(Design), nrow(Design), replace = FALSE)
    
    split.list <- chunk.fun(i.sample, n = n.cv) 
    
    
    cv.pred <- lapply(split.list, function(x){
      
      test.data <- Design[x,]
      train.data <- Design[-x,]
      
      y.test <- y[x]
      #colnames(y.test) <- colnames(protein.vector)
      y.train <- y[-x]
      #colnames(y.train) <- colnames(y)
      
      
      DTrf <- data.frame(train.data, y = y.train)
      
      rf.task <- makeRegrTask(data = DTrf, target = "y")
      
      if(tune.ranger == TRUE){
        
        set.seed(seed.set)
        rf.tune = tuneRanger(rf.task, measure = list(mse), num.trees = 1000, 
                             num.threads = 2, iters = 70, save.file.path = NULL,
                             show.info = TRUE)
        
        # ranger with tuned params
        fit.ranger <- ranger(x = train.data, y = y.train, data = DTrf, 
                             importance = "impurity_corrected", 
                             probability = FALSE, classification = FALSE,
                             num.trees = 500, oob.error = TRUE,
                             mtry = rf.tune$recommended.pars$mtry, 
                             min.node.size = rf.tune$recommended.pars$min.node.size, 
                             sample.fraction = rf.tune$recommended.pars$sample.fraction)
        
      }else{
        
        # ranger without tuned params
        fit.ranger <- ranger(x = train.data, y = y.train, data = DTrf, 
                             importance = "impurity_corrected", 
                             probability = FALSE, classification = FALSE,
                             num.trees = 500, oob.error = TRUE,
                             mtry = ncol(Design)/3, 
                             min.node.size = 5#, 
                             #sample.fraction = rf.tune$recommended.pars$sample.fraction
        )
        
      }
      
      
      #print(fit.ranger)
      #summary(fit.ranger)
      
      prediction <- predict(fit.ranger, test.data)
      
      prediction.value <- prediction$predictions
      
      return(prediction.value)
      
    })
    
    # full prediction vector (reordered to original order)
    prediction.vector <- unlist(cv.pred)[order(i.sample)]
    
    
  }else if(prediction.method == "random"){
    
    #TODO
    
  }
  
  
  # what to return --------------------------------------------------------------
  
  if(prediction.method == "full"){
    
    #importance_pvalues(fit.ranger)
    variable.importance <- importance(fit.ranger)
    
    #fit.ranger$prediction.error
    
    r2 <- fit.ranger$r.squared
    
    ## Adjustiertes R² berechnen
    n <- length(y)  # Anzahl der Beobachtungen
    p <- ncol(Design) - 1  # Anzahl der Prädiktoren
    
    adj.r2 <- 1 - ((1 - r2) * (n - 1) / (n - p - 1))

    
    # NOT CORRECT: 
    #tLL <- final.lasso$nulldev - deviance(final.lasso)[opt.glmnet]
    #lasso.aic <- -tLL+2*k#+2*k*(k+1)/(n-k-1) #correction term (AICc)
    #lasso.bic <-log(n)*k - tLL
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    # likelihood for the LASSO model
    ll_rf <- log_likelihood(y, prediction.vector)
    #p <- ncol(Design)
    #lasso.aic <- 2* p - 2 * ll_rf
    #lasso.bic <- log(n)* p - 2*ll_rf 
    
    # Varianz der Residuen schätzen
    sigma2 <- mean((y - prediction.vector)^2)
    
    # Anzahl der Parameter (p + 1 für Intercept)
    k <- p + 1
    
    # AIC und BIC berechnen
    aic <- n * log(sigma2) + 2 * k
    bic <- n * log(sigma2) + log(n) * k
    
    
    rf_list <- list(rf.model = fit.ranger,
                    gene = gene,
                       protein = protein,
                       #Design = Design,
                       variable.importance = variable.importance, 
                       #all.covar = coef.glmnet.opt,
                       k = k,
                       predictions.insample = prediction.vector,
                       correlation.pearson = correlation.pearson, 
                       correlation.spearman = correlation.spearman, 
                       y = y,
                       logLikelihood = ll_rf,
                       #deviance.explained = final.lasso$dev.ratio[opt.glmnet],
                       adj.R.squared = adj.r2,
                       AIC = aic,
                       BIC = bic)
    
    #names(lasso_list) <- protein 
    
    return(rf_list)
    
    
  }else if(prediction.method == "loo"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    rf_list <- list(rmse = rmse,
                       mae = mae,
                       rmse.scaled = rmse/sd.y,
                       mae.scaled = mae/sd.y,
                       correlation.pearson = correlation.pearson, 
                       correlation.spearman = correlation.spearman, 
                       y = y,
                       gene = gene,
                       protein = protein,
                       prediction = prediction.vector)
    
    #names(lasso_list) <- protein 
    
    return(rf_list)
    
    
  }else if(prediction.method == "cv"){
    
    # compute the prediction error (RMSE)
    rmse <- sqrt(mean((prediction.vector-y)^2))
    
    mae <- 1/length(y) * sum(abs(prediction.vector-y))
    
    sd.y <- sd(y)
    
    correlation.pearson <- cor(prediction.vector, y, method = "pearson", use = "complete.obs")
    correlation.spearman <- cor(prediction.vector, y, method = "spearman", use = "complete.obs")
    
    #print(geneprotein)
    
    rf_list <- list(rmse = rmse,
                       mae = mae,
                       rmse.scaled = rmse/sd.y,
                       mae.scaled = mae/sd.y,
                       correlation.pearson = correlation.pearson, 
                       correlation.spearman = correlation.spearman, 
                       y = y,
                       gene = gene,
                       protein = protein,
                       prediction = prediction.vector)
    
    #names(lasso_list) <- protein 
    
    return(rf_list)
    
    
  }else if(prediction.method == "random"){
    #TODO
  }
 
  
}


# FUNCTION TO PRESELECT COVARIATES FROM RF FEATURE IMPORTANCE ------------------
rf.preselect <- function(Design,
                         y,
                         w,
                         n,
                         gene = gene,
                         protein,
                         prediction.method = "loo", # loo, cv, full, random
                         seed.set = 23,
                         
                         #always.split = TRUE, always split for the unpenalized covariates
                         tune.ranger = FALSE,
                         #num.trees,
                         #mtry,
                         #min.node.size,
                         #sample.fraction,
                         #fraction of top importances or flat amount or percentage threshold?
                         top.n = 30){
  
  DTrf <- data.frame(Design, y = y)
  
  rf.task <- makeRegrTask(data = DTrf, target = "y")
  
  # ranger without tuned params
  set.seed(seed.set)
  fit.ranger <- ranger(x = Design, y = y, data = DTrf,
                       importance = "impurity_corrected", 
                       probability = FALSE, classification = FALSE,
                       num.trees = 2500, oob.error = TRUE,
                       mtry = floor(ncol(Design)/3), 
                       #always.split.variables = names(w)[which(w == 0)],
                       min.node.size = 5#, 
                       #sample.fraction = rf.tune$recommended.pars$sample.fraction
  )
  
  top.important.features <- sort(importance(fit.ranger), decreasing = TRUE)[1:top.n]
  
  return(top.important.features)
  
}



# ------------------------------------------------------------------------------
# FUNCTION THAT COVERS ALL ABOVE MODEL APPROACHES WITH ALL POSSIBLE DESIGN APPROACHES

protein.regression.function <- function(
  # Arguments for the Design.object:
  dataset,
  geneprotein,
  #groups, #vector of group membership
  weighted = FALSE,
  weight.method = NULL, # either "gca" similarity or "dipa" distance
  weight.matrix = NULL, # associated matrix of weights (from co-expression analysis or from dipa analysis)
  weight.matrix.p = NULL, # needed for GCA weights and combined design
  grouping = "all", # "all" genes as covariates, or "dipa", "coexpression"
  scaled = TRUE,
  RNA.log.scale = TRUE,
  manual.covar = NULL, 
  design.m = "genes", 
  na.process = "max.n", # complete cases ("cc"), maximum n ("max.n")
  include.treatment = FALSE,
  treatment.info = "full", #full, treatment, treatment.duration, separate
  duration.scale = "factor", #duration as a factor or numeric covariate
  treatment.penalty = FALSE, 
  PE.RNA.penalty = FALSE, # penalty for the protein encoding RNA
  baseline.interactions = FALSE,
  
  # Collinearity reduction:
  collinearity.reduce = FALSE,
  rho.cutoff = 0.75,
  
  # Arguments for the models:
  manual.weights = FALSE,
  weights.vector = NULL,
  model.method = "lasso", # lasso, svr, rf, lasso.rf.pre, rf.rf.pre, adaptive.lasso OR BASELINE
  prediction.method = "loo", # loo, cv, full, random
  
  # lasso specific
  lasso.fam = "gaussian", # family argument of glmnet
  type.measure = "deviance", # type.measure argument of glmnet
  alpha = 1,
  
  # lasso with rf preselection
  top.n = 30,
  enforce.treatment = FALSE,
  
  # lasso and svr
  nfolds = 10,
  
  # random forest specific
  tune.ranger = FALSE,
  
  n.cv = 10,
  
  #
  # prediction.sample = "random", #prediction measures based on which mice (random or timepoint)
  # prediction.group = "month12_ccl4", 
  # test.size.strata = 1,
  # test.size.random = 6, 
  seed.set = 23,
  sample.seed = NULL,
  intercept = TRUE){
  
  # create Design object 
  Design.obj <- design.function(dataset = dataset,
                                geneprotein = geneprotein,
                                #geneprotein_list = pairs.list[1:2]
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
                                baseline.interactions = baseline.interactions)
  
  # extract information from the Design object
  gene <- Design.obj$gene
  protein <- Design.obj$protein
  y <- Design.obj$y
  Design <- Design.obj$Design
  n <- Design.obj$n
  protein.vector <- Design.obj$protein.vector
  mice <- Design.obj$mice
  experiments <- Design.obj$experiments
  #attach(Design.obj)
  
  # manual weights given as weights vector as function argument or extracted from the Design object
  if(manual.weights == FALSE){
    
    w <- Design.obj$w
    
  }else{
    
    w <- weights.vector
    
  }
  
  
  #TODO work-in-progress: attempt to reduce multicollinearity in the variable selection
  #     by reducing collinearity in the design matrix beforehand
  if(collinearity.reduce  == TRUE){
    
    # Multicollinearity reduction
    Design_temp <- collinearityReduce(Design, rho.cutoff)
    w_temp <- w[match(colnames(Design_temp), names(w))]
    
    #make sure that the associated gene is in the design matrix
    i.g <- which(colnames(Design) == paste0(gene, "_G"))
    Design <- cbind(Design[,i.g], Design_temp)
    w <- c(w[i.g], w_temp)
    
    names(Design)[1] <- paste0(gene, "_G")
    
  }
  
  rf_design_cols <- ncol(Design)
  
  # BASELINE MODEL -------------------------------------------------------------
  
  
  baseline.model <- protein.regression.baseline(dataset = dataset,
                                                geneprotein = geneprotein,
                                                #form.obj.manual = NULL, #manual formula for model
                                                include.treatment = include.treatment,
                                                intercept = intercept,
                                                prediction.method = prediction.method# loo, cv, full, random
                                                )
  
  
  # MODEL ----------------------------------------------------------------------
  # cover different model approaches
  
  if(model.method == "baseline.only"){
    
    model <- baseline.model
    
  }else if(model.method == "lasso"){
    
    model <- lasso.mod(Design = Design,
                       y = y,
                       w = w,
                          n = n,
                       gene = gene,
                          protein = protein,
                          prediction.method = prediction.method, # loo, cv, full, random
                       alpha = alpha,
                          lasso.fam = lasso.fam,
                          nfolds = nfolds,
                          type.measure = type.measure,
                          n.cv = n.cv,
                          intercept = intercept,
                          seed.set = seed.set)
    
    
  }else if(model.method  == "adaptive.lasso"){
    
    model <- adaptive.lasso.mod(Design = Design,
                       y = y,
                       w = w,
                       n = n,
                       gene = gene,
                       protein = protein,
                       prediction.method = prediction.method, # loo, cv, full, random
                       #alpha = alpha,
                       lasso.fam = lasso.fam,
                       nfolds = nfolds,
                       type.measure = type.measure,
                       n.cv = n.cv,
                       intercept = intercept,
                       seed.set = seed.set)
    
    
  }else if(model.method == "lasso.rf.pre"){
    if (design.m == "all" & !is.null(manual.covar)) {
      rna_covars <- colnames(Design)[grep("_G$", colnames(Design))]
      variables_to_keep <- c(manual.covar, rna_covars)
      cols_to_keep <- intersect(variables_to_keep, colnames(Design))
      Design.rf.input <- Design[, cols_to_keep, drop = FALSE]
      w.rf.input <- w[cols_to_keep]
      rf_design_cols <- ncol(Design.rf.input)
    } else {
      Design.rf.input <- Design
      w.rf.input <- w
    }

    rf.pre.model <- rf.preselect(Design = Design.rf.input,
                                 y = y,
                                 w = w.rf.input,
                                 #num.trees,
                                 #mtry,
                                 #min.node.size,
                                 #sample.fraction,
                                 #n.cv = n.cv,
                                 top.n = top.n,
                                 seed.set = seed.set)
    
    
    # enforce treatment = include treatment variable although not pre-selected by RF
    if(enforce.treatment ==TRUE & include.treatment == TRUE){
      
      # add the treatment variables to the lasso covariates
      unpenalized.variables <- names(w)[which(w == 0)]
      
      # adapt the input for the lasso with rf preselection 
      Design.rf <- Design[,colnames(Design) %in% names(rf.pre.model) | colnames(Design) %in% unpenalized.variables]
      w.rf <- w[names(w) %in% names(rf.pre.model) | names(w) %in% unpenalized.variables]
      
      model <- lasso.mod(Design = Design.rf,
                         y = y,
                         w = w.rf,
                         n = n,
                         gene = gene,
                         protein = protein,
                         prediction.method = prediction.method, # loo, cv, full, random
                         alpha = alpha,
                         lasso.fam = lasso.fam,
                         nfolds = nfolds,
                         type.measure = type.measure,
                         n.cv = n.cv,
                         intercept = intercept,
                         seed.set = seed.set)
      
    }else{
      
      # adapt the input for the lasso with rf preselection 
      Design.rf <- Design[,colnames(Design) %in% names(rf.pre.model) | colnames(Design) %in% paste0(gene, "_G")]
      w.rf <- w[names(w) %in% names(rf.pre.model) | names(w) %in% paste0(gene, "_G")]
      
      model <- lasso.mod(Design = Design.rf,
                         y = y,
                         w = w.rf,
                         n = n,
                         gene = gene,
                         protein = protein,
                         prediction.method = prediction.method, # loo, cv, full, random
                         lasso.fam = lasso.fam,
                         nfolds = nfolds,
                         type.measure = type.measure,
                         n.cv = n.cv,
                         intercept = intercept,
                         seed.set = seed.set)
    }
    
  }else if(model.method == "lasso.predefined"){
    
    if (is.null(manual.covar)) {
      stop("Pour utiliser le modele 'lasso.predefined', vous devez fournir une liste de variables dans l'argument 'manual.covar'.")
    }

    # 1. On s'assure que le gène propre de la protéine est forcé dans les variables sélectionnées
    target_gene_name <- paste0(gene, "_G")
    if (design.m == "all") {
      rna_covars <- colnames(Design)[grep("_G$", colnames(Design))]
      variables_to_keep <- c(manual.covar, rna_covars)
    } else {
      variables_to_keep <- c(manual.covar, target_gene_name)
    }

    # enforce treatment = on inclut aussi la variable de traitement si demandé (où penalité w == 0)
    if (enforce.treatment == TRUE & include.treatment == TRUE) {
      unpenalized.variables <- names(w)[which(w == 0)]
      variables_to_keep <- c(variables_to_keep, unpenalized.variables)
    }

    # 2. On trouve l'intersection stricte entre les colonnes de notre design matrix et notre liste
    cols_to_keep <- intersect(variables_to_keep, colnames(Design))
    
    # glmnet crash prevention: it requires >= 2 columns
    if (length(cols_to_keep) < 2) {
      model <- "Aucune Mastery Protein trouvee dans ce cluster."
    } else {
      Design.predef <- Design[, cols_to_keep, drop = FALSE]
      w.predef <- w[cols_to_keep]

      # 3. On passe ces données au modèle Lasso classique
      model <- lasso.mod(Design = Design.predef,
                         y = y,
                         w = w.predef,
                         n = n,
                         gene = gene,
                         protein = protein,
                         prediction.method = prediction.method, # loo, cv, full, random
                         alpha = alpha,
                         lasso.fam = lasso.fam,
                         nfolds = nfolds,
                         type.measure = type.measure,
                         n.cv = n.cv,
                         intercept = intercept,
                         seed.set = seed.set)
    }
    
  }else if(model.method == "svr"){
    
    model <- svr.mod(Design = Design,
                       y = y,
                       #w = w,
                       n = n,
                     gene = gene,
                       protein = protein,
                       prediction.method = prediction.method, # loo, cv, full, random
                       nfolds = nfolds,
                       n.cv = n.cv,
                       seed.set = seed.set)
    
  }else if(model.method == "rf"){
  
    model <- rf.mod(Design = Design,
                    y = y,
                    #w = w,
                    n = n,
                    gene = gene,
                    protein = protein,
                    prediction.method = prediction.method, # loo, cv, full, random
                    tune.ranger = tune.ranger,
                    #num.trees,
                    #mtry,
                    #min.node.size,
                    #sample.fraction,
                    n.cv = n.cv,
                    seed.set = seed.set)
    
  }else if(model.method == "rf.rf.pre"){
    
    if (design.m == "all" & !is.null(manual.covar)) {
      rna_covars <- colnames(Design)[grep("_G$", colnames(Design))]
      variables_to_keep <- c(manual.covar, rna_covars)
      cols_to_keep <- intersect(variables_to_keep, colnames(Design))
      Design.rf.input <- Design[, cols_to_keep, drop = FALSE]
      w.rf.input <- w[cols_to_keep]
      rf_design_cols <- ncol(Design.rf.input)
    } else {
      Design.rf.input <- Design
      w.rf.input <- w
    }

    rf.pre.model <- rf.preselect(Design = Design.rf.input,
                                 y = y,
                                 w = w.rf.input,
                                 top.n = top.n,
                                 seed.set = seed.set)
    
    # enforce treatment = include treatment variable although not pre-selected by RF
    if(enforce.treatment == TRUE & include.treatment == TRUE){
      unpenalized.variables <- names(w)[which(w == 0)]
      selected_cols <- colnames(Design)[colnames(Design) %in% names(rf.pre.model) | colnames(Design) %in% unpenalized.variables]
    } else {
      selected_cols <- colnames(Design)[colnames(Design) %in% names(rf.pre.model) | colnames(Design) %in% paste0(gene, "_G")]
    }
    
    Design.rf <- Design[, selected_cols, drop = FALSE]
    
    model <- rf.mod(Design = Design.rf,
                    y = y,
                    w = w[selected_cols],
                    n = n,
                    gene = gene,
                    protein = protein,
                    prediction.method = prediction.method, # loo, cv, full, random
                    tune.ranger = tune.ranger,
                    n.cv = n.cv,
                    seed.set = seed.set)
    
  }
  
  
  # COMPARISON TO BASELINE -----------------------------------------------------
   # ----
  
  

  return(list(model = model,
              baseline.model = baseline.model,
              experiments = experiments,
              rf_mtry = if (model.method %in% c("lasso.rf.pre", "rf.rf.pre", "rf", "lasso.predefined")) floor(rf_design_cols / 3) else NULL))

  
}





# COMPLETE MODEL ---------------------------------------------------------------

# (prediction or only full data model)

protein.regression.complete <- function(
  # Arguments for the Design.object:
  dataset,
  geneprotein,
  #groups, #vector of group membership
  weighted = FALSE,
  weight.method = NULL, # either "gca" similarity or "dipa" distance
  weight.matrix = NULL, # associated matrix of weights (from co-expression analysis or from dipa analysis)
  weight.matrix.p = NULL, # needed for GCA weights and combined design
  grouping = "all", # "all" genes as covariates, or "dipa", "coexpression"
  scaled = TRUE,
  RNA.log.scale = TRUE,
  manual.covar = NULL, 
  design.m = "genes", 
  na.process = "max.n", # complete cases ("cc"), maximum n ("max.n")
  include.treatment = FALSE,
  treatment.info = "full", #full, treatment, treatment.duration, separate
  duration.scale = "factor", #duration as a factor or numeric covariate
  treatment.penalty = FALSE, 
  PE.RNA.penalty = FALSE, # penalty for the protein encoding RNA
  baseline.interactions = FALSE,
  
  prediction = TRUE,
  
  # Collinearity reduction:
  collinearity.reduce = FALSE,
  rho.cutoff = 0.75,
  
  # Arguments for the models:
  manual.weights = FALSE,
  weights.vector = NULL,
  model.method = "lasso", # lasso, svr, rf
  prediction.method = "loo", # loo, cv, full, random
  
  # lasso specific
  lasso.fam = "gaussian", # family argument of glmnet
  type.measure = "deviance", # type.measure argument of glmnet
  alpha = 1,
  
  # lasso with rf preselection
  top.n = 30,
  
  # lasso and svr
  nfolds = 10,
  enforce.treatment = FALSE,
  
  # random forest specific
  tune.ranger = FALSE,
  
  n.cv = 10,
  #prediction.sample = "random", #prediction measures based on which mice (random or timepoint)
  #prediction.group = "month12_ccl4", 
  #test.size.strata = 1,
  #test.size.random = 6, 
  seed.set = 23,
  #sample.seed = NULL,
  intercept = TRUE){
  
  model.obj <- protein.regression.function(dataset = dataset,
                                           geneprotein = geneprotein,
                                           #geneprotein_list = pairs.list[1:2]
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
                                           
                                           # which Model method (lasso, svr, rf)
                                           model.method = model.method,
                                           
                                           collinearity.reduce = collinearity.reduce,
                                           rho.cutoff = rho.cutoff,
                                           
                                           manual.weights = manual.weights,
                                           weights.vector = weights.vector,
                                           
                                           lasso.fam = lasso.fam, # family argument of glmnet
                                           alpha = alpha,
                                           type.measure = type.measure, # type.measure argument of glmnet
                                           
                                           prediction.method = "full", # loo, cv, full, random
                                           
                                           top.n = top.n,
                                           enforce.treatment = enforce.treatment,
                                           
                                           nfolds = nfolds,
                                           
                                           tune.ranger = tune.ranger,
                                           
                                           #prediction.sample = prediction.sample, #prediction measures based on which mice (random or timepoint)
                                           #prediction.group = prediction.group, 
                                           #test.size.strata = test.size.strata,
                                           #test.size.random = test.size.random, 
                                           seed.set = seed.set,
                                           #sample.seed = sample.seed,
                                           intercept = intercept)
  
  
  if(prediction == TRUE){
    
    # loo prediction
    prediction.obj <- protein.regression.function(dataset = dataset,
                                                  geneprotein = geneprotein,
                                                  #geneprotein_list = pairs.list[1:2]
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
                                                  
                                                  # which Model method (lasso, svr, rf)
                                                  model.method = model.method,
                                                  
                                                  prediction.method = "loo", # loo, cv, full, random
                                                  
                                                  top.n = top.n,
                                                
                                                  manual.weights = manual.weights,
                                                  weights.vector = weights.vector,
                                                  
                                                  lasso.fam = lasso.fam, # family argument of glmnet
                                                  alpha = alpha,
                                                  type.measure = type.measure, # type.measure argument of glmnet
                                                  
                                                  nfolds = nfolds,
                                                  
                                                  tune.ranger = tune.ranger,
                                                  
                                                  #prediction.sample = prediction.sample, #prediction measures based on which mice (random or timepoint)
                                                  #prediction.group = prediction.group, 
                                                  #test.size.strata = test.size.strata,
                                                  #test.size.random = test.size.random, 
                                                  seed.set = seed.set,
                                                  #sample.seed = sample.seed,
                                                  intercept = intercept)
    
    print(geneprotein)
    
    
    return(list(model.obj = model.obj,
                prediction.obj = prediction.obj))
    
  }else{
    
    return(list(model.obj = model.obj))
    
  }
}





