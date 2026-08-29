# GOODNESS- OF FIT MEASURES

# function for mean of the explained deviance for a group of genes
dev.mean.complete <- function(x, regression = "gam"){ # x output of one of the regression functions with full data
  
  devs <- sapply(x, function(model.list) model.list$model.obj$model$deviance.explained)
  
  return(mean(devs, na.rm = TRUE))
}


rsq.mean.complete <- function(x, method = "LASSO"){
  
  rsqs <- c()
  for(i in 1:length(x)){
    
    if(x[[i]]$model.obj$model$final.lasso$nobs <= x[[i]]$k && method == "LASSO"){
      rsqs[i] <- NA
    }else if(length(x[[i]]$model.obj$model$regression.object$y) <= x[[i]]$k && method == "LM"){
      rsqs[i] <- NA
    }else{
      rsqs[i] <- x[[i]]$model.obj$model$adj.R.squared
    }
    
  }
  
  rsq.mean <- mean(rsqs, na.rm = TRUE)
  rsq.rate <- sum(!is.na(rsqs)) / length(rsqs)
  
  return(list(r.squares = rsqs,
              r.squares.mean = rsq.mean,
              r.squares.rate = rsq.rate)
  )
}




aic.mean.complete <- function(x){ # x output of one of the regression functions with full data
  
  aics <- sapply(x, function(model.list) model.list$model.obj$model$AIC)
  
  return(mean(aics, na.rm = TRUE))
}


bic.mean.complete <- function(x){ # x output of one of the regression functions with full data
  
  bics <- sapply(x, function(model.list) model.list$model.obj$model$BIC)
  
  return(mean(bics, na.rm = TRUE))
}


adjR2.mean.complete <- function(x){
  R2s <- sapply(x, function(model.list) model.list$model.obj$model$adj.R.squared)
  
  R2s[which(R2s == -Inf)] <- NA
  
  return(mean(R2s, na.rm = TRUE))
}


rmse.mean.complete <- function(x){
  rmses <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse)
  
  return(mean(rmses, na.rm = TRUE))
}

rmse.scaled.mean.complete <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse.scaled)
  
  return(mean(rmses.scaled, na.rm = TRUE))
}

rmse.scaled.median.complete <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse.scaled)
  
  return(median(rmses.scaled, na.rm = TRUE))
}


mae.mean.complete <- function(x){
  maes <- sapply(x, function(model.list) model.list$prediction.obj$model$mae)
  
  return(mean(maes, na.rm = TRUE))
}

mae.scaled.mean.complete <- function(x){
  maes.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$mae.scaled)
  
  return(mean(maes.scaled, na.rm = TRUE))
}

mae.scaled.median.complete <- function(x){
  maes.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$mae.scaled)
  
  return(median(maes.scaled, na.rm = TRUE))
}


# Spearman correlation between y and predictions (as in Prohabar et al)
spearman.cor.y.prediction.complete <- function(x){
  predictions <- lapply(x, function(model.list) model.list$prediction.obj$model$prediction)
  ys <- lapply(x, function(model.list) model.list$prediction.obj$model$y)
  
  correlations <- mapply(predictions, ys, FUN = cor, method = "spearman", use = "complete.obs", SIMPLIFY = TRUE)
  
  return(correlations)
}

pearson.cor.y.prediction.complete <- function(x){
  predictions <- lapply(x, function(model.list) model.list$prediction.obj$model$prediction)
  ys <- lapply(x, function(model.list) model.list$prediction.obj$model$y)
  
  correlations <- mapply(predictions, ys, FUN = cor, method = "pearson", use = "complete.obs", SIMPLIFY = TRUE)
  
  return(correlations)
}


# FUNCTION TO RETURN ALL MEASURES AT ONCE --------------------------------------

gof.all.complete <- function(model,
                             prediction.method = "lasso" #lasso, svr, rf
                             ){ 
  
  
  res.aic <- aic.mean.complete(model)
  res.bic <- bic.mean.complete(model)
  res.rsq <- adjR2.mean.complete(model)
  
  #res.mse <- rmse.mean(pred)
  #res.mae <- mae.mean(pred)
  
  #res.mse.scaled <- rmse.scaled.mean.complete(model)
  #res.mae.scaled <- mae.scaled.mean.complete(model)
  
  res.mse <- rmse.mean.complete(model)
  res.mae <- mae.mean.complete(model)
  
  
  res.pearson.cor <- mean(pearson.cor.y.prediction.complete(model))
  res.spearman.cor <- mean(spearman.cor.y.prediction.complete(model))

  
  if(prediction.method == "lasso"){
    
    res.dev <- dev.mean.complete(model)
    
    res.all <- data.frame(mean.deviance.explained = res.dev,
                          mean.adj.R.squared = res.rsq,
                          mean.aic = res.aic,
                          mean.bic = res.bic,
                          mean.rmse = res.mse,
                          mean.mae = res.mae,
                          mean.pearson.cor = res.pearson.cor,
                          mean.spearman.cor = res.spearman.cor
    )
  }else if(prediction.method == "baseline"){
    
    
  }else{
    
    res.all <- data.frame(mean.adj.R.squared = res.rsq,
                          mean.aic = res.aic,
                          mean.bic = res.bic,
                          mean.rmse = res.mse.scaled,
                          mean.mae = res.mae.scaled,
                          mean.pearson.cor = res.pearson.cor,
                          mean.spearman.cor = res.spearman.cor
    )
    
  }
  
  
  return(res.all)
  
}



# COMPARISON TO BASELINE -------------------------------------------------------

rmse.scaled.mean.baseline <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$baseline.model$rmse.scaled)
  
  return(mean(rmses.scaled, na.rm = TRUE))
}

rmse.scaled.mean.model <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse.scaled)
  
  return(mean(rmses.scaled, na.rm = TRUE))
}


all.rmses.baseline <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$baseline.model$rmse.scaled)
  
  return(rmses.scaled)
}


all.rmses.model <- function(x){
  rmses.scaled <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse.scaled)
  
  return(rmses.scaled)
}


diffs.baseline <- function(x){
  diffs <- sapply(x, function(model.list) abs(as.numeric(model.list$prediction.obj$baseline.model$prediction) - model.list$prediction.obj$baseline.model$y))
  diffs.scaled <- sapply(x, function(model.list) abs(as.numeric(model.list$prediction.obj$baseline.model$prediction) - model.list$prediction.obj$baseline.model$y) / sd(model.list$prediction.obj$baseline.model$y))
  
  return(diffs.scaled)
}

diffs.model <- function(x){
  diffs <- sapply(x, function(model.list) abs(as.numeric(model.list$prediction.obj$model$prediction) - model.list$prediction.obj$model$y))
  diffs.scaled <- sapply(x, function(model.list) abs(as.numeric(model.list$prediction.obj$model$prediction) - model.list$prediction.obj$model$y) / sd(model.list$prediction.obj$baseline.model$y))
  
  return(diffs.scaled)
}


function.comparison.to.baseline <- function(x,
                                            thresh = 0.7){
  
  proteins <- sapply(x, function(model.list) model.list$prediction.obj$model$protein)
  
  # RMSEs
  rmses.baseline <- sapply(x, function(model.list) model.list$prediction.obj$baseline.model$rmse.scaled)
  
  rmses.model <- sapply(x, function(model.list) model.list$prediction.obj$model$rmse.scaled)
  
  # Correlations
  cors.baseline <- sapply(x, function(model.list) model.list$prediction.obj$baseline.model$correlation.pearson)
  
  cors.model <- sapply(x, function(model.list) model.list$prediction.obj$model$correlation.pearson)
  
  
  DTplot <- rbind(data.frame(protein = proteins, rmse = rmses.baseline, correlation = cors.baseline, type = "baseline"),
                  data.frame(protein = proteins, rmse = rmses.model, correlation = cors.model, type = "model"))
  
  DTcomp <- data.frame(protein = proteins, 
                       rmse.diff = rmses.model - rmses.baseline,
                       correlation.diff = cors.model - cors.baseline)
  
  #individual-wise difference (absolute)
  baseline.diffs <- diffs.baseline(x)
  model.diffs <- diffs.model(x)
  
  #TODO ADD PROTEINS AND MICE INFO (latter is missing in output)
  #proteins.stretched <- 
  #mice.info <-
  
  DTindividual <- rbind(data.frame(diffs = unlist(baseline.diffs), type = "baseline"),
                        data.frame(diffs = unlist(model.diffs), type = "model"))
  
  
  DTindividual.comp <- data.frame(unlist(baseline.diffs) - unlist(model.diffs))
  colnames(DTindividual.comp) <- c("comparison.to.baseline")
  
  # actual plots ---------------------------------------------------------------
  
  # colorblind-friendly palette
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # RMSE boxplots baseline and model
  ggRMSE <- ggplot(data = DTplot, aes(y = rmse, x = type, fill = type)) +
    geom_boxplot() +
    #ylim(c(-2,2)) + 
    ylab("RMSE") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(name = "Model type",
                      values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    #coord_flip() +
    #scale_x_continuous(breaks = NULL) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(alpha = "none")
  
  # RMSE boxplot model minus baseline
  ggRMSE_diff <- ggplot(data = DTcomp, aes(y = rmse.diff, x = "")) +
    geom_boxplot() +
    ylab("RMSE difference between mRNA Design Model and Baseline model") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    coord_flip() 
  
  # Correlation boxplots baseline and model
  ggpearson <- ggplot(data = DTplot, aes(y = correlation, x = type, fill = type)) +
    geom_boxplot() +
    ylab("Pearson correlation") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(name = "Model type",
                      values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    #coord_flip() +
    #scale_x_continuous(breaks = NULL) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(alpha = "none")
  
  # Correlation boxplot model minus baseline
  ggpearson_diff <- ggplot(data = DTcomp, aes(y = correlation.diff, x = "")) +
    geom_boxplot() +
    ylab("Correlation difference between mRNA Design Model and Baseline model") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    coord_flip() 
  
  
  # INDIVIDUAL-WISE PLOTS (abs(yhat-y))
  ggindividual <- ggplot(data = DTindividual, aes(y = diffs, x = type, fill = type)) +
    geom_boxplot() +
    ylab(expression(abs(hat(y) - y))) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(name = "Model type",
                      values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    #coord_flip() +
    #scale_x_continuous(breaks = NULL) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(alpha = "none")
  
  # INDIVIDUAL-WISE PLOTS (abs(yhat-y))
  ggindividual_diff <- ggplot(data = DTindividual.comp, aes(y = comparison.to.baseline, x = "")) +
    geom_boxplot() +
    ylab(paste0("|yhat - y|: mRNA Design minus baseline")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    #scale_fill_brewer(palette="Set2") +
    scale_fill_manual(name = "Design type",
                      values = cbp1) +
    theme_bw() + 
    #labs(fill = "Design type")
    coord_flip() +
    guides(alpha = "none")
  
  # TABLES ---------------------------------------------------------------------
  
  number.of.proteins.corr.0.7_baseline <- sum(cors.baseline > thresh)
  number.of.proteins.corr.0.7_model <- sum(cors.model > thresh)
  
  tab.number.of.proteins.corr.0.7 <- data.frame(number.of.proteins.corr.0.7_baseline, 
                                                number.of.proteins.corr.0.7_model)
  colnames(tab.number.of.proteins.corr.0.7) <- c("Proteins with corr. > 0.7 (baseline)", 
                                                 "Proteins with corr. > 0.7 (mRNA model)")
  
  
  tab.model.vs.baseline <- data.frame(rbind(c(sum(DTcomp$rmse.diff < 0)/nrow(DTcomp), sum(DTcomp$rmse.diff > 0)/nrow(DTcomp)),
                                      c(sum(DTcomp$correlation.diff > 0)/nrow(DTcomp), sum(DTcomp$correlation.diff < 0)/nrow(DTcomp)),
                                      c(sum(DTindividual.comp < 0)/nrow(DTindividual.comp), sum(DTindividual.comp > 0)/nrow(DTindividual.comp))))
  colnames(tab.model.vs.baseline) <- c("mRNA model better", "mRNA model worse")
  rownames(tab.model.vs.baseline) <- c("RMSE", "Correlation", "|yhat-y|")
  
  
  return(list(ggRMSE = ggRMSE,
              ggRMSE_diff = ggRMSE_diff,
              ggpearson = ggpearson,
              ggpearson_diff = ggpearson_diff,
              ggindividual = ggindividual,
              ggindividual_diff = ggindividual_diff,
              tab.number.of.proteins.corr.0.7 = tab.number.of.proteins.corr.0.7,
              tab.model.vs.baseline = tab.model.vs.baseline))
}



# PREDICTION AGAINST TRUE VALUES PLOTS -----------------------------------------

plot.prediction <- function(obj,
                            type = "model",
                            i.pair,
                            title = "Prediction vs true values"){
  
  if(type == "model"){
    
    DTplot <- data.frame(prediction = obj[[i.pair]]$prediction.obj$model$prediction,
                         y = obj[[i.pair]]$prediction.obj$model$y)
    
    correlation <-  obj[[i.pair]]$prediction.obj$model$correlation.pearson
    rmse <- obj[[i.pair]]$prediction.obj$model$rmse.scaled
    
  }else if(type == "baseline"){
    
    DTplot <- data.frame(prediction = obj[[i.pair]]$prediction.obj$baseline.model$prediction,
                         y = obj[[i.pair]]$prediction.obj$baseline.model$y)
    
    correlation <-  obj[[i.pair]]$prediction.obj$baseline.model$correlation.pearson
    rmse <- obj[[i.pair]]$prediction.obj$baseline.model$rmse.scaled
    
  }else if(type == "post.lasso"){
    
    DTplot <- data.frame(prediction = obj[[i.pair]]$prediction.obj$model$prediction_relaxed,
                         y = obj[[i.pair]]$prediction.obj$model$y)
    
    correlation <-  obj[[i.pair]]$prediction.obj$model$correlation.pearson_relaxed
    rmse <- obj[[i.pair]]$prediction.obj$model$rmse.scaled_relaxed
    
  }
  
  #for alt data
  
  # protein <- obj[[i.pair]]$model.obj$model$protein
  # DTplot$treatment <- DTccl4[ProteinID  == protein, Treatment][match(DTplot$y, DTccl4[ProteinID  == protein, ProteinIntensity])]
  # 
  # 
  # plot.pred <- ggplot(data = DTplot, aes(x = y, y = prediction, colour = as.factor(treatment))) + xlab(paste0("Observed Protein intensity of Protein ", protein)) + ylab("Predicted Protein intensity")+  
  #   geom_point()  + ggtitle(paste0(title, " for " ,protein, ", ","\ncorrelation =", round(correlation, 3), ", RMSE = ", round(rmse, 3))) +
  #   theme_bw() + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
  #                      axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  #                      plot.title = element_text(size = 25,  hjust = 0.5)) +
  #   labs(colour = "Treatment") + 
  #   ylim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE)) + xlim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE))
  # 
  # return(plot.pred)
  
  #for new data
  # protein <- obj[[i.pair]]$model.obj$model$protein
  # DTplot$treatment <- DT[ProteinID  == protein, Treatment][match(DTplot$y, DT[ProteinID  == protein, ProteinIntensity])]
  # 
  # 
  # plot.pred <- ggplot(data = DTplot, aes(x = y, y = prediction, colour = as.factor(treatment))) + xlab(paste0("Observed Protein intensity of Protein ", protein)) + ylab("Predicted Protein intensity")+  
  #   geom_point()  + ggtitle(paste0(title, " for " ,protein, ", ","\ncorrelation =", round(correlation, 3), ", RMSE = ", round(rmse, 3))) +
  #   theme_bw() + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
  #                      axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  #                      plot.title = element_text(size = 25,  hjust = 0.5)) +
  #   labs(colour = "Treatment") + 
  #   ylim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE)) + xlim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE))
  # 
  # return(plot.pred)
  
  #ohne treatment
  
  protein <- obj[[i.pair]]$model.obj$model$protein
  
  # Ici DTplot contient déjà y et prediction
  # On ne rajoute plus de colonne "treatment"
  
  plot.pred <- ggplot(data = DTplot, aes(x = y, y = prediction)) +
    xlab(paste0("Observed Protein intensity of Protein ", protein)) +
    ylab("Predicted Protein intensity") +
    geom_point() +
    ggtitle(paste0(title, " for ", protein, ", ",
                   "\ncorrelation = ", round(correlation, 3),
                   ", RMSE = ", round(rmse, 3))) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
                               plot.title = element_text(size = 25,  hjust = 0.5)) +
    ylim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE)) +
    xlim(min(DTplot[,1:2], na.rm = TRUE), max(DTplot[,1:2], na.rm = TRUE))
  
  return(plot.pred)
  
  
}



plot.y.vs.x <- function(protein,
                        covariate,
                        covariate.type = "protein",
                        title  = ""){
  
  #y <- DTccl4[ProteinID == protein, ProteinIntensity]
  #y <- DT_filtered[ProteinID == protein, ProteinIntensity]
  #y <- DTccl4_filtered[ProteinID == protein, ProteinIntensity]
  y <- DTccl4_DT_LCPM[ProteinID == protein, ProteinIntensity]
  
  if(covariate.type == "protein"){
    
    #x <- DTccl4[ProteinID == covariate, ProteinIntensity]
    
    #x <- DT_filtered[ProteinID == covariate, ProteinIntensity]
    #x <- DTccl4_filtered[ProteinID == covariate, ProteinIntensity]
    x <- DTccl4_DT_LCPM[ProteinID == covariate, ProteinIntensity]
    
  }else if(covariate.type == "rna"){
    
    #x <- DTccl4[GeneSyn == covariate, GeneCount]
    
    #x <- DT_filtered[GeneSyn == covariate, GeneCount]
    #x <- DTccl4_filtered[GeneSyn == covariate, GeneCount]
    x <- DTccl4_DT_LCPM[GeneSyn == covariate, GeneCount]
  }
  
  #treatment <- DTccl4[ProteinID == protein, Treatment]
  
  # treatment <- DT_filtered[ProteinID == protein, Treatment]
  # 
  # DTplot <- data.frame(y = y,
  #                      x = x, 
  #                      treatment = treatment)[complete.cases(data.frame(y=y, x=x)),]
  # 
  # correlation <- cor(DTplot$x, DTplot$y, method = "pearson")
  # 
  # 
  # plot.obj <- ggplot(data = DTplot, aes(x = x, y = y, colour = as.factor(treatment))) + xlab(paste0("Covariate expression: ", covariate)) + ylab(paste0("Target protein level: ", protein))+  
  #   geom_point()  + ggtitle(paste0(title, " ", covariate, " vs target protein " ,protein, ", ","\ncorrelation =", round(correlation, 3))) +
  #   theme_bw() + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
  #                      axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #                      axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  #                      plot.title = element_text(size = 25,  hjust = 0.5)) +
  #   labs(colour = "Treatment") + 
  #   ylim(min(DTplot$y, na.rm = TRUE), max(DTplot$y, na.rm = TRUE)) + xlim(min(DTplot$x, na.rm = TRUE), max(DTplot$x, na.rm = TRUE))
  # 
  # return(plot.obj)
  # Construire le data.frame
  DTplot <- data.frame(y = y, x = x)[complete.cases(data.frame(y=y, x=x)),]
  
  # Corrélation
  correlation <- cor(DTplot$x, DTplot$y, method = "pearson")
  
  # Plot
  plot.obj <- ggplot(data = DTplot, aes(x = x, y = y)) +
    xlab(paste0("Covariate expression: ", covariate)) +
    ylab(paste0("Target protein level: ", protein)) +
    geom_point() +
    ggtitle(paste0(title, " ", covariate, " vs target protein " ,protein,
                   ", \ncorrelation =", round(correlation, 3))) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                               axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
                               plot.title = element_text(size = 25,  hjust = 0.5))
  
  return(plot.obj)
  
}



plot.prediction.compare <- function(i,
                                    obj.rna,
                                    obj.protein){
  return(list(plot.prediction(obj.rna, type = "baseline", i, title = "Baseline model"),
              plot.prediction(obj.rna, type = "post.lasso", i, title = "mRNA model"),
              plot.prediction(obj.protein, type = "post.lasso", i, title = "Protein model")))
}

####### die Plots zu den Modellbeispielen ######################################
plot.example.full <- function(i,
                              obj.rna,
                              obj.protein){
  
  protein <- obj.rna[[i]]$model.obj$model$protein
  rna <- obj.rna[[i]]$model.obj$model$gene
  
  p_mainrna <- plot.y.vs.x(protein, rna, covariate.type = "rna", title = "Main mRNA")
  
  p_base <- plot.prediction(obj.rna, type = "baseline", i, title = "Baseline model")
  p_rna <- plot.prediction(obj.rna, type = "post.lasso", i, title = "mRNA model")
  p_protein <- plot.prediction(obj.protein, type = "post.lasso", i, title = "Protein model")
  
  models <- c("baseline", "rna", "protein")
  results <- c(obj.rna[[i]]$prediction.obj$baseline.model$correlation.pearson,
               obj.rna[[i]]$prediction.obj$model$correlation.pearson_relaxed,
               obj.protein[[i]]$prediction.obj$model$correlation.pearson_relaxed)
  
  i.best <- which(results == max(results))
  
  if(models[i.best] == "baseline"){
    
    covariates <- rna
    
    p_covar1 <- plot.y.vs.x(protein, covariates[1], covariate.type = "rna", title = "Main mRNA")
    
    
  }else if(models[i.best] == "rna"){
    
    covariates <- str_split(names(obj.rna[[i]]$model.obj$model$variable.selection[-1]), "_", simplify = TRUE)[,1]
    
    p_covar1 <- plot.y.vs.x(protein, covariates[1], covariate.type = "rna", title = "Main mRNA")
    
    if(length(covariates) > 1){
      
      for(i in 1:(length(covariates)-1)){
        assign(paste0("p_covar", i+1), plot.y.vs.x(protein, covariates[-1][i], covariate.type = "rna", title = "Covariate"))
      }
      
    }
    
    
  }else if(models[i.best] == "protein"){
    
    covariates <- str_split(names(obj.protein[[i]]$model.obj$model$variable.selection[-1]), "_", simplify = TRUE)[,1]
    
    p_covar1 <- plot.y.vs.x(protein, covariates[1], covariate.type = "rna", title = "Main mRNA")
    
    if(length(covariates) > 1){
      
       for(i in 1:(length(covariates)-1)){
        assign(paste0("p_covar", i+1), plot.y.vs.x(protein, covariates[-1][i], covariate.type = "protein", title = "Covariate"))
      }
      
    }
   
    
  }
  
  
  make_list <- function(n) {
    base_part <- list(
      mainrna = p_mainrna,
      Baseline = p_base,
      mRNA = p_rna,
      Protein = p_protein
    )
    
    numbered_names <- paste0("p_covar", 1:n)
    numbered_part <- mget(numbered_names, inherits = TRUE)
    
    return(c(base_part, numbered_part))
  }
  
  list.return <- make_list(length(covariates))
  
  
  return(list.return)
  
  
}



