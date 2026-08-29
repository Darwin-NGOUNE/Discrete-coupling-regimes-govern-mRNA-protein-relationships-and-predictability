# DESIGN FUNCTION --------------------------------------------------------------

library("data.table")
library("stringr" )

design.function <- function(dataset,
                            geneprotein,
                            #groups, #vector of group membership
                            weighted = FALSE,
                            weight.method = NULL, # either "gca" similarity or "dipa" distance
                            weight.matrix = NULL, # associated matrix of weights (from co-expression analysis or from dipa analysis)
                            weight.matrix.p = NULL, # needed for GCA weights and combined design
                            grouping = "all", # "all" genes as covariates, or "dipa", "coexpression", "manual"
                            scaled = TRUE,
                            RNA.log.scale = TRUE, # log-transformation of the mRNA data
                            manual.covar = NULL, 
                            design.m = "genes", 
                            na.process = "max.n", # complete cases ("cc"), maximum n ("max.n")
                            include.treatment = FALSE,
                            treatment.info = "full", #full, treatment, treatment.duration, separate
                            duration.scale = "factor", #duration as a factor or numeric covariate
                            treatment.penalty = FALSE, # penalty for the treatment covariate
                            PE.RNA.penalty = FALSE, # penalty for the protein encoding RNA
                            baseline.interactions = FALSE){
  
  # extract gene name and protein ID from geneprotein
  gene <- str_split(geneprotein, "_")[[1]][1]
  protein <- str_split(geneprotein, "_")[[1]][2]
  
  # all matching gene protein pairs in the dataset
  gene.protein.list <- data.table(str_split(unique(dataset$GeneProtein), "_", simplify = TRUE))
  colnames(gene.protein.list) <- c("GeneSyn", "ProteinID")
  
  gene.protein.list$GeneSyn <- str_replace(gene.protein.list$GeneSyn, "-",  ".")
  
  if(grouping == "dipa"){
    # sub-data set including the investigated gene and its cluster members
    genegroup <- unique(dataset[GeneProtein == geneprotein, ClusterDiPa])
    sub_data <- dataset[ClusterDiPa == genegroup,]
  }else if(grouping == "coexpression"){
    genegroup <- unique(dataset[GeneProtein == geneprotein, ClusterGCA])
    sub_data <- dataset[ClusterGCA == genegroup,]
  }else if(grouping == "km"){
    genegroup <- unique(dataset[GeneProtein == geneprotein, ClusterKmProtein])
    sub_data <- dataset[ClusterKmProtein == genegroup,]
  }else if(grouping == "all"){ # all genes as covariates
    sub_data <- dataset
  }else if(grouping == "manual"){
    
    if(design.m == "proteins"){
      sub_data <- dataset[ProteinID %in% c(protein, manual.covar),]
    }else if(design.m == "genes"){
      sub_data <- dataset[GeneSyn %in% c(gene, manual.covar),]
    }else if(design.m == "all"){
      sub_data <- dataset[GeneSyn %in% c(gene, manual.covar) | ProteinID %in% c(protein, manual.covar),]
    }
    
  }
  
  # for weight method DiPa, the weight matrix for proteins is the same as for RNA,
  # for weight method GCA, the weight matrix has to be supplied by the user
  if(weighted == TRUE){
    
    #if(weight.method == "dipa"){
    #weight.matrix.p <- weight.matrix
    #}
    
    #change "-" in Gene names to "." in order to make formulas work
    if(design.m == "genes"){
      
      colnames(weight.matrix) <- str_replace(colnames(weight.matrix), "-",  ".")
      rownames(weight.matrix) <- str_replace(colnames(weight.matrix), "-",  ".")
      
    }else if(design.m == "proteins"){
      
      #colnames(weight.matrix.p) <- str_replace(colnames(weight.matrix.p), "-",  ".")
      #rownames(weight.matrix.p) <- str_replace(colnames(weight.matrix.p), "-",  ".")
      
      if(weight.method == "dipa"){
        #weight.matrix.p <- weight.matrix
        colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        
      }else if(weight.method == "gca"){ #TODO TEMPORARY SOLUTION TO MAKE IT WORK WITH GENE CLUSTERING
        
        colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        
        #}else if(weight.method == "km"){
        
        #colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$ProteinID)]
        #rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$ProteinID)]
        
      }
      
    }else if(design.m == "all"){
      
      colnames(weight.matrix) <- str_replace(colnames(weight.matrix), "-",  ".")
      rownames(weight.matrix) <- str_replace(colnames(weight.matrix), "-",  ".")
      
      colnames(weight.matrix.p) <- str_replace(colnames(weight.matrix.p), "-",  ".")
      rownames(weight.matrix.p) <- str_replace(colnames(weight.matrix.p), "-",  ".")
      
      
      if(weight.method == "dipa"){
        #weight.matrix.p <- weight.matrix
        colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        
      }else if(weight.method == "gca"){ #TODO TEMPORARY SOLUTION TO MAKE IT WORK WITH GENE CLUSTERING
        
        colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$GeneSyn)]
        
        #}else if(weight.method == "km"){
        
        #colnames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$ProteinID)]
        #rownames(weight.matrix.p) <- gene.protein.list$ProteinID[match(colnames(weight.matrix.p), gene.protein.list$ProteinID)]
        
      }
      
    }
    
    
  }
  
  
  #change "-" in Gene names to "." in order to make formulas work
  gene <- str_replace(gene, "-",  ".")
  sub_data$GeneSyn <- str_replace(sub_data$GeneSyn, "-",  ".")
  
  # genes and proteins of the sub-dataset
  genes <- unique(sub_data$GeneSyn)
  proteins <- unique(sub_data$ProteinID)
  
  
  # wide format for both gene count data and protein data
  if(RNA.log.scale == TRUE){
    glong <- melt(sub_data[, .(MiceInfo, Treatment, GeneProtein, log2(GeneCount + 1))], id.vars = c("MiceInfo", "Treatment", "GeneProtein"))
  }else{
    glong <- melt(sub_data[, .(MiceInfo, Treatment, GeneProtein, GeneCount)], id.vars = c("MiceInfo", "Treatment", "GeneProtein"))
    
  }
  
  gwide <- dcast(glong, MiceInfo + Treatment  ~ GeneProtein)
  
  plong <- melt(sub_data[, .(MiceInfo, Treatment, GeneProtein, ProteinIntensity)], id.vars = c("MiceInfo", "Treatment", "GeneProtein"))
  pwide <- dcast(plong, MiceInfo + Treatment ~ GeneProtein)
  
  # position of the investigated protein within the protein sub-data
  ppos <- which(geneprotein == colnames(pwide))
  # position of the corresponding gene within the gene sub-data
  gpos <- which(geneprotein == colnames(gwide))
  
  
  # adjust colnames in order to comprehend, whether the data is gene counts or protein intensities
  # (keep the complete name including gene and protein in order to comprehend the pair affiliation)
  #colnames(gwide) <- paste0(colnames(gwide),"_g")
  #colnames(pwide) <- paste0(colnames(pwide),"_p")
  
  
  # only actual genes as columns
  colsg <- !(colnames(gwide) %in% c("MiceInfo","Treatment"))
  gdesign <- subset(gwide,,colsg)
  
  colsp <- !(colnames(pwide) %in% c("MiceInfo","Treatment"))
  pdesign <- subset(pwide,,colsp)
  
  # adjust colnames in order to comprehend, whether the data is gene counts or protein intensities
  # (keep the complete name including gene and protein in order to comprehend the pair affiliation)
  colnames(gdesign) <- paste0(str_split(colnames(gdesign), "_", simplify = TRUE)[,1], "_G")
  colnames(pdesign) <- paste0(str_split(colnames(pdesign), "_", simplify = TRUE)[,2], "_P")
  
  #change "-" in Gene names to "." in order to make formulas work
  colnames(gdesign) <- str_replace(colnames(gdesign), "-",  ".")
  
  #remove potential duplicated columns in the protein matrix
  #i.duplicated <- which(duplicated(colnames(pdesign)))
  #pdesign <- pdesign[,which(!duplicated(colnames(pdesign))), with = FALSE]
  
  
  i.maingene <- which(colnames(gdesign) == paste0(gene, "_G"))
  
  
  
  ## DESIGN MATRIX ---------------------------------------------------------------
  
  # (a) only gene counts as covariates
  # (b) gene counts and protein intensities (except the one to be estimated) as covariates
  
  # scale() the covariates for LASSO regression
  if(design.m == "genes"){ # RNA as covariates
    
    if(scaled == TRUE){
      design.raw <- scale(gdesign)
    }else{
      design.raw <- gdesign
    }
    
  }else if(design.m == "proteins"){ # Proteins as covariates
    
    #remove investigated protein from protein data:
    pdesign.raw <- copy(pdesign)[,paste0(protein, "_P"):=NULL]
    
    if(scaled == TRUE){
      design.raw <- scale(cbind(gdesign[, ..i.maingene],
                                pdesign.raw))
    }else{
      design.raw <- cbind(gdesign[, ..i.maingene],
                          pdesign.raw)
    }
    
  }else if(design.m == "all"){ # both RNA and proteins as covariates
    
    #remove investigated protein from protein data:
    pdesign.raw <- copy(pdesign)[,paste0(protein,"_P"):=NULL]
    
    if(scaled == TRUE){
      design.raw <- scale(cbind(gdesign, 
                                pdesign.raw))
    }else{
      design.raw <- cbind(gdesign, pdesign.raw)
    }
    
  }else if(design.m == "experiment"){
    design.raw <- NULL #TODO FUNKTIONIERT SO NOCH NICHT
  }
  
  # DEALING WITH NAs
  # option 1: remove rows (observations) containing NAs
  # option 2: remove columns (proteins) containing NAs, resulting in larger n (only relevant for Design approach (b))
  
  y_temp <- pwide[[ppos]] # full y vector (might contain NAs)
  mice_temp <- pwide$MiceInfo
  treatment_temp <- factor(pwide$Treatment,
                           levels = c("ccl4",
                                      "oil",
                                      "control", 
                                      "BDL",  
                                      "BDL_ASBTi"))
  
  tempdat <- data.frame(mice_temp, treatment_temp, y_temp, design.raw)
  colnames(tempdat)[-c(1:3)] <- colnames(design.raw)
  
  
  # if(na.process == "cc"){
  #   # omit NAs from the target variable
  #   i_rel <- which(complete.cases(tempdat))
  #   
  #   tempdat_final <- tempdat[complete.cases(tempdat), ]
  #   
  # }else if(na.process == "max.n"){
  #   # omit NAs from the target variable
  #   i_rel <- which(complete.cases(y_temp))
  #   
  #   # remove columns of tempdat_rel that contain NAs
  #   tempdat_rel <- tempdat[i_rel,]
  #   tempdat_final <- tempdat_rel[, colSums(is.na(tempdat_rel))==0]
  # }
  
  if (na.process == "cc") {
    # supprime toutes les lignes avec NA
    i_rel <- which(complete.cases(tempdat))
    tempdat_final <- tempdat[complete.cases(tempdat), ]
    
  } else if (na.process == "max.n") {
    # supprime les lignes avec NA dans y_temp et les colonnes avec NA
    i_rel <- which(complete.cases(y_temp))
    tempdat_rel <- tempdat[i_rel, ]
    tempdat_final <- tempdat_rel[, colSums(is.na(tempdat_rel)) == 0]
    
  } else if (na.process == "none") {
    # ne supprime rien, conserve les NA
    i_rel <- seq_len(nrow(tempdat))   # toutes les lignes
    tempdat_final <- tempdat          # toutes les colonnes inchangées
  }
  
  # associated vectors of y, mice information and experiment settings
  y <- y_temp[i_rel]
  n <- length(y)
  mice <- mice_temp[i_rel]
  
  
  experiments <- treatment_temp[i_rel]
  
  treatment <- treatment_temp[i_rel]
  
  # TREATMENT INFORMATION ------------------------------------------------------
  if(include.treatment == TRUE){
    
    if(treatment.info == "full"){
      
      t.info <- experiments
      design.experiments <- model.matrix(~ .-1, data.frame(Treatment.info = t.info))[,-1]
      
    }else if(treatment.info == "treatment"){
      
      t.info <- treatment
      design.experiments <- as.matrix(model.matrix(~ .-1, data.frame(Treatment.info = t.info))[,-1])
      colnames(design.experiments) <- "Treatmentccl4"
      
    }else if(treatment.info == "treatment.duration"){
      
      if(duration.scale == "factor"){
        treatment.duration <- factor(str_split(treatment_temp[i_rel], "_", simplify = TRUE)[,1],
                                     levels = c("month0", "month2", "month6", "month12"))
        
        t.info <- treatment.duration
        design.experiments <- model.matrix(~ .-1, data.frame(Treatment.info = t.info))[,-1]
        
      }else if(duration.scale == "numeric"){
        treatment.duration <- as.numeric(str_split(str_split(treatment_temp[i_rel], "_", simplify = TRUE)[,1], "h", simplify = TRUE)[,2]) 
        
        t.info <- treatment.duration
        design.experiments <- model.matrix(~ .-1, data.frame(Treatment.info = t.info))
      }
      
    }else if(treatment.info == "separate"){
      
      if(duration.scale == "factor"){
        treatment.duration <- factor(str_split(treatment_temp[i_rel], "_", simplify = TRUE)[,1],
                                     levels = c("month0", "month2", "month6", "month12"))
        
        
      }else if(duration.scale == "numeric"){
        treatment.duration <- as.numeric(str_split(str_split(treatment_temp[i_rel], "_", simplify = TRUE)[,1], "h", simplify = TRUE)[,2]) 
        
      }
      
      t.info <- data.frame(treatment, treatment.duration)
      design.experiments <- model.matrix(~ .-1, data.frame(Treatment.info = t.info))[,-1]
      
      
    }else( print("Treatment.info argument is missing or incorrect"))
    
    #TODO BASELINE INTERACTIONS UPDATEN (SIND AKTUELL NUR FÜR treatment.info == "full")
    if(baseline.interactions == TRUE){
      form.int <- as.formula(paste0("~", gene, "_G", ":", "treatment_temp"))
      design.gene.treatment.interaction <- model.matrix(form.int, tempdat_final)[,-c(1,2)]
      
      n.design.treatment <- ncol(design.experiments) + ncol(design.gene.treatment.interaction)
      
      Design <-  cbind(design.experiments, 
                       design.gene.treatment.interaction,
                       tempdat_final[, -c(1:3)])
    }else{
      Design <-  cbind(design.experiments, 
                       tempdat_final[, -c(1:3)])
      
      n.design.treatment <- ncol(design.experiments)
    }
  }else{ # no treatment covariates
    if(baseline.interactions == TRUE){
      
      #dummy variables for experiment factor variables
      design.experiments <- model.matrix(~ .-1, data.frame(Treatment = factor(experiments,
                                                                              levels = c("month0_oil", 
                                                                                         "month2_oil",  
                                                                                         "month12_oil", 
                                                                                         "month2_ccl4",
                                                                                         "month6_ccl4",
                                                                                         "month12_ccl4"))))[,-1]
      
      form.int <- as.formula(paste0("~", gene, "_G", ":", "treatment_temp"))
      design.gene.treatment.interaction <- model.matrix(form.int, tempdat_final)[,-c(1,2)]
      
      n.design.treatment <- ncol(design.experiments) + ncol(design.gene.treatment.interaction)
      
      Design <-  cbind(design.gene.treatment.interaction,
                       tempdat_final[, -c(1:3)])
    }else{
      Design <- tempdat_final[, -c(1:3)] 
    }
  }
  
  
  ## WEIGHTS FOR THE PENALTY ---------------------------------------------------
  
  
  # position of the "main gene" within the design matrix
  igene <- which(paste0(gene, "_G") == colnames(Design))
  
  if(weighted == TRUE){ 
    
    #if(weights == "co-expression"){ #input weight.matrix as GCA weight matrix or DiPa-distance
    
    
    if(design.m == "genes"){
      
      # positon of the "main gene" within the distance matrix
      itemp.g <- which(paste0(gene) == colnames(weight.matrix))
      
      gnames <- colnames(Design)[which(str_detect(colnames(Design), "_G"))]
      # match gene names from design matrix with weight matrix names
      match.g <- match(gnames, paste0(colnames(weight.matrix), "_G"))
      
      wtemp.g <- weight.matrix[, itemp.g]
      
      w.g <- wtemp.g[match.g]
      names(w.g) <- gnames
      
      # actual vector of weights
      if(include.treatment == TRUE){
        
        if(treatment.penalty == FALSE){ # should the treatment information be included with or without a penalty?
          
          w <- c(rep(0, sum(!(colnames(Design) %in% names(w.g)))), # weights for the experiment settings (=0)
                 w.g) # weights for the RNA covariates
          
        }else{
          
          w <- c(rep(1, sum(!(colnames(Design) %in% names(w.g)))), # weights for the experiment settings (=0)
                 w.g) # weights for the RNA covariates
          
        }
        
        
      }else{
        w <- w.g
      }
      
      names(w) <- colnames(Design)
      
      # should already be the case but make sure that the main gene is unweighted
      w[igene] <- 0
      
    }else if(design.m == "proteins"){
      
      
      # positon of the "main gene" within the distance matrix
      #itemp.g <- which(paste0(gene) == colnames(weight.matrix.p))
      itemp.p <- which(protein == colnames(weight.matrix.p))
      wtemp.p <- weight.matrix.p[, itemp.p]
      
      pnames <- colnames(Design)[which(str_detect(colnames(Design), "_P"))]
      
      # match protein names from design matrix with weight matrix names 
      match.p <- match(pnames, paste0(colnames(weight.matrix.p), "_P"))
      
      w.p <- wtemp.p[match.p]
      #print(match.p)
      names(w.p) <- pnames
      
      # actual vector of weights
      if(include.treatment == TRUE){
        
        if(treatment.penalty == FALSE){ # should the treatment information be included with or without a penalty?
          
          w <- c(rep(0, sum(!(colnames(Design) %in% names(w.p)))), # weights for the experiment settings (=0)
                 w.p) # weights for the protein covariates
          
        }else{
          
          w <- c(rep(1, sum(!(colnames(Design) %in% names(w.p)))), # weights for the experiment settings (=1)
                 w.p) # weights for the protein covariates
          
        }
        
      }else{
        w <- c(0, w.p)
      }
      
      names(w) <- colnames(Design)
      
    }else if(design.m == "all"){
      
      # positon of the "main gene" within the distance matrix
      itemp.g <- which(paste0(gene) == colnames(weight.matrix))
      itemp.p <- which(paste0(protein) == colnames(weight.matrix.p))
      
      gnames <- colnames(Design)[which(str_detect(colnames(Design), "_G"))]
      
      pnames <- colnames(Design)[which(str_detect(colnames(Design), "_P"))]
      
      # match gene names from design matrix with weight matrix names
      match.g <- match(gnames, paste0(colnames(weight.matrix), "_G"))
      
      wtemp.g <- weight.matrix[, itemp.g]
      wtemp.p <- weight.matrix.p[, itemp.p]
      
      w.g <- wtemp.g[match.g]
      names(w.g) <- gnames
      
      # match protein names from design matrix with weight matrix names 
      match.p <- match(pnames, paste0(colnames(weight.matrix.p), "_P"))
      
      w.p <- wtemp.p[match.p]
      names(w.p) <- pnames
      
      # actual vector of weights
      if(include.treatment == TRUE){
        
        if(treatment.penalty == FALSE){ # should the treatment information be included with or without a penalty?
          
          w <- c(rep(0, sum(!(colnames(Design) %in% c(names(w.g), names(w.p))))), # weights for the experiment settings (=0)
                 w.g, # weights for the RNA covariates
                 w.p) # weights for the protein covariates
          
        }else{
          
          w <- c(rep(1, sum(!(colnames(Design) %in% c(names(w.g), names(w.p))))), # weights for the experiment settings (=0)
                 w.g, # weights for the RNA covariates
                 w.p) # weights for the protein covariates
          
        }
        
        
      }else{
        w <- c(w.g, w.p)
      }
      
      names(w) <- colnames(Design)
      
      # should already be the case but make sure that the main gene is unweighted
      w[igene] <- 0
      
    }
    
    
    #}else if(weights == "dipa"){
    #  print("TODO")
    #}
    
  }else{ # no weights
    
    # ensure the "main gene" being selected in the lasso model
    w <-  rep(1, ncol(Design))
    
    if(include.treatment == TRUE){
      
      if(treatment.penalty == FALSE){ # should the treatment information be included with or without a penalty?
        
        w[c(1:n.design.treatment,igene)] <- 0
        
      }else{
        
        w[igene] <- 0
        
      }
      
    }else{
      w[igene] <- 0
    }
    
    if(PE.RNA.penalty == TRUE){
      w[igene] <- 1
    }
    
    names(w) <- colnames(Design)
    
  }
  
  # final steps ----------------------------------------------------------------
  
  protein.vector <- data.frame(y)
  colnames(protein.vector) <- protein
  
  
  return(list(gene = gene,
              protein = protein,
              y = y,
              protein.vector = protein.vector,
              Design = Design,
              w = w,
              n = n,
              mice = mice,
              experiments = experiments))
  
}

