### CCl4 - Omics Integration
#############################

# working directory
setwd(choose.dir())

## load packages
library(data.table)
library(mgcv)
library("stringr")
library("xtable")
library("ggplot2")
library(UniProt.ws)


## load data

load("Data/Data_CCl4_trancriptomics_raw.RData")
load("Data/Data_CCl4_proteomics_raw.RData")


################################################################################

### DATA PREPARATION
####################

## same notation for mice for RNA and protein data
##################################################

# rename mice of protein data set 
mice_p <- paste0(ccl4_pd$GROUP,"_", paste0("rep",ccl4_pd$SUBJECT))
#mice_p[grepl("0moil", mice_p)] <- sub("oil","control",mice_p[grepl("0moil", mice_p)])

mice_p <- sub("m","_", mice_p)
mice_p <- sub("CCl", "ccl", mice_p)
mice_p <- paste0("month", mice_p)


# add information on treatment time and treatment
ccl4_pd$treatment_time <- as.numeric(str_split(ccl4_pd$GROUP, "m", simplify = TRUE)[,1])
ccl4_pd$treatment <- str_split(ccl4_pd$GROUP, "m", simplify = TRUE)[,2]


# add mice info including treatment time, treatment type and replicate mouse
ccl4_pd$mice_info <- mice_p
head(ccl4_pd)

# mice from gene count data set
mice_g <- ccl4_gc$mice_info

# also add information on treatment time and treatment for the RNA data
ccl4_gc$treatment_time <- as.numeric(str_split(sub("month", "", ccl4_gc$mice_info), "_", simplify = TRUE)[,1])
ccl4_gc$treatment <- str_split(sub("month", "", ccl4_gc$mice_info), "_", simplify = TRUE)[,2]
ccl4_gc$treatment <- sub("control", "oil", ccl4_gc$treatment)

ccl4_gc$mice_info <- sub("control", "oil", ccl4_gc$mice_info)



# change protein Names

ccl4_pd$Protein <- str_split(ccl4_pd$Protein, "\\|", simplify = TRUE)[,2]


protein.list <- unique(ccl4_pd$Protein)

pairs.list.full <- mapUniProt("UniProtKB_AC-ID", "Gene_Name", query = protein.list)

length(protein.list)
nrow(pairs.list.full)
#ATTENTION: There might be duplicated gene names for several proteins

# which genes from the full transcriptomics data can be found in the pairs list
pairs.list.present <- pairs.list.full[pairs.list.full$To %in% unique(ccl4_gc$gene_syn),]

length(gene.list)


# only consider pairs with at least one observation in every treatment/time combination
list_reduce <- rep(FALSE, nrow(pairs.list.present))
for(i in 1:nrow(pairs.list.present)){
  j_temp <- which(ccl4_pd$Protein %in% pairs.list.present$From[i])
  list_reduce[i] <- all(unique(ccl4_pd$GROUP) %in% ccl4_pd[j_temp,]$GROUP)
}

length(which(list_reduce == TRUE))
#1255 gene-protein pairs are considered

# list with the 1255 considered gene-protein pairs
pairs.list_final <- pairs.list.present[which(list_reduce == TRUE),]

colnames(pairs.list_final) <- c("ProteinID", "GeneSyn")

pairs.list_final$GeneProtein <- paste0(pairs.list_final$GeneSyn, "_", pairs.list_final$ProteinID)



### MATCH GENE AND PROTEIN DATA TO CREATE JOINT DATA SET -----------------------

# Protein data 
DTprot <- ccl4_pd[ccl4_pd$Protein %in% pairs.list_final$ProteinID,]

# Gene count data
DTgene <- ccl4_gc[ccl4_gc$gene_syn %in% pairs.list_final$GeneSyn,]


# match the data by the mice information for genes
DTgene$gp_pair <- pairs.list_final[match(DTgene$gene_syn, pairs.list_final$GeneSyn),]$GeneProtein

#for proteins
DTprot$gp_pair <- pairs.list_final[match(DTprot$Protein, pairs.list_final$ProteinID),]$GeneProtein


# add unique column for pair name plus mice info
DTgene$PairInfo <- paste0(DTgene$gp_pair, "_", DTgene$mice_info)
DTprot$PairInfo <- paste0(DTprot$gp_pair, "_", DTprot$mice_info)

# match ProteinIntensities to RNA data's mice info 
#(as there are 36 mice in the RNA data and 31 mice in the protein data)
ProteinLevels <- DTprot$LogIntensities[match(DTgene$PairInfo, DTprot$PairInfo)]




# dataset containing all relevant information
DT <- data.frame(GeneProtein = DTgene$gp_pair,
                 GeneSyn = DTgene$gene_syn,
                 ProteinID = str_split(DTgene$gp_pair, "_", simplify = TRUE)[,2],
                 MiceInfo = DTgene$mice_info,
                 Treatment = relevel(as.factor(DTgene$treatment), ref = "oil"),
                 TreatmentTime = relevel(as.factor(DTgene$treatment_time), ref = "0"),
                 GeneCount = DTgene$counts,
                 ProteinIntensity = ProteinLevels
                 )

DT <- data.table(DT)

DT[, TreatmentTotal := paste0("month",TreatmentTime, "_", Treatment)]

#TODO temporary solution for the problem of the same protein connected to multiple genes
#for(i in 1:length(unique(DT$GeneProtein))){
#  temp <- unique(DT$GeneProtein)[i]
#  tempp <- str_split(temp, "_", simplify = TRUE)[,2]
#  DTtemp_p <- DTprot[DTprot$ProteinID == tempp,]
#  
#  i_temp <- match(DT[GeneProtein == unique(DT$GeneProtein)[i],]$MiceInfo,
#                  DTtemp_p$mice_info)
#  DT[GeneProtein == unique(DT$GeneProtein)[i], ProteinIntensity := DTtemp_p[i_temp,]$LogIntensities]
#}



## final processing steps

nrow(DT)

# remove lowly expressed genes
DT <- DT[-which(DT$GeneSyn == unique(DT$GeneSyn)[which(DT[, mean(GeneCount), by = GeneSyn]$V1 == 0)]),]
nrow(DT)

save(DT, file = "Data/Final_Data/Data_CCl4_Gene_Protein_full.RData")


# DEAL WITH MISSING DATA IN THE PROTEIN LEVELS ---------------------------------


# only look at complete cases (omit NAs from the ProteinIntensities)
DTcc <- DT[complete.cases(DT),]

# number of pairs with the respective amount of observations
table(DTcc[, .N, by = GeneSyn]$N)

# how many Pairs would be dropped when considering only pairs with x observations?
cumsum(table(DTcc[, .N, by = GeneSyn]$N))


## Option 1: only consider Gene-Protein-Pairs with 24 or more observations (6 pairs dropped)
Genedrop_n <- DTcc[, .N, by = GeneSyn]$GeneSyn[which(DTcc[, .N, by = GeneSyn]$N < 24)]

## Option 2: make sure that there are at least 4 observations in each experiment setting
DT_per_setting <- DT[, sum(!is.na(ProteinIntensity)), by = .(GeneSyn, Treatment, TreatmentTime)]

# Gene names of the pairs that are removed when considering only pairs with 
# a) at least 4 observations per setting: (234 Pairs removed)
Genedrop_setting4 <- unique(DT_per_setting[which(DT_per_setting$V1 < 4),]$GeneSyn)
# b) at least 3 observations per setting: (114 Pairs removed)
Genedrop_setting3 <- unique(DT_per_setting[which(DT_per_setting$V1 < 3),]$GeneSyn)


# List of genes with at least 24 corresponding protein level observations 
genelist_n <- unique(DT$GeneSyn)[!(unique(DT$GeneSyn) %in% Genedrop_n)] 
length(genelist_n) # list of 1168 genes

genelist_3_per_setting <- unique(DT$GeneSyn)[!(unique(DT$GeneSyn) %in% Genedrop_setting3)]
length(genelist_3_per_setting)
# 1135 genes

DTccl4 <- DT[GeneSyn %in% genelist_3_per_setting,]

save(DTccl4, file = "Data/Final_Data/Data_CCl4_Gene_Protein_full_final.RData")



################################################################################

