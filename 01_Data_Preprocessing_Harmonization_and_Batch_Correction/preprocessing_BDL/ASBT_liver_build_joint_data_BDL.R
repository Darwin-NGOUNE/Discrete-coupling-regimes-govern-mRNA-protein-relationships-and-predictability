## load packages
library(data.table)
library(mgcv)
library("stringr")
library("xtable")
library("ggplot2")
library(UniProt.ws)

#setwd(choose.dir())


#load("Data/Myllys_Data/Myllys_protein_data.RDATA")
#load("Data/Myllys_Data/Myllys_rna_data.RDATA")

load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_protein_data.RDATA")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_rna_data.RDATA")

load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_rna_LCPM_data.RDATA")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
#BiocManager::install(c("ensembldb", "org.Mm.eg.db", "biomaRt"))


library(ensembldb)
library(org.Mm.eg.db)
library(biomaRt)

### DATA PREPARATION
####################

## same notation for mice for RNA and protein data
##################################################

# mice of protein data set 
mice_p <- colnames(DTasbt_proteins)[-c(1,26)]


# mice from gene count data set
mice_g <- colnames(DTasbt_rna)



# update Protein data to only include mice available in both datasets
rownames(DTasbt_proteins) <- DTasbt_proteins$ProteinID
DTasbt_proteins <- DTasbt_proteins[,!colnames(DTasbt_proteins) %in% c("Protein", "ProteinID")]
DTasbt_proteins_updated <- DTasbt_proteins[,mice_p %in% mice_g]



# Match GeneSyn and Protein IDs from protein data
protein.list <- unique(DTasbt_proteins$ProteinID)

pairs.list.full <- mapUniProt("UniProtKB_AC-ID", "Gene_Name", query = protein.list) # unique 1527

pairs.list.ensembl <- mapUniProt("UniProtKB_AC-ID", "Ensembl", query = protein.list)
pairs.list.ensembl$To <- str_split(pairs.list.ensembl$To, "\\.", simplify = TRUE)[,1] #UNIQUE 1441

head(pairs.list.full)
head(pairs.list.ensembl)

# subset Protein data to only include matched protein
DTasbt_proteins_final <- DTasbt_proteins_updated[rownames(DTasbt_proteins_updated) %in% pairs.list.ensembl$From, ]

# subset RNA data to only include matched mRNA and simultaneously adjust column order to match protein data
DTasbt_rna_final <- DTasbt_rna[rownames(DTasbt_rna) %in% pairs.list.ensembl$To, 
                               intersect(colnames(DTasbt_proteins_final) , colnames(DTasbt_rna)), drop = FALSE]

################################################################################

# TEMP WORKAROUND WITHOUT DUPLICATES
map_unique <- pairs.list.ensembl[!duplicated(pairs.list.ensembl$From),]


mice.info <- colnames(DTasbt_proteins_final)[!colnames(DTasbt_proteins_final) %in% c("Protein", "ProteinID")]

# should be same order:
colnames(DTasbt_rna_final)


#  all combinations (gene, protein, mouse)
DT <- do.call(rbind, lapply(seq_len(nrow(map_unique)), function(i) {
  data.frame(
    GeneSyn = map_unique$To[i], # eigentlich EnsemblIDs, hier aber noch GeneSyn um mit DTccl4 übereinzustimmen
    ProteinID = map_unique$From[i],
    MiceInfo = mice.info,
    stringsAsFactors = FALSE
  )
}))


DT$Treatment <- str_split(DT$MiceInfo, "M", simplify = TRUE)[,1]
DT$Treatment[DT$Treatment == "Shamvehicle"] <- "control"
DT$Treatment[DT$Treatment == "BDLvehicle"] <- "BDL"
DT$Treatment[DT$Treatment == "BDLASBTi"] <- "BDL_ASBTi"

# get mouse values for the respective mRNA and protein
row.rna <- match(DT$GeneSyn, rownames(DTasbt_rna_final)) #Si DT$GeneSyn contient des gènes qui n’existent pas dans rownames(DTasbt_rna_final), alors match() retourne NA.
col.mice <- match(DT$MiceInfo, colnames(DTasbt_rna_final))
row.protein <- match(DT$ProteinID, rownames(DTasbt_proteins_final))

DT$GeneCount <- DTasbt_rna_final[cbind(row.rna, col.mice)] # 90 Missing values
DT$ProteinIntensity <- DTasbt_proteins_final[cbind(row.protein, col.mice)]

DT$GeneProtein <- paste0(DT$GeneSyn, "_", DT$ProteinID)

#data.table
DT <- data.table(DT)

###
# DT$GeneProtein <- paste(DT$GeneSyn, DT$ProteinID, sep = "_")
# DT <- DT[, c("GeneProtein", setdiff(names(DT), "GeneProtein"))]
# 
# DT$TreatmentTotal <- sub("M\\d+$", "", DT$MiceInfo)
# DT$TreatmentTotal <- paste(DT$TreatmentTotal,DT$Treatment, sep = "_")
# DT$Treatment <- as.factor(DT$Treatment)
# 
# DT$TreatmentTime <- 1
###

#DT_count <- DT

save(DT, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_asbt_Gene_Protein_full.RData")

save(DT_count, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_count_asbt_Gene_Protein_full.RData")

################################################################################
## Same Pipeline like for DTccl4, Log2 with TMM
################################################################################

mice_p <- colnames(DTasbt_proteins)[-c(1,26)]

mice_g_LCPM <- colnames(DTasbt_rna_LCPM)


# update Protein data to only include mice available in both datasets
rownames(DTasbt_proteins) <- DTasbt_proteins$ProteinID
DTasbt_proteins <- DTasbt_proteins[,!colnames(DTasbt_proteins) %in% c("Protein", "ProteinID")]
DTasbt_proteins_updated <- DTasbt_proteins[,mice_p %in% mice_g_LCPM]

# Match GeneSyn and Protein IDs from protein data
protein.list <- unique(DTasbt_proteins$ProteinID)

pairs.list.full <- mapUniProt("UniProtKB_AC-ID", "Gene_Name", query = protein.list) # unique 1527

pairs.list.ensembl <- mapUniProt("UniProtKB_AC-ID", "Ensembl", query = protein.list)
pairs.list.ensembl$To <- str_split(pairs.list.ensembl$To, "\\.", simplify = TRUE)[,1] #UNIQUE 1441

head(pairs.list.full)
head(pairs.list.ensembl)


# subset Protein data to only include matched protein
DTasbt_proteins_final <- DTasbt_proteins_updated[rownames(DTasbt_proteins_updated) %in% pairs.list.ensembl$From, ]

# subset RNA data to only include matched mRNA and simultaneously adjust column order to match protein data
DTasbt_rna_final <- DTasbt_rna_LCPM[rownames(DTasbt_rna_LCPM) %in% pairs.list.ensembl$To, 
                                    intersect(colnames(DTasbt_proteins_final) , colnames(DTasbt_rna_LCPM)), drop = FALSE]


# TEMP WORKAROUND WITHOUT DUPLICATES
map_unique <- pairs.list.ensembl[!duplicated(pairs.list.ensembl$From),]


mice.info <- colnames(DTasbt_proteins_final)[!colnames(DTasbt_proteins_final) %in% c("Protein", "ProteinID")]

# should be same order:
colnames(DTasbt_rna_final)

#  all combinations (gene, protein, mouse)
DT_LCPM <- do.call(rbind, lapply(seq_len(nrow(map_unique)), function(i) {
  data.frame(
    GeneSyn = map_unique$To[i], # eigentlich EnsemblIDs, hier aber noch GeneSyn um mit DTccl4 übereinzustimmen
    ProteinID = map_unique$From[i],
    MiceInfo = mice.info,
    stringsAsFactors = FALSE
  )
}))


DT_LCPM$Treatment <- str_split(DT_LCPM$MiceInfo, "M", simplify = TRUE)[,1]
DT_LCPM$Treatment[DT_LCPM$Treatment == "Shamvehicle"] <- "control"
DT_LCPM$Treatment[DT_LCPM$Treatment == "BDLvehicle"] <- "BDL"
DT_LCPM$Treatment[DT_LCPM$Treatment == "BDLASBTi"] <- "BDL_ASBTi"

# get mouse values for the respective mRNA and protein
row.rna <- match(DT_LCPM$GeneSyn, rownames(DTasbt_rna_final)) #Si DT$GeneSyn contient des gènes qui n’existent pas dans rownames(DTasbt_rna_final), alors match() retourne NA.
col.mice <- match(DT_LCPM$MiceInfo, colnames(DTasbt_rna_final))
row.protein <- match(DT_LCPM$ProteinID, rownames(DTasbt_proteins_final))

DT_LCPM$GeneCount <- DTasbt_rna_final[cbind(row.rna, col.mice)] # 90 Missing values
DT_LCPM$ProteinIntensity <- DTasbt_proteins_final[cbind(row.protein, col.mice)]

DT_LCPM$GeneProtein <- paste0(DT_LCPM$GeneSyn, "_", DT_LCPM$ProteinID)

#data.table
DT_LCPM <- data.table(DT_LCPM)

save(DT_LCPM, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_asbt_Gene_Protein_full_LCPM.RData")

colSums(is.na(DT_LCPM))
length(unique(DT_LCPM$GeneProtein))

################################################################################
################################################################################
################################################################################

