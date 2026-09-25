################################################################################
############### GENES AND PROTEIN THAT ARE IN THE BOTH DATASET #################
################################################################################

library(dplyr)

# Load the two datasets: one for CCl4 condition and one for the new dataset
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Data/Data_CCl4_Gene_Protein_full_final_with_KMclustering.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_asbt_Gene_Protein_full.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_asbt_Gene_Protein_full_LCPM.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/DT/New_Data/Data_count_asbt_Gene_Protein_full_with_DiPa.RData")
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_asbt_Gene_Protein_full_with_DiPa.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_LCPM_asbt_Gene_Protein_full_with_DiPa.RData")

################################################################################
# Echanger clusterdipa

DT_LCPM$ClusterDiPa <- DT_count$ClusterDiPa
DT_LCPM$ClusterDiPa_heal <- DT_count$ClusterDiPa_heal
DT_LCPM$ClusterDiPa_heal2 <- DT_count$ClusterDiPa_heal2
save(DT_LCPM, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data_LCPM_asbt_Gene_Protein_full_with_DiPa.RData")


################################################################################

# Check for missing values in each dataset
colSums(is.na(DTccl4))
#colSums(is.na(DT))
colSums(is.na(DT_LCPM))

# Count how many unique ProteinIDs in DTccl4 are also found in DT
#sum(unique(DTccl4$ProteinID) %in% DT$ProteinID)  # Expected: 950
sum(unique(DTccl4$ProteinID) %in% DT_LCPM$ProteinID)

# Retrieve gene names from Ensembl using biomaRt
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_ids <- unique(DT_LCPM$GeneSyn)  # Ensembl gene IDs from DT

# Query Ensembl to get external gene names
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = gene_ids,
                   mart = ensembl)

# Add gene names to DT without changing row order
#DT$GeneName <- gene_info$external_gene_name[match(DT$GeneSyn, gene_info$ensembl_gene_id)]

DT_LCPM$GeneName <- gene_info$external_gene_name[match(DT_LCPM$GeneSyn, gene_info$ensembl_gene_id)]

# Count how many GeneSyn values in DTccl4 match GeneName values in DT
sum(unique(DTccl4$GeneSyn) %in% DT_LCPM$GeneName)  # Expected: 943

# Count how many ProteinIDs are shared between both datasets
#sum(unique(DTccl4$ProteinID) %in% DT$ProteinID)  # Expected: 950
sum(unique(DTccl4$ProteinID) %in% DT_LCPM$ProteinID)

################################################################################

# Identify common gene names and protein IDs between both datasets
#common_genes <- intersect(unique(DT$GeneName), unique(DTccl4$GeneSyn))       # 943 genes
#common_proteins <- intersect(unique(DT$ProteinID), unique(DTccl4$ProteinID)) # 950 proteins

common_genes_LCPM <- intersect(unique(DT_LCPM$GeneName), unique(DTccl4$GeneSyn))       # 943 genes
common_proteins_LCPM <- intersect(unique(DT_LCPM$ProteinID), unique(DTccl4$ProteinID)) # 950 proteins

# Filter both datasets to keep only rows with common genes and proteins
#DT_filtered <- DT[DT$GeneName %in% common_genes & DT$ProteinID %in% common_proteins, ]
DT_LCPM_filtered <- DT_LCPM[DT_LCPM$GeneName %in% common_genes_LCPM & DT_LCPM$ProteinID %in% common_proteins_LCPM, ]
DTccl4_filtered <- DTccl4[DTccl4$GeneSyn %in% common_genes_LCPM & DTccl4$ProteinID %in% common_proteins_LCPM, ]

# Create a unique identifier for each gene–protein pair
#DT_filtered$GeneProteinPair <- paste(DT_filtered$GeneName, DT_filtered$ProteinID, sep = "_")
DT_LCPM_filtered$GeneProteinPair <- paste(DT_LCPM_filtered$GeneName, DT_LCPM_filtered$ProteinID, sep = "_")
DTccl4_filtered$GeneProteinPair <- paste(DTccl4_filtered$GeneSyn, DTccl4_filtered$ProteinID, sep = "_")

# Find the intersection of gene–protein pairs between both datasets
#common_pairs <- intersect(unique(DT_filtered$GeneProteinPair), unique(DTccl4_filtered$GeneProteinPair))
#length(common_pairs)  # Should return 943 — the number of fully matched gene–protein pairs

common_pairs <- intersect(unique(DT_LCPM_filtered$GeneProteinPair), unique(DTccl4_filtered$GeneProteinPair))
length(common_pairs) 

# Check for missing values in each dataset
colSums(is.na(DTccl4_filtered))
#colSums(is.na(DT_filtered))
colSums(is.na(DT_LCPM_filtered))

################################################################################

DTccl4_filtered$GeneProteinPair <- DTccl4_filtered$GeneProtein

# Supprimer plusieurs colonnes
DTccl4_filtered <- subset(DTccl4_filtered, select = -c(ClusterGCA, ClusterKmRNA, ClusterKmProtein))


#DT_filtered <- subset(DT_filtered, select = -c(ClusterDiPa_heal, ClusterDiPa_heal2,GeneName))

DT_LCPM_filtered <- subset(DT_LCPM_filtered, select = -c(ClusterDiPa_heal, ClusterDiPa_heal2,GeneName))

#DT_filtered$GeneName <- NULL

#DT_LCPM_filtered$GeneName <- NULL

save(DTccl4_filtered, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_filtered_Gene_Protein_full.RData")
#save(DT_filtered, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_filtered_Gene_Protein_full.RData")
save(DT_LCPM_filtered, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_LCPM_filtered_Gene_Protein_full.RData")

################################################################################
################################################################################
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_filtered_Gene_Protein_full.RData")
#load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_filtered_Gene_Protein_full.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_LCPM_filtered_Gene_Protein_full.RData")


#colnames(DT_filtered)[colnames(DT_filtered) == "GeneSyn"] <- "GeneSynOld"
#colnames(DT_filtered)[colnames(DT_filtered) == "GeneProtein"] <- "GeneProteinOld"
#colnames(DT_filtered)[colnames(DT_filtered) == "GeneName"] <- "GeneSyn"
#colnames(DT_filtered)[colnames(DT_filtered) == "GeneProteinPair"] <- "GeneProtein"


colnames(DT_LCPM_filtered)[colnames(DT_LCPM_filtered) == "GeneSyn"] <- "GeneSynOld"
colnames(DT_LCPM_filtered)[colnames(DT_LCPM_filtered) == "GeneProtein"] <- "GeneProteinOld"
colnames(DT_LCPM_filtered)[colnames(DT_LCPM_filtered) == "GeneName"] <- "GeneSyn"
colnames(DT_LCPM_filtered)[colnames(DT_LCPM_filtered) == "GeneProteinPair"] <- "GeneProtein"


#DT_filtered <- subset(DT_filtered, select = -c(GeneSynOld,GeneProteinOld))

DT_LCPM_filtered <- subset(DT_LCPM_filtered, select = -c(GeneSynOld,GeneProteinOld))

#DT_filtered$GeneSyn <- sapply(strsplit(DT_filtered$GeneProtein, "_"), `[`, 1)

DT_LCPM_filtered$GeneSyn <- sapply(strsplit(DT_LCPM_filtered$GeneProtein, "_"), `[`, 1)


DTccl4_filtered$GeneProteinPair <- NULL

DTccl4_filtered_final <- DTccl4_filtered
#DT_filtered_final <- DT_filtered
DT_LCPM_filtered_final <- DT_LCPM_filtered

save(DTccl4_filtered_final, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_filtered_final_Gene_Protein_full.RData")
#save(DT_filtered_final, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_filtered_final_Gene_Protein_full.RData")
save(DT_LCPM_filtered_final, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DT_LCPM_filtered_final_Gene_Protein_full.RData")

################################################################################


library(dplyr)

# Premier jeu : combiner DTccl4_filtered et DT_filtered
#DTccl4_DT <- bind_rows(DTccl4_filtered_final, DT_filtered_final)

# Deuxième jeu : combiner DTccl4_filtered et DT_LCPM_filtered
DTccl4_DT_LCPM <- bind_rows(DTccl4_filtered_final, DT_LCPM_filtered_final)

#save(DTccl4_DT, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_Gene_Protein_full.RData")
save(DTccl4_DT_LCPM, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData")


length(unique(DTccl4_DT_LCPM$GeneProtein))
colSums(is.na(DTccl4_DT_LCPM))
length(unique(DTccl4_DT_LCPM$MiceInfo))


###
#setequal(DT_filtered_final$GeneProteinPair, DTccl4_filtered_final$GeneProteinPair)

# Valeurs présentes dans df1 mais pas dans df2
setdiff(DT_LCPM_filtered_final$GeneProteinPair, DTccl4_filtered_final$GeneProteinPair)



