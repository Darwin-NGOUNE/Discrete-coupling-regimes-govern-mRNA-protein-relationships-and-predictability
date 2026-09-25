# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SWATH2stats")
# 
# 
# library(SWATH2stats)
# library(data.table)
# data('Spyogenes', package = 'SWATH2stats')

library(stringr)

# load dataset
#DTasbt_proteins <- read.delim2("Data/Myllys_Data/Myllys_BDL_curated.txt")
DTasbt_proteins <- read.delim2("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_BDL_curated.txt")

# log2 transformation
DTasbt_proteins[,-1] <- apply(DTasbt_proteins[,-1], 2, function(x) log2(as.numeric(x)))

# Protein IDs
DTasbt_proteins$ProteinID <- str_split(DTasbt_proteins[,1], "\\|", simplify = TRUE)[,2]

# change sample mice names to match with the mRNA data
colnames(DTasbt_proteins) <- str_replace_all(colnames(DTasbt_proteins), "\\." ,"")
colnames(DTasbt_proteins) <- str_replace_all(colnames(DTasbt_proteins), "sham" ,"Sham")

# DTasbt_proteins_long <- melt(DTasbt_proteins,
#                              id.vars = c("Protein", "ProteinID"),
#                              variable.name = "mice_info",
#                              value.name = "expression")
# 
# DTasbt_proteins <- DTasbt_proteins_long

#save(DTasbt_proteins, file = "Data/Myllys_Data/Myllys_protein_data.RDATA")
save(DTasbt_proteins, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_protein_data.RDATA")




