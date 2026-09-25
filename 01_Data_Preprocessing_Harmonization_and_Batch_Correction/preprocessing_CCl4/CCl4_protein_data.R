#################
### CCl4 mice ###
################


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("MSstats")

library("MSstats")

require(data.table)


setwd(choose.dir())

#ccl4 <- read.csv2("Data/Raw/Protein_data/CCl4-mice.csv", header = TRUE, sep = ",")
ccl4_new <- read.table("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_BDL_curated.txt", header = TRUE, sep = "\t")
data <- read.table("C:/Users/ngoune/Documents/Projet II/Protein/Myllys_WD_Proteomics_Results.txt", header = TRUE, sep = "\t")
###########################
### MSstats processing
###########################

DTsky <- SkylinetoMSstatsFormat(ccl4_new)

sky_p <- dataProcess(DTsky)
sky_p1 <- dataProcess(DTsky, normalization = "equalizeMedians")
sky_p2 <- dataProcess(DTsky, normalization = "QUANTILE")

ccl4_new_pd <- sky_p$ProteinLevelData

dataProcessPlots(sky_p2, type = "ProfilePlot")
#TODO umsortieren nach Monaten/Treatment

dataProcessPlots(sky_p, type = "QCPlot")


save(ccl4_new_pd, file = "ccl4_new_protein_data.RData")


################################################################################


