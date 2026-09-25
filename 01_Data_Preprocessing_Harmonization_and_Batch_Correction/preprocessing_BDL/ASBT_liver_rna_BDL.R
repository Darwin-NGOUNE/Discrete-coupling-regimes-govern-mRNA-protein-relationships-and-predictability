############################################################
###             RNA-Seq analysis pipeline                ###
###              Manual Step 10 : DESeq2                 ###
############################################################

# Contact: Julia Duda (duda@statistik.tu-dortmund.de)

# See here for what is happening (check if this is the most current version!)
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

################################################################################
# Setup
################################################################################

# Exchange this path accordingly so that it points at your
# 08_DESeq2_analysis directory:

setwd(choose.dir())

# If you do not have the DESeq2-package installed yet, run once:
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("DESeq2", force = TRUE)

# If you do not have the apeglm-package installed yet, run once:
# if (!requireNamespace('BiocManager', quietly = TRUE))
# install.packages('BiocManager')
# BiocManager::install("ashr")

# If you do not have the Enhanced-Volcano package installed yet, run once:
# if (!requireNamespace('BiocManager', quietly = TRUE))
# install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')

# If you do not have the topGO package installed yet, run once:
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("topGO")


library(stringr) 
library(DESeq2)
library(ashr) # For shrinkage of "only due to noise" large log2FCs
library(dplyr)
library(EnhancedVolcano)
library(ggplot2); theme_set(theme_bw())
library(patchwork) # to easily place ggplots next to each other
library(ggrepel) # to label samples nicely in PCA plot
library(writexl)
library(topGO) # for gene set enrichment / GO analysis
source("./00_functions.R")
#load("Data/Myllys_Data/tximeta_result_gse.RData")
load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/tximeta_result_gse.RData")
################################################################################
# Defining control and treatment
################################################################################

# Often, it is clear what is the control and what is the treatment.
# For example, we consider "WT" as the control genotype and "KO" as
# the treatment genotype.
# However, R does not now and simply checks which name comes first alphabetically.
# This might confuse later interpretations of the Fold-Change-calculations.

# We therefore define the control and treatment levels

# Take a look at the experiment:
colData(gse)

################################################################################
# Defining the model + Minimal pre-filtering of genes + Fit the model
################################################################################

# For this experiment, the following comparisons were of interest:

# We generate small groups of combinations for treatment and genotype
# to work around using interaction models

# gse$group <- paste0(gse$treatment,"_",gse$genotype) # not needed for these samples

# Tell R which is the reference or control level.
# If not specified, R will later just order the levels alphabetically.
# We now have the following names of levels or groups:

unique(gse$treatment)
# And we want "CHOW_WT" as reference or negative control level:
gse$treatment <- factor(gse$treatment, levels = c("BDL_vehicle", "BDL_ASBTi", "sham_vehicle"))

# We define the model now:
# (Put factor of interest at the end)
dds <- DESeqDataSet(gse, design = ~ treatment)

# Minimal pre-filtering
# We have this many genes:
nrow(dds)
# But we keep only genes that have (across all mice) at least 10 counts.
# This threshold is arbitrary and could be made larger (=stricter),
# e.g. if the data still seems noisy or if there are many samples.

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
# After filtering, we have this many genes left:
nrow(dds)
# [1] 16855

# Finally, fit the model (All p-values and fold changes are now
# calculated and must only be extracted correctly)

# Takes a bit
dds <- DESeq(dds)


# Counts
DTasbt_rna_count <- data.frame(assays(dds)[[1]])


################################################################################

colnames(DTasbt_rna) <- str_replace(colnames(DTasbt_rna), ".liver", "")

colnames(DTasbt_rna) <- str_replace_all(colnames(DTasbt_rna), "-", "")

colnames(DTasbt_rna) <- str_replace_all(colnames(DTasbt_rna), "\\.", "")

#save(DTasbt_rna, file = "Data/Myllys_Data/Myllys_rna_data.RDATA")
save(DTasbt_rna, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_rna_data.RDATA")




library("edgeR")
rna <- DGEList(counts = counts, group = rep(1:3, each=6))

dge <- calcNormFactors(rna, method = "TMM")
DTasbt_rna <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 3)

#counts(dds, normalized = TRUE)[1,1]


# adjust the sample names to match the protein data
colnames(DTasbt_rna) <- str_replace(colnames(DTasbt_rna), ".liver", "")

colnames(DTasbt_rna) <- str_replace_all(colnames(DTasbt_rna), "-", "")

colnames(DTasbt_rna) <- str_replace_all(colnames(DTasbt_rna), "\\.", "")

save(DTasbt_rna, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_rna_LCPM_data.RData")














