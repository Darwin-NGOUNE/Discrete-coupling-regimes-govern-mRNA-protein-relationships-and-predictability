setwd(choose.dir())


library("edgeR")
library("limma")

# read raw gene counts matrix
#rna_raw <- read.table(gzfile("Data/Raw/RNA_seq_data/GSE167216_Raw_gene_counts_matrix.txt.gz"))

rna_raw <- read.table(gzfile("GSE167216_Raw_gene_counts_matrix.txt.gz"))

#rna_raw_new <- read.table("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Myllys_BDL_curated.txt", header = TRUE, sep = "\t")


# SAMPLE INFORMATION -----------------------------------------------------------
#rna.info <- read.table(gzfile("Data/Raw/RNA_seq_data/GSE167216_series_matrix.txt.gz"), skip = 29, nrows = 20)

rna.info <- read.table(gzfile("GSE167216_series_matrix.txt.gz"), skip = 29, nrows = 20)

mice <- rna.info[which(rna.info$V1 == "!Sample_title"),2:ncol(rna.info)]

# rename the control group to match the names of the later DTccl4 dataset
mice[mice == c("month0_control_rep1")] <- "month0_oil_rep1"
mice[mice == c("month0_control_rep2")] <- "month0_oil_rep2"
mice[mice == c("month0_control_rep3")] <- "month0_oil_rep3"
mice[mice == c("month0_control_rep4")] <- "month0_oil_rep4"
mice[mice == c("month0_control_rep5")] <- "month0_oil_rep5"
mice[mice == c("month0_control_rep6")] <- "month0_oil_rep6"

sample.info <- rna.info[which(rna.info$V1 == "!Sample_description"),2:ncol(rna.info)]

length(mice) #36 mice
colnames(mice) <- colnames(rna_raw)


# RNA COUNT NORMALIZATION ------------------------------------------------------

# with experimental groups
rna <- DGEList(counts = rna_raw, group = rep(c(6,3,4,2,5,1),6))

#design <- model.matrix(~Group)
#TODO design matrix erstellen

#v <- voom(dge, design, plot=TRUE)

#TMM normalization and logCPM transformation
#keep <- filterByExpr(rna, rna$samples)
#norm_count_matrix = count_matrix[keep,,keep.lib.sizes=F] %>%
#  calcNormFactors() %>%
#  voom(design = design)

dge <- calcNormFactors(rna, method = "TMM")

logCPM <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 3)  # Log2 Counts per Million

#not normalized
#logcpmtest <- cpm(rna, log = TRUE)



# RESHAPING FOR LATER JOINT DATASET --------------------------------------------

# add treatment and time points information to the gene count matrix
rna_dat <- logCPM
colnames(rna_dat) <- as.character(mice)


# re-structure the data 
counts <- as.vector(t(rna_dat))
gt <- rownames(rna_dat)

temp1 <- sapply(gt, FUN = rep, ncol(rna_dat))
gene_syn <- as.vector(temp1)
length(gene_syn)

mice_info <- rep(as.character(mice), nrow(logCPM))

ccl4_gc <- data.frame(gene_syn, mice_info, counts)

#check:
#ccl4_gc[which(ccl4_gc$gene_syn == "0610007P14Rik"),3]
#as.numeric(rna_dat[1,])

save(ccl4_gc, file = "Data/ccl4_gene_count_data.RData")

################################################################################
################################################################################

# add treatment and time points information to the gene count matrix
rna_dat <- rna$counts
colnames(rna_dat) <- as.character(mice)


# re-structure the data 
counts <- as.vector(t(rna_dat))
gt <- rownames(rna_dat)

temp1 <- sapply(gt, FUN = rep, ncol(rna_dat))
gene_syn <- as.vector(temp1)
length(gene_syn)

mice_info <- rep(as.character(mice), nrow(rna$counts))

ccl4_gc <- data.frame(gene_syn, mice_info, counts)

#check:
#ccl4_gc[which(ccl4_gc$gene_syn == "0610007P14Rik"),3]
#as.numeric(rna_dat[1,])

save(ccl4_gc, file = "Data/ccl4_gene_count_data.RData")













