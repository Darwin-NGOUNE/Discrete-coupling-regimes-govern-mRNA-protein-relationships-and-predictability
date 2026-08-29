library(fmsb)
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggforce)


load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data/Data_count_filtered_asbt_Gene_Protein_full.RData")

load("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Data/Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData")


DTccl4_count_filtered <- DTccl4_count_filtered[!is.na(DTccl4_count_filtered$ProteinIntensity), ]

DT <- DT_count_filtered
DT <- DTccl4_count_filtered

colSums(is.na(DT))


################################################################################
## DT_count_filtered
################################################################################
# DIPA DATA PREPARATION --------------------------------------------------------
## (a) not normalized ----------------------------------------------------------

# build dataset for differentiation pattern analysis starting with mean gene expression 
# and mean protein intensity (raw gene counts and back-transformed log2 protein intensities)
DT_dipa <- DT[, .(mean(GeneCount, na.rm=TRUE),
                  mean(2^ProteinIntensity, na.rm=TRUE)), by = .(GeneProtein)]

DT_dipa[, PearsonCor := DT[,cor(GeneCount, 2^ProteinIntensity, method = "pearson", use="pairwise.complete.obs"), by=GeneProtein]$V1]

colnames(DT_dipa) <- c("GeneProtein", "meanGE", "meanPI", "PearsonCor")

# means for the respective treatments
DT_dipa[, `:=`(meanG_control = DT[Treatment == "control", mean(GeneCount), by = GeneProtein]$V1,
               meanG_BDL = DT[Treatment == "BDL", mean(GeneCount), by = GeneProtein]$V1,
               meanG_BDL_ASBTi = DT[Treatment == "BDL_ASBTi", mean(GeneCount), by = GeneProtein]$V1,
               meanP_control = DT[Treatment == "control", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1,
               meanP_BDL = DT[Treatment == "BDL", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1,
               meanP_BDL_ASBTi = DT[Treatment == "BDL_ASBTi", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]


#colSums(is.na(DT_dipa))

# overall mean ratio BDL to control
DT_dipa[, `:=`(meanRatioG = DT[Treatment == "BDL", mean(GeneCount), by = GeneProtein]$V1 / DT[Treatment == "control", mean(GeneCount), by = GeneProtein]$V1,
               meanRatioP = DT[Treatment == "BDL", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 / DT[Treatment == "control", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]


# overall mean ratio BDL_ASBTi to control
DT_dipa[, `:=`(meanRatioG_heal = DT[Treatment == "BDL_ASBTi", mean(GeneCount), by = GeneProtein]$V1 / DT[Treatment == "control", mean(GeneCount), by = GeneProtein]$V1,
               meanRatioP_heal = DT[Treatment == "BDL_ASBTi", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 / DT[Treatment == "control", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]


# overall mean ratio BDL_ASBTi to BDL
DT_dipa[, `:=`(meanRatioG_heal2 = DT[Treatment == "BDL_ASBTi", mean(GeneCount), by = GeneProtein]$V1 / DT[Treatment == "BDL", mean(GeneCount), by = GeneProtein]$V1,
               meanRatioP_heal2 = DT[Treatment == "BDL_ASBTi", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 / DT[Treatment == "BDL", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]

################################################################################
## DTccl4_count_filtered
################################################################################

DT_dipa <- DT[, .(mean(GeneCount, na.rm=TRUE),
                  mean(2^ProteinIntensity, na.rm=TRUE)), by = .(GeneProtein)]

DT_dipa[, PearsonCor := DT[,cor(GeneCount, 2^ProteinIntensity, method = "pearson", use="pairwise.complete.obs"), by=GeneProtein]$V1]

colnames(DT_dipa) <- c("GeneProtein", "meanGE", "meanPI", "PearsonCor_old")

# Median between treatments for Gene Counts and Protein Intensities
# and binary variables for correlations and differences

# means for the respective treatments
DT_dipa[, `:=`(meanG_control = DT[Treatment == "oil", mean(GeneCount), by = GeneProtein]$V1,
               meanG_ccl4 = DT[Treatment == "ccl4", mean(GeneCount), by = GeneProtein]$V1,
               meanP_control = DT[Treatment == "oil", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1,
               meanP_ccl4 = DT[Treatment == "ccl4", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]



DT_dipa[, `:=`(meanRatioG = DT[Treatment == "ccl4", mean(GeneCount, na.rm=TRUE), by = GeneProtein]$V1 / DT[Treatment == "oil", mean(GeneCount, na.rm=TRUE), by = GeneProtein]$V1,
               meanRatioP = DT[Treatment == "ccl4", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 / DT[Treatment == "oil", mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]

# for the durations of 2 and 12 months
DT_dipa[, `:=`(meanRatioG_2 = DT[Treatment == "ccl4" & TreatmentTime == 2, mean(GeneCount), by = GeneProtein]$V1 /
                 DT[Treatment == "oil" & TreatmentTime == 2, mean(GeneCount), by = GeneProtein]$V1,
               meanRatioP_2 = DT[Treatment == "ccl4" & TreatmentTime == 2, mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 /
                 DT[Treatment == "oil" & TreatmentTime == 2, mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]

DT_dipa[, `:=`(meanRatioG_12 = DT[Treatment == "ccl4" & TreatmentTime == 12, mean(GeneCount), by = GeneProtein]$V1 /
                 DT[Treatment == "oil" & TreatmentTime == 12, mean(GeneCount), by = GeneProtein]$V1,
               meanRatioP_12 = DT[Treatment == "ccl4" & TreatmentTime == 12, mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1 /
                 DT[Treatment == "oil" & TreatmentTime == 12, mean(2^ProteinIntensity, na.rm=TRUE), by = GeneProtein]$V1
)]



################################################################################
################################################################################


#TODO tune thresholds, for now chose them manually 
thresh_x <- 0.5
thresh_y <- 0.5
thresh_d <- 0.5

dipagrouping <- function(G,P,thresh_x,thresh_y,thresh_d, type = "full"){
  dipa.group <- c()
  if(type == "full"){
    # 0: constant
    dipa.group[which(abs(log2(G)) < thresh_x & abs(log2(P)) < thresh_y)] <- 0
    
    # I: RNA and Protein upregulation
    dipa.group[which(log2(G) > thresh_x 
                     & log2(P) > thresh_y)] <- 1
    
    # II: RNA and Protein downregulation
    dipa.group[which(log2(G) < -thresh_x 
                     & log2(P) < -thresh_y)] <- 2
    
    # III: Only Protein upregulation
    dipa.group[which(abs(log2(G)) <= thresh_x 
                     & log2(P) > thresh_y)] <- 3
    
    # IV: Only Protein downregulation
    dipa.group[which(abs(log2(G)) <= thresh_x 
                     & log2(P) < -thresh_y)] <- 4
    
    # V: Only RNA upregulation
    dipa.group[which(log2(G) > thresh_x 
                     & abs(log2(P)) < thresh_y)] <- 5
    
    # VI: Only RNA downregulation
    dipa.group[which(log2(G) < -thresh_x 
                     & abs(log2(P)) < thresh_y)] <- 6
    
    # VII: "crazy" group
    dipa.group[which(log2(G) <= -thresh_x 
                     & log2(P) >= thresh_y)] <- 7
    dipa.group[which(log2(G) >= thresh_x 
                     & log2(P) <= -thresh_y)] <- 7
  }else if(type == "belt"){
    dipa.group[which(log2(G)-log2(P) >= thresh_d)] <- 1
    dipa.group[which(log2(G)-log2(P) <= -thresh_d)] <- 2 
    dipa.group[which(abs(log2(G)-log2(P)) < thresh_d)] <- 0
  }
  
  return(dipa.group)
}

################################################################################
################################################################################
meanRatioGroups <- dipagrouping(G=DT_dipa$meanRatioG_12, P=DT_dipa$meanRatioP_12, thresh_x=thresh_x,
                                thresh_y=thresh_y, thresh_d=thresh_d)


DT_dipa[, DiPaGroups := meanRatioGroups]

DT_dipa_count_ccl4 <- DT_dipa

colSums(is.na(DT_dipa))
################################################################################
################################################################################

meanRatioGroups <- dipagrouping(G=DT_dipa$meanRatioG, P=DT_dipa$meanRatioP, thresh_x=thresh_x,
                                thresh_y=thresh_y, thresh_d=thresh_d)

DT_dipa[, DiPaGroups := meanRatioGroups]

DT_dipa_count_bdl <- DT_dipa

colSums(is.na(DT_dipa_count_bdl))
################################################################################
################################################################################
# add DiPa groups for the pairwise treatment comparisons to the full data as DiPa cluster groups
DT_count_filtered[, ClusterDiPa := DT_dipa$DiPaGroups[match(DT_count_filtered$GeneProtein, DT_dipa$GeneProtein)]]

DTccl4_count_filtered[, ClusterDiPa := DT_dipa$DiPaGroups[match(DTccl4_count_filtered$GeneProtein, DT_dipa$GeneProtein)]]

save(DTccl4_count_filtered, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/Data/Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData")
save(DT_count_filtered, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Data/Data_count_filtered_asbt_Gene_Protein_full.RData")


save(DT_dipa_count_bdl, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData")
save(DT_dipa_count_ccl4, file = "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/Data_CCL4_filtered_Gene_Protein_full_final_with_dipa.RData")

################################################################################
################################################################################


p1 <- ggplot(data = DT_dipa, 
             aes(x = log2(meanRatioG), y = log2(meanRatioP), color = as.factor(DiPaGroups))) +
  #geom_abline(intercept = thresh_y, slope = 1, color = "black", linewidth = 0.5) +
  #geom_abline(intercept = -thresh_y, slope = 1, color = "black", linewidth = 0.5) +
  geom_point() +
  geom_hline(yintercept=thresh_y, color = "black", linewidth = 0.5) +
  geom_hline(yintercept=-thresh_y, color = "black", linewidth = 0.5) +
  geom_vline(xintercept=thresh_x, color = "black", linewidth = 0.5) +
  geom_vline(xintercept=-thresh_x, color = "black", linewidth = 0.5) +
  xlab("Genes' mean ratio BDL/control (log2, 943 G_P)") + #
  ylab("Proteins' mean ratio BDL/control (log2, 943 G_P)") + #
  xlim(-4.5, 6.8) +
  ylim(-3,5) +
  annotate(geom="text", x=0, y=0, label="0", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=5.5, y=4.5, label="1", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=-3, y=-2.5, label="2", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = 0, y = 4.5, label="3", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=0, y=-2.5, label="4", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = 5.5, y = 0, label="5", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=-3, y=0, label="6", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = -3, y = 4.5, label="7", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=5.5, y=-2.5, label="7", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  labs(color='DiPa-Group')+
  theme(text = element_text(size = 20))



pdf("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/DiPa_newdata_BDL_control_count_943.pdf", height = 8,width=12)
p1
dev.off()

################################################################################
################################################################################

p2 <- ggplot(data = DT_dipa, 
             aes(x = log2(meanRatioG_12), y = log2(meanRatioP_12), color = as.factor(DiPaGroups))) +
  #geom_abline(intercept = thresh_y, slope = 1, color = "black", linewidth = 0.5) +
  #geom_abline(intercept = -thresh_y, slope = 1, color = "black", linewidth = 0.5) +
  geom_point() +
  geom_hline(yintercept=thresh_y, color = "black", linewidth = 0.5) +
  geom_hline(yintercept=-thresh_y, color = "black", linewidth = 0.5) +
  geom_vline(xintercept=thresh_x, color = "black", linewidth = 0.5) +
  geom_vline(xintercept=-thresh_x, color = "black", linewidth = 0.5) +
  xlab("Genes' mean ratio ccl4/control for 12 months duration (log2, 943 G_P)") + #
  ylab("Proteins' mean ratio ccl4/control for 12 months duration (log2, 943 G_P)") + #
  xlim(-4.5, 6.8) +
  ylim(-3,5) +
  annotate(geom="text", x=0, y=0, label="0", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=5.5, y=4.5, label="1", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=-3, y=-2.5, label="2", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = 0, y = 4.5, label="3", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=0, y=-2.5, label="4", color= 1, cex = 5) +
  geom_circle(aes(x0 = 0, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = 5.5, y = 0, label="5", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=-3, y=0, label="6", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = 0, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x = -3, y = 4.5, label="7", color= 1, cex = 5) +
  geom_circle(aes(x0 = -3, y0 = 4.5, r = 0.25), inherit.aes = FALSE) +
  annotate(geom="text", x=5.5, y=-2.5, label="7", color= 1, cex = 5) +
  geom_circle(aes(x0 = 5.5, y0 = -2.5, r = 0.25), inherit.aes = FALSE) +
  labs(color='DiPa-Group')+
  theme(text = element_text(size = 20))


pdf("C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Dipa/DiPa_ccl4_oil_12_count_943.pdf", height = 10,width=12)
p2
dev.off()



