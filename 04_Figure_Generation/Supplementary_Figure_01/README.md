# Supplementary Figure 1: Full Cohort PCA Trajectories & Global Expression Heatmaps (All 18 BDL & 36 $\text{CCl}_4$ Mice)
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_1.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_1-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 1** provides the complete, uncompressed multi-omics characterization across the entire sample cohorts:
* **All 18 BDL Mice:** Sham-operated controls + vehicle ($n = 6$), BDL + vehicle ($n = 6$), and BDL + AS0369 ASBT inhibitor ($n = 6$).
* **All 36 $\text{CCl}_4$ Mice:** Mineral oil vehicle controls ($n = 18$) and $\text{CCl}_4$ intoxication ($n = 18$) across 2, 6, and 12 months.

It presents:
1. **Full-Cohort Principal Component Analyses (Panel A):** Demonstrating reproducible disease separation in both models without sample exclusion.
2. **Full-Cohort Expression Landscapes (Panel B):** Complete $2 \times 2$ heatmaps displaying transcriptomic and proteomic expression across all 54 animals.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Left Column): Full Cohort PCA Scores (BDL & $\text{CCl}_4$)
* **BDL PCA (Top Left):** PCA on all 18 BDL mice showing clear clustering of Sham vs. BDL vs. BDL + ASBTi.
  * *Component File:* `Isolated_PCA_BDL_18Mice_normal_data.pdf` (26 KB)
  * *Generating Script:* [`Generate_A4_MASTER_BDL_Schwerpunkt_ALL_18Mice.R`](Generate_A4_MASTER_BDL_Schwerpunkt_ALL_18Mice.R)
* **$\text{CCl}_4$ PCA (Bottom Left):** PCA on all 36 $\text{CCl}_4$ mice illustrating progressive trajectory along PC1 across 2, 6, and 12 months.
  * *Component File:* `Isolated_PCA_CCL4_36Mice.pdf` (26 KB)
  * *Generating Script:* [`Generate_A4_MASTER_CCL4_Schwerpunkt_ALL_36Mice.R`](Generate_A4_MASTER_CCL4_Schwerpunkt_ALL_36Mice.R)

---

### 🔹 Panel B (Right Column): Full Cohort Expression Heatmaps ($2 \times 2$)
* **BDL Landscapes (Top Right):** RNA expression heatmap (Page 1) and Protein abundance heatmap (Page 1) across all 18 BDL mice.
* **$\text{CCl}_4$ Landscapes (Bottom Right):** RNA expression heatmap (Page 2) and Protein abundance heatmap (Page 2) across all 36 $\text{CCl}_4$ mice.
* *Component Files:* `Global_RNA_Landscape_Heatmaps_ALL_Mice.pdf` (389 KB), `Global_Protein_Landscape_Heatmaps_ALL_Mice.pdf` (343 KB).
* *Generating Scripts:*
  * [`Global_RNA_Landscape_Heatmaps_ALL_Mice.R`](Global_RNA_Landscape_Heatmaps_ALL_Mice.R)
  * [`Global_Protein_Landscape_Heatmaps_ALL_Mice.R`](Global_Protein_Landscape_Heatmaps_ALL_Mice.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_1_Master_Plate.R`](Generate_Supplementary_Figure_1_Master_Plate.R)
* **Execution:** Combines Panel A (left column) and Panel B (right column) at 300 DPI into a 16:9 widescreen canvas using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_1_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_1.pdf` (1.61 MB).

---

## 4. Directory Layout

```text
Supplementary_Figure_01/
│
├── 📄 Supplementary_Figure_1.pdf                        # Master Publication Plate (1.61 MB)
├── 📄 Generate_Supplementary_Figure_1_Master_Plate.R    # Master assembly script (cowplot)
│
├── 📄 Isolated_PCA_BDL_18Mice_normal_data.pdf           # Panel A1 (BDL PCA)
├── 📄 Isolated_PCA_CCL4_36Mice.pdf                      # Panel A2 (CCl4 PCA)
├── 📄 Global_RNA_Landscape_Heatmaps_ALL_Mice.pdf        # Panel B (RNA Heatmaps)
├── 📄 Global_Protein_Landscape_Heatmaps_ALL_Mice.pdf    # Panel B (Protein Heatmaps)
│
├── 📄 Generate_A4_MASTER_BDL_Schwerpunkt_ALL_18Mice.R   # Generates BDL PCA
├── 📄 Generate_A4_MASTER_CCL4_Schwerpunkt_ALL_36Mice.R  # Generates CCl4 PCA
├── 📄 Global_RNA_Landscape_Heatmaps_ALL_Mice.R          # Generates RNA heatmaps
├── 📄 Global_Protein_Landscape_Heatmaps_ALL_Mice.R      # Generates Protein heatmaps
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 1: Full-cohort principal component analyses and global transcriptomic and proteomic heatmaps across all experimental animals.**  
> **(A)** Principal component analyses (PCA) of RNA expression and protein abundance across all 18 BDL mice (top, Sham + vehicle, BDL + vehicle, BDL + AS0369) and all 36 $\text{CCl}_4$ mice (bottom, mineral oil controls and $\text{CCl}_4$ treated animals at 2, 6, and 12 months).  
> **(B)** Hierarchically clustered expression heatmaps displaying all 943 shared transcripts (left) and proteins (right) for BDL (top row) and $\text{CCl}_4$ (bottom row) across all individual animals.
