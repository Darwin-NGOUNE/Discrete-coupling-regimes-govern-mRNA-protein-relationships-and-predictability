# Supplementary Figures 1A & 1B: Full Cohort PCA Trajectories & Global Expression Heatmaps
## (Nature Communications Supplementary Information)

[![Vector PDF - 1A](https://img.shields.io/badge/Format-Supp_Fig_1A_PDF-red.svg)](Supplementary_Figure_1a.pdf)
[![Vector PDF - 1B](https://img.shields.io/badge/Format-Supp_Fig_1B_PDF-blue.svg)](Supplementary_Figure_1b.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Figs-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figures 1A and 1B** provide the complete, uncompressed multi-omics characterization across the entire sample cohorts ($N = 54$ animals without exclusion):
* **All 18 BDL Mice:** Sham-operated controls + vehicle ($n = 6$), BDL + vehicle ($n = 6$), and BDL + AS0369 ASBT inhibitor ($n = 6$).
* **All 36 $\text{CCl}_4$ Mice:** Mineral oil vehicle controls ($n = 18$) and $\text{CCl}_4$ intoxication ($n = 18$) across 2, 6, and 12 months.

Following editorial feedback to maximize graphical visibility, the data is partitioned into two standalone publication-grade figures:
1. **Supplementary Figure 1A (Full-Cohort PCAs):** Large-scale Principal Component Analyses demonstrating disease trajectories and therapeutic protection in BDL and $\text{CCl}_4$.
2. **Supplementary Figure 1B (Full-Cohort Expression Landscapes):** Full-page $2 \times 2$ heatmaps displaying transcriptomic (RNA Log2 CPM) and proteomic (LFQ intensity) patterns across all 943 cognate pairs.

---

## 2. Component Figures & Descriptions

### 🔹 Supplementary Figure 1A: Full Cohort PCA Scores (BDL & $\text{CCl}_4$)
* **BDL PCA (Top):** PCA on all 18 BDL mice showing clear clustering of Sham vs. BDL vs. BDL + ASBTi.
* **$\text{CCl}_4$ PCA (Bottom):** PCA on all 36 $\text{CCl}_4$ mice illustrating progressive trajectory along PC1 across 2, 6, and 12 months.
* **Output File:** [`Supplementary_Figure_1a.pdf`](Supplementary_Figure_1a.pdf) (320 KB, 300 DPI vector).
* **Generating Script:** [`Generate_Supplementary_Figure_1a_and_1b.R`](Generate_Supplementary_Figure_1a_and_1b.R)

---

### 🔹 Supplementary Figure 1B: Full Cohort Expression Heatmaps ($2 \times 2$)
* **BDL Landscapes (Top Row):** RNA expression heatmap (left) and Protein abundance heatmap (right) across all 18 BDL mice.
* **$\text{CCl}_4$ Landscapes (Bottom Row):** RNA expression heatmap (left) and Protein abundance heatmap (right) across all 36 $\text{CCl}_4$ mice.
* **Output File:** [`Supplementary_Figure_1b.pdf`](Supplementary_Figure_1b.pdf) (1.36 MB, 300 DPI vector).
* **Generating Script:** [`Generate_Supplementary_Figure_1a_and_1b.R`](Generate_Supplementary_Figure_1a_and_1b.R)

---

## 3. How to Reproduce & Generate the Figures

To generate both standalone figures (`Supplementary_Figure_1a.pdf` and `Supplementary_Figure_1b.pdf`):
```bash
Rscript Generate_Supplementary_Figure_1a_and_1b.R
```

To generate the combined single-plate version (`Supplementary_Figure_1.pdf`):
```bash
Rscript Generate_Supplementary_Figure_1_Master_Plate.R
```

---

## 4. Directory Layout

```text
Supplementary_Figure_01/
│
├── 📄 Supplementary_Figure_1a.pdf                       # Standalone Figure 1A: Full PCAs (BDL & CCl4) (320 KB)
├── 📄 Supplementary_Figure_1b.pdf                       # Standalone Figure 1B: Full 2x2 Heatmaps (1.36 MB)
├── 📄 Supplementary_Figure_1.pdf                        # Combined Master Plate (1.67 MB)
│
├── 📄 Generate_Supplementary_Figure_1a_and_1b.R         # Standalone generation script (High Visibility)
├── 📄 Generate_Supplementary_Figure_1_Master_Plate.R    # Combined assembly script
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

## 5. Nature Communications Supplementary Figure Legends

### Legend for Supplementary Figure 1A:
> **Supplementary Figure 1A: Full-cohort principal component analyses across all experimental animals.**  
> Principal component analyses (PCA) of RNA expression and protein abundance across all 18 BDL mice (top, Sham + vehicle, BDL + vehicle, BDL + AS0369) and all 36 $\text{CCl}_4$ mice (bottom, mineral oil controls and $\text{CCl}_4$ treated animals at 2, 6, and 12 months).

### Legend for Supplementary Figure 1B:
> **Supplementary Figure 1B: Global transcriptomic and proteomic expression landscapes across all experimental animals.**  
> Hierarchically clustered expression heatmaps displaying all 943 shared transcripts (left) and proteins (right) for BDL (top row, $n = 18$) and $\text{CCl}_4$ (bottom row, $n = 36$) across all individual animals without sample exclusion.
