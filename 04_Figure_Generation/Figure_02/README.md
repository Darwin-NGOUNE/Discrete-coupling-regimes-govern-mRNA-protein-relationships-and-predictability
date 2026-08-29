# Figure 2: Global Omics Landscapes, Trajectory PCAs, Volcano Distributions, and Coordinated Protein Correlation Networks
## (Nature Communications Manuscript)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Figure_2.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_2-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 2** provides the comprehensive multi-omics baseline characterizing the global transcriptomic and proteomic remodeling across both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis. It demonstrates three foundational biological insights:

1. **Coordinated Block Expression:** Transcriptomic and proteomic abundance shifts follow tightly coordinated modular patterns rather than independent fluctuations.
2. **Clear Disease Trajectory Separation:** Principal Component Analyses (PCA) and Volcano distributions confirm massive, reproducible molecular divergence between control and fibrotic disease states.
3. **Universal Correlation Networks ("Best Friends" Principle):** Proteome-wide correlation analyses reveal that essentially every detected protein possesses at least one highly correlated partner ($\rho_{\text{BP}} > 0.95$), demonstrating that proteins operate within densely interconnected co-regulatory networks.

---

## 2. Panel Breakdown & Generating Scripts

### 🔹 Panels A & B: Global Expression Landscapes (RNA & Protein Heatmaps)
* **Biological Meaning:** Hierarchically clustered expression heatmaps of all 943 shared genes and proteins across all mice in both disease cohorts (Controls, BDL, BDL + ASBTi, $\text{CCl}_4$ 2M, 6M, 12M).
* **Generating Scripts:**
  * `Global_RNA_Landscape_Analysis.R` (Transcriptomic landscape)
  * `Global_Protein_Landscape_Analysis.R` (Proteomic landscape)
* **Input Dataset:** `DTccl4_DT_LCPM_Gene_Protein_full.RData` / `DTccl4_DT_LCPM_BatchCorrected.RData`.

---

### 🔹 Panels C & D: Isolated PCA Trajectories & Volcano Plots (BDL & $\text{CCl}_4$)
* **Biological Meaning:**
  * **PCA Scores (Panels C1, D1):** 2D projection displaying clear separation between healthy controls and disease progression along PC1 and PC2.
  * **Volcano Plots (Panels C2, D2):** Statistical distribution of $\log_2$ Fold-Change versus $-\log_{10}(p\text{-value})$, highlighting significant differentially expressed genes and proteins.
* **Generating Scripts:**
  * `Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R` (BDL PCA & Volcano generation)
  * `Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R` ($\text{CCl}_4$ PCA & Volcano generation)
* **Input Datasets:** `Data_count_filtered_asbt_Gene_Protein_full.RData`, `Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData`.

---

### 🔹 Panel E: Shared-Order Correlation Heatmaps Triptych (BDL $\leftrightarrow$ $\text{CCl}_4$ $\leftrightarrow$ $\Delta r$)
* **Biological Meaning:** Pairwise Pearson correlation matrices ($943 \times 943$) computed in BDL and projected onto $\text{CCl}_4$ using identical hierarchical clustering order, alongside the difference matrix $\Delta r = r_{\text{BDL}} - r_{\text{CCl}_4}$. This demonstrates preserved modular correlation structures across distinct fibrosis etiologies.
* **Generating Scripts:**
  * `jan_correlation_heatmaps_shared_order.R` (Proteome correlation triptych)
  * `jan_correlation_heatmaps_mrna_shared_order.R` (Transcriptome correlation triptych)
* **Output Component:** `Correlation_Heatmaps_Triptych_BDL_CCl4_DeltaR.pdf`.

---

### 🔹 Panel F: Protein Correlation Density ("Best Friends" vs. Random Background)
* **Biological Meaning:** Density distribution of the single maximum correlation coefficient ($\text{Top-1 Partner}$, $\rho_{\text{BP}} > 0.95$) for each protein compared to the broad background distribution of random pairwise correlations centered around 0.
* **Generating Scripts:**
  * `Protein_Correlation_Density_BDL.R` (BDL cohort density analysis, Page 1)
  * `Protein_Correlation_Density_CCL4.R` ($\text{CCl}_4$ cohort density analysis, Page 1)
* **Output Components:** `Protein_Correlation_Density_BDL.pdf`, `Protein_Correlation_Density_CCL4.pdf`.

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_2_Master_Plate.R`](Generate_Figure_2_Master_Plate.R)
* **Execution:** Combines all individual vector and rasterized grob components using `cowplot` and `pdftools` into the publication-ready multi-panel vector file:
  ```bash
  Rscript Generate_Figure_2_Master_Plate.R
  ```
* **Output File:** `Figure_2.pdf` (4.29 MB).

---

## 4. Directory Layout

```text
Figure_02/
│
├── 📄 Figure_2.pdf                                      # Master Publication Plate (4.29 MB)
├── 📄 Generate_Figure_2_Master_Plate.R                  # Master assembly script (cowplot)
│
├── 📄 Global_RNA_Landscape_Analysis.R                  # Panel A: Global RNA expression heatmap
├── 📄 Global_Protein_Landscape_Analysis.R              # Panel B: Global Protein expression heatmap
├── 📄 Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R # Panel C: BDL PCA and Volcano plots
├── 📄 Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R# Panel D: CCl4 PCA and Volcano plots
├── 📄 jan_correlation_heatmaps_shared_order.R          # Panel E: Protein correlation triptych
├── 📄 jan_correlation_heatmaps_mrna_shared_order.R     # Panel E (supp): RNA correlation triptych
├── 📄 Protein_Correlation_Density_BDL.R                # Panel F: BDL Best Friends density (Page 1)
├── 📄 Protein_Correlation_Density_CCL4.R               # Panel F: CCl4 Best Friends density (Page 1)
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Figure Legend

> **Figure 2: Global multi-omics landscapes, disease trajectories, and coordinated protein correlation networks across cholestatic and toxic liver fibrosis.**  
> **(A, B)** Hierarchically clustered expression heatmaps of 943 shared transcripts **(A)** and proteins **(B)** across BDL ($n = 18$) and $\text{CCl}_4$ ($n = 36$) murine cohorts.  
> **(C, D)** Principal component analyses (PCA) and volcano plots ($\log_2\text{FC}$ vs $-\log_{10}p$) illustrating disease separation in BDL **(C)** and $\text{CCl}_4$ **(D)**.  
> **(E)** Proteome-wide pairwise Pearson correlation matrices for BDL (left) and $\text{CCl}_4$ (middle) displayed in identical hierarchical order, with the pairwise difference matrix ($\Delta r = r_{\text{BDL}} - r_{\text{CCl}_4}$, right).  
> **(F)** Density distributions of maximum pairwise correlation coefficients (Top-1 "Best Friend" partner, red) versus random pairwise background correlations (gray) across BDL and $\text{CCl}_4$.
