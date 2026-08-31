# Figure 2: Global Omics Landscapes, Disease Trajectories, Volcano Distributions, and Coordinated Correlation Networks
## *(Nature Communications Manuscript)*

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF_300DPI-red.svg)](Figure_2.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_2-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 2** provides the comprehensive multi-omics baseline characterizing the global transcriptomic and proteomic remodeling across both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis cohorts. It establishes four foundational biological insights:

1. **Clear Disease Trajectory Separation (Panel A):** Principal Component Analyses (PCA) and Volcano distributions confirm massive, highly reproducible molecular divergence between control and fibrotic disease states across both mRNA and protein layers.
2. **Coordinated Block Expression Landscapes (Panel B):** Transcriptomic and proteomic abundance shifts follow tightly coordinated modular block patterns across individual animals rather than isolated fluctuations.
3. **Preserved Modular Correlation Networks (Panel C):** Cross-cohort pairwise correlation matrices ($943 \times 943$) reveal shared modular architectures between cholestatic and toxic fibrosis, with conserved co-regulatory gene and protein communities.
4. **Universal Nearest-Neighbor Density ("Best Friends" Principle, Panel D):** Proteome-wide correlation density distributions demonstrate that essentially every quantified protein possesses at least one highly correlated partner ($|\rho_{\text{BP}}| > 0.90$), establishing dense interconnected co-regulation.

---

## 2. Panel Breakdown & Generating Scripts

### 🔹 Panel A: PCA Trajectories & Volcano Distributions (BDL & $\text{CCl}_4$)
* **Biological Content:**
  * **PCA Projections (Top row):** 2D score plots along PC1 and PC2 for transcriptomic and proteomic landscapes in BDL ($N = 18$) and $\text{CCl}_4$ ($N = 31$).
  * **Volcano Distributions (Bottom row):** $\log_2(\text{fold change})$ versus $-\log_{10}(P\text{ value})$, quantifying significantly up- and down-regulated genes and proteins ($P < 0.05, |\log_2\text{FC}| > 1$).
* **Generating Scripts:**
  * [`Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R)
  * [`Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R)
* **Intermediate Vector Files:** `Isolated_PCA_BDL.pdf`, `Isolated_PCA_CCL4.pdf`, `Isolated_Volcano_BDL.pdf`, `Isolated_Volcano_CCL4.pdf`.

---

### 🔹 Panel B: Global Expression Landscapes (RNA & Protein Heatmaps)
* **Biological Content:**
  * Hierarchically clustered expression heatmaps of all 943 shared core transcripts (top row) and proteins (bottom row) across individual mice in BDL (left) and $\text{CCl}_4$ (right).
  * Clear sample-level annotation bars indicating individual animal IDs (`No. 1` to `No. 12`) and disease conditions.
* **Generating Scripts:**
  * [`Global_RNA_Landscape_Analysis.R`](Global_RNA_Landscape_Analysis.R) (Transcriptomic landscape heatmaps)
  * [`Global_Protein_Landscape_Analysis.R`](Global_Protein_Landscape_Analysis.R) (Proteomic landscape heatmaps)
* **Intermediate Vector Files:** `Global_RNA_Landscape_Summary.pdf` (Pages 2 & 3), `Global_Protein_Landscape_Summary.pdf` (Pages 2 & 3).

---

### 🔹 Panel C: Shared-Order Correlation Matrices Triptych (RNA & Protein)
* **Biological Content:**
  * Unified $2 \times 3$ matrix showing pairwise Pearson correlation matrices ($943 \times 943$) for mRNA (top) and Protein (bottom).
  * Columns: **BDL data**, **$\text{CCl}_4$ data** (projected onto identical hierarchical clustering order), and **Difference of BDL on $\text{CCl}_4$ data** ($\Delta \rho_{\text{BP}} = \rho_{\text{BDL}} - \rho_{\text{CCl}_4}$).
  * 3 dedicated colorbars centered beneath each column for $\rho_{\text{BP}}$ and $\Delta\rho_{\text{BP}}$.
* **Generating Scripts:**
  * [`Generate_Panel_C_Combined_Correlation_Heatmaps.R`](Generate_Panel_C_Combined_Correlation_Heatmaps.R) (Unified $2 \times 3$ generator)
  * [`jan_correlation_heatmaps_mrna_shared_order.R`](jan_correlation_heatmaps_mrna_shared_order.R)
  * [`jan_correlation_heatmaps_shared_order.R`](jan_correlation_heatmaps_shared_order.R)
* **Intermediate Vector File:** `Panel_C_Combined_Correlation_Heatmaps.pdf`.

---

### 🔹 Panel D: Protein Correlation Density Distributions (Top 1 to Top 100 Ranks)
* **Biological Content:**
  * Empirical probability density distributions of absolute Pearson correlation coefficients ($|\rho_{\text{BP}}|$) across hierarchical nearest-neighbor ranks ($\text{Top 1}, \text{Top 2}, \dots, \text{Top 100}$) for BDL (left) and $\text{CCl}_4$ (right).
  * A single, centered shared legend for correlation rank tiers.
* **Generating Scripts:**
  * [`Generate_Panel_D_Combined_Density_Plots.R`](Generate_Panel_D_Combined_Density_Plots.R) (Unified 1x2 generator with shared legend)
  * [`Protein_Correlation_Density_BDL.R`](Protein_Correlation_Density_BDL.R)
  * [`Protein_Correlation_Density_CCL4.R`](Protein_Correlation_Density_CCL4.R)
* **Intermediate Vector File:** `Panel_D_Combined_Density_Plots.pdf`.

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_2_Master_Plate.R`](Generate_Figure_2_Master_Plate.R)
* **Execution:** Combines all panels at high resolution (300 DPI) using `cowplot` and `magick`:
  ```bash
  Rscript Generate_Figure_2_Master_Plate.R
  ```
* **Final Publication Output:** [`Figure_2.pdf`](Figure_2.pdf)

---

## 4. Directory Layout

```text
Figure_02/
│
├── 📄 Figure_2.pdf                                      # Master Publication Plate (300 DPI)
├── 📄 Generate_Figure_2_Master_Plate.R                  # Master cowplot assembly script
│
├── 📄 Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R  # Panel A: BDL PCA & Volcano
├── 📄 Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R # Panel A: CCl4 PCA & Volcano
├── 📄 Global_RNA_Landscape_Analysis.R                   # Panel B: Global RNA heatmaps
├── 📄 Global_Protein_Landscape_Analysis.R               # Panel B: Global Protein heatmaps
├── 📄 Generate_Panel_C_Combined_Correlation_Heatmaps.R  # Panel C: Unified 2x3 Correlation Triptych
├── 📄 Generate_Panel_D_Combined_Density_Plots.R         # Panel D: Unified Correlation Density curves
├── 📄 jan_correlation_heatmaps_mrna_shared_order.R      # Panel C: RNA correlation analysis
├── 📄 jan_correlation_heatmaps_shared_order.R           # Panel C: Protein correlation analysis
├── 📄 Protein_Correlation_Density_BDL.R                 # Panel D: BDL Density & Centrality
├── 📄 Protein_Correlation_Density_CCL4.R                # Panel D: CCl4 Density & Centrality
│
└── 📄 README.md                                         # Comprehensive documentation
```

---

## 5. Nature Communications Figure Legend

> **Figure 2: Global transcriptomic and proteomic landscapes, disease trajectory segregation, and modular correlation networks in cholestatic and toxic liver injury.**  
> **(A)** Principal component analyses (top) and volcano distributions (bottom, $\log_2(\text{fold change})$ versus $-\log_{10}(P\text{ value})$) demonstrating robust molecular separation between control and diseased animals across transcriptomes and proteomes in BDL ($N = 18$) and $\text{CCl}_4$ ($N = 31$) cohorts. Red points denote significantly regulated features ($P < 0.05, |\log_2\text{FC}| > 1$).  
> **(B)** Hierarchically clustered expression heatmaps of all 943 shared core transcripts (top) and proteins (bottom) across individual animals in BDL (left) and $\text{CCl}_4$ (right). Column annotation bars display individual mouse identifiers (`No. 1`–`No. 12`) and disease states.  
> **(C)** Pairwise Pearson correlation matrices ($943 \times 943$) for mRNA (top) and Protein (bottom) across BDL and $\text{CCl}_4$, sorted according to consensus hierarchical clustering order, alongside differential correlation matrices ($\Delta \rho_{\text{BP}} = \rho_{\text{BDL}} - \rho_{\text{CCl}_4}$). Dedicated colorbars denote correlation coefficient scale ($-1.0$ to $+1.0$) and differential correlation scale ($-1.5$ to $+1.5$).  
> **(D)** Empirical probability density distributions of absolute pairwise Pearson correlation coefficients ($|\rho_{\text{BP}}|$) across hierarchical nearest-neighbor ranks ($\text{Top 1}$ to $\text{Top 100}$) for BDL (left) and $\text{CCl}_4$ (right), highlighting the universal high-correlation "Best Friends" principle in the proteome.
