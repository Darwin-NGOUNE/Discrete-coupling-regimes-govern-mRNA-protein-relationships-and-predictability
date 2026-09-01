# Figure 5: Cross-Cohort Machine Learning Predictability & 75th-Percentile Prototypical Pair Validations (Procedure 1)
## *(Nature Communications Manuscript)*

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF_300DPI-red.svg)](Figure_5.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_5-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 5** establishes the fundamental predictability of protein abundance from transcriptomics under independent cross-disease validation (**Procedure 1: Independent Testing**). Models trained exclusively on cholestatic injury (BDL) are evaluated strictly on toxic chemical injury ($\text{CCl}_4$), and reciprocally ($\text{CCl}_4 \rightarrow \text{BDL}$).

Key insights:
1. **Regime-Governed Predictability (Panels A & B):** Highly concordant pairs (DiPa 1 & 2) yield superior cross-cohort prediction performance ($\rho_{\text{BP}} \ge 0.5$ in $\sim 70\%$ of pairs, $\rho_{\text{BP}} \ge 0.8$ in $\sim 28\%$). In contrast, discordant or shifted regimes (DiPa 3 & 4, DiPa 5 & 6, DiPa 8) fail to generalize across etiologies.
2. **Prototypical Pair Trajectory Validations (Panels C & D):** Representative gene–protein pairs at the 75th performance percentile illustrate exact quantitative agreement between predicted and observed protein levels across disease severity stages (Control, Intermediate, Disease).

---

## 2. Component Panels & Generating Scripts

### 🔹 Panels A & B (Top Row): Baseline Model Cross-Cohort Pearson Correlations
* **Panel A:** Direction 1 — **$\mathbf{\text{Train BDL, test } \text{CCl}_4}$**
* **Panel B:** Direction 2 — **$\mathbf{\text{Train } \text{CCl}_4 \text{, test BDL}}$**
* **Visual Specifications:**
  * Clean boxplots overlaid with individual protein points.
  * Dashed benchmark thresholds at $\rho_{\text{BP}} = 0.5$ (red) and $\rho_{\text{BP}} = 0.8$ (green).
  * Separation grid segments stopping strictly at $y = -1.00$.
  * Summary percentage annotations ($\ge 0.8$, $\ge 0.5$) clearly positioned below $-1.00$.
* **Generating Script:** [`generate_baseline_pearson_boxplots_Proc1.R`](generate_baseline_pearson_boxplots_Proc1.R)
* **Intermediate Vector Files:**
  * `Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf`
  * `Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf`

---

### 🔹 Panel C (Middle Row): 75th-Percentile Prototypical Scatters (Subset Animals)
* **Title:** **$\mathbf{\text{Baseline model: train BDL, test } \text{CCl}_4}$**
* **Layout:** 4 panels corresponding to DiPa groups (DiPa 1 & 2, DiPa 3 & 4, DiPa 5 & 6, DiPa 8).
* **Metrics:** Embedded box displaying $\boldsymbol{\rho_{\text{BP}}}$, $\mathbf{\text{RMSE}}$, $\mathbf{\text{NRMSE}}$, and $\mathbf{R^2}$ for both group median and the individual prototypical pair.
* **Generating Script:** [`generate_75th_percentile_scatterplots_Proc1.R`](generate_75th_percentile_scatterplots_Proc1.R)
* **Intermediate Vector File:** `Proc1_3Pages_1x4_Scatterplot_75thPercentile_Richtung1_Train_BDL_Test_CCl4.pdf` (Page 1)

---

### 🔹 Panel D (Bottom Row): 75th-Percentile Prototypical Scatters (All Animals)
* **Title:** **$\mathbf{\text{Baseline model: train BDL, test } \text{CCl}_4 \text{ (all animals)}}$**
* **Layout:** 4 panels across all animals in the full cohort, confirming stability of high-percentile pairs across all biological replicates.
* **Generating Script:** [`generate_full_dataset_75th_percentile_scatterplots_Proc1.R`](generate_full_dataset_75th_percentile_scatterplots_Proc1.R)
* **Intermediate Vector File:** `Proc1_Full_3Pages_1x4_Scatterplot_75thPercentile_BDL_CCL4.pdf` (Page 1)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_5_Master_Plate.R`](Generate_Figure_5_Master_Plate.R)
* **Execution:** Combines all 4 components (Panels A, B, C, D) into a unified 3-tier master plate at 300 DPI:
  ```bash
  Rscript Generate_Figure_5_Master_Plate.R
  ```
* **Output Publication File:** [`Figure_5.pdf`](Figure_5.pdf) (High-resolution vector PDF).

---

## 4. Directory Layout

```text
Figure_05/
│
├── 📄 Figure_5.pdf                                                              # Master Publication Plate (300 DPI Vector PDF)
├── 📄 Generate_Figure_5_Master_Plate.R                                          # Master cowplot assembly script
│
├── 📄 Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf           # Panel A: Direction 1 Boxplots
├── 📄 Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf           # Panel B: Direction 2 Boxplots
├── 📄 Proc1_3Pages_1x4_Scatterplot_75thPercentile_Richtung1_Train_BDL_Test_CCl4.pdf # Panel C: 75th-pct Scatters (Subset)
├── 📄 Proc1_Full_3Pages_1x4_Scatterplot_75thPercentile_BDL_CCL4.pdf            # Panel D: 75th-pct Scatters (All animals)
│
├── 📄 generate_baseline_pearson_boxplots_Proc1.R                                # Generates Panels A & B
├── 📄 generate_75th_percentile_scatterplots_Proc1.R                             # Generates Panel C
├── 📄 generate_full_dataset_75th_percentile_scatterplots_Proc1.R                # Generates Panel D
│
└── 📄 README.md                                                                 # Comprehensive documentation
```

---

## 5. Nature Communications Figure Legend

> **Figure 5: Discrete mRNA–protein coupling regimes govern cross-cohort machine learning predictability.**  
> **(A, B)** Cross-cohort evaluation of baseline machine learning models predicting protein abundance from transcriptomics under Procedure 1 (independent cross-disease testing). Models trained on cholestatic fibrosis (BDL) and tested on toxic liver injury ($\text{CCl}_4$, **A**), or trained on $\text{CCl}_4$ and tested on BDL (**B**). Boxplots display Bravais–Pearson correlation coefficients ($\rho_{\text{BP}}$) stratified by DiPa co-regulation groups (DiPa 1 & 2, DiPa 3 & 4, DiPa 5 & 6, DiPa 8). Dashed lines indicate moderate ($\rho_{\text{BP}} = 0.5$, red) and high ($\rho_{\text{BP}} = 0.8$, green) predictive thresholds; bottom values indicate percentages of pairs achieving these criteria.  
> **(C, D)** Representative prototypical gene–protein pairs selected at the 75th performance percentile across DiPa groups for subset animals (**C**, $\text{Baseline model: train BDL, test } \text{CCl}_4$) and all animals (**D**, $\text{Baseline model: train BDL, test } \text{CCl}_4 \text{ (all animals)}$). Scatter points represent individual murine samples colored by disease progression stage (Control, Intermediate, Disease). Inset annotations indicate median group metrics and single-pair validation statistics ($\rho_{\text{BP}}$, RMSE, NRMSE, $R^2$).
