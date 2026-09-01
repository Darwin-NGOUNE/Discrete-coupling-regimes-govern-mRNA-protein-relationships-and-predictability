# Figure 6: Procedure 3 Merged-Cohort Validation, Multi-Model Benchmarking & 75th-Percentile Predictions
## (Nature Communications Manuscript)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Figure_6.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_6-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 6** represents the culminating machine learning synthesis of the manuscript. It evaluates model accuracy under **Procedure 3 (Merged Cohort Batch-Corrected Cross-Validation, $n = 54$ mice)** and presents a systematic cross-cohort **multi-architecture benchmark (Procedure 1, $\text{BDL} \rightarrow \text{CCl}_4$)**:

1. **Procedure 3 4-Models Comparison (Panel A, $\mathbf{\text{BDL} + \text{CCl}_4}$):** Evaluates Baseline ($B_j$), Mastery 50, RF + LASSO, and RF + RF on the unified batch-corrected dataset across all DiPa quadrants with vertical dividers bounded at $-1.00$ and bold $\mathbf{\rho_{\text{BP}}}$.
2. **Procedure 3 75th-Percentile Trajectory Scatters (Panels B & C):** $1 \times 4$ scatterplots displaying measured vs. predicted protein levels at the $75^{\text{th}}$ performance percentile on both the **conserved subset** ($\mathbf{\text{Protein model: BDL} + \text{CCl}_4}$) and the **full cohort** ($\mathbf{\text{Protein model: BDL} + \text{CCl}_4\text{ (all animals)}}$).
3. **Procedure 1 Multi-Model Benchmark (Panel D, $\mathbf{\text{Train BDL, test } \text{CCl}_4}$):** Side-by-side comparison of 4 machine learning architectures:
   * **Baseline ($B_j$):** Cognate mRNA only.
   * **Mastery 50:** Direct linear Lasso with top-50 correlated protein covariates.
   * **RF + LASSO:** Two-stage non-linear RF pre-selection + sparse Lasso regression.
   * **RF + RF:** Two-stage non-linear RF pre-selection + Random Forest regression.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Tier 1): Procedure 3 4-Models Performance Comparison ($\text{BDL} + \text{CCl}_4$)
* **Content:** Pearson correlation ($\rho_{\text{BP}}$) boxplots across DiPa quadrants comparing Baseline, Mastery 50, RF + LASSO, and RF + RF on the merged batch-corrected cohort.
* **Component File:** `Procedure_3_4Models_Pearson_Merged_Batch.pdf`
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc3.R`](generate_4models_pearson_boxplots_Proc3.R)

---

### 🔹 Panel B (Tier 2): Procedure 3 75th-Percentile Scatters ($\text{Protein model: BDL} + \text{CCl}_4$)
* **Content:** $1 \times 4$ multi-panel scatterplots at the $75^{\text{th}}$ accuracy percentile across DiPa quadrants 1--2, 3--4, 5--6, and 8 evaluated on the conserved subset (Control vs Disease).
* **Component File:** `Proc3_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf` (Page 3)
* **Generating Script:** [`generate_75th_percentile_scatterplots_Proc3.R`](generate_75th_percentile_scatterplots_Proc3.R)

---

### 🔹 Panel C (Tier 3): Procedure 3 75th-Percentile Scatters ($\text{Protein model: BDL} + \text{CCl}_4\text{ (all animals)}$)
* **Content:** $1 \times 4$ multi-panel scatterplots at the $75^{\text{th}}$ accuracy percentile across DiPa quadrants evaluated on the full merged cohort (Control, Intermediate, Disease).
* **Component File:** `Proc3_Full_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf` (Page 3)
* **Generating Script:** [`generate_full_dataset_75th_percentile_scatterplots_Proc3.R`](generate_full_dataset_75th_percentile_scatterplots_Proc3.R)

---

### 🔹 Panel D (Tier 4): Procedure 1 Multi-Model Benchmark ($\text{Train BDL, test } \text{CCl}_4$)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing Baseline, Mastery 50, RF + LASSO, and RF + RF under Procedure 1 ($\text{Train BDL} \rightarrow \text{Test }\text{CCl}_4$).
* **Component File:** `Procedure_1_4Models_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf`
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc1.R`](generate_4models_pearson_boxplots_Proc1.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_6_Master_Plate.R`](Generate_Figure_6_Master_Plate.R)
* **Execution:** Stitches Tier 1 (Panel A), Tier 2 (Panel B), Tier 3 (Panel C), and Tier 4 (Panel D) at 300 DPI into a 4-tier master plate canvas using `cowplot`:
  ```bash
  Rscript Generate_Figure_6_Master_Plate.R
  ```
* **Output File:** `Figure_6.pdf`.

---

## 4. Directory Layout

```text
Figure_06/
│
├── 📄 Figure_6.pdf                                      # Master Publication Plate (4-Tier)
├── 📄 Generate_Figure_6_Master_Plate.R                  # Master assembly script (cowplot)
│
├── 📄 Procedure_3_4Models_Pearson_Merged_Batch.pdf     # Panel A (4 Models Proc 3: BDL + CCl4)
├── 📄 Proc3_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf      # Panel B (Subset Scatters)
├── 📄 Proc3_Full_3Pages_1x4_Scatterplot_75thPercentile_Merged_Batch.pdf # Panel C (Full Scatters)
├── 📄 Procedure_1_4Models_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf     # Panel D (4 Models Proc 1)
│
├── 📄 generate_4models_pearson_boxplots_Proc3.R        # Generates Panel A
├── 📄 generate_75th_percentile_scatterplots_Proc3.R     # Generates Panel B
├── 📄 generate_full_dataset_75th_percentile_scatterplots_Proc3.R # Generates Panel C
├── 📄 generate_4models_pearson_boxplots_Proc1.R        # Generates Panel D
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Figure Legend

> **Figure 6: Merged cohort modeling performance, representative predictions, and multi-model cross-cohort benchmarking.**  
> **(A)** Test prediction accuracy (Bravais-Pearson correlation $\mathbf{\rho_{\text{BP}}}$) across DiPa quadrants (DiPa 1--2, DiPa 3--4, DiPa 5--6, and DiPa 8) under Procedure 3 ($\mathbf{\text{BDL} + \text{CCl}_4}$, $n = 54$) comparing Baseline, Mastery 50, RF + LASSO, and RF + RF architectures. Percentages indicate the proportion of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) accuracy below the $-1.00$ baseline.  
> **(B, C)** Representative $1 \times 4$ scatterplots displaying measured versus model-predicted protein abundances at the $75^{\text{th}}$ performance percentile across DiPa quadrants under Procedure 3 subset ($\mathbf{\text{Protein model: BDL} + \text{CCl}_4}$, **B**) and full cohort ($\mathbf{\text{Protein model: BDL} + \text{CCl}_4\text{ (all animals)}}$, **C**) evaluations. Inset metrics show non-bold group median and pair-specific Bravais-Pearson correlation $\rho_{\text{BP}}$, $\text{RMSE}$, $\text{nRMSE}$, and $R^2$.  
> **(D)** Cross-cohort benchmark under Procedure 1 ($\mathbf{\text{Train BDL, test } \text{CCl}_4}$) comparing 4 machine learning architectures: Baseline, Mastery 50, RF + LASSO, and RF + RF.
