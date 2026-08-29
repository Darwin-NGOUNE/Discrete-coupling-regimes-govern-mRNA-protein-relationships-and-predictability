# Supplementary Figure 3: Multi-Model Machine Learning Benchmarking across Procedure 3 & Procedure 5
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_3.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_3-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 3** delivers the comprehensive multi-model performance benchmarking across all four machine learning architectures evaluated under:
1. **Procedure 3 (Merged Cohort Batch-Corrected Cross-Validation, Panel A):** Cross-validation across all 54 animals in the unified dataset.
2. **Procedure 5 (Intra-Cohort Cross-Validation, Panels B & C):** Leave-One-Out Cross-Validation (LOOCV) within individual cohorts:
   * **Panel B:** Intra-BDL cohort ($n = 18$ mice).
   * **Panel C:** Intra-$\text{CCl}_4$ cohort ($n = 36$ mice).

### The 4 Evaluated Machine Learning Architectures:
* **1. Baseline Model ($B_j$):** Cognate mRNA predictor only ($\hat{P}_{E,ji} = \beta_{0j} + \beta_{1j} G_{Eji}$).
* **2. Mastery 50 Model:** Direct linear penalized Lasso regression using top-50 correlated protein covariates.
* **3. Protein Model (RF + Lasso):** Two-stage modeling using Random Forest feature importance pre-selection followed by sparse Lasso regression.
* **4. Protein Model (RF + RF):** Two-stage non-linear ensemble modeling using Random Forest for both feature pre-selection and final regression.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Tier 1, Top): Procedure 3 4-Models Benchmark (Merged Batch)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing Baseline, Mastery 50, Protein (RF + Lasso), and Protein (RF + RF) on the merged batch-corrected cohort.
* **Component File:** `Procedure_3_4Models_Pearson_Merged_Batch.pdf` (102 KB)
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc3.R`](generate_4models_pearson_boxplots_Proc3.R)

---

### 🔹 Panel B (Tier 2, Middle): Procedure 5 4-Models Benchmark (Intra-BDL)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing the 4 architectures within the BDL cohort ($n = 18$ mice).
* **Component File:** `Procedure_5_4Models_Pearson_Intra_BDL.pdf` (221 KB)
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc5.R`](generate_4models_pearson_boxplots_Proc5.R)

---

### 🔹 Panel C (Tier 3, Bottom): Procedure 5 4-Models Benchmark (Intra-$\text{CCl}_4$)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing the 4 architectures within the $\text{CCl}_4$ cohort ($n = 36$ mice).
* **Component File:** `Procedure_5_4Models_Pearson_Intra_CCl4.pdf` (218 KB)
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc5.R`](generate_4models_pearson_boxplots_Proc5.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_3_Master_Plate.R`](Generate_Supplementary_Figure_3_Master_Plate.R)
* **Execution:** Arranges Tier 1 (Panel A), Tier 2 (Panel B), and Tier 3 (Panel C) in a 3-tier vertical layout at 300 DPI using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_3_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_3.pdf` (2.83 MB).

---

## 4. Directory Layout

```text
Supplementary_Figure_03/
│
├── 📄 Supplementary_Figure_3.pdf                        # Master Publication Plate (2.83 MB)
├── 📄 Generate_Supplementary_Figure_3_Master_Plate.R    # Master assembly script (cowplot)
│
├── 📄 Procedure_3_4Models_Pearson_Merged_Batch.pdf     # Panel A (Proc 3 Merged Batch 4-Models)
├── 📄 Procedure_5_4Models_Pearson_Intra_BDL.pdf        # Panel B (Proc 5 Intra BDL 4-Models)
├── 📄 Procedure_5_4Models_Pearson_Intra_CCl4.pdf       # Panel C (Proc 5 Intra CCl4 4-Models)
│
├── 📄 generate_4models_pearson_boxplots_Proc3.R        # Generates Panel A
├── 📄 generate_4models_pearson_boxplots_Proc5.R        # Generates Panels B & C
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 3: Comprehensive multi-model benchmarking across merged and intra-cohort evaluation procedures.**  
> **(A)** Test prediction accuracy (Bravais-Pearson correlation $\rho_{\text{BP}}$) across DiPa quadrants under Procedure 3 (merged batch-corrected cross-validation, $n = 54$) comparing 4 machine learning architectures: Baseline ($B_j$), Mastery 50, Protein (RF + Lasso), and Protein (RF + RF).  
> **(B, C)** Test prediction accuracy under Procedure 5 (intra-cohort cross-validation) comparing the 4 architectures within BDL ($n = 18$, **B**) and $\text{CCl}_4$ ($n = 36$, **C**).
