# Supplementary Figure 3: Multi-Model Machine Learning Benchmarking across Procedure 5 Intra-Cohort Cross-Validation
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_3.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_3-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 3** delivers the comprehensive multi-model performance benchmarking across all four machine learning architectures evaluated under **Procedure 5 (Intra-Cohort Cross-Validation)** in a 2-Tier vertical presentation:
1. **Intra-BDL Cohort (Panel A, $\mathbf{\text{Intra BDL}}$, $n = 18$ mice):** Leave-One-Out Cross-Validation within the BDL cohort across all four DiPa quadrants.
2. **Intra-$\text{CCl}_4$ Cohort (Panel B, $\mathbf{\text{Intra } \text{CCl}_4}$, $n = 36$ mice):** Leave-One-Out Cross-Validation within the $\text{CCl}_4$ cohort across all four DiPa quadrants.

### The 4 Evaluated Machine Learning Architectures:
* **1. Baseline Model ($B_j$):** Cognate mRNA predictor only ($\hat{P}_{E,ji} = \beta_{0j} + \beta_{1j} G_{Eji}$).
* **2. Mastery 50 Model:** Direct linear penalized Lasso regression using top-50 correlated protein covariates.
* **3. RF + LASSO:** Two-stage modeling using Random Forest feature importance pre-selection followed by sparse Lasso regression.
* **4. RF + RF:** Two-stage non-linear ensemble modeling using Random Forest for both feature pre-selection and final regression.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Tier 1, Top): Procedure 5 4-Models Benchmark ($\text{Intra BDL}$)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing Baseline, Mastery 50, RF + LASSO, and RF + RF within the BDL cohort ($n = 18$ mice). Vertical divider lines terminate strictly at $-1.00$ with bold $\mathbf{\rho_{\text{BP}}}$.
* **Component File:** `Procedure_5_4Models_Pearson_Intra_BDL.pdf`
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc5.R`](generate_4models_pearson_boxplots_Proc5.R)

---

### 🔹 Panel B (Tier 2, Bottom): Procedure 5 4-Models Benchmark ($\text{Intra } \text{CCl}_4$)
* **Content:** Pearson correlation boxplots across DiPa quadrants comparing the 4 architectures within the $\text{CCl}_4$ cohort ($n = 36$ mice). Vertical divider lines terminate strictly at $-1.00$ with bold $\mathbf{\rho_{\text{BP}}}$.
* **Component File:** `Procedure_5_4Models_Pearson_Intra_CCl4.pdf`
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc5.R`](generate_4models_pearson_boxplots_Proc5.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_3_Master_Plate.R`](Generate_Supplementary_Figure_3_Master_Plate.R)
* **Execution:** Arranges Tier 1 (Panel A) and Tier 2 (Panel B) in a 2-tier vertical layout at 300 DPI using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_3_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_3.pdf`.

---

## 4. Directory Layout

```text
Supplementary_Figure_03/
│
├── 📄 Supplementary_Figure_3.pdf                        # Master Publication Plate (2-Tier Vertical Stack)
├── 📄 Generate_Supplementary_Figure_3_Master_Plate.R    # Master assembly script (cowplot)
│
├── 📄 Procedure_5_4Models_Pearson_Intra_BDL.pdf        # Panel A (Proc 5 Intra BDL 4-Models)
├── 📄 Procedure_5_4Models_Pearson_Intra_CCl4.pdf       # Panel B (Proc 5 Intra CCl4 4-Models)
│
├── 📄 generate_4models_pearson_boxplots_Proc5.R        # Generates Panels A & B
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 3: Comprehensive multi-model benchmarking across intra-cohort cross-validation procedures.**  
> **(A, B)** Test prediction accuracy (Bravais-Pearson correlation $\mathbf{\rho_{\text{BP}}}$) across DiPa quadrants (DiPa 1--2, DiPa 3--4, DiPa 5--6, and DiPa 8) under Procedure 5 comparing 4 machine learning architectures: Baseline ($B_j$), Mastery 50, RF + LASSO, and RF + RF within BDL ($n = 18$, **A**) and $\text{CCl}_4$ ($n = 36$, **B**). Percentages indicate the proportion of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) correlation sitting below the bounded $-1.00$ baseline.
