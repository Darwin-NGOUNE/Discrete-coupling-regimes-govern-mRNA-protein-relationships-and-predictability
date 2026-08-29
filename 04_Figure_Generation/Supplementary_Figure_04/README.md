# Supplementary Figure 4: Cross-Cohort Multi-Model Benchmarking (Procedure 1 Direction 2: $\text{Train }\text{CCl}_4 \rightarrow \text{Test BDL}$)
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_4.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_4-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 4** presents the reciprocal cross-cohort multi-model benchmark under **Procedure 1 Direction 2**:
$$\text{Train on Toxic Liver Fibrosis } (\text{CCl}_4, n = 36) \longrightarrow \text{Test Blindly on Cholestatic Liver Fibrosis } (\text{BDL}, n = 18)$$

It complements Figure 6D (Direction 1: $\text{BDL} \rightarrow \text{CCl}_4$) by demonstrating that model transferability is strictly bidirectional and robust across four machine learning architectures:
1. **Baseline Model ($B_j$):** Cognate mRNA only.
2. **Mastery 50 Model:** Direct linear penalized Lasso regression with top-50 correlated protein covariates.
3. **Protein Model (RF + Lasso):** Two-stage modeling using Random Forest feature importance pre-selection followed by sparse Lasso regression.
4. **Protein Model (RF + RF):** Two-stage non-linear ensemble modeling using Random Forest for both feature pre-selection and final regression.

---

## 2. Component Panels & Generating Scripts

* **Content:** Bravais-Pearson correlation ($\rho_{\text{BP}}$) boxplots across DiPa quadrants (DiPa 1--2, DiPa 3--4, DiPa 5--6, DiPa 8, and All pairs) comparing all 4 machine learning models. Percentages indicate the proportion of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) accuracy.
* **Component / Output File:** `Supplementary_Figure_4.pdf` / `Procedure_1_4Models_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf` (102 KB).
* **Generating Script:** [`generate_4models_pearson_boxplots_Proc1.R`](generate_4models_pearson_boxplots_Proc1.R).
* **Input Dataset:** `Models_mastery50_CCL4_BDL_*.RData`, `Models_rf_lasso_full_testing_new_Batch_prime_*.RData` (from `Batch_corrected_prime_analysis/`).

---

## 3. Directory Layout

```text
Supplementary_Figure_04/
│
├── 📄 Supplementary_Figure_4.pdf                        # Master Publication Plate (102 KB)
├── 📄 Procedure_1_4Models_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf # Component Boxplot PDF
│
├── 📄 generate_4models_pearson_boxplots_Proc1.R        # Generating R script
│
└── 📄 README.md                                         # This documentation file
```

---

## 4. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 4: Cross-cohort multi-model benchmarking under Procedure 1 Direction 2 ($\text{Train }\text{CCl}_4 \rightarrow \text{Test BDL}$).**  
> Test prediction accuracy (Bravais-Pearson correlation $\rho_{\text{BP}}$) across DiPa quadrants under Procedure 1 Direction 2 comparing 4 machine learning architectures: Baseline ($B_j$), Mastery 50, Protein (RF + Lasso), and Protein (RF + RF). Percentages indicate the fraction of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) test performance when trained on toxic necrosis and tested on cholestatic biliary fibrosis.
