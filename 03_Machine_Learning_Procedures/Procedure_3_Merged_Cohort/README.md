# Procedure 3: Merged-Cohort Batch-Corrected Cross-Validation
## (Unified Consortium Dataset Modeling, $N = 54$ Mice & $N = 24$ Subset)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Validation](https://img.shields.io/badge/Validation-10--Fold_CV-orange.svg)](https://en.wikipedia.org/wiki/Cross-validation_(statistics))
[![Nature Portfolio](https://img.shields.io/badge/Target-Nature_Communications-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Procedure 3** evaluates machine learning model performance using **10-Folds Cross-Validation** on the **unified, batch-corrected consortium dataset** combining both cholestatic (BDL) and toxic ($\text{CCl}_4$) murine cohorts:
* **Full Merged Cohort ($N = 54$ mice):** All 18 BDL animals + all 36 $\text{CCl}_4$ animals.
* **Merged Subset Cohort (`_24`, $N = 24$ mice):** Exactly 12 BDL subset mice (`Treatment %in% c("BDL", "control")`) + 12 $\text{CCl}_4$ subset mice (`(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)`).

---

## 2. Target Response Pairs vs. 943 Proteome Feature Selection Pool

### 🔹 Stratification of Target Proteins to Predict:
Models predict target protein abundance for specific gene--protein pairs stratified by their **DiPa co-regulation quadrants**:
* **DiPa 1--2:** **118** gene--protein pairs.
* **DiPa 3--4:** **67** gene--protein pairs.
* **DiPa 5--6:** **23** gene--protein pairs.
* **DiPa 8 / 0:** **202** gene--protein pairs.
* **Total Conserved Matched Pairs:** $118 + 67 + 23 + 202 = \mathbf{410}$ pairs maintaining identical quadrant assignments across both diseases.

### 🔹 Candidate Covariate Search Space (Full 943 Proteome):
While the target response variables belong to the specific DiPa quadrants above, the feature selection mechanisms in **RF + Lasso**, **RF + RF**, and **Mastery 50** search across the **entire proteome of 943 gene--protein pairs** to select optimal predictive covariates.

---

## 3. Cohort Definitions: Full Merged Cohort vs. `_24` Subset Cohort

### 🔹 Full Cohort (All 54 Animals):
* Contains all experimental groups (Sham, BDL, BDL + ASBTi, Oil controls, and $\text{CCl}_4$ 2M, 6M, 12M).

### 🔹 Subset Cohort (`_24`, Exactly 24 Animals):
* **12 BDL Subset Mice:** `Treatment %in% c("BDL", "control")` (Sham $n=6$, untreated BDL $n=6$).
* **12 $\text{CCl}_4$ Subset Mice:** `(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)` (12-month chronic $\text{CCl}_4$ $n=6$, baseline oil controls $n=6$).

---

## 4. Evaluated Machine Learning Architectures

### 🔹 1. Mastery 50 Model (Direct Penalized Lasso)
* **Location:** `Mastery50/`
* Direct penalized Lasso regression (`glmnet`, $\alpha = 1$) utilizing the top 50 Master Hub proteins (`Top50_Mastery_Hubs_Protein_Merged_Batch.RData`).

### 🔹 2. Two-Stage Protein Model (RF + Lasso)
* **Location:** `RF_Lasso/`
* **Stage 1:** Unpruned Random Forest ranks all 943 candidate proteins by impurity_corrected Feature Importance.
* **Stage 2:** Sparse Lasso regression fits the top pre-selected covariates from stage 1 + cognate mRNA.

### 🔹 3. Two-Stage Non-Linear Model (RF + RF) -- Conserved Subset (`_24`)
* **Location:** `RF_RF/`
* **Stage 1:** Random Forest feature importance pre-selection.
* **Stage 2:** Non-linear Random Forest regression evaluated on the **conserved subset cohort (`_24`)**.

---

## 5. Directory Layout & Script Nomenclature

```text
Procedure_3_Merged_Cohort/
│
├── 📄 Function_protein_modelling_full_new_version_Meg.R # Merged Modeling Engine
├── 📄 Function_design_filtered.R                        # Design Construction
├── 📄 Function_model_measures.R                         # Evaluation Metrics
├── 📄 Top50_Mastery_Hubs_Protein_Merged_Batch.RData     # Merged Hubs (Batch)
│
├── 📁 Mastery50/                                        # Direct Lasso with Top-50 Master Hubs
│   ├── 📄 Models_mastery50_mergedata_batch_over_0.R     # DiPa 8 (202 pairs, Full cohort N=54)
│   ├── 📄 Models_mastery50_mergedata_batch_over_0_24.R  # DiPa 8 (202 pairs, Subset N=24)
│   ├── 📄 Models_mastery50_mergedata_batch_over_1_2.R   # DiPa 1-2 (118 pairs, Full cohort N=54)
│   ├── 📄 Models_mastery50_mergedata_batch_over_1_2_24.R# DiPa 1-2 (118 pairs, Subset N=24)
│   ├── 📄 Models_mastery50_mergedata_batch_over_3_4.R   # DiPa 3-4 (67 pairs, Full cohort N=54)
│   ├── 📄 Models_mastery50_mergedata_batch_over_3_4_24.R# DiPa 3-4 (67 pairs, Subset N=24)
│   ├── 📄 Models_mastery50_mergedata_batch_over_5_6.R   # DiPa 5-6 (23 pairs, Full cohort N=54)
│   └── 📄 Models_mastery50_mergedata_batch_over_5_6_24.R# DiPa 5-6 (23 pairs, Subset N=24)
│
├── 📁 RF_Lasso/                                         # Two-Stage RF + Lasso
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_0.R     # DiPa 8 (Full cohort N=54)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_0_24.R  # DiPa 8 (Subset N=24)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_1_2.R   # DiPa 1-2 (Full cohort N=54)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_1_2_24.R# DiPa 1-2 (Subset N=24)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_3_4.R   # DiPa 3-4 (Full cohort N=54)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_3_4_24.R# DiPa 3-4 (Subset N=24)
│   ├── 📄 Models_rf_lasso_mergedata_batch_prime_over_5_6.R   # DiPa 5-6 (Full cohort N=54)
│   └── 📄 Models_rf_lasso_mergedata_batch_prime_over_5_6_24.R# DiPa 5-6 (Subset N=24)
│
├── 📁 RF_RF/                                            # Two-Stage RF + RF (Conserved Subset N=24)
│   ├── 📄 Models_rf_rf_mergedata_batch_prime_over_0_24.R    # DiPa 8 (202 pairs, Subset N=24)
│   ├── 📄 Models_rf_rf_mergedata_batch_prime_over_1_2_24.R  # DiPa 1-2 (118 pairs, Subset N=24)
│   ├── 📄 Models_rf_rf_mergedata_batch_prime_over_3_4_24.R  # DiPa 3-4 (67 pairs, Subset N=24)
│   └── 📄 Models_rf_rf_mergedata_batch_prime_over_5_6_24.R  # DiPa 5-6 (23 pairs, Subset N=24)
│
└── 📄 README.md                                         # This documentation file
```
