# Procedure 1: Cross-Cohort Generalizability & Bidirectional Testing
## (Train $\text{BDL} \leftrightarrow \text{Test }\text{CCl}_4$ under ComBat Batch Correction)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Validation](https://img.shields.io/badge/Validation-Cross--Cohort_Blind_Testing-success.svg)](https://en.wikipedia.org/wiki/Cross-validation_(statistics))
[![Nature Portfolio](https://img.shields.io/badge/Target-Nature_Communications-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Procedure 1** evaluates the cross-disease generalizability and predictive transferability of mRNA-to-protein models across distinct fibrosis etiologies:
* **Direction 1 ($\text{BDL} \rightarrow \text{CCl}_4$):** Train on Cholestatic BDL $\longrightarrow$ Blind external validation on Toxic $\text{CCl}_4$.
* **Direction 2 ($\text{CCl}_4 \rightarrow \text{BDL}$):** Train on Toxic $\text{CCl}_4$ $\longrightarrow$ Blind external validation on Cholestatic BDL.

---

## 2. Target Response Pairs vs. 943 Proteome Feature Selection Pool

### 🔹 Stratification of Target Proteins to Predict:
Models predict target protein abundance for specific gene--protein pairs stratified by their **DiPa co-regulation quadrants**:
* **DiPa 1--2:** **118** gene--protein pairs.
* **DiPa 3--4:** **67** gene--protein pairs.
* **DiPa 5--6:** **23** gene--protein pairs.
* **DiPa 8 / 0 :** **202** gene--protein pairs .
* **Total Conserved Matched Pairs:** $118 + 67 + 23 + 202 = \mathbf{410}$ pairs maintaining identical quadrant assignments across both diseases.

### 🔹 Candidate Covariate Search Space (Full 943 Proteome):
While the target response variables belong to the specific DiPa quadrants above, the feature selection mechanisms in **RF + Lasso**, **RF + RF**, and **Mastery 50** search across the **entire proteome of 943 gene--protein pairs** to identify optimal predictive covariates.

---

## 3. Cohort Definitions: Full Cohort vs. Subset Cohort (`sub`)

### 🔹 Full Cohort (All Animals):
* **Full BDL:** All $n = 18$ mice (Sham $n=6$, BDL $n=6$, BDL + ASBTi $n=6$).
* **Full $\text{CCl}_4$:** All $n = 36$ mice (Oil controls and $\text{CCl}_4$ across 2, 6, 12 months).

### 🔹 Subset Cohort (`sub`, $n = 12$ mice each):
* **BDL Subset ($n = 12$):** Filtered for `Treatment %in% c("BDL", "control")` (Sham control $n=6$, BDL untreated $n=6$, excluding ASBTi).
* **$\text{CCl}_4$ Subset ($n = 12$):** Filtered for `(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)` (12-month chronic $\text{CCl}_4$ fibrosis $n=6$, baseline oil controls $n=6$).

---

## 4. Evaluated Machine Learning Architectures

### 🔹 1. Mastery 50 Model (Direct Penalized Lasso)
* **Location:** `Mastery50/`
* Direct penalized Lasso regression (`glmnet`, $\alpha = 1$) utilizing the top 50 Master Hub proteins (`Top50_Mastery_Hubs_Protein_*_Batch.RData`) selected by proteome-wide covariance coverage.

### 🔹 2. Two-Stage Protein Model (RF + Lasso)
* **Location:** `RF_Lasso/`
* **Stage 1:** Unpruned Random Forest ranks all 943 candidate proteins by impurity_corrected Importance.
* **Stage 2:** Sparse Lasso regression fits the top pre-selected covariates from stage 1 + cognate mRNA.

### 🔹 3. Two-Stage Non-Linear Model (RF + RF) -- Conserved Subset
* **Location:** `RF_RF/`
* **Stage 1:** Random Forest feature importance pre-selection.
* **Stage 2:** Non-linear Random Forest regression for final abundance prediction evaluated on the **conserved subset cohort** (`_sub`).

---

## 5. Directory Layout & Script Nomenclature

```text
Procedure_1_Cross_Cohort/
│
├── 📄 Function_protein_modelling_testing_full_new.R     # Testing Engine
├── 📄 Function_protein_modelling_testing_full.R         # Testing Engine (Legacy/Full)
├── 📄 Function_design.R                                 # Design Construction
├── 📄 Function_design_newdata_filtered.R                # Filtered Newdata Design
├── 📄 Function_model_measures.R                         # Metrics (Pearson, RMSE, R2)
├── 📄 Top50_Mastery_Hubs_Protein_BDL_Batch.RData        # BDL Hubs (Batch)
├── 📄 Top50_Mastery_Hubs_Protein_CCL4_Batch.RData       # CCl4 Hubs (Batch)
│
├── 📁 Mastery50/                                        # Direct Lasso with Top-50 Master Hubs
│   ├── 📄 BDL_CCL4_Mastery50_0_prime.R                 # BDL -> CCl4: DiPa 8 (202 pairs, Full cohort)
│   ├── 📄 BDL_CCL4_Mastery50_0_prime_sub.R             # BDL -> CCl4: DiPa 8 (202 pairs, Subset 12 mice)
│   ├── 📄 BDL_CCL4_Mastery50_1_2_prime.R               # BDL -> CCl4: DiPa 1-2 (118 pairs, Full cohort)
│   ├── 📄 BDL_CCL4_Mastery50_1_2_prime_sub.R           # BDL -> CCl4: DiPa 1-2 (118 pairs, Subset)
│   ├── 📄 BDL_CCL4_Mastery50_3_4_prime.R               # BDL -> CCl4: DiPa 3-4 (67 pairs, Full cohort)
│   ├── 📄 BDL_CCL4_Mastery50_3_4_prime_sub.R           # BDL -> CCl4: DiPa 3-4 (67 pairs, Subset)
│   ├── 📄 BDL_CCL4_Mastery50_5_6_prime.R               # BDL -> CCl4: DiPa 5-6 (23 pairs, Full cohort)
│   ├── 📄 BDL_CCL4_Mastery50_5_6_prime_sub.R           # BDL -> CCl4: DiPa 5-6 (23 pairs, Subset)
│   │
│   ├── 📄 CCL4_BDL_Mastery50_0_prime.R                 # CCl4 -> BDL: DiPa 8 (202 pairs, Full cohort)
│   ├── 📄 CCL4_BDL_Mastery50_0_prime_sub.R             # CCl4 -> BDL: DiPa 8 (202 pairs, Subset)
│   ├── 📄 CCL4_BDL_Mastery50_1_2_prime.R               # CCl4 -> BDL: DiPa 1-2 (118 pairs, Full cohort)
│   ├── 📄 CCL4_BDL_Mastery50_1_2_prime_sub.R           # CCl4 -> BDL: DiPa 1-2 (118 pairs, Subset)
│   ├── 📄 CCL4_BDL_Mastery50_3_4_prime.R               # CCl4 -> BDL: DiPa 3-4 (67 pairs, Full cohort)
│   ├── 📄 CCL4_BDL_Mastery50_3_4_prime_sub.R           # CCl4 -> BDL: DiPa 3-4 (67 pairs, Subset)
│   ├── 📄 CCL4_BDL_Mastery50_5_6_prime.R               # CCl4 -> BDL: DiPa 5-6 (23 pairs, Full cohort)
│   └── 📄 CCL4_BDL_Mastery50_5_6_prime_sub.R           # CCl4 -> BDL: DiPa 5-6 (23 pairs, Subset)
│
├── 📁 RF_Lasso/                                         # Two-Stage RF + Lasso (16 scripts)
│   ├── 📄 BDL_CCL4_Batch_0_prime.R / _sub.R            # DiPa 8 (202 pairs)
│   ├── 📄 BDL_CCL4_Batch_1_2_prime.R / _sub.R          # DiPa 1-2 (118 pairs)
│   ├── 📄 BDL_CCL4_Batch_3_4_prime.R / _sub.R          # DiPa 3-4 (67 pairs)
│   ├── 📄 BDL_CCL4_Batch_5_6_prime.R / _sub.R          # DiPa 5-6 (23 pairs)
│   └── 📄 CCL4_BDL_Batch_*_prime*.R                    # CCl4 -> BDL scripts
│
├── 📁 RF_RF/                                            # Two-Stage RF + RF (Conserved Subset)
│   ├── 📄 RF_BDL_CCL4_Batch_*_prime_sub.R              # BDL -> CCl4 RF + RF scripts
│   └── 📄 RF_CCL4_BDL_Batch_*_prime_sub.R              # CCl4 -> BDL RF + RF scripts
│
└── 📄 README.md                                         # This documentation file
```
