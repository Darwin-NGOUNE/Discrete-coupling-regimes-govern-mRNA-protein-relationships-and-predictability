# Procedure 5: Intra-Cohort Leave-One-Out Cross-Validation (LOOCV)
## (Intra-BDL & Intra-$\text{CCl}_4$ Machine Learning Modeling)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Validation](https://img.shields.io/badge/Validation-LOOCV-orange.svg)](https://en.wikipedia.org/wiki/Cross-validation_(statistics))
[![Nature Portfolio](https://img.shields.io/badge/Target-Nature_Communications-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Procedure 5** assesses the internal predictability of protein abundance using **Leave-One-Out Cross-Validation (LOOCV)** within each disease cohort separately:
* **Intra-BDL Cohort:** $n = 18$ mice (Sham, BDL, BDL + ASBTi).
* **Intra-$\text{CCl}_4$ Cohort:** $n = 36$ mice (Oil controls and $\text{CCl}_4$ across 2, 6, 12 months).

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
While the target response variables belong to the specific DiPa quadrants above, the feature selection mechanisms in **RF + Lasso**, **RF + RF**, and **Mastery 50** search across the **entire proteome of 943 gene--protein pairs** to identify optimal predictive covariates.

---

## 3. Cohort Definitions: Full Cohort vs. Subset Cohort (`sub`)

### 🔹 Full Cohort (All Animals):
* **Full BDL ($n = 18$):** Sham $n=6$, BDL $n=6$, BDL + ASBTi $n=6$.
* **Full $\text{CCl}_4$ ($n = 36$):** Oil controls and $\text{CCl}_4$ across 2, 6, 12 months.

### 🔹 Subset Cohort (`sub`, $n = 12$ mice each):
* **BDL Subset ($n = 12$):** Filtered for `Treatment %in% c("BDL", "control")` (Sham control $n=6$, untreated BDL $n=6$, excluding ASBTi).
* **$\text{CCl}_4$ Subset ($n = 12$):** Filtered for `(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)` (12-month chronic $\text{CCl}_4$ fibrosis $n=6$, baseline oil controls $n=6$).

---

## 4. Evaluated Machine Learning Architectures

### 🔹 1. Mastery 50 Model (Direct Penalized Lasso)
* **Location:** `Mastery50/`
* Uses the top 50 "Master Hub" proteins (`Top50_Mastery_Hubs_Protein_*.RData`) selected by proteome-wide covariance coverage.
* Fits penalized Lasso regression via `glmnet` ($\alpha = 1$, `prediction.method = "loo"`).

### 🔹 2. Two-Stage Protein Model (RF + Lasso)
* **Location:** `RF_Models/`
* **Stage 1:** Unpruned Random Forest ranks all 943 candidate proteins by Permutation Feature Importance.
* **Stage 2:** Cross-validated Lasso regression fits the top pre-selected covariates + cognate mRNA.

### 🔹 3. Two-Stage Non-Linear Model (RF + RF) -- Conserved Subset
* **Location:** `RF_RF_Subset/`
* **Stage 1:** Random Forest feature importance pre-selection.
* **Stage 2:** Non-linear Random Forest regression for final abundance prediction evaluated on the **conserved subset cohort** (`_sub`).

---

## 5. Directory Layout & Script Nomenclature

```text
Procedure_5_Intra_Cohort/
│
├── 📄 Function_protein_modelling_full_new_version_CD.R  # Modeling Engine
├── 📄 Function_design_newdata.R                         # Design Construction
├── 📄 Function_model_measures.R                         # Evaluation Metrics
├── 📄 Top50_Mastery_Hubs_Protein_BDL_Raw.RData          # BDL Hubs (Raw)
├── 📄 Top50_Mastery_Hubs_Protein_CCL4_Raw.RData         # CCl4 Hubs (Raw)
│
├── 📁 Mastery50/                                        # Direct Lasso with Top-50 Master Hubs
│   ├── 📄 Models_mastery50_cluster0.R                  # BDL: DiPa 8 (202 pairs, Full cohort n=18)
│   ├── 📄 Models_mastery50_cluster0sub.R               # BDL: DiPa 8 (202 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_cluster12.R                 # BDL: DiPa 1-2 (118 pairs, Full cohort n=18)
│   ├── 📄 Models_mastery50_cluster12sub.R              # BDL: DiPa 1-2 (118 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_cluster34.R                 # BDL: DiPa 3-4 (67 pairs, Full cohort n=18)
│   ├── 📄 Models_mastery50_cluster34sub.R              # BDL: DiPa 3-4 (67 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_cluster56.R                 # BDL: DiPa 5-6 (23 pairs, Full cohort n=18)
│   ├── 📄 Models_mastery50_cluster56sub.R              # BDL: DiPa 5-6 (23 pairs, Subset n=12)
│   │
│   ├── 📄 Models_mastery50_DTccl4_cluster0.R           # CCl4: DiPa 8 (202 pairs, Full cohort n=36)
│   ├── 📄 Models_mastery50_DTccl4_cluster0sub.R        # CCl4: DiPa 8 (202 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_DTccl4_cluster12.R          # CCl4: DiPa 1-2 (118 pairs, Full cohort n=36)
│   ├── 📄 Models_mastery50_DTccl4_cluster12sub.R       # CCl4: DiPa 1-2 (118 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_DTccl4_cluster34.R          # CCl4: DiPa 3-4 (67 pairs, Full cohort n=36)
│   ├── 📄 Models_mastery50_DTccl4_cluster34sub.R       # CCl4: DiPa 3-4 (67 pairs, Subset n=12)
│   ├── 📄 Models_mastery50_DTccl4_cluster56.R          # CCl4: DiPa 5-6 (23 pairs, Full cohort n=36)
│   └── 📄 Models_mastery50_DTccl4_cluster56sub.R       # CCl4: DiPa 5-6 (23 pairs, Subset n=12)
│
├── 📁 RF_Models/                                       # Two-Stage RF + Lasso
│   ├── 📄 Models_rf_lasso_prime_cluster*.R             # BDL RF + Lasso scripts
│   └── 📄 Models_rf_lasso_DTccl4_prime_*.R             # CCl4 RF + Lasso scripts
│
├── 📁 RF_RF_Subset/                                    # Two-Stage RF + RF (Conserved Subset n=12)
│   ├── 📄 Models_rf_rf_prime_cluster0sub.R             # BDL: DiPa 8 (202 pairs, Subset n=12)
│   ├── 📄 Models_rf_rf_prime_cluster12sub.R            # BDL: DiPa 1-2 (118 pairs, Subset n=12)
│   ├── 📄 Models_rf_rf_prime_cluster34sub.R            # BDL: DiPa 3-4 (67 pairs, Subset n=12)
│   ├── 📄 Models_rf_rf_prime_cluster56sub.R            # BDL: DiPa 5-6 (23 pairs, Subset n=12)
│   │
│   ├── 📄 Models_rf_rf_DTccl4_prime_cluster0sub.R      # CCl4: DiPa 8 (202 pairs, Subset n=12)
│   ├── 📄 Models_rf_rf_DTccl4_prime_cluster12sub.R     # CCl4: DiPa 1-2 (118 pairs, Subset n=12)
│   ├── 📄 Models_rf_rf_DTccl4_prime_cluster34sub.R     # CCl4: DiPa 3-4 (67 pairs, Subset n=12)
│   └── 📄 Models_rf_rf_DTccl4_prime_cluster56sub.R     # CCl4: DiPa 5-6 (23 pairs, Subset n=12)
│
└── 📄 README.md                                        # This documentation file
```
