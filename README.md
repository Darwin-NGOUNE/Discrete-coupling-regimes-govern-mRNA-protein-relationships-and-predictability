# Discrete-coupling-regimes-govern-mRNA-protein-relationships-and-predictability

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R: >= 4.1.0](https://img.shields.io/badge/R-%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Batch Correction: ComBat](https://img.shields.io/badge/Batch_Correction-ComBat-green.svg)](https://bioconductor.org/packages/release/bioc/html/sva.html)
[![DiPa Classification](https://img.shields.io/badge/Methodology-DiPa_Quadrants-purple.svg)](02_DiPa_Classification_and_Analysis/)
[![Nature Communications](https://img.shields.io/badge/Article-In_Submission-teal.svg)](https://www.nature.com/ncomms/)

---

## 📖 Scientific Overview

This repository provides the complete, self-contained, reproducible codebase and analysis pipelines for our study on **machine learning-driven proteome modeling across cholestatic and toxic liver fibrosis**.

By integrating parallel transcriptomic (RNA-Seq) and quantitative proteomic (LC-MS/MS) measurements from:
1. **Biliary Obstruction Model (BDL):** $n = 18$ mice (Sham, BDL, BDL + AS0369 ASBT inhibitor).
2. **Toxic Necrosis Model ($\text{CCl}_4$):** $n = 36$ mice (Mineral oil controls and $\text{CCl}_4$ across 2, 6, and 12 months).
3. **Conserved Harmonized Manifold:** **943** core matched gene--protein pairs analyzed across $N = 54$ animals.

We establish the **Differentiation Pattern (DiPa)** co-regulation taxonomy and benchmark four distinct machine learning architectures across **intra-cohort**, **cross-cohort blind transfer**, and **merged-cohort** evaluation frameworks.

---

## 🎯 Target Protein Stratification vs. Candidate Covariate Pool

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

## 🐁 Cohort Definitions: Full Cohort vs. Subset Cohorts

### 🔹 Full Cohorts (All Experimental Animals):
* **Full BDL Cohort ($n = 18$):** Sham ($n=6$), BDL ($n=6$), BDL + ASBTi ($n=6$).
* **Full $\text{CCl}_4$ Cohort ($n = 36$):** Mineral oil vehicle controls ($n=18$) and $\text{CCl}_4$ ($n=18$) across 2, 6, and 12 months.
* **Full Merged Cohort ($N = 54$):** Combined ComBat batch-corrected manifold of all 54 animals.

### 🔹 Subset Cohorts (`sub` in Proc 1 & 5, `_24` in Proc 3):
* **BDL Subset ($n = 12$):** Filtered for `Treatment %in% c("BDL", "control")` (Sham control $n=6$, BDL untreated $n=6$, excluding ASBTi).
* **$\text{CCl}_4$ Subset ($n = 12$):** Filtered for `(Treatment == "ccl4" & TreatmentTime == 12) | (Treatment == "oil" & TreatmentTime == 0)` (12-month chronic fibrotic mice $n=6$, baseline oil controls $n=6$).
* **Merged Subset (`_24`, $N = 24$):** Exactly combining the 12 BDL subset mice + 12 $\text{CCl}_4$ subset mice.

---

## 📂 Repository Architecture

```text
Protein_Modeling_Consortium_GitHub/
│
├── 📁 01_Data_Harmonization_and_Batch_Correction/ ──► Cross-cohort mapping & ComBat batch correction
├── 📁 02_DiPa_Classification_and_Analysis/        ──► Raw count 2D DiPa quadrant classification (thresh = 0.5)
│
├── 📁 03_Machine_Learning_Procedures/
│   ├── 📁 Procedure_1_Cross_Cohort/               ──► Cross-cohort blind testing (Train BDL <-> Test CCl4)
│   ├── 📁 Procedure_3_Merged_Cohort/              ──► Merged cohort cross-validation (N = 54 & N = 24 subset)
│   ├── 📁 Procedure_5_Intra_Cohort/               ──► Intra-cohort LOOCV (Intra-BDL & Intra-CCl4)
│   └── 📁 Shared_Modeling_Functions/              ──► Core modeling engine functions
│
└── 📁 04_Figure_Generation/
    ├── 📁 Figure_01/ to Figure_06/                ──► Main Figures (Vector PDFs, Master Scripts & READMEs)
    └── 📁 Supplementary_Figure_01/ to 05/         ──► Supplementary Figures (Vector PDFs, Master Scripts)
```

---

## 🛠️ Required R Packages & Installation

To install all required dependencies:

```r
install.packages(c(
  "data.table", "ggplot2", "cowplot", "magick", "pdftools",
  "glmnet", "randomForest", "mgcv", "stringr", "pbapply", "pheatmap"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("sva", "limma", "edgeR"))
```

---

## 📜 Citation & Author Information

**Author:**  
**Darwin Brunel Ngoune Domo**  
*Department of Statistics, TU Dortmund University, Dortmund, Germany*  

When using the code or datasets from this repository, please cite:
> Ngoune, D. et al. *Discrete-coupling-regimes-govern-mRNA-protein-relationships-and-predictability.* **Nature Communications** (2026).

---

## 📄 License
This repository is licensed under the [MIT License](LICENSE).
