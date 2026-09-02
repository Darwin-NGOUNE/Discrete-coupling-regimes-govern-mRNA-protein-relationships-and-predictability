# Stage 01: Data Harmonization & Batch Effect Correction (ComBat)
## (Liver Fibrosis Multi-Omics Consortium)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioc-sva_%26_biomaRt-green.svg)](https://bioconductor.org/packages/sva/)

---

## 1. Overview & Master Datasets

This stage establishes the foundational dataset harmonization and batch effect correction pipeline for the **Consortium Liver Fibrosis Project**. It bridges the two independent disease cohorts into the **two master datasets** used across all subsequent downstream analyses (Procedures 1, 3, 5, and Figures 1--6):

1. **Master Dataset 1 (Uncorrected / Raw Merged):** `DTccl4_DT_LCPM_Gene_Protein_full.RData`
2. **Master Dataset 2 (ComBat Batch-Corrected Merged):** `DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData`

---

## 2. Directory Structure & Key Files

```text
01_Data_Harmonization_and_Batch_Correction/
│
├── 📁 raw_data/                               # Initial cohort datasets
│   ├── Data_LCPM_asbt_Gene_Protein_full_with_DiPa.RData      (BDL cohort, n = 18 mice)
│   └── Data_CCl4_Gene_Protein_full_final_with_KMclustering.RData (CCl4 cohort, n = 36 mice)
│
├── 📄 01_harmonize_and_merge_datasets.R        # Exact script: DTccl4_DT_Gene_Protein_full_final.R
│                                              # Maps UniProt/Ensembl IDs via biomaRt to isolate 943 shared proteins
│
├── 📄 02_combat_batch_correction.R             # Exact script: Batch_Effect_Correction_3.R
│                                              # Applies ComBat (sva) preserving biological disease groups (~ DiseaseGroup)
│
├── 📁 processed_data/                         # The essential master datasets
│   ├── DT_LCPM_filtered_final_Gene_Protein_full.RData        (BDL 943 features, raw)
│   ├── DTccl4_filtered_final_Gene_Protein_full.RData        (CCl4 943 features, raw)
│   ├── DTccl4_DT_LCPM_Gene_Protein_full.RData               (MASTER 1: Merged Raw)
│   ├── DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData (BDL 943 features, ComBat)
│   ├── DTccl4_filtered_final_Gene_Protein_full_Batch.RData  (CCl4 943 features, ComBat)
│   └── DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData                  (MASTER 2: Merged ComBat Corrected)
│
└── 📄 README.md                               # This documentation file
```

---

## 3. Workflow Steps

### Step 1: Identifier Harmonization & Dataset Fusion
* **Script:** `01_harmonize_and_merge_datasets.R` (`DTccl4_DT_Gene_Protein_full_final.R`)
* **Operation:** Matches Ensembl gene IDs to external gene names and UniProt accessions, aligns sample identifiers across BDL ($n = 18$) and $\text{CCl}_4$ ($n = 36$), and extracts the **943 core mRNA--protein pairs**.
* **Output:** `processed_data/DTccl4_DT_LCPM_Gene_Protein_full.RData`.

### Step 2: Biological-Preserving ComBat Batch Correction
* **Script:** `02_combat_batch_correction.R` (`Batch_Effect_Correction_3.R`)
* **Operation:** Constructs full expression matrices and applies Empirical Bayes ComBat adjustment from the `sva` package using `mod <- model.matrix(~ DiseaseGroup)` to protect the true biological disease trajectory (`Control`, `BDL`, `BDL_ASBTi`, `CCl4`, `Oil`) while eliminating inter-batch platform shifts.
* **Output:** `processed_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData`.


## Final Combined Dataset

The final combined dataset is provided as an Excel file:

`DTccl4_DT_LCPM_Raw_BatchCorrected.xlsx`

This file contains the combined datasets from `BDL` and `CCL4`, including both the raw and batch-corrected data.