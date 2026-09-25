# Stage 01: Raw Data Preprocessing, Harmonization & Batch Effect Correction (ComBat)
## (Liver Fibrosis Multi-Omics Consortium)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioc-sva_%26_edgeR_%26_DESeq2_%26_MSstats-green.svg)](https://bioconductor.org/)

---

## 1. Overview & Pipeline Scope

This stage establishes the foundational, end-to-end multi-omics data processing pipeline for the **Consortium Liver Fibrosis Study**. It covers the entire trajectory from raw transcriptomic (RNA-Seq counts) and quantitative proteomic (LC-MS/MS) measurements down to the **two master datasets** used across all downstream modeling procedures (Procedures 1, 3, 5, and Figures 1--6):

1. **Master Dataset 1 (Uncorrected / Raw Merged):** `DTccl4_DT_LCPM_Gene_Protein_full.RData`
2. **Master Dataset 2 (ComBat Batch-Corrected Merged):** `DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData`

---

## 2. Directory Architecture & Key Files

```text
01_Data_Preprocessing_Harmonization_and_Batch_Correction/
│
├── 📁 preprocessing_BDL/                       # Upstream preprocessing for Biliary Obstruction cohort (N = 18)
│   ├── 📄 ASBT_liver_rna_BDL.R                # DESeq2 & edgeR TMM log-CPM normalization for RNA-Seq
│   ├── 📄 ASBT_liver_proteins_BDL.R           # Log2 transformation and sample ID alignment for SWATH proteomics
│   └── 📄 ASBT_liver_build_joint_data_BDL.R   # UniProt/Ensembl mapping to create joint BDL multi-omics dataset
│
├── 📁 preprocessing_CCl4/                      # Upstream preprocessing for Toxic Fibrosis cohort (N = 36)
│   ├── 📄 CCl4_gene_counts_normalized.R       # edgeR TMM log-CPM normalization for RNA-Seq
│   ├── 📄 CCl4_protein_data.R                 # MSstats processing of Skyline quantitative proteomics
│   └── 📄 CCl4_gene_protein_build_joint_data.R# UniProt mapping & complete-case filtering for CCl4 joint dataset
│
├── 📁 raw_data/                                # Cohort joint datasets prior to cross-cohort harmonization
│   ├── Data_LCPM_asbt_Gene_Protein_full_with_DiPa.RData          (BDL cohort, N = 18 mice)
│   └── Data_CCl4_Gene_Protein_full_final_with_KMclustering.RData (CCl4 cohort, N = 36 mice)
│
├── 📄 01_harmonize_and_merge_datasets.R        # Cross-cohort Ensembl/UniProt harmonization (943 core pairs)
├── 📄 02_combat_batch_correction.R             # Empirical Bayes ComBat correction preserving ~ DiseaseGroup
├── 📄 generate_supplementary_table_1_wide.R    # Journal-compliant wide Excel matrix generation
│
├── 📁 processed_data/                         # Final master datasets
│   ├── DT_LCPM_filtered_final_Gene_Protein_full.RData        (BDL 943 features, raw)
│   ├── DTccl4_filtered_final_Gene_Protein_full.RData        (CCl4 943 features, raw)
│   ├── DTccl4_DT_LCPM_Gene_Protein_full.RData               (MASTER 1: Merged Raw)
│   ├── DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData (BDL 943 features, ComBat)
│   ├── DTccl4_filtered_final_Gene_Protein_full_Batch.RData  (CCl4 943 features, ComBat)
│   └── DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData         (MASTER 2: Merged ComBat Corrected)
│
├── 📊 Supplementary_Table_1_RNA_Protein_Pairs.xlsx           # OFFICIAL PUBLICATION SUPPLEMENTARY TABLE:
│                                                            # Wide matrix (54 mice x 3,775 columns, 0 NAs)
├── 📊 DTccl4_DT_LCPM_Raw_BatchCorrected.xlsx                 # Comprehensive long-format dataset (50,922 rows x 21 cols)
│
└── 📄 README.md                               # This documentation file
```

---

## 3. Detailed Workflow Steps

### Phase 0: Individual Cohort Preprocessing

#### A. BDL Cohort Preprocessing (`preprocessing_BDL/`)
* **RNA-Seq Normalization (`ASBT_liver_rna_BDL.R`):**  
  Applies `DESeq2` filtering (`rowSums(counts) > 10`) followed by `edgeR` Trimmed Mean of M-values (**TMM**) normalization and $\log_2$-Counts Per Million ($\log\text{CPM}$, `prior.count = 3`) to yield normalized transcript abundances across the 18 experimental mice (Sham $n=6$, BDL $n=6$, BDL + ASBTi $n=6$).
* **Proteomics Processing (`ASBT_liver_proteins_BDL.R`):**  
  Extracts curated protein spectral intensities, applies $\log_2$-transformation, parses UniProt accessions, and harmonizes sample identifiers to match the RNA dataset.
* **Joint Dataset Construction (`ASBT_liver_build_joint_data_BDL.R`):**  
  Integrates parallel RNA and protein expression matrices via `mapUniProt` (mapping UniProtKB accessions to Ensembl gene identifiers) into `Data_LCPM_asbt_Gene_Protein_full_with_DiPa.RData`.

#### B. $\text{CCl}_4$ Cohort Preprocessing (`preprocessing_CCl4/`)
* **RNA-Seq Normalization (`CCl4_gene_counts_normalized.R`):**  
  Performs `edgeR` TMM normalization and $\log\text{CPM}$ transformation on the raw gene count matrix across 36 mice spanning mineral oil controls and $\text{CCl}_4$ treatments across 0, 2, 6, and 12 months.
* **Proteomics Processing (`CCl4_protein_data.R`):**  
  Processes Skyline export data through `MSstats` (`dataProcess` with median equalization) to compute sample-level protein log-intensities.
* **Joint Dataset Construction (`CCl4_gene_protein_build_joint_data.R`):**  
  Aligns transcriptomic and proteomic profiles across the time-course and retains gene--protein pairs with complete observation profiles across treatment conditions, outputting `Data_CCl4_Gene_Protein_full_final_with_KMclustering.RData`.

---

### Phase 1: Cross-Cohort Harmonization & Mapping
* **Script:** `01_harmonize_and_merge_datasets.R`
* **Operation:** Matches Ensembl gene IDs to external gene names and UniProt accessions across BDL ($N = 18$) and $\text{CCl}_4$ ($N = 36$), isolating the **943 conserved mRNA--protein pairs** present in both cohorts.
* **Output:** `processed_data/DTccl4_DT_LCPM_Gene_Protein_full.RData` (Master Dataset 1: Merged Raw).

---

### Phase 2: Biological-Preserving ComBat Batch Effect Correction
* **Script:** `02_combat_batch_correction.R`
* **Operation:** Applies Empirical Bayes ComBat adjustment (`sva` package) using `mod <- model.matrix(~ DiseaseGroup)` to protect the true biological disease trajectories (`Control`, `BDL`, `BDL_ASBTi`, `CCl4`, `Oil`) while eliminating inter-batch platform shifts between the two independent experimental series.
* **Output:** `processed_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData` (Master Dataset 2: Merged ComBat Corrected).

---

### Phase 3: Generation of Journal Supplementary Tables
* **Script:** `generate_supplementary_table_1_wide.R`
* **Operation:** Reshapes the 50,922 multi-omics observations into the journal-compliant wide data matrix.
* **Structure of `Supplementary_Table_1_RNA_Protein_Pairs.xlsx`:**
  - **Rows (54 mice):** Numbered 1 to 54, with `Mouse number`, `Intervention group` (`BDL`, `CCl4`), and `Treatment` (Sham, BDL, BDL+ASBTi, 0/2/6/12 months oil/CCl4).
  - **Columns (3,775 columns):** 3 sample metadata columns followed by 4 quantitative columns per cognate gene-protein pair across all 943 pairs:
    1. `Protein [GeneName], without batch correction`
    2. `Protein [GeneName], with batch correction`
    3. `mRNA [GeneName], without batch correction`
    4. `mRNA [GeneName], with batch correction`
  - **Data Integrity:** Exactly 0 missing values (0 NAs) across all 203,850 quantitative cells.

---

## 4. Summary of Master Datasets

| Dataset File | Description | Cohort | Dimensions |
| :--- | :--- | :--- | :--- |
| `processed_data/DT_LCPM_filtered_final_Gene_Protein_full.RData` | Raw Uncorrected | BDL ($n=18$) | 943 pairs $\times$ 18 mice |
| `processed_data/DTccl4_filtered_final_Gene_Protein_full.RData` | Raw Uncorrected | $\text{CCl}_4$ ($n=36$) | 943 pairs $\times$ 36 mice |
| `processed_data/DTccl4_DT_LCPM_Gene_Protein_full.RData` | **Master 1 (Merged Raw)** | Both ($N=54$) | 50,922 observations |
| `processed_data/DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData` | ComBat Corrected | BDL ($n=18$) | 943 pairs $\times$ 18 mice |
| `processed_data/DTccl4_filtered_final_Gene_Protein_full_Batch.RData` | ComBat Corrected | $\text{CCl}_4$ ($n=36$) | 943 pairs $\times$ 36 mice |
| `processed_data/DTccl4_DT_LCPM_BatchCorrected_3_Groups.RData` | **Master 2 (Merged ComBat)** | Both ($N=54$) | 50,922 observations |
| `Supplementary_Table_1_RNA_Protein_Pairs.xlsx` | Journal Wide Matrix | Both ($N=54$) | 54 rows $\times$ 3,775 cols (0 NAs) |
| `DTccl4_DT_LCPM_Raw_BatchCorrected.xlsx` | Long Data Table | Both ($N=54$) | 50,922 rows $\times$ 21 cols |