# Supplementary Figure 2: Master Proteins & Central Proteomic Correlation Hubs
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_2.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_2-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 2** (Slide 2 of the manuscript package) identifies and characterizes the **"Master Proteins" (Central Hubs)** within the murine proteomic landscape across both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis.

A "Master Protein" is defined as a highly connected hub protein that maintains strong, widespread correlation with the vast majority of other detected proteins across the entire proteome ($\rho_{\text{BP}} > 0.85$ across dozens of partners), acting as a global molecular barometer for pathological disease progression.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Left): BDL Master Hub Identification (Page 4 Analysis)
* **Biological Meaning:** Identification and ranking of the top central hub proteins in the BDL cohort ($n = 18$ mice), illustrating their high average inter-correlation across the proteome.
* **Component File:** `Protein_Correlation_Density_BDL.pdf` (Page 4, 86 KB)
* **Generating Script:** [`Protein_Correlation_Density_BDL.R`](Protein_Correlation_Density_BDL.R)
* **Input Dataset:** `Data_count_filtered_asbt_Gene_Protein_full.RData` / `DT_LCPM_filtered_final_Gene_Protein_full_Batch.RData`.

---

### 🔹 Panel B (Right): $\text{CCl}_4$ Master Hub Identification (Page 4 Analysis)
* **Biological Meaning:** Identification and ranking of the top central hub proteins in the $\text{CCl}_4$ cohort ($n = 36$ mice), confirming that master hub architecture is conserved despite distinct injury mechanisms.
* **Component File:** `Protein_Correlation_Density_CCL4.pdf` (Page 4, 91 KB)
* **Generating Script:** [`Protein_Correlation_Density_CCL4.R`](Protein_Correlation_Density_CCL4.R)
* **Input Dataset:** `Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData` / `DTccl4_filtered_final_Gene_Protein_full_Batch.RData`.

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_2_Master_Plate.R`](Generate_Supplementary_Figure_2_Master_Plate.R)
* **Execution:** Extracts Page 4 from both PDF density analyses at 300 DPI and arranges them side-by-side using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_2_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_2.pdf` (502 KB).

---

## 4. Directory Layout

```text
Supplementary_Figure_02/
│
├── 📄 Supplementary_Figure_2.pdf                        # Master Publication Plate (502 KB)
├── 📄 Generate_Supplementary_Figure_2_Master_Plate.R    # Master assembly script (cowplot)
│
├── 📄 Protein_Correlation_Density_BDL.pdf               # Full 4-page density analysis for BDL
├── 📄 Protein_Correlation_Density_CCL4.pdf              # Full 4-page density analysis for CCl4
│
├── 📄 Protein_Correlation_Density_BDL.R                 # Generates BDL density PDF (including Page 4)
├── 📄 Protein_Correlation_Density_CCL4.R                # Generates CCl4 density PDF (including Page 4)
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 2: Identification of Master Proteins and central proteomic hubs across liver injury models.**  
> **(A, B)** Ranking and correlation profiles of central "Master Proteins" identified in BDL ($n = 18$, **A**) and $\text{CCl}_4$ ($n = 36$, **B**). Master proteins maintain widespread high Pearson correlation coefficients ($\rho_{\text{BP}} > 0.85$) with large subsets of the proteome, serving as representative hub biomarkers of global fibrotic remodeling.
