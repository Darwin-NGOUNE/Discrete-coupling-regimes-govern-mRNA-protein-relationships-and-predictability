# Figure 5: Procedure 1 Cross-Cohort Validation & 75th-Percentile Representative Predictions
## (Nature Communications Manuscript)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Figure_5.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_5-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 5** evaluates the **cross-cohort generalizability (Procedure 1)** of baseline mRNA-to-protein modeling when transferring predictive algorithms across fundamentally different disease mechanisms:

1. **Direction 1 ($\text{BDL} \rightarrow \text{CCl}_4$):** Train strictly on Cholestatic Liver Fibrosis ($n = 18$ mice) $\longrightarrow$ Blind test on Toxic Liver Fibrosis ($n = 36$ mice).
2. **Direction 2 ($\text{CCl}_4 \rightarrow \text{BDL}$):** Train strictly on Toxic Liver Fibrosis ($n = 36$ mice) $\longrightarrow$ Blind test on Cholestatic Liver Fibrosis ($n = 18$ mice).

It provides:
* **Quantified Transferability (Panels A & B):** Bravais-Pearson correlation ($\rho_{\text{BP}}$) boxplots across all DiPa quadrants, indicating the percentages of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) cross-cohort test accuracy.
* **Representative Trajectory Predictions (Panels C & D):** $1 \times 4$ scatterplots displaying measured versus model-predicted protein abundances at the $75^{\text{th}}$ performance percentile across DiPa quadrants 1--2, 3--4, 5--6, and 8 on both the **conserved subset** and the **full cohort**.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Top Left): Procedure 1 Baseline Accuracy (Direction 1: BDL $\rightarrow$ $\text{CCl}_4$)
* **Content:** Pearson correlation boxplots across DiPa quadrants under Procedure 1 Direction 1.
* **Component File:** `Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf` (30 KB)
* **Generating Script:** [`generate_baseline_pearson_boxplots_Proc1.R`](generate_baseline_pearson_boxplots_Proc1.R)

---

### 🔹 Panel B (Top Right): Procedure 1 Baseline Accuracy (Direction 2: $\text{CCl}_4$ $\rightarrow$ BDL)
* **Content:** Pearson correlation boxplots across DiPa quadrants under Procedure 1 Direction 2.
* **Component File:** `Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf` (30 KB)
* **Generating Script:** [`generate_baseline_pearson_boxplots_Proc1.R`](generate_baseline_pearson_boxplots_Proc1.R)

---

### 🔹 Panel C (Middle Row): 75th-Percentile Representative Predictions (Conserved Subset)
* **Content:** $1 \times 4$ multi-panel scatterplots displaying measured vs. predicted protein abundances for representative proteins at the $75^{\text{th}}$ percentile of test accuracy across DiPa quadrants (DiPa 1--2, 3--4, 5--6, 8) evaluated on the strictly conserved subset.
* **Component File:** `Proc1_3Pages_1x4_Scatterplot_75thPercentile_Richtung1_Train_BDL_Test_CCl4.pdf` (35 KB)
* **Generating Script:** [`generate_75th_percentile_scatterplots_Proc1.R`](generate_75th_percentile_scatterplots_Proc1.R)

---

### 🔹 Panel D (Bottom Row): 75th-Percentile Representative Predictions (Full Cohort)
* **Content:** $1 \times 4$ multi-panel scatterplots displaying measured vs. predicted protein abundances at the $75^{\text{th}}$ percentile evaluated across the entire full cohort.
* **Component File:** `Proc1_Full_3Pages_1x4_Scatterplot_75thPercentile_BDL_CCL4.pdf` (51 KB)
* **Generating Script:** [`generate_full_dataset_75th_percentile_scatterplots_Proc1.R`](generate_full_dataset_75th_percentile_scatterplots_Proc1.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_5_Master_Plate.R`](Generate_Figure_5_Master_Plate.R)
* **Execution:** Stitches Tier 1 (Panels A & B), Tier 2 (Panel C), and Tier 3 (Panel D) at 300 DPI into a 3-tier master plate canvas using `cowplot`:
  ```bash
  Rscript Generate_Figure_5_Master_Plate.R
  ```
* **Output File:** `Figure_5.pdf` (1.60 MB).

---

## 4. Directory Layout

```text
Figure_05/
│
├── 📄 Figure_5.pdf                                      # Master Publication Plate (1.60 MB)
├── 📄 Generate_Figure_5_Master_Plate.R                  # Master assembly script (cowplot)
│
├── 📄 Procedure_1_Baseline_Pearson_Richtung1_Train_BDL_Test_CCl4.pdf # Panel A
├── 📄 Procedure_1_Baseline_Pearson_Richtung2_Train_CCl4_Test_BDL.pdf # Panel B
├── 📄 Proc1_3Pages_1x4_Scatterplot_75thPercentile_Richtung1_Train_BDL_Test_CCl4.pdf # Panel C
├── 📄 Proc1_Full_3Pages_1x4_Scatterplot_75thPercentile_BDL_CCL4.pdf                 # Panel D
│
├── 📄 generate_baseline_pearson_boxplots_Proc1.R        # Generates Panels A & B
├── 📄 generate_75th_percentile_scatterplots_Proc1.R     # Generates Panel C
├── 📄 generate_full_dataset_75th_percentile_scatterplots_Proc1.R # Generates Panel D
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Figure Legend

> **Figure 5: Cross-cohort validation and representative 75th-percentile predictions under Procedure 1.**  
> **(A, B)** Test prediction accuracy (Bravais-Pearson correlation $\rho_{\text{BP}}$) across DiPa quadrants under Procedure 1 Direction 1 ($\text{Train BDL} \rightarrow \text{Test }\text{CCl}_4$, **A**) and Direction 2 ($\text{Train }\text{CCl}_4 \rightarrow \text{Test BDL}$, **B**). Percentages indicate the fraction of proteins achieving high ($\rho_{\text{BP}} \ge 0.8$) and moderate ($\rho_{\text{BP}} \ge 0.5$) accuracy.  
> **(C, D)** Representative $1 \times 4$ scatterplots displaying measured versus model-predicted protein abundances at the $75^{\text{th}}$ performance percentile across DiPa quadrants under Procedure 1 subset **(C)** and full cohort **(D)** evaluations.
