# Supplementary Figure 5: Proportions of Conserved DiPa Pairs Achieving High Predictive Accuracy
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_5.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_5-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 5** provides a granular quantitative evaluation of the **Conserved DiPa Pairs** across both cross-cohort prediction directions:
* **Direction 1:** $\text{Train on BDL } (n = 18) \longrightarrow \text{Test on }\text{CCl}_4 (n = 36)$.
* **Direction 2:** $\text{Train on }\text{CCl}_4 (n = 36) \longrightarrow \text{Test on BDL } (n = 18)$.

It quantifies the exact proportion and absolute count of conserved mRNA--protein pairs that achieve:
1. **High Prediction Accuracy:** $\rho_{\text{BP}} \ge 0.8$ (strong predictive transferability).
2. **Moderate Prediction Accuracy:** $0.5 \le \rho_{\text{BP}} < 0.8$.
3. **Poor / Uncoupled Prediction:** $\rho_{\text{BP}} < 0.5$.

This breakdown demonstrates that conserved pairs belonging to synergistic quadrants (DiPa 1 & 2) overwhelmingly retain superior cross-disease predictability compared to buffered or discordant pairs.

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Tier 1, Top): Conserved Pairs Accuracy in Direction 1 ($\text{BDL} \rightarrow \text{CCl}_4$)
* **Content:** Stratified bar plots showing the proportion and count of conserved pairs reaching $\rho_{\text{BP}} \ge 0.8$ vs. $\rho_{\text{BP}} \ge 0.5$ across DiPa quadrants under Direction 1.
* **Component File:** `Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_Subset.pdf` (46 KB)
* **Generating Script:** [`jan_dipa_conserved_pairs_analysis.R`](jan_dipa_conserved_pairs_analysis.R)

---

### 🔹 Panel B (Tier 2, Bottom): Conserved Pairs Accuracy in Direction 2 ($\text{CCl}_4$ $\rightarrow$ BDL)
* **Content:** Stratified bar plots showing the proportion and count of conserved pairs reaching $\rho_{\text{BP}} \ge 0.8$ vs. $\rho_{\text{BP}} \ge 0.5$ across DiPa quadrants under Direction 2.
* **Component File:** `Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_Subset.pdf` (46 KB)
* **Generating Script:** [`jan_dipa_conserved_pairs_analysis.R`](jan_dipa_conserved_pairs_analysis.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_5_Master_Plate.R`](Generate_Supplementary_Figure_5_Master_Plate.R)
* **Execution:** Stitches Panel A (top) and Panel B (bottom) at 300 DPI into a vertical 2-tier plate using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_5_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_5.pdf` (643 KB).

---

## 4. Directory Layout

```text
Supplementary_Figure_05/
│
├── 📄 Supplementary_Figure_5.pdf                        # Master Publication Plate (643 KB)
├── 📄 Generate_Supplementary_Figure_5_Master_Plate.R    # Master assembly script (cowplot)
│
├── 📄 Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_Subset.pdf # Panel A
├── 📄 Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_Subset.pdf # Panel B
│
├── 📄 jan_dipa_conserved_pairs_analysis.R              # Generates Panels A & B
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Supplementary Figure Legend

> **Supplementary Figure 5: Proportions of conserved DiPa pairs achieving high cross-cohort prediction accuracy.**  
> **(A, B)** Bar graphs displaying the proportion and number of conserved transcript--protein pairs achieving high ($\rho_{\text{BP}} \ge 0.8$, dark blue) and moderate ($\rho_{\text{BP}} \ge 0.5$, light blue) cross-cohort prediction accuracy across DiPa quadrants for Direction 1 ($\text{BDL} \rightarrow \text{CCl}_4$, **A**) and Direction 2 ($\text{CCl}_4 \rightarrow \text{BDL}$, **B**).
