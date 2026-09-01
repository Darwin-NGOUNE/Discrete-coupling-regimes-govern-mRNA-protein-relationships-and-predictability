# Supplementary Figure 5: Conservation of protein best-partner pairs ($\rho_{\text{BP}}$)
## (Nature Communications Supplementary Information)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Supplementary_Figure_5.pdf)
[![Supplementary Information](https://img.shields.io/badge/Article-Supp_Fig_5-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 5** evaluates the cross-cohort conservation of mRNA--protein best-partner abundance relationships ($\rho_{\text{BP}}$) across DiPa categories:
* **Panel A (Direction 1):** Best-partner pairs identified in BDL ($n = 12$) and evaluated for conservation in $\text{CCl}_4$ ($n = 12$).
* **Panel B (Direction 2):** Best-partner pairs identified in $\text{CCl}_4$ ($n = 12$) and evaluated for conservation in BDL ($n = 12$).

Conserved partner relationships are stratified into four rigorous conservation tiers:
1. **High:** $\rho_{\text{BP}} \ge 0.8$ (strong cross-disease conservation).
2. **Moderate:** $0.5 \le \rho_{\text{BP}} < 0.8$ (retained correlation).
3. **Low:** $0.2 \le \rho_{\text{BP}} < 0.5$ (weak correlation).
4. **Lost:** $\rho_{\text{BP}} < 0.2$ (decoupled / noise).

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Tier 1, Top): Best-partner BDL and conservation in $\text{CCl}_4$
* **Content:** Stacked bar plots depicting the proportion of conserved pairs across DiPa groups when trained on BDL and evaluated on $\text{CCl}_4$.
* **Component File:** `Proportion_Conserved_Pairs_per_DiPa_Group_Direction1_BDL_to_CCl4_Pearson_Subset.pdf`
* **Generating Script:** [`jan_dipa_conserved_pairs_analysis.R`](jan_dipa_conserved_pairs_analysis.R)

---

### 🔹 Panel B (Tier 2, Bottom): Best-partner $\text{CCl}_4$ and conservation in BDL
* **Content:** Stacked bar plots depicting the proportion of conserved pairs across DiPa groups when trained on $\text{CCl}_4$ and evaluated on BDL.
* **Component File:** `Proportion_Conserved_Pairs_per_DiPa_Group_Direction2_CCl4_to_BDL_Pearson_Subset.pdf`
* **Generating Script:** [`jan_dipa_conserved_pairs_analysis.R`](jan_dipa_conserved_pairs_analysis.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Supplementary_Figure_5_Master_Plate.R`](Generate_Supplementary_Figure_5_Master_Plate.R)
* **Execution:** Stitches Panel A (top) and Panel B (bottom) at 300 DPI into a vertical 2-tier plate using `cowplot`:
  ```bash
  Rscript Generate_Supplementary_Figure_5_Master_Plate.R
  ```
* **Output File:** `Supplementary_Figure_5.pdf`

---

## 4. Directory Layout

```text
Supplementary_Figure_05/
│
├── 📄 Supplementary_Figure_5.pdf                        # Master Publication Plate
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

> **Supplementary Figure 5: Conservation of protein best-partner pairs ($\rho_{\text{BP}}$).**  
> **(A, B)** Stacked bar charts depicting the proportion of conserved best-partner pairs across DiPa groups for **(A)** Best-partner BDL and conservation in $\text{CCl}_4$ and **(B)** Best-partner $\text{CCl}_4$ and conservation in BDL. Slices denote High ($\rho_{\text{BP}} \ge 0.8$), Moderate ($0.5 \le \rho_{\text{BP}} < 0.8$), Low ($0.2 \le \rho_{\text{BP}} < 0.5$), and Lost ($\rho_{\text{BP}} < 0.2$) conservation tiers.
