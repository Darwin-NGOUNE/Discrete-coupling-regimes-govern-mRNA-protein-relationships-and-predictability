# Figure 3: Differentiation Pattern (DiPa) 2D Coordinate Clouds & Co-Regulation Regimes
## (Nature Communications Manuscript)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Figure_3.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_3-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 3** visualizes the systematic mapping of all **943 core mRNA--protein pairs** within the two-dimensional **Differentiation Pattern (DiPa)** space across both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis. 

Instead of distributing continuously, transcript--protein relationships segregate into **discrete coupling quadrants**, illustrating:
1. **Synergistic Co-regulation (DiPa 1 & 2):** Concordant Up-Up and Down-Down regulation where changes in transcript abundance directly translate into proportional protein changes.
2. **Transcriptomic Buffering (DiPa 3--6):** Regimes where mRNA alterations occur without matching protein changes (or vice-versa), reflecting translational control or post-transcriptional buffering.
3. **Cross-Cohort Structural Consistency:** The geometric topology of the DiPa cloud and quadrant distributions remains strikingly preserved between acute/subacute biliary obstruction (BDL) and chronic toxic necrosis ($\text{CCl}_4$).

---

## 2. Component Panels & Generating Scripts

### 🔹 Left Panel: BDL DiPa Cloud & 8 Exemplary Quadrant Scatterplots
* **Content:** 2D scatter cloud of $\log_2(\text{Ratio}_G)$ vs $\log_2(\text{Ratio}_P)$ for BDL ($n = 18$ mice), color-coded by DiPa quadrants, accompanied by representative $1 \times 8$ trajectory scatterplots showing observed vs. predicted trends for archetypal genes in each quadrant.
* **Component File:** `Isolated_DiPa_Wolken_BDL.pdf` (100 KB)
* **Generating Script:** [`Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R)
* **Input Dataset:** `Data_count_filtered_asbt_Gene_Protein_full.RData`

---

### 🔹 Right Panel: $\text{CCl}_4$ DiPa Cloud & 8 Exemplary Quadrant Scatterplots
* **Content:** 2D scatter cloud of $\log_2(\text{Ratio}_G)$ vs $\log_2(\text{Ratio}_P)$ for $\text{CCl}_4$ ($n = 36$ mice), color-coded by DiPa quadrants, accompanied by representative $1 \times 8$ trajectory scatterplots for each quadrant.
* **Component File:** `Isolated_DiPa_Wolken_CCL4.pdf` (98 KB)
* **Generating Script:** [`Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R)
* **Input Dataset:** `Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData`

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_3_Master_Plate.R`](Generate_Figure_3_Master_Plate.R)
* **Execution:** Combines the left (BDL) and right ($\text{CCl}_4$) isolated DiPa clouds at 300 DPI into a widescreen vector plate using `cowplot`:
  ```bash
  Rscript Generate_Figure_3_Master_Plate.R
  ```
* **Output File:** `Figure_3.pdf` (2.32 MB).

---

## 4. Directory Layout

```text
Figure_03/
│
├── 📄 Figure_3.pdf                                      # Master Publication Plate (2.32 MB)
├── 📄 Generate_Figure_3_Master_Plate.R                  # Master assembly script (cowplot)
│
├── 📄 Isolated_DiPa_Wolken_BDL.pdf                     # Left Component: BDL DiPa Cloud + 8 Quadrants
├── 📄 Isolated_DiPa_Wolken_CCL4.pdf                    # Right Component: CCl4 DiPa Cloud + 8 Quadrants
│
├── 📄 Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R # Generates Isolated_DiPa_Wolken_BDL.pdf
├── 📄 Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R# Generates Isolated_DiPa_Wolken_CCL4.pdf
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Figure Legend

> **Figure 3: Differentiation Pattern (DiPa) classification reveals discrete mRNA--protein co-regulation regimes in cholestatic and toxic liver fibrosis.**  
> **(A, B)** Two-dimensional DiPa coordinate clouds mapping transcript fold-change ($\log_2 \text{Ratio}_G$, x-axis) versus protein fold-change ($\log_2 \text{Ratio}_P$, y-axis) for BDL ($n = 18$, left) and $\text{CCl}_4$ ($n = 36$, right) models across 943 shared pairs. Quadrants are delineated by $\pm 0.5$ thresholds ($\sim 1.41$-fold change). Surrounding multi-panel plots display representative measured abundance trajectories for prototypical genes belonging to each of the 8 DiPa quadrants.
