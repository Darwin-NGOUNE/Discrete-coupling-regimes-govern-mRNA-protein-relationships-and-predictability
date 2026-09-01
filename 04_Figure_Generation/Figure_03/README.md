# Figure 3: Differentiation Pattern (DiPa) 2D Coordinate Clouds & Discrete Coupling Regimes
## *(Nature Communications Manuscript)*

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF_300DPI-red.svg)](Figure_3.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_3-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 3** presents the systematic two-dimensional mapping of all **943 matched mRNA–protein pairs** within the **Differentiation Pattern (DiPa)** manifold across both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis cohorts.

Rather than adhering to a continuous uniform distribution, transcript–protein dynamics segregate into **discrete coupling regimes**, illustrating four fundamental post-transcriptional phenomena:

1. **Concordant Synergistic Coupling (DiPa 1 & 2):** Directionally aligned Up-Up ($\text{DiPa 1}$) and Down-Down ($\text{DiPa 2}$) regulation, where transcriptional activation or repression translates directly into proportional protein changes.
2. **Transcriptomic Buffering / Post-Transcriptional Attenuation (DiPa 5 & 6):** Substantial mRNA shifts with minimal protein change ($\text{DiPa 5}$: mRNA Up / Protein unchanged; $\text{DiPa 6}$: mRNA Down / Protein unchanged), indicative of translational buffering or protein half-life stabilization.
3. **Protein-Dominant Remodeling (DiPa 3 & 4):** Prominent protein-level accumulation or depletion in the absence of significant transcriptomic changes ($\text{DiPa 3}$: Protein Up / mRNA unchanged; $\text{DiPa 4}$: Protein Down / mRNA unchanged), reflecting altered translation rates or selective proteasomal degradation.
4. **Discordant Opposing Coupling (DiPa 7) & Non-Responsive Homeostasis (DiPa 8):** Inversion regimes where transcript and protein levels move in opposite directions, and unperturbed baseline manifolds.
5. **Cross-Etiology Geometric Conservation:** The topological distribution across DiPa quadrants is strikingly preserved between acute cholestasis (BDL) and chronic chemical poisoning ($\text{CCl}_4$).

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A: BDL 2D DiPa Coordinate Cloud
* **Biological Content:** 2D scatter manifold of $\log_2(\text{fold change}_{\text{RNA}})$ ($x$-axis) versus $\log_2(\text{fold change}_{\text{Protein}})$ ($y$-axis) for the BDL cohort ($N = 18$).
* **Thresholds:** $\pm 0.5$ $\log_2\text{FC}$ boundaries dividing the plane into 8 discrete DiPa coupling quadrants (1: Concordant Up, 2: Concordant Down, 3: Prot Up, 4: Prot Down, 5: RNA Up, 6: RNA Down, 7: Discordant, 8: Center).
* **Generating Script:** [`Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R)
* **Intermediate Vector File:** `Isolated_DiPa_Wolken_BDL.pdf` (Page 1, Left side)

---

### 🔹 Panel B: BDL Prototypical Pair Trajectories (Pairs 1 to 8)
* **Biological Content:** Individual $2 \times 4$ multi-panel scatter plots showing sample-level mRNA vs. protein abundance for archetypal genes chosen from each of the 8 DiPa quadrants in BDL.
* **Features:** 
  * Individual mouse points colored by condition (Blue: Control, Red: BDL).
  * Trajectory vectors connecting Control and Disease group centroids (*Schwerpunkte*, marked by `+` crosses).
  * Bold Bravais-Pearson correlation statistics: overall ($\rho_{\text{BP}}$), control-specific ($\rho_{\text{C}}$), and disease-specific ($\rho_{\text{D}}$).
* **Generating Script:** [`Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R)

---

### 🔹 Panel C: $\text{CCl}_4$ 2D DiPa Coordinate Cloud
* **Biological Content:** 2D scatter manifold of $\log_2(\text{fold change}_{\text{RNA}})$ vs $\log_2(\text{fold change}_{\text{Protein}})$ for the $\text{CCl}_4$ cohort ($N = 31$, Month 12 $\text{CCl}_4$ vs. Month 0 Oil).
* **Generating Script:** [`Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R)
* **Intermediate Vector File:** `Isolated_DiPa_Wolken_CCL4.pdf` (Page 1, Left side)

---

### 🔹 Panel D: $\text{CCl}_4$ Prototypical Pair Trajectories (Pairs 1 to 8)
* **Biological Content:** Individual $2 \times 4$ multi-panel scatter plots showing sample-level mRNA vs. protein abundance for archetypal genes across all 8 DiPa quadrants in $\text{CCl}_4$.
* **Features:** Individual mouse points colored by condition (Blue: Oil Control, Red: $\text{CCl}_4$), centroid vectors, and stratified correlation metrics ($\rho_{\text{BP}}, \rho_{\text{C}}, \rho_{\text{D}}$).
* **Generating Script:** [`Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R`](Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R)

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_3_Master_Plate.R`](Generate_Figure_3_Master_Plate.R)
* **Execution:** Combines the left (BDL, Panels A & B) and right ($\text{CCl}_4$, Panels C & D) components at high resolution (300 DPI) using `cowplot` and `magick`:
  ```bash
  Rscript Generate_Figure_3_Master_Plate.R
  ```
* **Output File:** [`Figure_3.pdf`](Figure_3.pdf) (High-resolution vector graphic).

---

## 4. Directory Layout

```text
Figure_03/
│
├── 📄 Figure_3.pdf                                      # Master Publication Plate (300 DPI Vector PDF)
├── 📄 Generate_Figure_3_Master_Plate.R                  # Master cowplot assembly script
│
├── 📄 Isolated_DiPa_Wolken_BDL.pdf                     # Left Plate: BDL DiPa Cloud + 8 Quadrant Scatters
├── 📄 Isolated_DiPa_Wolken_CCL4.pdf                    # Right Plate: CCl4 DiPa Cloud + 8 Quadrant Scatters
│
├── 📄 Generate_A4_MASTER_BDL_Schwerpunkt_normal_data.R # Generates Isolated_DiPa_Wolken_BDL.pdf
├── 📄 Generate_A4_MASTER_CCL4_Schwerpunkt_normal_data.R# Generates Isolated_DiPa_Wolken_CCL4.pdf
│
└── 📄 README.md                                         # Comprehensive documentation
```

---

## 5. Nature Communications Figure Legend

> **Figure 3: Differentiation Pattern (DiPa) classification resolves discrete coupling and decoupling regimes between mRNA and protein dynamics in liver injury.**  
> **(A, C)** Two-dimensional DiPa coordinate manifolds mapping transcript fold-change ($\log_2(\text{fold change}_{\text{RNA}})$, $x$-axis) against protein fold-change ($\log_2(\text{fold change}_{\text{Protein}})$, $y$-axis) for BDL ($N = 18$, **A**) and $\text{CCl}_4$ ($N = 31$, **C**) cohorts across all 943 matched core pairs. Dashed threshold lines ($\pm 0.5$, corresponding to $\sim 1.41$-fold change) delineate the 8 discrete DiPa coupling quadrants: Concordant Up/Down (1 & 2), Protein-dominant (3 & 4), mRNA-dominant / Buffered (5 & 6), Discordant (7), and Central Non-regulated (8). Archetypal representative pairs from each quadrant are highlighted with black diamond markers.  
> **(B, D)** Measured sample-level abundance scatter plots for prototypical transcript–protein pairs (Pairs 1–8) representing each DiPa quadrant in BDL (**B**) and $\text{CCl}_4$ (**D**). Blue points indicate healthy control animals, red points denote fibrotic disease animals, and colored crosses (`+`) mark respective subgroup centroids (*Schwerpunkte*) linked by trajectory vectors. Bold annotations summarize total Bravais-Pearson correlation ($\rho_{\text{BP}}$), intra-control correlation ($\rho_{\text{C}}$), and intra-disease correlation ($\rho_{\text{D}}$).
