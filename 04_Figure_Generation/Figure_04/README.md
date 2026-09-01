# Figure 4: Cross-Cohort DiPa Group Overlap, Statistical Conservation, and Centroid Slope Distributions
## *(Nature Communications Manuscript)*

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF_300DPI-red.svg)](Figure_4.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_4-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 4** investigates the cross-disease conservation and quantitative trajectory scaling of transcript–protein coupling across cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis. 

It addresses two central mechanistic questions:
1. **Statistical Conservation of Regulatory Regimes (Panel A):** To what extent do specific gene–protein pairs maintain identical co-regulation behaviors (DiPa groups 1 to 8) when transitioning from cholestasis to toxic necroinflammation?
2. **Quantitative Scaling of Centroid Slopes (Panels B & C):** Do the directional transition vectors connecting control and diseased centroids (*Schwerpunkte*) exhibit identical slope distributions ($\Delta\text{Protein} / \Delta\text{RNA} \approx 1$) across concordant groups in both etiologies?

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Top Row): DiPa Group Cross-Cohort Overlap & Fisher Exact Significance
* **Biological Content:** Mini-Venn diagrams for each of the 8 DiPa groups comparing gene–protein assignments between BDL ($D$, yellow) and $\text{CCl}_4$ ($T$, blue) against the background universe of all core pairs ($A = 943$).
* **Quantitative Metrics:**
  * Number of overlapping pairs ($O$) and category-specific counts ($D, T$).
  * Overlap Index ($D_i = \frac{O \cdot A}{D \cdot T}$), where $D_i > 1$ indicates significant enrichment.
  * Fisher's Exact Test $p$-values evaluating whether overlap exceeds random expectation.
* **Generating Script:** [`DiPa_Overlap_Analysis_BDL_CCL4.R`](DiPa_Overlap_Analysis_BDL_CCL4.R)
* **Intermediate Vector File:** `DiPa_Overlap_Analysis_Results.pdf`

---

### 🔹 Panel B (Bottom Left): Centroid Slope Distributions in BDL
* **Biological Content:** Empirical probability density distributions of centroid trajectory slopes ($m = \Delta\text{Protein} / \Delta\text{RNA}$) for concordant DiPa group 1 (Up-Up, blue, $N_1 = 146$) and DiPa group 2 (Down-Down, red, $N_2 = 137$) in BDL.
* **Key Finding:** Centroid slopes peak tightly around $m = 1.0$ (dashed bounds at $[0.5, 1.5]$), demonstrating proportional transcript-to-protein synthesis scaling during cholestatic progression.
* **Generating Script:** [`Schwerpunkt_Slope_Distribution_BDL.R`](Schwerpunkt_Slope_Distribution_BDL.R)
* **Intermediate Vector File:** `Centroid_Slope_Frequency_BDL.pdf`

---

### 🔹 Panel C (Bottom Right): Centroid Slope Distributions in $\text{CCl}_4$
* **Biological Content:** Empirical probability density distributions of centroid trajectory slopes for DiPa group 1 (blue, $N_1 = 95$) and DiPa group 2 (red, $N_2 = 76$) in the chronic $\text{CCl}_4$ cohort.
* **Key Finding:** The density profiles closely mirror the BDL distribution, validating that linear stoichiometric scaling ($m \approx 1$) is an etiology-independent property of concordant genes.
* **Generating Script:** [`Schwerpunkt_Slope_Distribution_CCL4.R`](Schwerpunkt_Slope_Distribution_CCL4.R)
* **Intermediate Vector File:** `Centroid_Slope_Frequency_CCL4.pdf`

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_4_Master_Plate.R`](Generate_Figure_4_Master_Plate.R)
* **Execution:** Combines the top row (Panel A: Overlap Analysis) and bottom row (Panels B & C: Centroid Slope Densities) at high resolution (300 DPI) using `cowplot` and `magick`:
  ```bash
  Rscript Generate_Figure_4_Master_Plate.R
  ```
* **Output Publication File:** [`Figure_4.pdf`](Figure_4.pdf) (High-resolution vector plate).

---

## 4. Directory Layout

```text
Figure_04/
│
├── 📄 Figure_4.pdf                                      # Master Publication Plate (300 DPI Vector PDF)
├── 📄 Generate_Figure_4_Master_Plate.R                  # Master cowplot assembly script
│
├── 📄 DiPa_Overlap_Analysis_Results.pdf                # Panel A: DiPa Overlap Venn Diagrams
├── 📄 Centroid_Slope_Frequency_BDL.pdf                 # Panel B: BDL Centroid Slope Densities
├── 📄 Centroid_Slope_Frequency_CCL4.pdf                # Panel C: CCl4 Centroid Slope Densities
│
├── 📄 DiPa_Overlap_Analysis_BDL_CCL4.R                 # Generates Panel A
├── 📄 Schwerpunkt_Slope_Distribution_BDL.R             # Generates Panel B
├── 📄 Schwerpunkt_Slope_Distribution_CCL4.R            # Generates Panel C
│
└── 📄 README.md                                         # Comprehensive documentation
```

---

## 5. Nature Communications Figure Legend

> **Figure 4: Cross-cohort DiPa group overlap, statistical conservation, and linear centroid slope dynamics across liver injury etiologies.**  
> **(A)** Cross-cohort overlap analysis comparing DiPa group assignments of 943 shared mRNA–protein pairs between cholestatic (BDL, $D$, yellow) and toxic ($\text{CCl}_4$, $T$, blue) murine models. For each DiPa group (groups 1–8), inner counts display group-specific and shared pairs ($O$), alongside the total universe ($A = 943$). Bottom annotations report the Overlap Index ($D_i = (O \cdot A)/(D \cdot T)$) and Fisher’s exact test $P$-values ($^{**}P < 0.01, ^{***}P < 0.001$, ns: not significant).  
> **(B, C)** Empirical probability density distributions of centroid trajectory slopes ($m = \Delta\text{Protein}/\Delta\text{RNA}$) for concordant DiPa group 1 (Up-Up, blue) and DiPa group 2 (Down-Down, red) in BDL ($N_1 = 146, N_2 = 137$, **B**) and $\text{CCl}_4$ ($N_1 = 95, N_2 = 76$, **C**). Dashed vertical lines indicate the proportional stoichiometric bounds ($0.5 \le m \le 1.5$).
