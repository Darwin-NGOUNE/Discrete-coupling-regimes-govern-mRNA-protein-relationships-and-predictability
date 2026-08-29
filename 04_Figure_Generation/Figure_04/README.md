# Figure 4: Cross-Cohort DiPa Quadrant Overlap & Centroid Slope Dynamics
## (Nature Communications Manuscript)

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF-red.svg)](Figure_4.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Figure_4-teal.svg)](https://www.nature.com/ncomms/)

---

## 1. Scientific Overview & Objective

**Figure 4** investigates the cross-disease conservation and geometric trajectory dynamics of transcript--protein coupling across cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis. 

It addresses two fundamental questions:
1. **Quadrant Overlap & Conserved Pairs (Panel A):** To what extent do specific gene--protein pairs retain identical co-regulation behaviors (DiPa quadrants) when transitioning from biliary cholestasis to chronic chemical necrosis?
2. **Centroid Trajectory Slopes (Panels B & C):** How do the directional slopes ($\Delta \text{Protein} / \Delta \text{mRNA}$) within individual DiPa quadrants distribute across disease stages, characterizing the quantitative gain of protein synthesis per unit of transcriptional change?

---

## 2. Component Panels & Generating Scripts

### 🔹 Panel A (Top Row): DiPa Overlap Matrix & Conserved Pair Mapping
* **Biological Meaning:** Quantitative transition matrix and alluvial overlap mapping of the 943 core mRNA--protein pairs between BDL and $\text{CCl}_4$, identifying strictly conserved pairs (e.g. DiPa 1--2 concordant pairs maintaining up/downregulation across both diseases).
* **Component File:** `DiPa_Overlap_Analysis_Results.pdf` (7 KB)
* **Generating Script:** [`DiPa_Overlap_Analysis_BDL_CCL4.R`](DiPa_Overlap_Analysis_BDL_CCL4.R)
* **Input Datasets:** `Data_count_filtered_asbt_Gene_Protein_full.RData`, `Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData`.

---

### 🔹 Panel B (Bottom Left): Centroid Slope Distributions in BDL
* **Biological Meaning:** Frequency histograms of the centroid regression slopes ($S_j = \Delta P_j / \Delta G_j$) across DiPa quadrants in the BDL cohort ($n = 18$ mice), highlighting the strong non-zero positive slope distribution in synergistic quadrants.
* **Component File:** `Centroid_Slope_Frequency_BDL.pdf` (34 KB)
* **Generating Script:** [`Schwerpunkt_Slope_Distribution_BDL.R`](Schwerpunkt_Slope_Distribution_BDL.R)
* **Input Dataset:** `Data_count_filtered_asbt_Gene_Protein_full.RData`.

---

### 🔹 Panel C (Bottom Right): Centroid Slope Distributions in $\text{CCl}_4$
* **Biological Meaning:** Frequency histograms of the centroid regression slopes across DiPa quadrants in the $\text{CCl}_4$ cohort ($n = 36$ mice), confirming quantitative conservation of response gains across both etiologies.
* **Component File:** `Centroid_Slope_Frequency_CCL4.pdf` (36 KB)
* **Generating Script:** [`Schwerpunkt_Slope_Distribution_CCL4.R`](Schwerpunkt_Slope_Distribution_CCL4.R)
* **Input Dataset:** `Data_CCL4_filtered_Gene_Protein_full_final_with_DiPa.RData`.

---

## 3. Master Plate Assembly

* **Master Assembly Script:** [`Generate_Figure_4_Master_Plate.R`](Generate_Figure_4_Master_Plate.R)
* **Execution:** Combines the top row (Panel A: Overlap Analysis) and bottom row (Panels B & C: Centroid Slope Histograms) into a 2-row layout canvas using `cowplot`:
  ```bash
  Rscript Generate_Figure_4_Master_Plate.R
  ```
* **Output File:** `Figure_4.pdf` (714 KB).

---

## 4. Directory Layout

```text
Figure_04/
│
├── 📄 Figure_4.pdf                                      # Master Publication Plate (714 KB)
├── 📄 Generate_Figure_4_Master_Plate.R                  # Master assembly script (cowplot)
│
├── 📄 DiPa_Overlap_Analysis_Results.pdf                # Panel A: DiPa Overlap Matrix
├── 📄 Centroid_Slope_Frequency_BDL.pdf                 # Panel B: BDL Centroid Slopes
├── 📄 Centroid_Slope_Frequency_CCL4.pdf                # Panel C: CCl4 Centroid Slopes
│
├── 📄 DiPa_Overlap_Analysis_BDL_CCL4.R                 # Generates Panel A
├── 📄 Schwerpunkt_Slope_Distribution_BDL.R             # Generates Panel B
├── 📄 Schwerpunkt_Slope_Distribution_CCL4.R            # Generates Panel C
│
└── 📄 README.md                                         # This documentation file
```

---

## 5. Nature Communications Figure Legend

> **Figure 4: Cross-cohort DiPa quadrant overlap and centroid slope dynamics across liver injury etiologies.**  
> **(A)** Cross-cohort transition analysis comparing DiPa quadrant assignments of 943 shared mRNA--protein pairs between cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis. Numbers indicate the count and percentage of pairs maintaining conserved versus shifted regulatory coupling regimes.  
> **(B, C)** Frequency distributions of centroid trajectory slopes ($\Delta \text{Protein} / \Delta \text{mRNA}$) for BDL **(B)** and $\text{CCl}_4$ **(C)** across DiPa co-regulation quadrants.
