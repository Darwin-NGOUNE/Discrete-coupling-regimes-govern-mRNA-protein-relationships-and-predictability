# Stage 02: Differentiation Pattern (DiPa) Classification
## (Mapping mRNA--Protein Co-Regulation Regimes)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Method: DiPa](https://img.shields.io/badge/Framework-DiPa_2D_Mapping-purple.svg)](https://pubmed.ncbi.nlm.nih.gov/)

---

## 1. Scientific Overview & Core Principle

The **Differentiation Pattern (DiPa)** framework classifies transcript--protein pairs according to their joint response to disease intervention compared to baseline controls.

### ⚠️ Important Note on Input Data (Raw Counts vs. Log-CPM):
To compute mathematically unbiased Fold-Changes ($\text{FC}$):
* **RNA:** Raw gene counts (`GeneCount`), **NOT log-CPM data**, are used to calculate the mean expression in treatment and control conditions.
* **Protein:** Mass spectrometry $\log_2$ intensities are back-transformed to the natural linear scale ($2^{\text{ProteinIntensity}}$) prior to computing treatment-to-control ratios.

---

## 2. Mathematical Definition of Fold-Change Ratios

For each gene--protein pair $j$:

$$\text{Ratio}_{G,j} = \frac{\frac{1}{n_{\text{treat}}} \sum_{i \in \text{treat}} \text{GeneCount}_{ji}}{\frac{1}{n_{\text{ctrl}}} \sum_{k \in \text{ctrl}} \text{GeneCount}_{jk}}$$

$$\text{Ratio}_{P,j} = \frac{\frac{1}{n_{\text{treat}}} \sum_{i \in \text{treat}} 2^{\text{ProteinIntensity}_{ji}}}{\frac{1}{n_{\text{ctrl}}} \sum_{k \in \text{ctrl}} 2^{\text{ProteinIntensity}_{jk}}}$$

The 2D coordinate for each pair is then defined by:
$$(x_j, y_j) = \left( \log_2(\text{Ratio}_{G,j}), \; \log_2(\text{Ratio}_{P,j}) \right)$$

---

## 3. The 8 DiPa Quadrants ($\text{threshold} = 0.5$)

Using thresholds $\theta_x = 0.5$ ($\sim 1.41$-fold on RNA) and $\theta_y = 0.5$ ($\sim 1.41$-fold on Protein), pairs are assigned to discrete co-regulation regimes:

| `ClusterDiPa` | DiPa Group | Biological Regime | Coordinate Definition | Biological Description |
| :---: | :---: | :---: | :---: | :--- |
| **`0` / `8`** | DiPa 0 / 8 | **Neutral / Unchanged** | $\lvert x \rvert < 0.5 \text{ and } \lvert y \rvert < 0.5$ | Baseline central cloud; non-responsive |
| **`1`** | DiPa 1 | **Synergistic Upregulation** | $x > 0.5 \text{ and } y > 0.5$ | Concordant Up-Up co-induction |
| **`2`** | DiPa 2 | **Synergistic Downregulation** | $x < -0.5 \text{ and } y < -0.5$ | Concordant Down-Down co-repression |
| **`3`** | DiPa 3 | **Protein Upregulation Only** | $\lvert x \rvert \le 0.5 \text{ and } y > 0.5$ | Post-transcriptional protein accumulation |
| **`4`** | DiPa 4 | **Protein Downregulation Only** | $\lvert x \rvert \le 0.5 \text{ and } y < -0.5$ | Post-transcriptional protein degradation |
| **`5`** | DiPa 5 | **RNA Upregulation Only** | $x > 0.5 \text{ and } \lvert y \rvert < 0.5$ | Transcriptional buffering (translation block) |
| **`6`** | DiPa 6 | **RNA Downregulation Only** | $x < -0.5 \text{ and } \lvert y \rvert < 0.5$ | Transcriptional buffering (protein persistence) |
| **`7`** | DiPa 7 | **Antagonistic / Inverted** | $(x \le -0.5, y \ge 0.5) \text{ or } (x \ge 0.5, y \le -0.5)$ | Inverted regulation / discordant response |

---

## 4. Script & File Layout

```text
02_DiPa_Classification_and_Analysis/
│
├── 📄 01_calculate_dipa_quadrants.R           # Exact script: DiPa_new_old_data_filtered.R
│                                              # Computes ratios on raw counts and assigns DiPa groups
└── 📄 README.md                               # This documentation file
```

---

## 5. Output Annotations in Master Datasets

The resulting classification is mapped onto all observations as the categorical column:
```r
dataset$ClusterDiPa
```
This column is then directly used across all machine learning procedures (Procedures 1, 3, 5) to stratify models by co-regulation behavior.
