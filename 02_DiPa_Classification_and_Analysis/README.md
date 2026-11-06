# Stage 02: Differentiation Pattern (DiPa) & Multi-Omics Cloud Morphology Classification
## (Mapping mRNA–Protein Co-Regulation Regimes & Geometric Dispersion Shapes)

[![R](https://img.shields.io/badge/Language-R_%3E%3D_4.1.0-blue.svg)](https://www.r-project.org/)
[![Method: DiPa](https://img.shields.io/badge/Framework-DiPa_2D_Mapping-purple.svg)](https://pubmed.ncbi.nlm.nih.gov/)
[![Clustering: K--Means 2D](https://img.shields.io/badge/Morphology-K--Means_2D_Data--Driven-green.svg)](https://github.com/)

---

## 1. Scientific Overview & Core Architecture

Stage 02 establishes a **dual-layer characterization** of multi-omics regulation during liver injury across the **943** matched gene–protein pairs:

1. **Macro-Level Regulatory Shift (DiPa Groups 0 to 7):**  
   Spatial 2D coordinate positioning $(x_j, y_j) = (\log_2\text{FC}_{\text{RNA}}, \log_2\text{FC}_{\text{Protein}})$ partitioning the response into 8 discrete co-regulation quadrants.
2. **Micro-Level Sample Dispersion Geometry (Cloud Categories):**  
   Unsupervised data-driven classification of sample-level bivariate scatterplot morphologies (*Diagonal*, *Horizontal*, *Vertical*, *Round*) derived via 2D $K$-Means clustering.

---

## 2. Mathematical Definition of DiPa Quadrant Stratification

To compute mathematically unbiased Fold-Changes ($\text{FC}$):
* **RNA:** Raw gene counts (`GeneCount`), **NOT log-CPM data**, are used to calculate the mean expression in treatment and control conditions.
* **Protein:** Mass spectrometry $\log_2$ intensities are back-transformed to the natural linear scale ($2^{\text{ProteinIntensity}}$) prior to computing treatment-to-control ratios.

$$\text{Ratio}_{G,j} = \frac{\frac{1}{n_{\text{treat}}} \sum_{i \in \text{treat}} \text{GeneCount}_{ji}}{\frac{1}{n_{\text{ctrl}}} \sum_{k \in \text{ctrl}} \text{GeneCount}_{jk}}, \qquad \text{Ratio}_{P,j} = \frac{\frac{1}{n_{\text{treat}}} \sum_{i \in \text{treat}} 2^{\text{ProteinIntensity}_{ji}}}{\frac{1}{n_{\text{ctrl}}} \sum_{k \in \text{ctrl}} 2^{\text{ProteinIntensity}_{jk}}}$$

The 2D coordinate for each pair is: $(x_j, y_j) = \left( \log_2(\text{Ratio}_{G,j}), \; \log_2(\text{Ratio}_{P,j}) \right)$.

| `ClusterDiPa` | DiPa Group | Biological Regime | Coordinate Definition | $\text{CCl}_4$ Pairs ($N$) | BDL Pairs ($N$) | Conserved Target Evaluated Core ($N$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`0` / `8`** | DiPa 0 / 8 | **Baseline / Neutral (non-responsive)** | $\lvert x \rvert < 0.5 \text{ and } \lvert y \rvert < 0.5$ | 342 | 386 | **202** (DiPa 0/8) |
| **`1`** | DiPa 1 | **Synergistic Co-Upregulation (Up-Up)** | $x > 0.5 \text{ and } y > 0.5$ | 94 | 146 | \multirow{2}{*}{**118** (DiPa 1–2)} |
| **`2`** | DiPa 2 | **Synergistic Co-Downregulation (Down-Down)** | $x < -0.5 \text{ and } y < -0.5$ | 88 | 137 | |
| **`3`** | DiPa 3 | **Protein Upregulation Only** (Post-transcriptional accumulation) | $\lvert x \rvert \le 0.5 \text{ and } y > 0.5$ | 57 | 56 | \multirow{2}{*}{**67** (DiPa 3–4)} |
| **`4`** | DiPa 4 | **Protein Downregulation Only** (Post-transcriptional degradation) | $\lvert x \rvert \le 0.5 \text{ and } y < -0.5$ | 209 | 126 | |
| **`5`** | DiPa 5 | **RNA Upregulation Only** (Translational buffering) | $x > 0.5 \text{ and } \lvert y \rvert < 0.5$ | 113 | 68 | \multirow{2}{*}{**23** (DiPa 5–6)} |
| **`6`** | DiPa 6 | **RNA Downregulation Only** (Protein persistence) | $x < -0.5 \text{ and } \lvert y \rvert < 0.5$ | 14 | 17 | |
| **`7`** | DiPa 7 | **Antagonistic / Inverted Regulation** | $(x \le -0.5, y \ge 0.5) \lor (x \ge 0.5, y \le -0.5)$ | 26 | 7 | — |
| **Total** | | | | **943** | **943** | **410** |

> [!NOTE]
> **Cluster Harmonization Note:** Across all datasets, the canonical DiPa assignments are strictly synchronized: Cluster 7 contains exactly **26** pairs in $\text{CCl}_4$ and **7** pairs in BDL. The fold-change calculation evaluates all available biological replicates without premature omission of samples, ensuring exact alignment between raw count tables and normalized machine learning master matrices.

---

## 3. Unsupervised Data-Driven Cloud Morphology Derivation (`Cloud_Threshold_DataDriven.R`)

Rather than relying on arbitrary manual thresholds, the boundaries defining the 4 canonical scatterplot shapes were derived using **2D $K$-Means clustering ($K = 4$, $n_{\text{start}} = 50$, $\text{seed} = 42$)** on the normalized space of variance ratio $\text{Ratio} = \frac{\text{SD}(\text{RNA})}{\text{SD}(\text{Protein})}$ and absolute Pearson correlation $|R|$:

| Cloud Category | Color Code | Derived Decision Boundary | Biological Meaning |
| :--- | :---: | :--- | :--- |
| **Diagonal** | **Blue** (`#2196F3`) | $|R| > 0.467 \;\land\; 0.528 \le \text{Ratio} \le 1.014$ | Strong translational coupling; RNA directly drives Protein. |
| **Horizontal** | **Orange** (`#FF9800`) | $\text{Ratio} > 1.014 \;\land\; |R| < 0.673$ | High transcriptomic variance; post-transcriptional buffering. |
| **Vertical** | **Purple** (`#9C27B0`) | $\text{Ratio} < 0.528 \;\land\; |R| < 0.673$ | High proteomic variance; dominant post-translational regulation. |
| **Round** | **Green** (`#4CAF50`) | $0.528 \le \text{Ratio} \le 1.014 \;\land\; |R| < 0.070$ | Isotropic basal biological noise. |
| **Unclassified** | **Grey** (`#BDBDBD`) | Intermediate / unassigned | Non-distinct dispersion pattern. |

---

## 4. Script & File Layout in Repository

```text
02_DiPa_Classification_and_Analysis/
│
├── 📄 01_calculate_dipa_quadrants.R           # Computes fold-change ratios on raw counts & assigns DiPa 0 to 7
├── 📄 Cloud_Threshold_DataDriven.R            # Unsupervised 2D K-means clustering deriving optimal cloud boundaries
├── 📊 Cloud_Threshold_DataDriven.pdf          # 4-panel master diagnostic figure of cluster boundaries & distributions
└── 📄 README.md                               # Stage 02 documentation
```

### 🔹 Input & Output Datasets for Stage 02:
* **Inputs:**
  - BDL Raw Counts: `01_Data_Preprocessing.../raw_data/Data_count_filtered_asbt_Gene_Protein_full.RData`
  - $\text{CCl}_4$ Raw Counts: `01_Data_Preprocessing.../raw_data/Data_CCL4_filtered_Gene_Protein_full_final_with_dipa.RData`
* **Outputs:**
  - Vector PDF Quadrant Plots (BDL: `DiPa_newdata_BDL_control_count_943.pdf`; $\text{CCl}_4$: `DiPa_ccl4_oil_12_count_943.pdf`)
  - Updated RData files with exact canonical `ClusterDiPa` annotations propagated to Stage 01 Master Datasets.

---

## 5. Output Columns in Master Datasets

The resulting multi-omics classifications are stored as factors for downstream machine learning:
* `dataset$ClusterDiPa` $\in \{0, 1, 2, 3, 4, 5, 6, 7\}$ (Macro-level DiPa quadrant).
* `dataset$CloudCategory` $\in \{\text{"Diagonal"}, \text{"Horizontal"}, \text{"Vertical"}, \text{"Round"}, \text{"Unclassified"}\}$ (Micro-level shape).
