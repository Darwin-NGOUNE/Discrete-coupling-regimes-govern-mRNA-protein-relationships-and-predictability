# Supplementary Figure 6: Constrained Bivariate Gaussian Mixture Modeling of mRNA–Protein Fold Changes & Model Selection
## *(Nature Communications Manuscript)*

[![Vector PDF](https://img.shields.io/badge/Format-Vector_PDF_300DPI-red.svg)](Supplementary_Figure_06.pdf)
[![Nature Communications](https://img.shields.io/badge/Article-Supplementary_Figure_6-teal.svg)](https://www.nature.com/ncomms/)
[![Supplementary Table 5](https://img.shields.io/badge/Table-Supplementary_Table_5-blue.svg)](Supplementary_Table_5_Mixture_Models.csv)

---

## 1. Scientific Overview & Objective

**Supplementary Figure 6** provides rigorous, data-driven statistical validation for the **Differentiation Pattern (DiPa)** coupling regimes presented in **Figure 3**. 

A fundamental question arising from the 2D coordinate clouds of $\log_2(\text{fold change}_{\text{RNA}})$ vs $\log_2(\text{fold change}_{\text{Protein}})$ is whether the observed functional sub-regimes (canonical diagonal coupling, protein-dominant decoupling, and transcriptomic buffering) represent **bona fide biological subpopulations**, or whether they could be artifacts of the empirical $\pm 0.5$ $\log_2\text{FC}$ boundaries.

To answer this conclusively without subjective thresholding:
1. **Unbiased Latent Mixture Modeling:** We fitted constrained bivariate Gaussian Mixture Models ($K \in \{1, 2, 3\}$) using an Expectation-Maximization (EM) algorithm to the joint distributions of all 943 matched mRNA–protein pairs in both cholestatic (BDL) and toxic ($\text{CCl}_4$) liver fibrosis cohorts.
2. **Overfitting Protection via Elongation Constraints:** To prevent auxiliary components from capturing arbitrary local point density fluctuations, extra non-diagonal components were constrained to axis-aligned orientations with minimum elongation factors $c \in \{1.3, 1.5, 2.0, 2.5, 4.0\}$.
3. **Information-Theoretic Model Selection (BIC):** The optimal component complexity was selected strictly via the Bayesian Information Criterion (BIC), which severely penalizes model complexity.

---

## 2. Component Panels

### 🔹 Panel A: BDL Winning Constrained Mixture Model ($K = 2$)
* **Data:** 943 matched mRNA–protein pairs under cholestatic liver fibrosis ($N = 18$).
* **Winner Configuration:** $K = 2$ components, minimum elongation factor $c = 1.3$ (aspect ratio $\sqrt{c} = 1.14$).
* **Information Criteria:** $\text{BIC} = 3390.5$ ($\Delta\text{BIC} = -304.8$ compared to the single bivariate Gaussian null model $K = 1$).
* **Mixing Proportions ($\pi_k$):**
  * **Component 1 (Diagonal, Blue):** $\boldsymbol{\pi_1 = 90.9\%}$ ($\rho_1 \in [0.70, 0.90]$). The overwhelming majority of genes adhere to canonical, synchronous transcription–translation coupling.
  * **Component 2 (Vertical, Purple):** $\boldsymbol{\pi_2 = 9.1\%}$ ($\sigma_{y2}^2 \ge c \cdot \sigma_{x2}^2$). Captures a distinct subpopulation exhibiting selective protein-dominant remodeling in the absence of significant transcriptomic changes.
* **Contours:** Fitted bivariate density ellipses at the 50%, 90%, and 95% confidence intervals, overlaid on DiPa quadrants 1–8.

---

### 🔹 Panel B: $\text{CCl}_4$ Winning Constrained Mixture Model ($K = 3$)
* **Data:** 943 matched mRNA–protein pairs under toxic liver fibrosis ($N = 31$, Month 12 $\text{CCl}_4$ vs. Month 0 Oil).
* **Winner Configuration:** $K = 3$ components, minimum elongation factor $c = 1.3$ (aspect ratio $\sqrt{c} = 1.14$).
* **Information Criteria:** $\text{BIC} = 3490.6$ ($\Delta\text{BIC} = -222.3$ compared to $K = 1$).
* **Mixing Proportions ($\pi_k$):**
  * **Component 1 (Diagonal, Blue):** $\boldsymbol{\pi_1 = 61.6\%}$. Canonical synchronous coupling drops markedly relative to BDL (from $90.9\%$ down to $61.6\%$).
  * **Component 2 (Vertical, Purple):** $\boldsymbol{\pi_2 = 30.2\%}$. A dramatic threefold surge in protein-dominant deregulation, capturing acute toxic proteostatic stress.
  * **Component 3 (Horizontal, Orange):** $\boldsymbol{\pi_3 = 8.2\%}$ ($\sigma_{x3}^2 \ge c \cdot \sigma_{y3}^2$). Resolves post-transcriptional buffering (large mRNA changes where protein levels remain strictly stabilized).
* **Contours:** Fitted bivariate density ellipses at the 50%, 90%, and 95% confidence intervals.

---

### 🔹 Panel C: Systematic Model Selection via Bayesian Information Criterion (BIC)
* **Biological Content:** Comprehensive bar chart (*Balkendiagramm*) comparing all tested constrained Gaussian mixture models across both cohorts:
  * Single bivariate Gaussian null model ($K = 1, P = 5$)
  * 2-component constrained models ($K = 2, P = 10$) across $c \in \{1.3, 1.5, 2.0, 2.5, 4.0\}$
  * 3-component constrained models ($K = 3, P = 15$) across $c \in \{1.3, 1.5, 2.0, 2.5, 4.0\}$
* **Bar Chart Architecture & Visual Encoding:**
  * **Bars (`geom_col(width = 0.68)`):** Each bar represents a candidate mixture model specification. Bar height reflects the absolute BIC value (where lower BIC denotes superior model fit and parsimony).
  * **Color Fill Scheme (`ModelType`):** Discrete palette distinguishing component complexity:
    * Grey (`#90A4AE`): Single bivariate Gaussian ($K = 1$).
    * Purple (`#7E57C2`): 2-component mixture ($K = 2$).
    * Turquoise (`#26A69A`): 3-component mixture ($K = 3$).
  * **Numerical Values atop Bars:** Exact BIC printed directly above each bar (`geom_text(aes(label = sprintf("%.1f", BIC)))`) for transparent quantitative inspection.
  * **Winner Asterisk Indicator:** Prominent red asterisks (`shape = 8`, color `#D32F2F`) positioned 38 pt above the lowest-BIC bar in each facet, marking the optimal model without graphical overlap.
* **Faceted Sorting Mechanism (`scale_x_reordered`):**
  * Models within each etiology facet are ranked and displayed in strict ascending order of BIC (from best to worst fit).
  * In `ggplot2`, facet-specific reordering is implemented via composite level keys (`ModelInFacet = paste0(Cohort, "___", ModelClean)`), sorted via `setorder(full_tab, Cohort, BIC)`.
  * The custom scale transformation:
    ```r
    scale_x_reordered <- function(...) {
      scale_x_discrete(labels = function(x) gsub(".*___", "", x), ...)
    }
    ```
    strips the cohort prefix dynamically, ensuring clean, legible category labels (`K=2 (c=1.3)`, `K=3 (c=2.5)`, etc.) sorted independently within each facet without factor conflicts across cohorts.
  * **Mathematical Facet Strips:** Facet headers are rendered using `labeller = label_parsed`, producing bold, formatted labels **BDL** and $\mathbf{CCl_4}$ with mathematical subscript formatting.

---

## 3. Accompanying Supplementary Table 5

The numerical estimates and criteria for all 22 tested models are compiled into **Supplementary Table 5**:
* 📄 **CSV Spreadsheet:** [`Supplementary_Table_5_Mixture_Models.csv`](Supplementary_Table_5_Mixture_Models.csv)
* 📄 **LaTeX Table Snippet:** [`Supplementary_Table_5_Snippet.tex`](Supplementary_Table_5_Snippet.tex)
* 📄 **Standalone PDF Table:** [`Supplementary_Table_5_Standalone.pdf`](Supplementary_Table_5_Standalone.pdf)

### Summary of Key Quantitative Findings:
| Cohort | Optimal Model | Elongation ($c$) | Parameters ($P$) | Log-Likelihood ($\log L$) | BIC | $\Delta\text{BIC}$ vs $K=1$ | $\pi_1$ (Diag) | $\pi_2$ (Vert) | $\pi_3$ (Horiz) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **BDL** | **$K = 2$ (Winner)** | **1.3** | **10** | **-1660.9** | **3390.5** | **-304.8** | **90.9%** | **9.1%** | **—** |
| BDL | $K = 2$ | 2.5 | 10 | -1671.4 | 3411.6 | -283.7 | 92.0% | 8.0% | — |
| BDL | $K = 3$ | 1.3 | 15 | -1647.7 | 3398.5 | -296.8 | 85.8% | 8.8% | 5.4% |
| **$\text{CCl}_4$** | **$K = 3$ (Winner)** | **1.3** | **15** | **-1693.9** | **3490.6** | **-222.3** | **61.6%** | **30.2%** | **8.2%** |
| $\text{CCl}_4$ | $K = 3$ | 2.5 | 15 | -1703.4 | 3509.7 | -203.2 | 65.3% | 30.7% | 4.0% |
| $\text{CCl}_4$ | $K = 2$ | 1.3 | 10 | -1769.1 | 3606.4 | -106.5 | 69.5% | 30.5% | — |

---

## 4. Manuscript Integration Text

The following paragraph is incorporated into the main results section of the manuscript:

> *"To assess whether the apparent coupling regimes were supported by the underlying data independently of the predefined DiPa thresholds, we fitted constrained Gaussian mixture models to the joint distributions of mRNA and protein log$_2$ fold changes. Model selection using the Bayesian information criterion (BIC) favored a two-component model for BDL ($\Delta\mathrm{BIC} = -304.8$ relative to a single-component model), comprising a dominant diagonal component (90.9\%) and a smaller vertically oriented component (9.1\%), whereas a three-component model was favored for CCl$_4$ ($\Delta\mathrm{BIC} = -222.3$), resolving diagonal (61.6\%), vertical (30.2\%) and horizontal (8.2\%) components (\textbf{Supplementary Fig.~6 and Supplementary Table~5}). These results indicate that a single bivariate distribution was insufficient to describe the observed data and support the presence of distinct statistical subpopulations with different mRNA--protein response geometries. Thus, the observed patterns are not adequately explained by a single uniform relationship between mRNA and protein abundance."*

---

## 5. Execution & Reproduction

To reproduce Supplementary Figure 6 and Supplementary Table 5 from scratch, execute:

```bash
Rscript "Generate_Supplementary_Figure_6_Master_Plate.R"
pdflatex -interaction=nonstopmode "Supplementary_Table_5_Standalone.tex"
```
