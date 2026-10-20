# ==============================================================================
# SCRIPT: Generate_Supplementary_Figure_6_Master_Plate.R
# LOCATION: C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Supplementary_Figure_06/
# PURPOSE: Assemble Supplementary Figure 6:
#          - Panel A (Top Left):  BDL Winner Model (K = 2) - Clean, no subtitle
#          - Panel B (Top Right): CCl4 Winner Model (K = 3) - Clean, no subtitle
#          - Panel C (Bottom):    BIC Bar Chart sorted by increasing BIC
#                                 - Clean, no title/subtitle, only axes & legend
#          - Capital tags: A, B, C
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

cat("=== Building Streamlined Supplementary Figure 6 Master Plate ===\n")

data_dir    <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share"
out_fig_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_Consortium_GitHub/04_Figure_Generation/Supplementary_Figure_06"
share_graf_dir <- "C:/Users/ngoune/Documents/Projet I/Protein_Modeling_share/New_Data/Grafiken_Paper_1"

if (!dir.exists(out_fig_dir)) dir.create(out_fig_dir, recursive = TRUE)
if (!dir.exists(share_graf_dir)) dir.create(share_graf_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. LOAD EXPRESSION DATA
# ------------------------------------------------------------------------------
load(file.path(data_dir, "Cross_Data_Analysis/DTccl4_DT_LCPM_Gene_Protein_full.RData"))
Full_DT <- copy(DTccl4_DT_LCPM); setDT(Full_DT)
Full_DT[, Protein_Raw := ProteinIntensity]
Full_DT[, Dataset := ifelse(is.na(TreatmentTime), "BDL Dataset", "CCl4 Dataset")]

# BDL
env_dipa_bdl <- new.env()
load(file.path(data_dir, "New_Data/Dipa/Data_count_filtered_asbt_Gene_Protein_full_dipa.RData"), envir = env_dipa_bdl)
DT_bdl <- as.data.table(env_dipa_bdl$DT_dipa_count_bdl)
DT_bdl[, log2G := log2(meanRatioG)]; DT_bdl[, log2P := log2(meanRatioP)]
DT_bdl <- DT_bdl[is.finite(log2G) & is.finite(log2P) & !is.na(DiPaGroups)]

Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "control", DiseaseGroup := "Control"]
Full_DT[Dataset == "BDL Dataset" & as.character(Treatment) == "BDL",     DiseaseGroup := "Disease"]
bdl_DT <- Full_DT[Dataset == "BDL Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]
Cloud_BDL <- bdl_DT[, .(N = .N, SD_RNA = sd(GeneCount), SD_Prot = sd(Protein_Raw), Pearson_R = suppressWarnings(cor(GeneCount, Protein_Raw))), by = GeneProtein]
Cloud_BDL <- Cloud_BDL[N >= 10 & !is.na(Pearson_R) & SD_RNA > 0 & SD_Prot > 0]
Combined_BDL <- merge(DT_bdl[, .(GeneProtein, log2G, log2P)], Cloud_BDL[, .(GeneProtein)], by = "GeneProtein")
X_bdl <- as.matrix(Combined_BDL[, .(log2G, log2P)])

# CCl4
env_dipa_ccl4 <- new.env()
load(file.path(data_dir, "New_Data/Dipa/Data_CCL4_filtered_Gene_Protein_full_final_with_dipa.RData"), envir = env_dipa_ccl4)
dipa_var <- ls(env_dipa_ccl4)[sapply(ls(env_dipa_ccl4), function(x) inherits(env_dipa_ccl4[[x]], "data.frame"))][1]
DT_ccl4 <- copy(env_dipa_ccl4[[dipa_var]]); setDT(DT_ccl4)
DT_ccl4[, log2G := log2(meanRatioG_12)]; DT_ccl4[, log2P := log2(meanRatioP_12)]
DT_ccl4 <- DT_ccl4[is.finite(log2G) & is.finite(log2P) & !is.na(DiPaGroups)]

Full_DT[, DiseaseGroup := NA_character_]
Full_DT[Dataset == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) == 0, DiseaseGroup := "Control"]
Full_DT[Dataset == "CCl4 Dataset" & as.numeric(as.character(TreatmentTime)) == 12 & as.character(Treatment) == "ccl4", DiseaseGroup := "Disease"]
ccl4_DT <- Full_DT[Dataset == "CCl4 Dataset" & !is.na(GeneCount) & !is.na(Protein_Raw) & !is.na(DiseaseGroup)]
Cloud_CCL4 <- ccl4_DT[, .(N = .N, SD_RNA = sd(GeneCount), SD_Prot = sd(Protein_Raw), Pearson_R = suppressWarnings(cor(GeneCount, Protein_Raw))), by = GeneProtein]
Cloud_CCL4 <- Cloud_CCL4[N >= 10 & !is.na(Pearson_R) & SD_RNA > 0 & SD_Prot > 0]
Combined_CCL4 <- merge(DT_ccl4[, .(GeneProtein, log2G, log2P)], Cloud_CCL4[, .(GeneProtein)], by = "GeneProtein", all.x = TRUE)
X_ccl4 <- as.matrix(Combined_CCL4[, .(log2G, log2P)])

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS (MATHEMATICAL FOUNDATIONS)
# ------------------------------------------------------------------------------

# Log-density of bivariate normal distribution N_2(mu, Sigma)
# Implemented using pure Base R matrix algebra (bypassing external dependencies)
# with analytical 2x2 determinant, matrix inversion, and Mahalanobis distance.
log_dmvnorm2 <- function(X, mu, Sigma) {
  diff_m <- sweep(X, 2, mu, "-")
  detS   <- Sigma[1, 1] * Sigma[2, 2] - Sigma[1, 2]^2
  if (detS <= 1e-12) detS <- 1e-12 # Numerical guard against singularity
  invS   <- matrix(c(Sigma[2, 2], -Sigma[1, 2], -Sigma[2, 1], Sigma[1, 1]), nrow = 2) / detS
  quad   <- rowSums((diff_m %*% invS) * diff_m)
  -0.5 * (2 * log(2 * pi) + log(detS) + quad)
}

# Spectral decomposition confidence ellipses:
# E_k(alpha) = mu_k + sqrt(chi^2_2(alpha)) * V * Lambda^(1/2) * [cos(phi), sin(phi)]^T
# Computes exact Mahalanobis contours enclosing alpha in {0.50, 0.90, 0.95} mass.
get_ellipse_coords <- function(mu, Sigma, levels = c(0.50, 0.90, 0.95), n_points = 200, label_prefix = "") {
  theta_seq   <- seq(0, 2 * pi, length.out = n_points)
  unit_circle <- rbind(cos(theta_seq), sin(theta_seq))
  eig         <- eigen(Sigma) # Spectral decomposition: Sigma = V %*% Lambda %*% V^T
  scale_mat   <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8)), 2)
  res <- lapply(levels, function(lvl) {
    c_rad  <- sqrt(qchisq(lvl, df = 2)) # Radii: r_0.50 ~ 1.177, r_0.90 ~ 2.146, r_0.95 ~ 2.448
    coords <- mu + c_rad * (scale_mat %*% unit_circle)
    data.table(Component = label_prefix, Level = paste0(round(lvl * 100), "%"), x = coords[1, ], y = coords[2, ])
  })
  rbindlist(res)
}

quadrant_labels <- data.frame(
  x = c(0, 3.8, -3.8, 0, 0, 3.8, -3.8, -3.8, 3.8),
  y = c(0, 3.8, -3.8, 3.8, -3.8, 0, 0, 3.8, -3.8),
  label = c("8", "1", "2", "3", "4", "5", "6", "7", "7")
)

c_colors <- c(
  "Diagonal"                      = "#1976D2", # Blue (Component 1: Stoichiometric coupling)
  "Vertical (protein-dominant)"   = "#7B1FA2", # Purple (Component 2: Protein-dominant decoupling)
  "Horizontal (protein-buffered)" = "#E65100"  # Orange (Component 3: Post-transcriptional buffering)
)

# ------------------------------------------------------------------------------
# 3. FIT WINNING MODELS FOR VISUALIZATION (PANEL A & B)
# ------------------------------------------------------------------------------

# --- BDL Cohort Winner: K = 2, c = 1.3 ---
# Biophysical constraints enforced:
#   - Component 1: rho_1 in [0.70, 0.90]
#   - Component 2: rho_2 = 0, sigma_y2^2 >= 1.3 * sigma_x2^2
N_bdl <- nrow(X_bdl)

# Initialization: Deterministic coarse partition
is_vert_bdl <- abs(X_bdl[, 1]) < 0.55 & abs(X_bdl[, 2]) > 0.40
gamma1_b <- ifelse(!is_vert_bdl, 0.85, 0.15)
gamma2_b <- 1 - gamma1_b

ll_prev <- -Inf
for (iter in 1:400) { # max_iter = 400 safety cap; converges via |Delta ell| < 1e-6
  # --- M-Step: 
  N1 <- sum(gamma1_b); N2 <- sum(gamma2_b)
  pi_bdl <- c(N1/N_bdl, N2/N_bdl)
  mu1_k <- colSums(X_bdl * gamma1_b) / N1; mu2_k <- colSums(X_bdl * gamma2_b) / N2
  diff1 <- sweep(X_bdl, 2, mu1_k, "-"); diff2 <- sweep(X_bdl, 2, mu2_k, "-")
  
  # Component 1 (Diagonal): Update variances and clamp correlation rho_1 in [0.70, 0.90] 
  vx1 <- sum(gamma1_b * diff1[,1]^2)/N1; vy1 <- sum(gamma1_b * diff1[,2]^2)/N1
  cov1 <- sum(gamma1_b * diff1[,1]*diff1[,2])/N1
  rho1 <- min(max(cov1 / (sqrt(vx1)*sqrt(vy1)), 0.70), 0.90)
  Sig1_bdl <- matrix(c(vx1, rho1*sqrt(vx1*vy1), rho1*sqrt(vx1*vy1), vy1), 2)
  
  # Component 2 (Vertical): rho_2 = 0, enforce vertical elongation vy2 >= 1.3 * vx2 
  vx2 <- sum(gamma2_b * diff2[,1]^2)/N2; vy2 <- sum(gamma2_b * diff2[,2]^2)/N2
  if (vy2 < 1.3 * vx2) vy2 <- 1.3 * vx2
  Sig2_bdl <- diag(c(vx2, vy2))
  
  # --- E-Step: 
  lf1 <- log_dmvnorm2(X_bdl, mu1_k, Sig1_bdl); lf2 <- log_dmvnorm2(X_bdl, mu2_k, Sig2_bdl)
  ln1 <- log(pi_bdl[1]) + lf1; ln2 <- log(pi_bdl[2]) + lf2
  ml  <- pmax(ln1, ln2) # Log-sum-exp stabilization to avoid underflow
  lden <- ml + log(exp(ln1 - ml) + exp(ln2 - ml))
  ll   <- sum(lden) # Evaluated total log-likelihood
  gamma1_b <- exp(ln1 - lden); gamma2_b <- exp(ln2 - lden)
  
  # Convergence check: |Delta ell| < 1e-6 (McLachlan & Peel, 2000)
  if (abs(ll - ll_prev) < 1e-6) break
  ll_prev <- ll
}

el_bdl <- rbind(
  get_ellipse_coords(mu1_k, Sig1_bdl, levels = c(0.50, 0.90, 0.95), label_prefix = "Diagonal"),
  get_ellipse_coords(mu2_k, Sig2_bdl, levels = c(0.50, 0.90, 0.95), label_prefix = "Vertical (protein-dominant)")
)

# --- CCl4 Cohort Winner: K = 3, c = 1.3 ---
# Biophysical constraints enforced:
#   - Component 1: rho_1 in [0.70, 0.90]
#   - Component 2: rho_2 = 0, sigma_y2^2 >= 1.3 * sigma_x2^2
#   - Component 3: rho_3 = 0, sigma_x3^2 >= 1.3 * sigma_y3^2
N_ccl4 <- nrow(X_ccl4)

# Initialization: Deterministic coarse partition
is_vert_c  <- abs(X_ccl4[, 1]) < 0.55 & abs(X_ccl4[, 2]) > 0.40
is_horiz_c <- abs(X_ccl4[, 1]) > 0.50 & abs(X_ccl4[, 2]) < 0.50
is_diag_c  <- !is_vert_c & !is_horiz_c
gamma1_c <- ifelse(is_diag_c, 0.70, 0.15)
gamma2_c <- ifelse(is_vert_c, 0.70, 0.15)
gamma3_c <- ifelse(is_horiz_c, 0.70, 0.15)
norm_g   <- gamma1_c + gamma2_c + gamma3_c
gamma1_c <- gamma1_c / norm_g; gamma2_c <- gamma2_c / norm_g; gamma3_c <- gamma3_c / norm_g

ll_prev <- -Inf
for (iter in 1:400) { # max_iter = 400 safety cap; converges via |Delta ell| < 1e-6
  # --- M-Step: 
  N1 <- sum(gamma1_c); N2 <- sum(gamma2_c); N3 <- sum(gamma3_c)
  pi_ccl4 <- c(N1/N_ccl4, N2/N_ccl4, N3/N_ccl4)
  mu1_k_c <- colSums(X_ccl4 * gamma1_c) / N1
  mu2_k_c <- colSums(X_ccl4 * gamma2_c) / N2
  mu3_k_c <- colSums(X_ccl4 * gamma3_c) / N3
  diff1 <- sweep(X_ccl4, 2, mu1_k_c, "-"); diff2 <- sweep(X_ccl4, 2, mu2_k_c, "-"); diff3 <- sweep(X_ccl4, 2, mu3_k_c, "-")
  
  # Component 1 (Diagonal): rho_1 in [0.70, 0.90] 
  vx1 <- sum(gamma1_c * diff1[,1]^2)/N1; vy1 <- sum(gamma1_c * diff1[,2]^2)/N1
  cov1 <- sum(gamma1_c * diff1[,1]*diff1[,2])/N1
  rho1 <- min(max(cov1 / (sqrt(vx1)*sqrt(vy1)), 0.70), 0.90)
  Sig1_ccl4 <- matrix(c(vx1, rho1*sqrt(vx1*vy1), rho1*sqrt(vx1*vy1), vy1), 2)
  
  # Component 2 (Vertical): rho_2 = 0, vy2 >= 1.3 * vx2 
  vx2 <- sum(gamma2_c * diff2[,1]^2)/N2; vy2 <- sum(gamma2_c * diff2[,2]^2)/N2
  if (vy2 < 1.3 * vx2) vy2 <- 1.3 * vx2
  Sig2_ccl4 <- diag(c(vx2, vy2))
  
  # Component 3 (Horizontal): rho_3 = 0, vx3 >= 1.3 * vy3 
  vx3 <- sum(gamma3_c * diff3[,1]^2)/N3; vy3 <- sum(gamma3_c * diff3[,2]^2)/N3
  if (vx3 < 1.3 * vy3) vx3 <- 1.3 * vy3
  Sig3_ccl4 <- diag(c(vx3, vy3))
  
  # --- E-Step: 
  lf1 <- log_dmvnorm2(X_ccl4, mu1_k_c, Sig1_ccl4)
  lf2 <- log_dmvnorm2(X_ccl4, mu2_k_c, Sig2_ccl4)
  lf3 <- log_dmvnorm2(X_ccl4, mu3_k_c, Sig3_ccl4)
  ln1 <- log(pi_ccl4[1]) + lf1; ln2 <- log(pi_ccl4[2]) + lf2; ln3 <- log(pi_ccl4[3]) + lf3
  ml  <- pmax(ln1, ln2, ln3) # Log-sum-exp stabilization
  lden <- ml + log(exp(ln1 - ml) + exp(ln2 - ml) + exp(ln3 - ml))
  ll   <- sum(lden) # Evaluated total log-likelihood
  gamma1_c <- exp(ln1 - lden); gamma2_c <- exp(ln2 - lden); gamma3_c <- exp(ln3 - lden)
  
  # Convergence check: |Delta ell| < 1e-6 (McLachlan & Peel, 2000)
  if (abs(ll - ll_prev) < 1e-6) break
  ll_prev <- ll
}

el_ccl4 <- rbind(
  get_ellipse_coords(mu1_k_c, Sig1_ccl4, levels = c(0.50, 0.90, 0.95), label_prefix = "Diagonal"),
  get_ellipse_coords(mu2_k_c, Sig2_ccl4, levels = c(0.50, 0.90, 0.95), label_prefix = "Vertical (protein-dominant)"),
  get_ellipse_coords(mu3_k_c, Sig3_ccl4, levels = c(0.50, 0.90, 0.95), label_prefix = "Horizontal (protein-buffered)")
)

# ------------------------------------------------------------------------------
# 4. BUILD PANEL A (BDL) & PANEL B (CCl4) - STREAMLINED (NO SUBTITLE)
# ------------------------------------------------------------------------------
DT_bdl_plot <- as.data.table(X_bdl); names(DT_bdl_plot) <- c("log2G", "log2P")
DT_ccl4_plot <- as.data.table(X_ccl4); names(DT_ccl4_plot) <- c("log2G", "log2P")

theme_panel_clean <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_blank(),
    axis.title    = element_text(size = 14, face = "bold", color = "black"),
    axis.text     = element_text(size = 13, face = "bold", color = "black"),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text   = element_text(size = 13, face = "bold"),
    legend.title  = element_blank(),
    legend.margin = margin(t = 2),
    plot.margin   = margin(t = 12, r = 12, b = 8, l = 12)
  )

p_a <- ggplot(DT_bdl_plot, aes(log2G, log2P)) +
  geom_point(color = "#607D8B", alpha = 0.40, size = 1.3) +
  geom_path(data = el_bdl, aes(x, y, group = interaction(Component, Level), color = Component), linewidth = 1.2) +
  scale_color_manual(values = c_colors) +
  geom_hline(yintercept = c(-0.5, 0.5), color = "grey35", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), color = "grey35", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7) +
  geom_label(data = quadrant_labels, aes(x, y, label = label), inherit.aes = FALSE, color = "grey30", fill = "white", linewidth = 0.3, fontface = "bold", size = 4.2) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5)) +
  labs(
    title = "BDL (K = 2)",
    x = expression(bold(RNA*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)")),
    y = expression(bold(Protein*":"~BDL~vehicle~vs.~sham~vehicle~~log[2]~"(fold change)"))
  ) +
  theme_panel_clean +
  guides(color = guide_legend(nrow = 1))

p_b <- ggplot(DT_ccl4_plot, aes(log2G, log2P)) +
  geom_point(color = "#607D8B", alpha = 0.40, size = 1.3) +
  geom_path(data = el_ccl4, aes(x, y, group = interaction(Component, Level), color = Component), linewidth = 1.2) +
  scale_color_manual(values = c_colors) +
  geom_hline(yintercept = c(-0.5, 0.5), color = "grey35", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), color = "grey35", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.7) +
  geom_label(data = quadrant_labels, aes(x, y, label = label), inherit.aes = FALSE, color = "grey30", fill = "white", linewidth = 0.3, fontface = "bold", size = 4.2) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5)) +
  labs(
    title = expression(bold(~CCl[4]~ "(K = 3)")),
    x = expression(bold("RNA: "~CCl[4]~"(Month 12) vs. oil (Month 0)"~~log[2]~"(fold change)")),
    y = expression(bold("Protein: "~CCl[4]~"(Month 12) vs. oil (Month 0)"~~log[2]~"(fold change)"))
  ) +
  theme_panel_clean +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# ------------------------------------------------------------------------------
# 5. PANEL C: BAR CHART SORTED BY INCREASING BIC (NO TITLE / SUBTITLE)
# ------------------------------------------------------------------------------
# Tested elongation factor values: c in {1.3, 1.5, 2.0, 2.5, 4.0}
# Corresponding to standard deviation ratios sqrt(c) in [1.14, 2.00]
c_vals <- c(1.3, 1.5, 2.0, 2.5, 4.0)

# Function to fit all 11 model candidates per cohort:
#   - Model 1: Null single bivariate Gaussian (K = 1, p = 5)
#   - Models 2-6: Constrained 2-component mixture across 5 elongation factors (K = 2,  p = 10)
#   - Models 7-11: Constrained 3-component mixture across 5 elongation factors (K = 3, p = 15)
# Parameter counts p and BIC penalization are defined in Section 4,
fit_cohort_models <- function(X, cohort_name) {
  N <- nrow(X)
  
  # --- K = 1: Unconstrained Single Bivariate Gaussian ---
  # Parameters: 2 means (mu_x, mu_y) + 2 variances (s_x^2, s_y^2) + 1 covariance (s_xy) => P = 5
  mu1  <- colMeans(X)
  S1   <- cov(X)
  ll1  <- sum(log_dmvnorm2(X, mu1, S1))
  bic1 <- -2 * ll1 + 5 * log(N) #  with p = 5
  
  res <- list()
  res[[1]] <- data.table(
    Cohort = cohort_name, K = 1, c = NA_real_, Ratio = 1.0, P = 5,
    LogLik = ll1, BIC = bic1, dBIC = 0, pi1 = 1.0, pi2 = 0.0, pi3 = 0.0,
    Status = "Null Single Gaussian"
  )
  
  # --- K = 2: Constrained Mixture (Diagonal + Vertical) across elongation factors c ---
  # Parameters: 1 weight (pi_1) + 4 means (2 centroids) + 3 dispersion (Comp 1) + 2 dispersion (Comp 2) => P = 10
  # Constraints: rho_1 in [0.70, 0.90];  rho_2 = 0, sigma_y2^2 >= c * sigma_x2^2
  is_vert <- abs(X[, 1]) < 0.55 & abs(X[, 2]) > 0.40
  for (f in c_vals) {
    gamma1 <- ifelse(!is_vert, 0.85, 0.15)
    gamma2 <- 1 - gamma1
    ll_prev <- -Inf
    for (iter in 1:400) { # EM algorithm: converges via |Delta ell| < 1e-6
      N1 <- sum(gamma1); N2 <- sum(gamma2)
      pi_k2 <- c(N1/N, N2/N)
      mu1_k <- colSums(X * gamma1) / N1; mu2_k <- colSums(X * gamma2) / N2
      diff1 <- sweep(X, 2, mu1_k, "-"); diff2 <- sweep(X, 2, mu2_k, "-")
      
      # M-Step: with constraints 
      vx1 <- sum(gamma1 * diff1[,1]^2)/N1; vy1 <- sum(gamma1 * diff1[,2]^2)/N1
      cov1 <- sum(gamma1 * diff1[,1]*diff1[,2])/N1
      rho1 <- min(max(cov1 / (sqrt(vx1)*sqrt(vy1)), 0.70), 0.90)
      Sig1 <- matrix(c(vx1, rho1*sqrt(vx1*vy1), rho1*sqrt(vx1*vy1), vy1), 2)
      
      vx2 <- sum(gamma2 * diff2[,1]^2)/N2; vy2 <- sum(gamma2 * diff2[,2]^2)/N2
      if (vy2 < f * vx2) vy2 <- f * vx2
      Sig2 <- diag(c(vx2, vy2))
      
      # E-Step: with log-sum-exp stabilization
      lf1 <- log_dmvnorm2(X, mu1_k, Sig1); lf2 <- log_dmvnorm2(X, mu2_k, Sig2)
      ln1 <- log(pi_k2[1]) + lf1; ln2 <- log(pi_k2[2]) + lf2
      ml  <- pmax(ln1, ln2)
      lden <- ml + log(exp(ln1 - ml) + exp(ln2 - ml))
      ll   <- sum(lden)
      gamma1 <- exp(ln1 - lden); gamma2 <- exp(ln2 - lden)
      
      if (abs(ll - ll_prev) < 1e-6) break
      ll_prev <- ll
    }
    bic2 <- -2 * ll + 10 * log(N) #  with p = 10
    st <- if (cohort_name == "BDL" && f == 1.3) "Winning Model (Lowest BIC)" else "Constrained Model"
    res[[length(res) + 1]] <- data.table(
      Cohort = cohort_name, K = 2, c = f, Ratio = sqrt(f), P = 10,
      LogLik = ll, BIC = bic2, dBIC = bic2 - bic1, pi1 = pi_k2[1], pi2 = pi_k2[2], pi3 = 0.0,
      Status = st
    )
  }
  
  # --- K = 3: Constrained Mixture (Diagonal + Vertical + Horizontal) across elongation factors c ---
  # Parameters: 2 weights + 6 means (3 centroids) + 3 (Comp 1) + 2 (Comp 2) + 2 (Comp 3) => p = 15
  # Constraints:  rho_1 in [0.70, 0.90];  rho_2 = 0, sigma_y2^2 >= c*sigma_x2^2;
  #              rho_3 = 0, sigma_x3^2 >= c*sigma_y3^2
  is_vert  <- abs(X[, 1]) < 0.55 & abs(X[, 2]) > 0.40
  is_horiz <- abs(X[, 1]) > 0.50 & abs(X[, 2]) < 0.50
  is_diag  <- !is_vert & !is_horiz
  for (f in c_vals) {
    gamma1 <- ifelse(is_diag, 0.70, 0.15)
    gamma2 <- ifelse(is_vert, 0.70, 0.15)
    gamma3 <- ifelse(is_horiz, 0.70, 0.15)
    norm_g <- gamma1 + gamma2 + gamma3
    gamma1 <- gamma1 / norm_g; gamma2 <- gamma2 / norm_g; gamma3 <- gamma3 / norm_g
    ll_prev <- -Inf
    for (iter in 1:400) { # EM algorithm: converges via |Delta ell| < 1e-6
      N1 <- sum(gamma1); N2 <- sum(gamma2); N3 <- sum(gamma3)
      pi_k3 <- c(N1/N, N2/N, N3/N)
      mu1_k <- colSums(X * gamma1) / N1; mu2_k <- colSums(X * gamma2) / N2; mu3_k <- colSums(X * gamma3) / N3
      diff1 <- sweep(X, 2, mu1_k, "-"); diff2 <- sweep(X, 2, mu2_k, "-"); diff3 <- sweep(X, 2, mu3_k, "-")
      
      # M-Step:  with constraints 
      vx1 <- sum(gamma1 * diff1[,1]^2)/N1; vy1 <- sum(gamma1 * diff1[,2]^2)/N1
      cov1 <- sum(gamma1 * diff1[,1]*diff1[,2])/N1
      rho1 <- min(max(cov1 / (sqrt(vx1)*sqrt(vy1)), 0.70), 0.90)
      Sig1 <- matrix(c(vx1, rho1*sqrt(vx1*vy1), rho1*sqrt(vx1*vy1), vy1), 2)
      
      vx2 <- sum(gamma2 * diff2[,1]^2)/N2; vy2 <- sum(gamma2 * diff2[,2]^2)/N2
      if (vy2 < f * vx2) vy2 <- f * vx2
      Sig2 <- diag(c(vx2, vy2))
      
      vx3 <- sum(gamma3 * diff3[,1]^2)/N3; vy3 <- sum(gamma3 * diff3[,2]^2)/N3
      if (vx3 < f * vy3) vx3 <- f * vy3
      Sig3 <- diag(c(vx3, vy3))
      
      # E-Step:  with log-sum-exp stabilization
      lf1 <- log_dmvnorm2(X, mu1_k, Sig1); lf2 <- log_dmvnorm2(X, mu2_k, Sig2); lf3 <- log_dmvnorm2(X, mu3_k, Sig3)
      ln1 <- log(pi_k3[1]) + lf1; ln2 <- log(pi_k3[2]) + lf2; ln3 <- log(pi_k3[3]) + lf3
      ml  <- pmax(ln1, ln2, ln3)
      lden <- ml + log(exp(ln1 - ml) + exp(ln2 - ml) + exp(ln3 - ml))
      ll   <- sum(lden)
      gamma1 <- exp(ln1 - lden); gamma2 <- exp(ln2 - lden); gamma3 <- exp(ln3 - lden)
      
      if (abs(ll - ll_prev) < 1e-6) break
      ll_prev <- ll
    }
    bic3 <- -2 * ll + 15 * log(N) #  with p = 15
    st <- if (cohort_name == "CCl4" && f == 1.3) "Winning Model (Lowest BIC)" else "Constrained Model"
    res[[length(res) + 1]] <- data.table(
      Cohort = cohort_name, K = 3, c = f, Ratio = sqrt(f), P = 15,
      LogLik = ll, BIC = bic3, dBIC = bic3 - bic1, pi1 = pi_k3[1], pi2 = pi_k3[2], pi3 = pi_k3[3],
      Status = st
    )
  }
  rbindlist(res)
}

tab_bdl  <- fit_cohort_models(X_bdl, "BDL")
tab_ccl4 <- fit_cohort_models(X_ccl4, "CCl4")
full_tab <- rbind(tab_bdl, tab_ccl4)

# Sort by increasing BIC within each cohort
setorder(full_tab, Cohort, BIC)
full_tab[, Rank := 1:.N, by = Cohort]
full_tab[, ModelClean := ifelse(K == 1, "K=1\nSingle", paste0("K=", K, "\nc=", c))]
full_tab[, ModelType := factor(paste0("K = ", K), levels = c("K = 1", "K = 2", "K = 3"))]
full_tab[, ModelInFacet := factor(paste0(Cohort, "___", ModelClean), levels = paste0(Cohort, "___", ModelClean))]
full_tab[, FacetCohort := factor(ifelse(Cohort == "BDL", "bold(BDL)", "bold(CCl[4])"),
                                 levels = c("bold(BDL)", "bold(CCl[4])"))]

scale_x_reordered <- function(...) {
  scale_x_discrete(labels = function(x) gsub(".*___", "", x), ...)
}

# Panel C Plot - Clean: NO TITLE, NO SUBTITLE, ONLY AXES AND LEGEND
p_c <- ggplot(full_tab, aes(x = ModelInFacet, y = BIC, fill = ModelType)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.0f", BIC)), vjust = -0.4, size = 4.2, fontface = "bold") +
  # Red asterisk placed safely above the winning bar (Rank == 1) with comfortable clearance, hidden from fill legend
  geom_point(data = full_tab[Rank == 1], aes(x = ModelInFacet, y = BIC + 48), 
             shape = 8, size = 4.5, color = "#D32F2F", stroke = 1.8, show.legend = FALSE) +
  coord_cartesian(ylim = c(3320, 3840)) +
  facet_wrap(~ FacetCohort, scales = "free_x", ncol = 2, labeller = label_parsed) +
  scale_x_reordered() +
  scale_fill_manual(values = c("K = 1" = "#90A4AE", "K = 2" = "#7E57C2", "K = 3" = "#26A69A")) +
  labs(
    x = "Model specification",
    y = "BIC",
    fill = "Component structure:"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    axis.text.x   = element_text(size = 12, face = "bold", color = "black", lineheight = 0.9),
    axis.text.y   = element_text(size = 13, face = "bold", color = "black"),
    axis.title    = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_rect(fill = "#ECEFF1", color = "black", linewidth = 1.0),
    strip.text    = element_text(size = 14, face = "bold"),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.0),
    legend.position = "bottom",
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.margin   = margin(t = 6, r = 12, b = 8, l = 12)
  )

# ------------------------------------------------------------------------------
# 6. ASSEMBLE MASTER PLATE WITH CAPITAL LABELS: A, B, C
# ------------------------------------------------------------------------------
cat("Assembling Streamlined Master Plate (A, B, C)...\n")

top_row <- plot_grid(
  p_a, p_b,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 24,
  label_fontface = "bold"
)

bottom_row <- plot_grid(
  p_c,
  ncol = 1,
  labels = c("C"),
  label_size = 24,
  label_fontface = "bold"
)

master_plate <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1.30, 0.95)
)

out_fig_pdf <- file.path(out_fig_dir, "Supplementary_Figure_6.pdf")
out_fig_png <- file.path(out_fig_dir, "Supplementary_Figure_06.png")
ggsave(out_fig_pdf, master_plate, width = 18, height = 16, device = cairo_pdf)
ggsave(out_fig_png, master_plate, width = 14.5, height = 14.5, dpi = 300)
ggsave(file.path(out_fig_dir, "Supplementary_Figure_6.png"), master_plate, width = 14.5, height = 14.5, dpi = 300)
cat(sprintf("[SUCCESS] Saved %s and %s\n", out_fig_pdf, out_fig_png))

# Mirror copy to Grafiken_Paper_1
file.copy(out_fig_pdf, file.path(share_graf_dir, basename(out_fig_pdf)), overwrite = TRUE)
file.copy(out_fig_png, file.path(share_graf_dir, basename(out_fig_png)), overwrite = TRUE)
file.copy(file.path(out_fig_dir, "Supplementary_Figure_6.png"), file.path(share_graf_dir, "Supplementary_Figure_6.png"), overwrite = TRUE)

# Copy to artifact directory for preview
file.copy(out_fig_png, "C:/Users/ngoune/.gemini/antigravity/brain/880061ef-cc44-4c0a-a528-730c28eaaa33/Supplementary_Figure_06.png", overwrite = TRUE)

# ------------------------------------------------------------------------------
# 7. EXPORT SUPPLEMENTARY TABLE 5
# ------------------------------------------------------------------------------
export_tab <- copy(full_tab[, .(Cohort, K, c, Ratio, P, LogLik, BIC, dBIC, pi1, pi2, pi3, Status)])
csv_path <- file.path(out_fig_dir, "Supplementary_Table_5_Mixture_Models.csv")

safe_fwrite <- function(dt, path) {
  tryCatch({
    fwrite(dt, path)
    cat(sprintf("[SUCCESS] Exported table to %s\n", path))
  }, error = function(e) {
    cat(sprintf("[NOTE] Could not overwrite %s (file is open in another program such as Excel). Preserved existing file.\n", basename(path)))
  })
}

safe_fwrite(export_tab, csv_path)
safe_fwrite(export_tab, file.path(share_graf_dir, "Supplementary_Table_5_Mixture_Models.csv"))

cat("=== Streamlined Supplementary Figure 6 completed! ===\n")
