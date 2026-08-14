#!/usr/bin/env Rscript
# Model_Synthetic.R
# Run MCMC on a CSV dataset (e.g., synthetic or user-provided PPMI data)
#
# Input:
#   synthetic_Y.csv  -- longitudinal responses (PATNO, T_AGE, V1..VJ)
#   synthetic_X.csv  -- covariates (PATNO, covariate columns)
#
# Usage: Rscript Model_Synthetic.R
#   or use f_MCMC_parallel.cpp for within-chain parallelism (faster for single runs)

rm(list = ls())
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(RcppParallel)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

G          <- 5
IR         <- 10000
BI         <- 5000
seed_set   <- 31

set.seed(seed_set)

dataY_long <- read.csv("synthetic_Y.csv")
dataX      <- read.csv("synthetic_X.csv")

cat("Data loaded\n")
cat(sprintf("  Y: %d rows x %d cols\n", nrow(dataY_long), ncol(dataY_long)))
cat(sprintf("  X: %d rows x %d cols\n", nrow(dataX), ncol(dataX)))

{
  var_cols <- colnames(dataY_long)[!(colnames(dataY_long) %in% c("PATNO", "T_AGE"))]
  J        <- length(var_cols)
  dj <- sapply(var_cols, function(v) max(dataY_long[[v]], na.rm = TRUE))
  names(dj) <- NULL

  Ni <- sort(unique(dataY_long$PATNO))
  N  <- length(Ni)

  ri <- as.numeric(table(dataY_long$PATNO))
  TT <- max(ri)

  dataX <- dataX[order(dataX$PATNO), ]
  XDATA <- as.matrix(dataX)
  B     <- ncol(XDATA) - 1

  tig <- matrix(0, TT, N)
  dataY_long <- dataY_long[order(dataY_long$PATNO, dataY_long$T_AGE), ]
  idx <- 0
  for (i in 1:N) {
    for (t in 1:ri[i]) {
      idx <- idx + 1
      tig[t, i] <- dataY_long$T_AGE[idx]
    }
  }

  Y <- matrix(0, J, TT * N)
  idx <- 0
  for (i in 1:N) {
    for (t in 1:ri[i]) {
      idx <- idx + 1
      for (j in 1:J) {
        Y[j, t + TT * (i - 1)] <- dataY_long[idx, var_cols[j]]
      }
    }
  }

  max_ri <- max(ri)
  max_dj <- max(dj)
  Yijtc  <- matrix(0, max_ri * J, max_dj * N)
  for (i in 1:N) {
    for (j in 1:J) {
      for (t in 1:ri[i]) {
        cc <- Y[j, t + TT * (i - 1)]
        if (cc > 0 && cc <= dj[j]) {
          Yijtc[t + max_ri * (j - 1), cc + max_dj * (i - 1)] <- 1
        }
      }
    }
  }
}

cat(sprintf("Matrices built: N=%d, K=%d, B=%d, G=%d\n", N, J, B, G))

Rcpp::sourceCpp("../code/f_MCMC.cpp")
RcppParallel::setThreadOptions(numThreads = 1)

X <- matrix(0, B, N * G)
for (b in 1:B) for (i in 1:N) X[b, i + N * (1 - 1)] <- XDATA[i, b + 1]
for (g in 2:G) for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- X[b, i + N * (1 - 1)]

repeat {
  XI <- MCMCpack::rdirichlet(1, rep(1, G))
  s  <- sample(1:G, J, TRUE, XI)
  if (length(unique(s)) == G && all(table(s) >= 2)) break
}

Theta <- matrix(0, max(dj), 2 * J)
for (j in 1:J) {
  for (k in 1:2) {
    repeat {
      samp <- MCMCpack::rdirichlet(1, rep(1, dj[j]))
      if (any(samp < 0.1)) next
      Theta[1:dj[j], k + 2 * (j - 1)] <- samp
      if (k == 2)
        if (Theta[1, 1 + 2 * (j - 1)] < Theta[1, 2 + 2 * (j - 1)]) next
      break
    }
  }
}

order_s      <- 3
LengthVKnots <- 20
knots_s      <- seq(min(tig), max(tig), length = LengthVKnots)
Lknots       <- LengthVKnots - 2 + order_s

sig0            <- 2
coef_range_beta <- 10
coef_range_ui   <- 10
a_eta           <- 0.1
b_eta           <- 0.1

eta    <- rgamma(G, a_eta, rate = b_eta)
beta   <- matrix(0, B, G)
gammaL <- matrix(0.1, nrow = Lknots, ncol = G)

sigma_mean     <- 0.5
sigma_variance <- 1^2
sigma_alpha    <- sigma_mean^2 / sigma_variance + 2
sigma_beta_pr  <- sigma_mean * (sigma_alpha - 1)
SigmaG <- MCMCpack::rinvgamma(1, sigma_alpha, sigma_beta_pr)
Ui     <- rnorm(N, 0, sqrt(SigmaG))

cat(sprintf("\nRunning MCMC: IR=%d, BI=%d, G=%d\n", IR, BI, G))

res <- f_mcmc(IR, BI, N, J, G, B, order_s, Lknots, ri, dj, Y, s, knots_s, Theta,
              Yijtc, beta, gammaL, eta, X, tig, Ui,
              sig0, coef_range_beta, coef_range_ui,
              a_eta, b_eta, SigmaG, sigma_alpha, sigma_beta_pr,
              1L)

cat("\n")
cat(sprintf("  Model Results  (G = %d, IR = %d, BI = %d)\n", G, IR, BI))

dir.create("../output", showWarnings = FALSE)

# ── Variable name mapping: V1..V30 → MDS-UPDRS item names ──
item_names <- c("NP1SLPN", "NP1SLPD", "NP1PAIN", "NP1URIN", "NP1CNST", "NP1FATG",
                "NP2SPCH", "NP2SALV", "NP2HWRT", "NP2HOBB", "NP2TRMR",
                "NP3RIGN", "NP3FACXP", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL",
                "NP3FTAPR", "NP3FTAPL", "NP3HMOVR", "NP3HMOVL", "NP3PRSPR", "NP3PRSPL",
                "NP3TTAPR", "NP3TTAPL", "NP3LGAGL", "NP3POSTR", "NP3BRADY",
                "NP3RTARU", "NP3RTCON")
variable_description <- c("Sleep problems at night", "Daytime sleepiness",
                          "Pain and other sensations", "Urinary problems",
                          "Constipation problems", "Fatigue",
                          "Speech", "Saliva and drooling", "Handwriting",
                          "Doing hobbies and other activities", "Tremor",
                          "Rigidity - Neck", "Facial expression",
                          "Rigidity - RUE", "Rigidity - LUE",
                          "Rigidity - RLE", "Rigidity - LLE",
                          "Finger Tapping - R", "Finger Tapping - L",
                          "Hand movements - R", "Hand movements - L",
                          "Pronation-Supination - R", "Pronation-Supination - L",
                          "Toe tapping - R", "Toe tapping - L",
                          "Leg agility - L", "Posture",
                          "Global spontaneity of movement",
                          "Rest tremor amplitude", "Constancy of rest tremor")
if (J <= length(item_names)) {
  display_names <- item_names[1:J]
  display_desc  <- variable_description[1:J]
} else {
  display_names <- var_cols
  display_desc  <- var_cols
}

# ── Covariate names ──
cov_names <- colnames(dataX)[-1]  # exclude PATNO

# ── Table 1: Regression Coefficients (Mean, SD, 95% CI) + SigmaG ──
{
  n_post <- IR - BI
  beta_mean <- matrix(0, B, G)
  beta_sd   <- matrix(0, B, G)
  beta_q025 <- matrix(0, B, G)
  beta_q975 <- matrix(0, B, G)
  for (b in 1:B) {
    for (g in 1:G) {
      chain <- res$beta_e[b, g + G * (BI:(IR - 1))]
      beta_mean[b, g] <- mean(chain)
      beta_sd[b, g]   <- sd(chain)
      beta_q025[b, g] <- quantile(chain, 0.025)
      beta_q975[b, g] <- quantile(chain, 0.975)
    }
  }

  sig_chain <- res$SigmaG_e[(BI + 1):IR]
  sig_mean  <- mean(sig_chain)
  sig_sd    <- sd(sig_chain)
  sig_q025  <- quantile(sig_chain, 0.025)
  sig_q975  <- quantile(sig_chain, 0.975)

  tab_rows <- data.frame(
    Covariate  = character(0),
    Parameter  = character(0),
    Domain     = character(0),
    Mean       = numeric(0),
    SD         = numeric(0),
    CI_95      = character(0),
    stringsAsFactors = FALSE
  )
  for (b in 1:B) {
    for (g in 1:G) {
      tab_rows <- rbind(tab_rows, data.frame(
        Covariate = cov_names[b],
        Parameter = sprintf("beta_%d%d", b, g),
        Domain    = paste0("Domain ", g),
        Mean      = round(beta_mean[b, g], 4),
        SD        = round(beta_sd[b, g], 4),
        CI_95     = sprintf("[%.4f, %.4f]", beta_q025[b, g], beta_q975[b, g]),
        stringsAsFactors = FALSE
      ))
    }
  }
  tab_rows <- rbind(tab_rows, data.frame(
    Covariate = "",
    Parameter = "sigma^2",
    Domain    = "",
    Mean      = round(sig_mean, 4),
    SD        = round(sig_sd, 4),
    CI_95     = sprintf("[%.4f, %.4f]", sig_q025, sig_q975),
    stringsAsFactors = FALSE
  ))

  write.csv(tab_rows, "../output/Table1.csv", row.names = FALSE)
  cat("\n--- Table 1 saved to ../output/Table1.csv ---\n")
  print(tab_rows, row.names = FALSE)
}

# ── Domain Structure ──
s_est <- as.vector(res$mode_S)
L_est <- matrix(0, J, G)
for (j in 1:J) L_est[j, s_est[j]] <- 1
rownames(L_est) <- display_names
colnames(L_est) <- as.character(1:G)

# ── Figure 3: Domain Structure Heatmap ──
{
  L_new <- L_est
  for (i in 1:ncol(L_new)) {
    L_new[L_new[, i] == 1, i] <- i
    L_new[, i] <- as.character(L_new[, i])
  }
  data_long_fig3 <- melt(L_new)
  colnames(data_long_fig3) <- c("Variable", "Group", "value")
  data_long_fig3$Group <- as.character(data_long_fig3$Group)
  data_long_fig3$Variable <- as.character(data_long_fig3$Variable)
  data_long_fig3$Group <- factor(data_long_fig3$Group, levels = as.character(1:G))
  colors <- c('#237B9F', '#71BFB2', '#BF91D5', '#EC817E', '#F5B309')
  if (G > 5) colors <- rep(colors, length.out = G)
  color_mapping <- setNames(colors[1:G], as.character(1:G))
  data_long_fig3$fill_color <- ifelse(as.character(data_long_fig3$value) == "0", "white",
                                      as.character(data_long_fig3$Group))

  # Color for y-axis labels (colored by domain assignment)
  color_for_label <- sapply(rev(display_names), function(v) {
    idx <- which(display_names == v)
    color_mapping[as.character(s_est[idx])]
  })
  color_for_label <- unname(color_for_label)

  # Deterministic y mapping: NP1SLPN on top, NP3RTCON at bottom
  data_long_fig3$VariableNum <- J - match(data_long_fig3$Variable, display_names) + 1

  group_names <- paste0("Domain ", 1:G)

  p3 <- data_long_fig3 %>%
    ggplot(aes(x = Group, y = VariableNum, fill = fill_color)) +
    geom_tile(color = NA) +
    scale_y_continuous(
      breaks = seq_len(J),
      labels = rev(display_names),
      sec.axis = ggplot2::dup_axis(
        labels = rev(display_desc),
        name = ""
      )
    ) +
    scale_x_discrete(labels = group_names) +
    coord_cartesian(ylim = c(1.8, J - 0.8), xlim = c(1.1, G - 0.1)) +
    scale_fill_manual(values = c("white" = "white", color_mapping)) +
    theme_minimal() +
    labs(title = "", x = "Clinical domains", y = "MDS-UPDRS items") +
    theme(
      axis.text.x  = element_text(angle = 0, size = 12, face = "bold", margin = margin(t = 5),
                                   color = color_mapping[as.character(1:G)]),
      axis.text.y  = element_text(angle = 0, size = 12, color = color_for_label,
                                   margin = margin(r = 3), face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
      axis.text.y.right = element_text(hjust = 0, margin = margin(l = 10), size = 10),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks.y = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(-0.15, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", inherit.blank = FALSE, fill = NA, linewidth = 1),
      plot.margin  = margin(5, 5, 5, 5),
      legend.position = "none"
    )
  ggsave("../output/Figure_3.pdf", plot = p3, width = 11, height = 8.5, units = "in", dpi = 300)
  cat("--- Figure 3 saved to ../output/Figure_3.pdf ---\n")
}

# ── Figure 4: Emission Probabilities (Theta) ──
{
  K <- 2  # latent states: Normal, Severe
  theta_est <- res$mean_Theta
  max_dj <- max(dj)

  Theta_show1 <- matrix(0, J, max_dj * K)
  for (djcount in 1:max_dj) {
    for (k in 1:K) {
      for (j in 1:J) {
        Theta_show1[j, djcount + max_dj * (k - 1)] <- theta_est[djcount, k + K * (j - 1)]
      }
    }
  }
  rownames(Theta_show1) <- display_names

  n_panels  <- ceiling(J / 10)
  p_list    <- list()
  for (indexplot in 1:n_panels) {
    j_start <- 1 + 10 * (indexplot - 1)
    j_end   <- min(10 * indexplot, J)
    Theta_sub <- Theta_show1[j_start:j_end, , drop = FALSE]
    items_order <- rownames(Theta_sub)
    attributes_order <- c("Normal", "Severe")

    Theta_stack <- rbind(Theta_sub[, 1:max_dj], Theta_sub[, (max_dj + 1):(max_dj * K)])
    n_items <- nrow(Theta_sub)
    Theta_stack <- cbind(
      c(rownames(Theta_sub), rownames(Theta_sub)),
      c(rep("Normal", n_items), rep("Severe", n_items)),
      Theta_stack
    )
    col_vals <- paste0("value", 1:max_dj)
    colnames(Theta_stack) <- c("Items", "Attribute_Profiles", col_vals)
    Theta_stack <- as.data.frame(Theta_stack, stringsAsFactors = FALSE)
    Theta_stack$Items <- factor(Theta_stack$Items, levels = items_order)
    Theta_stack$Attribute_Profiles <- factor(Theta_stack$Attribute_Profiles, levels = attributes_order)

    dl <- Theta_stack %>%
      pivot_longer(cols = starts_with("value"), names_to = "value_type", values_to = "probability")
    dl$probability <- as.numeric(dl$probability)
    dl$value_type  <- gsub("value", "", dl$value_type)

    p_list[[indexplot]] <- ggplot(dl, aes(x = value_type, y = probability, fill = value_type)) +
      geom_bar(stat = "identity", position = "identity", width = 1, color = "black", linewidth = 0.8) +
      labs(x = "", y = "", title = "") +
      coord_polar(start = 0, direction = -1) +
      facet_grid(Items ~ Attribute_Profiles, switch = "y") +
      theme(
        legend.position  = "none",
        axis.ticks       = element_blank(),
        axis.text.x      = element_text(size = 8),
        axis.text.y      = element_blank(),
        panel.background = element_blank(),
        panel.spacing    = unit(-0.1, "cm"),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text.x     = element_text(size = 12, face = "bold", color = "black"),
        strip.text.y     = element_text(angle = 0, color = "black", size = 9, face = "bold",
                                        margin = margin(0, -0.05, 0, 0)),
        strip.placement  = "outside"
      ) +
      scale_fill_manual(values = colorRampPalette(c("#FFCCCC", "#FF3333"))(max_dj)) +
      geom_hline(yintercept = 1, color = "black", linewidth = 0.8) +
      ylim(0, 1)
  }
  pdf("../output/Figure_4.pdf", width = 3 * n_panels, height = 11)
  grid.arrange(grobs = p_list, ncol = n_panels, padding = unit(-3, "cm"))
  grid.text("MDS-UPDRS items", x = unit(0.015, "npc"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text("Latent states", y = unit(0.015, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  cat("--- Figure 4 saved to ../output/Figure_4.pdf ---\n")
}

# ── Save compact posterior result object ──
{
  res_save <- res
  res_save$Ui_e <- NULL
  res_save$gammaL_e <- NULL
  res_save$PRigt_e <- NULL
  res_save$beta_e <- NULL
  res_save$Theta_e <- NULL
  res_save$SigmaG_e <- NULL
  res_save$s_e <- NULL
  save(res_save, file = "../output/Fitting_Synthetic.RData")
  cat("--- Posterior summary saved to ../output/Fitting_Synthetic.RData ---\n")
}

cat("\n========================================================\n")
cat("Done.\n")
