#!/usr/bin/env Rscript
# pre-process.R
# Pre-process PPMI (Parkinson's Progression Markers Initiative) data for model input.
#
# ── Data Source ──
# PPMI data can be obtained from: https://www.ppmi-info.org/access-data-specimens/download-data
# After registering and obtaining access, download relevant data files from the following categories:
#   - "Motor Assessments" (MDS-UPDRS Part III items)
#   - "Non-motor Assessments" (MDS-UPDRS Part I/II items)
#   - "Subject Characteristics" (demographics, DAT scan, etc.)
#
# ── IMPORTANT: Data Preparation Before Running This Script ──
# 1. Download the relevant CSV files from PPMI.
# 2. Merge tables by PATNO (subject ID) and EVENT_ID (visit ID), and select variables to create:
#      - Y_data.csv: Longitudinal measurements. Columns include PATNO, visit/event info
#        or T_AGE (computed age at visit ), and MDS-UPDRS items (NP1*, NP2*, NP3* columns).
#      - X_data.csv: Subject-level covariates. Columns include PATNO, SEX, EDUCYRS
#        (education years), PUTAMEN_R (right putamen).

# ── PPMI Disclaimer ──
# The PPMI database is continuously evolving and updated over time.
# If the database content changes, the code output may differ from published results.
#
# ── Output ──
# After running this script, the following variables are ready for model input:
#   N, J, B, G, ri, dj, Y, X, tig, Yijtc, knots_s, etc.


# ── Number of latent groups (pre-specified) ──
G <- 5

{
  # ──────────────────────────────────────────────
  # Step 1: Read raw data
  # ──────────────────────────────────────────────
  # setwd('/path/to/your/data')  # <-- Set to the directory containing your data files
  dataY <- read.csv("Y_data.csv")
  dataX <- read.csv("X_data.csv")

  # ──────────────────────────────────────────────
  # Step 2: Select covariates and handle missing values
  # ──────────────────────────────────────────────
  dataX[dataX == ""] <- NA

  # Remove subjects with any missing covariates
  NA_PATNO <- dataX$PATNO[apply(dataX, 1, function(row) any(is.na(row)))]
  dataY <- dataY[!dataY$PATNO %in% NA_PATNO, ]

  # ──────────────────────────────────────────────
  # Step 3: Filter subjects with at least 5 visits
  # ──────────────────────────────────────────────
  SeriesLength_PATNO <- table(dataY$PATNO)
  SeriesLengthFilter <- names(SeriesLength_PATNO[SeriesLength_PATNO >= 5])
  dataY <- dataY[dataY$PATNO %in% SeriesLengthFilter, ]

  # ──────────────────────────────────────────────
  # Step 4: Clean response data
  # ──────────────────────────────────────────────
  dataY[, 2:ncol(dataY)] <- lapply(dataY[, 2:ncol(dataY)], as.numeric)

  # Shift categories: original 0-4 → 1-5, so 0 can represent missing
  dataY[, 3:ncol(dataY)] <- dataY[, 3:ncol(dataY)] + 1

  # Code abnormal/out-of-range entries as missing (0)
  dataY[, 3:ncol(dataY)][dataY[, 3:ncol(dataY)] == 102 | dataY[, 3:ncol(dataY)] == 10] <- NA
  dataY[, 3:ncol(dataY)][is.na(dataY[, 3:ncol(dataY)])] <- 0

  # ──────────────────────────────────────────────
  # Step 5: Select MDS-UPDRS variables of interest
  # ──────────────────────────────────────────────
  coll <- c("PATNO", "T_AGE",
            "NP1SLPN", "NP1SLPD", "NP1PAIN", "NP1URIN", "NP1CNST", "NP1FATG",
            "NP2SPCH", "NP2SALV", "NP2HWRT", "NP2HOBB", "NP2TRMR",
            "NP3RIGN", "NP3FACXP", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL",
            "NP3FTAPR", "NP3FTAPL", "NP3HMOVR", "NP3HMOVL", "NP3PRSPR", "NP3PRSPL",
            "NP3TTAPR", "NP3TTAPL", "NP3LGAGL", "NP3POSTR", "NP3BRADY",
            "NP3RTARU", "NP3RTCON")
  dataY <- dataY[, coll]

  # ──────────────────────────────────────────────
  # Step 6: Merge sparse response categories
  # ──────────────────────────────────────────────
  # Merge highest categories successively if they contain < 5% of observations
  for (coll in 3:ncol(dataY)) {
    dataY[, coll][dataY[, coll] == 5] <- 4  # Merge category 5 → 4
  }
  for (coll in 3:ncol(dataY)) {
    tableY <- as.numeric(table(dataY[, coll])) / length(dataY$PATNO)
    if (tableY[length(tableY)] < 0.05) dataY[, coll][dataY[, coll] == 4] <- 3
  }
  for (coll in 3:ncol(dataY)) {
    tableY <- as.numeric(table(dataY[, coll])) / length(dataY$PATNO)
    if (tableY[length(tableY)] < 0.05) dataY[, coll][dataY[, coll] == 3] <- 2
  }

  # Remove variables where highest remaining category still < 5%
  variables_delete <- c()
  for (coll in 3:ncol(dataY)) {
    tableY <- as.numeric(table(dataY[, coll])) / length(dataY$PATNO)
    if (tableY[length(tableY)] < 0.05) variables_delete <- c(variables_delete, names(dataY[, coll]))
  }
  dataY <- dataY[, !names(dataY) %in% variables_delete]
  variable_name <- colnames(dataY[, 3:ncol(dataY)])

  # ──────────────────────────────────────────────
  # Step 7: Finalize covariate matrix
  # ──────────────────────────────────────────────
  Ni <- unique(dataY$PATNO)
  dataX <- dataX[dataX$PATNO %in% Ni, ]
  if (sum(is.na(dataX)) != 0) stop("Error: NA values remain in covariate data.")
  dataX[2:ncol(dataX)] <- lapply(dataX[2:ncol(dataX)], as.numeric)
  XDATA <- as.matrix(dataX)

  # Standardize continuous covariates (columns 3-4: EDUCYRS, PUTAMEN_R)
  XDATA[, 3:4] <- apply(XDATA[, 3:4], 2, function(x) (x - mean(x)) / sd(x))

  # ──────────────────────────────────────────────
  # Step 8: Compute model dimensions
  # ──────────────────────────────────────────────
  N  <- length(unique(dataY$PATNO))  # Number of subjects
  ri <- as.numeric(table(dataY$PATNO))  # Visits per subject
  dj <- apply(dataY[, 3:ncol(dataY)], 2, function(x) length(unique(x))) - 1  # Categories per variable
  J  <- length(dj)   # Number of variables
  TT <- max(ri)      # Maximum number of visits
  B  <- ncol(XDATA) - 1  # Number of covariates

  # ──────────────────────────────────────────────
  # Step 9: Build visit-time matrix tig (max_ri × N)
  # ──────────────────────────────────────────────
  YDATA <- as.matrix(dataY[, -1])
  tig <- matrix(0, TT, N)
  for (t in 1:ri[1]) tig[t, 1] <- YDATA[t, 1]
  for (i in 2:N) {
    pre_index <- sum(ri[1:(i - 1)])
    for (t in 1:ri[i]) tig[t, i] <- YDATA[t + pre_index, 1]
  }

  # ──────────────────────────────────────────────
  # Step 10: Build response matrix Y (J × TT*N)
  # ──────────────────────────────────────────────
  Y <- matrix(0, J, TT * N)
  i <- 1
  for (j in 1:J) for (t in 1:ri[1]) Y[j, t + TT * (i - 1)] <- YDATA[t, j + 1]
  for (i in 2:N) {
    start_idx <- sum(ri[1:(i - 1)]) + 1
    end_idx   <- sum(ri[1:i])
    for (j in 1:J) Y[j, (1:ri[i]) + TT * (i - 1)] <- YDATA[start_idx:end_idx, j + 1]
  }

  # ──────────────────────────────────────────────
  # Step 11: Build covariate matrix X (B × N*G)
  # ──────────────────────────────────────────────
  X <- matrix(0, B, N * G)
  for (b in 1:B) for (i in 1:N) X[b, i + N * (1 - 1)] <- XDATA[i, b + 1]
  for (g in 2:G) for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- X[b, i + N * (1 - 1)]

  # ──────────────────────────────────────────────
  # Step 12: Build indicator tensor Yijtc (max_ri*J × max_dj*N)
  # ──────────────────────────────────────────────
  max_ri <- max(ri)
  max_dj <- max(dj)
  Yijtc <- matrix(0, max_ri * J, max_dj * N)
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

  cat(sprintf("Pre-processing complete: N=%d, J=%d, B=%d, G=%d\n", N, J, B, G))

  # ──────────────────────────────────────────────
  # Step 13: Run the MCMC with the pre-processed data
  # ──────────────────────────────────────────────
  Rcpp::sourceCpp("../code/f_MCMC.cpp")
  RcppParallel::setThreadOptions(numThreads = 1)

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

}