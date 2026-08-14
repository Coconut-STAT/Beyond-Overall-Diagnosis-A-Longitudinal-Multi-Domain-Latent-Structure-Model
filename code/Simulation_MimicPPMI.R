# Simulation_MimicPPMI.R
# Mimic-PPMI simulation with varying G (3 to 7)
# Output: ../output/MimicPPMI_G_{G}.Rdata
# Usage: Rscript Simulation_MimicPPMI.R <3-7>

rm(list = ls())
dir.create("../output", showWarnings = FALSE)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(RcppParallel)
library(dqrng)
library(doParallel)
library(doRNG)

n_core_set <- as.integer(Sys.getenv("N_CORES", unset = "102"))
seed_set <- 31
DS <- 100

set_para_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

dqset.seed(seed_set)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(min(n_core_set, n_cores))
registerDoParallel(cl)
registerDoRNG(seed_set)

IR <- 4000
BI <- 2000

load("Mimic_PPMI_data.RData")
load("Fitting_PPMI.RData")
XDATA_shared <- XDATA
dj_shared <- dj
ri_shared <- ri
tig_shared <- tig
missing_mask_shared <- missing_mask
Fitting_PPMI_shared <- Fitting_PPMI

result <- foreach(ds = 1:DS, .packages = c("Rcpp", "RcppArmadillo", "MASS", "RcppParallel", "MCMCpack", "dqrng")) %dorng% {
    Rcpp::sourceCpp("f_MCMC.cpp")
    RcppParallel::setThreadOptions(numThreads = 1)

    iter_seed <- seed_set * 1000 + ds
    set.seed(iter_seed)
    dqrng::dqset.seed(iter_seed)

    XDATA <- XDATA_shared
    dj <- dj_shared
    ri <- ri_shared
    tig <- tig_shared
    missing_mask <- missing_mask_shared
    Fitting_PPMI <- Fitting_PPMI_shared

    {
        G <- 5
        beta <- round(Fitting_PPMI$mean_beta, 2)
        B <- dim(beta)[1]
        N <- 849
        J <- length(dj)
        TT <- max(ri)

        X <- matrix(0, B, N * G)
        g <- 1
        for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- XDATA[i, b + 1]
        for (g in 2:G) for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- X[b, i + N * (1 - 1)]
        SigmaG <- round(Fitting_PPMI$mean_SigmaG, 2)
        Ui <- rep(0, N)
        Ui <- rnorm(N, 0, sqrt(SigmaG))
        Tig <- matrix(0, N, G)
        gammaL_est <- Fitting_PPMI$mean_gammaL

        order <- 3
        LengthVKnots <- 20
        knots <- seq(min(tig), max(tig), length = LengthVKnots)
        get_Haz_at_t <- function(t, g_idx, gammaL_est, knots, order) {
            bs_basis <- f_Ispline(order, t, knots)
            return(as.numeric(gammaL_est[, g_idx] %*% bs_basis))
        }

        for (i in 1:N) {
            for (g in 1:G) {
                U <- runif(1)
                V <- -log(U) / exp(sum(XDATA[i, -1] * beta[, g]) + Ui[i])
                res <- try(uniroot(function(t) get_Haz_at_t(t, g, gammaL_est, knots, order) - V,
                    interval = c(0, 100)
                ), silent = TRUE)
                if (class(res) == "try-error") {
                    Tig[i, g] <- 100 + 2
                } else {
                    Tig[i, g] <- res$root
                }
            }
        }
        Tig <- round(Tig, 3)

        InterL <- matrix(0, N, G)
        InterR <- matrix(0, N, G)
        status <- matrix(0, N, G)
        R <- matrix(0, N, G)
        for (i in 1:N) {
            for (g in 1:G) {
                if (findInterval(Tig[i, g], vec = tig[1:ri[i], i], rightmost.closed = TRUE) == 0) {
                    InterL[i, g] <- 0
                    InterR[i, g] <- tig[1, i]
                    status[i, g] <- 0
                    R[i, g] <- 1
                }
                if (findInterval(Tig[i, g], vec = tig[1:ri[i], i], rightmost.closed = TRUE) == ri[i]) {
                    InterL[i, g] <- tig[ri[i], i]
                    InterR[i, g] <- NA
                    status[i, g] <- 2
                    R[i, g] <- ri[i] + 1
                }
                if (findInterval(Tig[i, g], vec = tig[1:ri[i], i], rightmost.closed = TRUE) > 0 & findInterval(Tig[i, g], vec = tig[1:ri[i], i], rightmost.closed = TRUE) < ri[i]) {
                    Interval_index <- findInterval(Tig[i, g], vec = tig[1:ri[i], i], rightmost.closed = TRUE)
                    InterL[i, g] <- tig[Interval_index, i]
                    InterR[i, g] <- tig[Interval_index + 1, i]
                    status[i, g] <- 1
                    R[i, g] <- Interval_index + 1
                }
            }
        }

        CR <- rep(0, G)
        for (g in 1:G) {
            CR[g] <- sum(status[, g] == 2) / (N)
        }
        CR

        L <- matrix(0, J, G)
        for (j in 1:J) L[j, Fitting_PPMI$mode_S[j]] <- 1
        Theta <- Fitting_PPMI$mean_Theta

        Z <- matrix(0, G, max(ri) * N)
        for (i in 1:N) {
            for (g in 1:G) {
                Z[g, ((1 + TT * (i - 1)):(R[i, g] - 1 + TT * (i - 1)))] <- 1
                if (R[i, g] < ri[i] + 1) {
                    Z[g, ((R[i, g] + TT * (i - 1)):(ri[i] + TT * (i - 1)))] <- 2
                }
            }
        }

        Y <- matrix(0, J, TT * N)
        LL <- which(L == 1, arr.ind = TRUE)
        for (i in 1:N) {
            for (j in 1:J) {
                g <- unname(LL[which(LL[, 1] == j), 2])
                for (t in 1:ri[i]) {
                    Y[j, t + TT * (i - 1)] <- sample(1:dj[j], 1, TRUE, Theta[1:dj[j], Z[g, t + TT * (i - 1)] + 2 * (j - 1)])
                }
            }
        }
        Y[missing_mask] <- 0

        s <- numeric(J)
        for (j in 1:J) {
            for (g in 1:G) {
                if (L[j, g] == 1) {
                    s[j] <- g
                }
            }
        }

        Yijtc <- matrix(0, max(ri) * J, max(dj) * N)
        for (i in 1:N) {
            for (j in 1:J) {
                for (t in 1:ri[i]) {
                    for (c in 1:dj[j]) {
                        if (Y[j, t + TT * (i - 1)] == c) {
                            Yijtc[t + max(ri) * (j - 1), c + max(dj) * (i - 1)] <- 1
                        } else {
                            Yijtc[t + max(ri) * (j - 1), c + max(dj) * (i - 1)] <- 0
                        }
                    }
                }
            }
        }

        beta_real <- beta
        R_real <- R
        L_real <- L
        Theta_real <- Theta
        status_real <- status
        InterL_real <- InterL
        InterR_real <- InterR
        Ui_real <- Ui
        SigmaG_real <- SigmaG
        s_real <- s
    }
    {
        G <- set_para_num
        X <- matrix(0, B, N * G)
        g <- 1
        for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- XDATA[i, b + 1]
        for (g in 2:G) for (b in 1:B) for (i in 1:N) X[b, i + N * (g - 1)] <- X[b, i + N * (1 - 1)]
        repeat {
            XI <- rdirichlet(1, rep(1, G))
            s <- sample(1:G, J, TRUE, XI)
            if (length(unique(s)) == G & all(table(s) >= 2)) {
                break
            }
        }

        Theta <- matrix(0, max(dj), 2 * J)
        for (j in 1:J) {
            for (k in 1:2) {
                repeat {
                    sample <- (rdirichlet(1, rep(1, dj[j])))
                    if (any(sample < 0.1)) {
                        next
                    }
                    Theta[1:dj[j], k + 2 * (j - 1)] <- sample
                    if (k == 2) {
                        if (Theta[1, 1 + 2 * (j - 1)] < Theta[1, 2 + 2 * (j - 1)]) {
                            next
                        }
                    }
                    break
                }
            }
        }

        order <- 3
        LengthVKnots <- 20
        knots <- seq(min(tig), max(tig), length = LengthVKnots)
        Lknots <- LengthVKnots - 2 + order

        sig0 <- 2
        coef_range_beta <- 10
        coef_range_ui <- 10
        a_eta <- 0.1
        b_eta <- 0.1

        eta <- rgamma(G, a_eta, rate = b_eta)
        beta <- matrix(0, B, G)
        gammaL <- matrix(0.1, nrow = Lknots, ncol = G)

        sigma_mean <- 0.5
        sigma_variance <- 1^2
        sigma_alpha <- sigma_mean^2 / sigma_variance + 2
        sigma_beta <- sigma_mean * (sigma_alpha - 1)
        SigmaG <- rinvgamma(1, sigma_alpha, sigma_beta)
        Ui <- rnorm(N, 0, sqrt(SigmaG))
    }

    res <- f_mcmc(IR, BI, N, J, G, B, order, Lknots, ri, dj, Y, s, knots, Theta, Yijtc, beta, gammaL, eta, X, tig, Ui, sig0, coef_range_beta, coef_range_ui, a_eta, b_eta, SigmaG, sigma_alpha, sigma_beta, ds)
    res$Ui_e <- NULL
    res$gammaL_e <- NULL
    res$PRigt_e <- NULL
    res$beta_e <- NULL
    res$Theta_e <- NULL
    res$SigmaG_e <- NULL
    res$s_e <- NULL
    res$R_real <- R_real
    return(res)
}

stopCluster(cl)
save(result, file = paste0("../output/MimicPPMI_G_", set_para_num, ".Rdata"))

cat(sprintf(
    "Done. Saved ../output/%s\n",
    paste0("MimicPPMI_G_", set_para_num, ".Rdata")
))

