# Simulation_ModelSelection.R
# Model selection simulation: G candidates vary around true G
# Settings 1-6: (K,G,N) = (15,3,100), (15,3,200), (25,5,100), (25,5,200), (35,5,100), (35,5,200)
# Output: ../output/ModelSelection_{setting}_N{N}.Rdata
# Usage: Rscript Simulation_ModelSelection.R <1-6>

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

set_para_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

n_core_set <- as.integer(Sys.getenv("N_CORES", unset = "102"))
seed_set <- 31
DS <- 100

get_para_by_num <- function(num) {
    all_settings <- list(
        list(15, 3, 100),
        list(15, 3, 200),
        list(25, 5, 100),
        list(25, 5, 200),
        list(35, 5, 100),
        list(35, 5, 200)
    )
    return(all_settings[[num]])
}

para <- get_para_by_num(set_para_num)
J_true <- para[[1]]
G_true <- para[[2]]
N <- para[[3]]

if (G_true == 3) {
    G_candidates <- 2:6
} else if (G_true == 5) {
    G_candidates <- 3:7
} else {
    stop("Invalid G_true; only 3 or 5 supported.")
}

cat(sprintf(
    "Running Setting #%d: J=%d, G_true=%d, N=%d, CR=20%%\n",
    set_para_num, J_true, G_true, N
))
cat(sprintf("G candidates: %s\n", paste(G_candidates, collapse = ", ")))

dqset.seed(seed_set)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(min(n_core_set, n_cores))
registerDoParallel(cl)
registerDoRNG(seed_set)

IR <- 4000
BI <- 2000

get_true_setting <- function(J, G, N) {
    B <- 2
    if (G == 5) {
        beta <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
    } else if (G == 3) {
        beta <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
    } else {
        stop("Invalid G value; only 3 or 5 supported for true model.")
    }

    if (J / G == 7) {
        theta_vecs <- list(
            c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
            c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
            c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
            c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
            c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2),
            c(0.5, 0.4, 0.1, 0.2, 0.1, 0.7),
            c(0.2, 0.7, 0.1, 0.1, 0.05, 0.85)
        )
    } else if (J / G == 5) {
        theta_vecs <- list(
            c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
            c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
            c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
            c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
            c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2)
        )
    } else {
        stop("Invalid J/G ratio; only 5 or 7 supported.")
    }

    baseline_nu <- 1
    expRate <- 3

    list(
        N = N,
        G = G,
        J = J,
        beta = beta,
        theta_vecs = theta_vecs,
        baseline_nu = baseline_nu,
        expRate = expRate
    )
}

cfg_true <- get_true_setting(J_true, G_true, N)

result <- foreach(ds = 1:DS, .packages = c("Rcpp", "RcppArmadillo", "MASS", "RcppParallel", "MCMCpack", "dqrng")) %dorng% {
    Rcpp::sourceCpp("f_MCMC_ModelSelection.cpp")
    RcppParallel::setThreadOptions(numThreads = 1)

    iter_seed <- seed_set * 1000 + ds
    set.seed(iter_seed)
    dqrng::dqset.seed(iter_seed)

    {
        N <- cfg_true$N
        G <- cfg_true$G
        J <- cfg_true$J
        beta <- cfg_true$beta
        theta_vecs <- cfg_true$theta_vecs
        baseline_nu <- cfg_true$baseline_nu
        expRate <- cfg_true$expRate

        B <- 2
        X <- matrix(0, B, N * G)
        for (g in 1:G) {
            if (g == 1) {
                for (i in 1:N) {
                    X[, i + N * (g - 1)] <- rnorm(B, 0, 1)
                }
            } else {
                for (i in 1:N) {
                    X[, i + N * (g - 1)] <- X[, i + N * (1 - 1)]
                }
            }
        }

        X_real <- X[, 1:N]

        SigmaG <- 0.5
        Ui <- rnorm(N, 0, sqrt(SigmaG))
        Unif_gener <- matrix(runif(N * G, 0, 1), nrow = N, ncol = G)
        Tig <- matrix(0, N, G)
        lambda0 <- 1
        for (i in 1:N) {
            for (g in 1:G) {
                Tig[i, g] <- ((-log(Unif_gener[i, g]) * exp(-t(beta[, g]) %*% X[, i + N * (g - 1)] - Ui[i]))^(1 / baseline_nu)) / lambda0
            }
        }
        Tig <- round(Tig, 3)

        lower_len <- 7
        upper_len <- 12
        ri <- sample(lower_len:upper_len, N, replace = TRUE)
        TT <- max(ri)
        tig <- matrix(0, TT, N)
        for (i in 1:N) {
            tig[1, i] <- rexp(1, rate = expRate)
            for (r in 2:ri[i]) {
                tig[r, i] <- rexp(1, rate = expRate) + tig[r - 1, i]
            }
        }
        tig <- round(tig, 3)

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

        Z <- matrix(0, G, max(ri) * N)
        for (i in 1:N) {
            for (g in 1:G) {
                Z[g, ((1 + TT * (i - 1)):(R[i, g] - 1 + TT * (i - 1)))] <- 1
                if (R[i, g] < ri[i] + 1) {
                    Z[g, ((R[i, g] + TT * (i - 1)):(ri[i] + TT * (i - 1)))] <- 2
                }
            }
        }

        L <- matrix(0, J, G)
        for (g in 1:G) {
            L[((g - 1) * (J / G) + 1):(g * (J / G)), g] <- 1
        }
        dj <- rep(3, J)
        Theta <- matrix(0, max(dj), 2 * J)
        Theta_list <- lapply(seq_along(theta_vecs), function(k) {
            matrix(theta_vecs[[k]], nrow = dj[k], ncol = 2, byrow = FALSE)
        })
        Theta <- do.call(cbind, rep(Theta_list, times = G))

        Y <- matrix(0, J, TT * N)
        LL <- which(L == 1, arr.ind = TRUE)
        for (i in 1:N) {
            for (j in 1:J) {
                g <- unname(LL[which(LL[, 1] == j), 2])
                for (t in 1:ri[i]) {
                    Y[j, t + TT * (i - 1)] <- sample(1:dj[j], 1, TRUE, Theta[, Z[g, t + TT * (i - 1)] + 2 * (j - 1)])
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
        G_real <- G
        J_real <- J
    }

    result_all_G <- vector("list", length(G_candidates))
    names(result_all_G) <- paste0("G", G_candidates)

    for (g_idx in seq_along(G_candidates)) {
        G_fit <- G_candidates[g_idx]

        X <- matrix(0, B, N * G_fit)
        for (g in 1:G_fit) {
            for (i in 1:N) {
                X[, i + N * (g - 1)] <- X_real[, i]
            }
        }

        J_fit <- J_real

        best_res <- NULL
        best_loglik <- -Inf

        for (mcmc_run in 1:3) {{
            repeat {
                XI <- rdirichlet(1, rep(1, G_fit))
                s <- sample(1:G_fit, J_fit, TRUE, XI)
                if (length(unique(s)) == G_fit & all(table(s) >= 2)) {
                    break
                }
            }

            dj <- rep(3, J_fit)
            Theta <- matrix(0, max(dj), 2 * J_fit)
            for (j in 1:J_fit) {
                for (k in 1:2) {
                    repeat {
                        sample_theta <- (rdirichlet(1, rep(1, dj[j])))
                        if (any(sample_theta < 0.1)) {
                            next
                        }
                        Theta[1:dj[j], k + 2 * (j - 1)] <- sample_theta
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
            sig0 <- 2
            coef_range_beta <- 10
            coef_range_ui <- 10
            a_eta <- 0.1
            b_eta <- 0.1

            LengthVKnots <- 20
            knots <- seq(min(tig), max(tig), length = LengthVKnots)
            Lknots <- LengthVKnots - 2 + order

            eta <- rgamma(G_fit, a_eta, rate = b_eta)
            beta <- matrix(0, B, G_fit)
            gammaL <- matrix(rgamma(Lknots * G_fit, 1, 1), nrow = Lknots, ncol = G_fit)

            sigma_mean <- 1
            sigma_variance <- 0.5^2
            sigma_alpha <- sigma_mean^2 / sigma_variance + 2
            sigma_beta <- sigma_mean * (sigma_alpha - 1)
            sigma_beta <- sigma_alpha
            SigmaG <- rinvgamma(1, sigma_alpha, sigma_beta)
            Ui <- rnorm(N, 0, sqrt(SigmaG))

            Yijtc <- matrix(0, max(ri) * J_fit, max(dj) * N)
            for (i in 1:N) {
                for (j in 1:J_fit) {
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
        }

        res_run <- f_mcmc(IR, BI, N, J_fit, G_fit, B, order, Lknots, ri, dj, Y, s, knots, Theta, Yijtc, beta, gammaL, eta, X, tig, Ui, sig0, coef_range_beta, coef_range_ui, a_eta, b_eta, SigmaG, sigma_alpha, sigma_beta, ds)

        current_loglik <- ifelse(is.null(res_run$LogLikelihood), -Inf, res_run$LogLikelihood)
        if (current_loglik > best_loglik) {
            best_loglik <- current_loglik
            best_res <- res_run
        }}

        res <- best_res

        res$Ui_e <- NULL
        res$gammaL_e <- NULL
        res$PRigt_e <- NULL
        res$beta_e <- NULL
        res$Theta_e <- NULL
        res$SigmaG_e <- NULL
        res$s_e <- NULL

        res$G_fit <- G_fit
        res$G_true <- G_real
        res$J <- J_fit
        res$N <- N

        result_all_G[[g_idx]] <- res
    }

    result_all_G$R_real <- R_real
    result_all_G$G_true <- G_real
    result_all_G$J_true <- J_real
    result_all_G$beta_real <- beta_real
    result_all_G$L_real <- L_real

    return(result_all_G)
}

stopCluster(cl)

setting_names <- ceiling(set_para_num / 2)

save(result, file = paste0("../output/ModelSelection_", setting_names, "_N", N, ".Rdata"))

cat(sprintf(
    "Done. Saved ../output/%s\n",
    paste0("ModelSelection_", setting_names, "_N", N, ".Rdata")
))

