# Simulation_Main.R
# Settings 1-4: (G,K,lambda) = (3,15,1), (5,25,1), (5,35,1), (3,15,2t)
# Output: ../output/{s}+20%+{N}.Rdata, ../output/{s}+40%+{N}.Rdata
# Usage: Rscript Simulation_Main.R <1-8>

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
seed_set <- 31
n_core_set <- as.integer(Sys.getenv("N_CORES", unset = "102"))
DS <- 100

get_para_by_num <- function(num) {
    all_settings <- list(
        list(1, 200),
        list(1, 500),
        list(2, 200),
        list(2, 500),
        list(3, 200),
        list(3, 500),
        list(4, 200),
        list(4, 500)
    )
    return(all_settings[[num]])
}

para <- get_para_by_num(set_para_num)
set_para1 <- para[[1]]
N_set <- para[[2]]

get_setting_num <- function(id, N) {
    B <- 2

    theta5 <- list(
        c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
        c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
        c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
        c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
        c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2)
    )
    theta7 <- list(
        c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
        c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
        c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
        c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
        c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2),
        c(0.5, 0.4, 0.1, 0.2, 0.1, 0.7),
        c(0.2, 0.7, 0.1, 0.1, 0.05, 0.85)
    )

    if (id == 1) {
        G <- 3
        J <- 15
        baseline_nu <- 1
        beta <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
        theta_vecs <- theta5
    } else if (id == 2) {
        G <- 5
        J <- 25
        baseline_nu <- 1
        beta <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        theta_vecs <- theta5
    } else if (id == 3) {
        G <- 5
        J <- 35
        baseline_nu <- 1
        beta <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        theta_vecs <- theta7
    } else if (id == 4) {
        G <- 3
        J <- 15
        baseline_nu <- 2
        beta <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
        theta_vecs <- theta5
    } else {
        stop("Invalid id; supported: 1, 2, 3, 4.")
    }

    if (G == 3) {
        extra_cens_probs <- c(0.4, 0.4, 0.15, 0.05)
    } else if (G == 5) {
        extra_cens_probs <- c(0.01, 0.5, 0.3, 0.15, 0.03, 0.01)
    }

    list(
        N = N, G = G, J = J, B = B, beta = beta,
        theta_vecs = theta_vecs, baseline_nu = baseline_nu, lambda0 = 1,
        extra_cens_probs = extra_cens_probs
    )
}

cfg <- get_setting_num(set_para1, N_set)
if (set_para1 == 4) {
    cfg$r20 <- 0.2
} else {
    cfg$r20 <- 0.13
}

IR <- 10000
BI <- 5000

fix_strict_increasing <- function(x, digits = 4, lower = 0, upper = NULL) {
    step <- 10^(-digits)
    if (!is.null(upper)) x[x > upper] <- upper
    x[x < lower] <- lower
    x <- sort(x)
    for (r in 2:length(x)) if (x[r] <= x[r - 1]) x[r] <- x[r - 1] + step
    if (!is.null(upper) && x[length(x)] > upper) {
        shift <- x[length(x)] - upper
        x <- x - shift
        x[x < lower] <- lower
        x <- sort(x)
        for (r in 2:length(x)) if (x[r] <= x[r - 1]) x[r] <- x[r - 1] + step
    }
    x <- round(x, digits)
    for (r in 2:length(x)) if (x[r] <= x[r - 1]) x[r] <- x[r - 1] + step
    round(x, digits)
}

dqset.seed(seed_set)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(min(n_core_set, n_cores))
registerDoParallel(cl)
registerDoRNG(seed_set)

cat(sprintf(
    "Running Setting #%d: setting_id=%d, N=%d, (K, G, lambda_0g(t)) = (%d, %d, %s)\n",
    set_para_num, set_para1, N_set, cfg$J, cfg$G, if (cfg$baseline_nu == 1) "1" else "2t"
))

result <- foreach(
    ds = 1:DS,
    .packages = c(
        "Rcpp", "RcppArmadillo", "MASS",
        "RcppParallel", "MCMCpack", "dqrng"
    )
) %dorng% {
    Rcpp::sourceCpp("f_MCMC.cpp")
    RcppParallel::setThreadOptions(numThreads = 1)

    iter_seed <- seed_set * 1000 + ds
    set.seed(iter_seed)
    dqrng::dqset.seed(iter_seed)

    N <- cfg$N
    G <- cfg$G
    J <- cfg$J
    B <- cfg$B
    beta <- cfg$beta
    baseline_nu <- cfg$baseline_nu
    lambda0 <- cfg$lambda0
    theta_vecs <- cfg$theta_vecs
    r20 <- cfg$r20
    extra_cens_probs <- cfg$extra_cens_probs

    X <- matrix(0, B, N * G)
    for (g in 1:G) {
        if (g == 1) {
            for (i in 1:N) X[, i + N * (g - 1)] <- rnorm(B, 0, 1)
        } else {
            for (i in 1:N) X[, i + N * (g - 1)] <- X[, i + N * (1 - 1)]
        }
    }

    SigmaG <- 0.5
    Ui <- rnorm(N, 0, sqrt(SigmaG))
    Unif_gener <- matrix(runif(N * G, 0, 1), nrow = N, ncol = G)
    Tig <- matrix(0, N, G)
    for (i in 1:N) {
        for (g in 1:G) {
            Tig[i, g] <- ((-log(Unif_gener[i, g]) *
                exp(-t(beta[, g]) %*% X[, i + N * (g - 1)] - Ui[i]))^
                (1 / baseline_nu)) / lambda0
        }
    }
    Tig <- round(Tig, 4)

    lower_len <- 7
    upper_len <- 12
    ri <- sample(lower_len:upper_len, N, replace = TRUE)
    TT <- max(ri)

    tig20 <- matrix(0, TT, N)
    tig40 <- matrix(0, TT, N)
    eps <- 1e-4
    tiny <- 1e-6

    for (i in 1:N) {
        m <- ri[i]
        Tmax20 <- rexp(1, rate = r20)
        gaps <- rexp(m - 2, rate = 1)
        visits <- cumsum(sort(gaps))
        inner <- visits / max(visits) * (Tmax20 - tiny) * (m - 2) / (m - 1) + tiny
        t20_raw <- c(tiny, inner, Tmax20)
        t20 <- fix_strict_increasing(round(t20_raw, 4), digits = 4, lower = tiny)
        tig20[1:m, i] <- t20

        ncens20_i <- sum(Tig[i, ] > Tmax20)
        extra_raw <- sample(0:(length(extra_cens_probs) - 1), 1,
            prob = extra_cens_probs
        )
        ncens40_i <- min(ncens20_i + extra_raw, G)
        extra_cens_i <- ncens40_i - ncens20_i

        if (extra_cens_i <= 0) {
            tig40[1:m, i] <- t20
            next
        }

        T_sorted <- sort(Tig[i, ])
        k <- G - ncens40_i + 1
        thr <- T_sorted[k] - eps

        last_kept <- if (k >= 2) T_sorted[k - 1] else 0
        rb_idx <- which(t20 > last_kept)[1]

        if (is.na(rb_idx) || t20[rb_idx] > thr) {
            tig40[1:m, i] <- t20
            next
        }

        n_anchor <- rb_idx
        n_compress <- m - n_anchor
        last_anchor <- t20[n_anchor]

        if (n_compress <= 0) {
            t40 <- t20
        } else {
            lo40 <- last_anchor + eps
            hi40 <- thr
            if (hi40 > lo40 + n_compress * eps) {
                new_pts <- seq(lo40, hi40, length.out = n_compress)
            } else {
                new_pts <- last_anchor + seq_len(n_compress) * eps
            }
            t40 <- c(t20[1:n_anchor], round(new_pts, 4))
        }

        t40 <- fix_strict_increasing(round(t40, 4),
            digits = 4,
            lower = tiny, upper = thr
        )

        if (length(t40) < m) {
            pad <- t40[length(t40)] + seq_len(m - length(t40)) * tiny
            t40 <- fix_strict_increasing(round(c(t40, pad), 4),
                digits = 4, lower = tiny, upper = thr
            )
        }
        if (length(t40) > m) t40 <- t40[1:m]

        tig40[1:m, i] <- t40
    }

    tig20 <- round(tig20, 4)
    tig40 <- round(tig40, 4)

    L <- matrix(0, J, G)
    for (g in 1:G) L[((g - 1) * (J / G) + 1):(g * (J / G)), g] <- 1
    dj <- rep(3, J)
    Theta_list <- lapply(seq_along(theta_vecs), function(k) {
        matrix(theta_vecs[[k]], nrow = dj[k], ncol = 2, byrow = FALSE)
    })
    Theta_real <- do.call(cbind, rep(Theta_list, times = G))

    s_real <- numeric(J)
    for (j in 1:J) for (g in 1:G) if (L[j, g] == 1) s_real[j] <- g

    beta_real <- beta
    L_real <- L
    Ui_real <- Ui
    SigmaG_real <- SigmaG
    X_real <- X

    run_scenario <- function(tig_use) {
        InterL <- matrix(0, N, G)
        InterR <- matrix(0, N, G)
        status <- matrix(0, N, G)
        R <- matrix(0, N, G)
        for (i in 1:N) {
            for (g in 1:G) {
                idx <- findInterval(Tig[i, g],
                    vec = tig_use[1:ri[i], i],
                    rightmost.closed = TRUE
                )
                if (idx == 0) {
                    InterL[i, g] <- 0
                    InterR[i, g] <- tig_use[1, i]
                    status[i, g] <- 0
                    R[i, g] <- 1
                }
                if (idx == ri[i]) {
                    InterL[i, g] <- tig_use[ri[i], i]
                    InterR[i, g] <- NA
                    status[i, g] <- 2
                    R[i, g] <- ri[i] + 1
                }
                if (idx > 0 && idx < ri[i]) {
                    InterL[i, g] <- tig_use[idx, i]
                    InterR[i, g] <- tig_use[idx + 1, i]
                    status[i, g] <- 1
                    R[i, g] <- idx + 1
                }
            }
        }
        CR <- sapply(1:G, function(g) sum(status[, g] == 2) / N)

        Z <- matrix(0, G, TT * N)
        for (i in 1:N) {
            for (g in 1:G) {
                Z[g, (1 + TT * (i - 1)):(R[i, g] - 1 + TT * (i - 1))] <- 1
                if (R[i, g] < ri[i] + 1) {
                    Z[g, (R[i, g] + TT * (i - 1)):(ri[i] + TT * (i - 1))] <- 2
                }
            }
        }

        Theta <- Theta_real
        LL_mat <- which(L == 1, arr.ind = TRUE)
        Y <- matrix(0, J, TT * N)
        for (i in 1:N) {
            for (j in 1:J) {
                g <- unname(LL_mat[which(LL_mat[, 1] == j), 2])
                for (t in 1:ri[i]) {
                    Y[j, t + TT * (i - 1)] <- sample(
                        1:dj[j], 1, TRUE,
                        Theta[, Z[g, t + TT * (i - 1)] + 2 * (j - 1)]
                    )
                }
            }
        }

        Yijtc <- matrix(0, max(ri) * J, max(dj) * N)
        for (i in 1:N) {
            for (j in 1:J) {
                for (t in 1:ri[i]) {
                    for (cc in 1:dj[j]) {
                        Yijtc[t + max(ri) * (j - 1), cc + max(dj) * (i - 1)] <-
                            as.numeric(Y[j, t + TT * (i - 1)] == cc)
                    }
                }
            }
        }

        R_real_sc <- R
        status_real_sc <- status

        repeat {
            XI <- rdirichlet(1, rep(1, G))
            s <- sample(1:G, J, TRUE, XI)
            if (length(unique(s)) == G & all(table(s) >= 2)) break
        }

        Theta_init <- Theta
        for (j in 1:J) {
            for (k in 1:2) {
                repeat {
                    samp <- rdirichlet(1, rep(1, dj[j]))
                    if (any(samp < 0.1)) next
                    Theta_init[1:dj[j], k + 2 * (j - 1)] <- samp
                    if (k == 2) {
                        if (Theta_init[1, 1 + 2 * (j - 1)] < Theta_init[1, 2 + 2 * (j - 1)]) {
                            next
                        }
                    }
                    break
                }
            }
        }

        order_s <- 3
        sig0 <- 2
        coef_range_beta <- 10
        coef_range_ui <- 10
        a_eta <- 0.1
        b_eta <- 0.1
        LengthVKnots <- 10

        T_range_max <- quantile(Tig, 0.99)
        knots <- seq(0, T_range_max, length.out = LengthVKnots)
        Lknots <- LengthVKnots - 2 + order_s

        eta <- rgamma(G, a_eta, rate = b_eta)
        beta_init <- matrix(0, B, G)
        gammaL <- matrix(rgamma(Lknots * G, 1, 1), nrow = Lknots, ncol = G)

        sigma_mean <- 0.5
        sigma_variance <- 1^2
        sigma_alpha <- sigma_mean^2 / sigma_variance + 2
        sigma_beta_pr <- sigma_mean * (sigma_alpha - 1)
        SigmaG_init <- rinvgamma(1, sigma_alpha, sigma_beta_pr)
        Ui_init <- rnorm(N, 0, sqrt(SigmaG_init))

        res <- f_mcmc(
            IR, BI, N, J, G, B, order_s, Lknots, ri, dj,
            Y, s, knots, Theta_init, Yijtc, beta_init, gammaL, eta,
            X, tig_use, Ui_init, sig0, coef_range_beta, coef_range_ui,
            a_eta, b_eta, SigmaG_init, sigma_alpha, sigma_beta_pr,
            ds
        )

        res$Ui_e <- NULL
        res$gammaL_e <- NULL
        res$PRigt_e <- NULL
        res$beta_e <- NULL
        res$Theta_e <- NULL
        res$SigmaG_e <- NULL
        res$s_e <- NULL
        res$R_real <- R_real_sc
        return(res)
    }

    cat(sprintf("  ds=%d  running 20%% ...\n", ds))
    res20 <- run_scenario(tig20)

    cat(sprintf("  ds=%d  running 40%% ...\n", ds))
    res40 <- run_scenario(tig40)

    list(res20 = res20, res40 = res40)
}

stopCluster(cl)

result_20 <- lapply(result, function(x) x$res20)
result_40 <- lapply(result, function(x) x$res40)

result <- result_20
save(result, file = paste0("../output/", set_para1, "+20%+", N_set, ".Rdata"))

result <- result_40
save(result, file = paste0("../output/", set_para1, "+40%+", N_set, ".Rdata"))

cat(sprintf(
    "Done. Saved ../output/%s and ../output/%s\n",
    paste0(set_para1, "+20%+", N_set, ".Rdata"),
    paste0(set_para1, "+40%+", N_set, ".Rdata")
))

