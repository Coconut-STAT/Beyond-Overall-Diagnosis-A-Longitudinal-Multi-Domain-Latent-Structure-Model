# TableS4_AdditionalRobustness.R
# Produces Table S4: Additional simulation with different theta settings
# Input: ../output/Additional_{s}+{CR}+{N}+{theta}.Rdata (from Simulation_Addition.R)

library(clue)
library(mclust)

base_path <- "../output"
DS        <- 100
SigmaG_real <- 0.5
B           <- 2
G           <- 3
J           <- 15

beta_real <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
s_real    <- rep(1:G, each = J / G)

theta_vecs_type1 <- list(
    c(0.7, 0.2, 0.1, 0.5, 0.1, 0.4),
    c(0.3, 0.6, 0.1, 0.1, 0.4, 0.5),
    c(0.6, 0.2, 0.2, 0.2, 0.3, 0.5),
    c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
    c(0.8, 0.1, 0.1, 0.4, 0.1, 0.5)
)
theta_vecs_type2 <- list(
    c(0.7, 0.2, 0.1, 0.5, 0.3, 0.2),
    c(0.3, 0.6, 0.1, 0.2, 0.5, 0.3),
    c(0.6, 0.2, 0.2, 0.4, 0.3, 0.3),
    c(0.9, 0.05, 0.05, 0.7, 0.15, 0.15),
    c(0.5, 0.4, 0.1, 0.4, 0.3, 0.3)
)

make_Theta_real <- function(theta_vecs) {
    t_list <- lapply(theta_vecs, function(v) matrix(v, nrow = 3, ncol = 2))
    do.call(cbind, rep(t_list, times = G))
}
Theta_real_type1 <- make_Theta_real(theta_vecs_type1)
Theta_real_type2 <- make_Theta_real(theta_vecs_type2)

all_settings <- list(
    list(s = 1, cen = "20%", n = 200, theta = 1),
    list(s = 1, cen = "20%", n = 200, theta = 2),
    list(s = 1, cen = "20%", n = 500, theta = 1),
    list(s = 1, cen = "20%", n = 500, theta = 2),
    list(s = 1, cen = "40%", n = 200, theta = 1),
    list(s = 1, cen = "40%", n = 200, theta = 2),
    list(s = 1, cen = "40%", n = 500, theta = 1),
    list(s = 1, cen = "40%", n = 500, theta = 2),
    list(s = 4, cen = "20%", n = 200, theta = 1),
    list(s = 4, cen = "20%", n = 200, theta = 2),
    list(s = 4, cen = "20%", n = 500, theta = 1),
    list(s = 4, cen = "20%", n = 500, theta = 2),
    list(s = 4, cen = "40%", n = 200, theta = 1),
    list(s = 4, cen = "40%", n = 200, theta = 2),
    list(s = 4, cen = "40%", n = 500, theta = 1),
    list(s = 4, cen = "40%", n = 500, theta = 2)
)

col_names <- c(
    "ARI_M", "ARI_IQR",
    "Theta_M", "Theta_IQR",
    "Beta_M", "Beta_IQR",
    "Sigma_M", "Sigma_IQR"
)

process_settings <- function(settings_list, theta_type) {
    n <- length(settings_list)
    tbl      <- matrix(NA, nrow = n, ncol = 8)
    colnames(tbl) <- col_names
    rnames   <- character(n)
    rig_mae  <- setNames(rep(NA_real_, n), seq_len(n))
    ari_acc  <- setNames(rep(NA_real_, n), seq_len(n))

    Theta_real <- if (theta_type == 1) Theta_real_type1 else Theta_real_type2

    for (idx in seq_along(settings_list)) {
        cfg <- settings_list[[idx]]
        file_path <- file.path(
            base_path,
            paste0("Additional_", cfg$s, "+", cfg$cen, "+", cfg$n, "+", cfg$theta, ".Rdata")
        )
        rnames[idx] <- sprintf("S%d_%s_N%d", cfg$s, cfg$cen, cfg$n)

        if (!file.exists(file_path)) {
            message("Skipping (not found): ", file_path)
            next
        }

        load(file_path)  

        S_ARI      <- numeric(DS)
        RMSE_Theta <- numeric(DS)
        RMSE_beta  <- numeric(DS)
        RMSE_Sigma <- numeric(DS)
        MAE_Rig    <- numeric(DS)
        R_abs      <- matrix(0, G, DS)

        for (ds in 1:DS) {
            res_ds <- result[[ds]]

            S_ARI[ds] <- adjustedRandIndex(s_real, t(res_ds$mode_S))

            RMSE_Theta[ds] <- sqrt(mean((res_ds$mean_Theta - Theta_real)^2))

            inner_prod   <- t(beta_real) %*% res_ds$mean_beta
            cost_matrix  <- max(inner_prod) - inner_prod
            perm         <- as.vector(solve_LSAP(cost_matrix))
            beta_aligned <- res_ds$mean_beta[, perm]
            RMSE_beta[ds] <- sqrt(mean((beta_aligned - beta_real)^2))

            RMSE_Sigma[ds] <- abs(SigmaG_real - res_ds$mean_SigmaG)

            R_real <- res_ds$R_real
            R_est  <- res_ds$mode_R
            inner_product <- t(R_real) %*% R_est
            permutation   <- apply(inner_product, 1, which.max)
            for (g in 1:G) {
                R_abs[g, ds] <- mean(abs(R_real[, g] - R_est[, permutation[g]]))
            }
            MAE_Rig[ds] <- mean(R_abs[, ds])
        }

        tbl[idx, ] <- c(
            median(S_ARI),      IQR(S_ARI),
            median(RMSE_Theta), IQR(RMSE_Theta),
            median(RMSE_beta),  IQR(RMSE_beta),
            median(RMSE_Sigma), IQR(RMSE_Sigma)
        )
        rig_mae[idx] <- mean(MAE_Rig)
        ari_acc[idx] <- mean(S_ARI == 1)

        rm(result); gc()
    }

    rownames(tbl) <- rnames
    names(rig_mae) <- rnames
    names(ari_acc) <- rnames
    list(table = round(tbl, 3), rig_mae = round(rig_mae, 3), ari_acc = round(ari_acc, 3))
}

settings_t1 <- all_settings[sapply(all_settings, function(x) x$theta == 1)]
settings_t2 <- all_settings[sapply(all_settings, function(x) x$theta == 2)]

out1 <- process_settings(settings_t1, theta_type = 1)
out2 <- process_settings(settings_t2, theta_type = 2)

print(out1$table)
cat("\nR_ig MAE\n")
print(out1$rig_mae)

cat("Domain accuracy:\n")
print(out2$ari_acc)

cat(sprintf("\nARI   - mean: %.3f,  IQR mean: %.3f\n",
            mean(out2$table[, "ARI_M"],    na.rm = TRUE),
            mean(out2$table[, "ARI_IQR"],  na.rm = TRUE)))
cat(sprintf("Theta - mean: %.3f\n",
            mean(out2$table[, "Theta_M"],  na.rm = TRUE)))
cat(sprintf("Beta  - mean: %.3f\n",
            mean(out2$table[, "Beta_M"],   na.rm = TRUE)))

write.csv(out1$table, "../output/TableS4.csv")
message("Saved: ../output/TableS4.csv")
