# TableS1_MainSimulation_fast.R
# Quick verification: produce Table S1 output for available settings only
# Unlike TableS1_MainSimulation.R (which requires all 8 settings), this script processes
# whichever result files exist in ../output/ and skips missing ones.
# Usage: Rscript TableS1_MainSimulation_fast.R

library(clue)
library(mclust)

base_path <- "../output"
SigmaG_real <- 0.5
all_theta_vecs <- list(
    c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
    c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
    c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
    c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
    c(0.8, 0.1, 0.1, 0.2, 0.6, 0.2),
    c(0.5, 0.4, 0.1, 0.2, 0.1, 0.7),
    c(0.2, 0.7, 0.1, 0.1, 0.05, 0.85)
)
scenarios <- c(1, 2, 3, 4)
censoring <- c("20%", "40%")
samples <- c(200, 500)

results_list <- list()

for (s in scenarios) {
    if (s %in% c(1, 4)) {
        G <- 3; J <- 15; B <- 2
        beta_real <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:5], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    } else if (s == 2) {
        G <- 5; J <- 25; B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:5], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    } else if (s == 3) {
        G <- 5; J <- 35; B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:7], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    }
    s_real <- rep(1:G, each = J / G)

    for (cen in censoring) {
        for (n in samples) {
            file_path <- file.path(base_path, paste0(s, "+", cen, "+", n, ".Rdata"))
            if (!file.exists(file_path)) {
                message("Skipping (not found): ", file_path)
                next
            }
            load(file_path)
            DS <- length(result)
            label <- paste0("S", s, "_", cen, "_N", n)

            S_ARI <- rep(0, DS)
            RMSE_Theta <- rep(0, DS)
            RMSE_beta <- rep(0, DS)
            RMSE_SigmaG <- rep(0, DS)

            for (ds in 1:DS) {
                res_ds <- result[[ds]]
                S_ARI[ds] <- adjustedRandIndex(s_real, t(res_ds$mode_S))
                RMSE_Theta[ds] <- sqrt(mean((res_ds$mean_Theta - Theta_real)^2))

                inner_prod <- t(beta_real) %*% res_ds$mean_beta
                cost_matrix <- max(inner_prod) - inner_prod
                perm <- as.vector(solve_LSAP(cost_matrix))
                beta_mean_aligned <- res_ds$mean_beta[, perm]
                RMSE_beta[ds] <- sqrt(mean((beta_mean_aligned - beta_real)^2))
                RMSE_SigmaG[ds] <- abs(SigmaG_real - res_ds$mean_SigmaG)
            }

            results_list[[label]] <- c(
                ARI_M = median(S_ARI), ARI_IQR = IQR(S_ARI),
                Theta_M = median(RMSE_Theta), Theta_IQR = IQR(RMSE_Theta),
                Beta_M = median(RMSE_beta), Beta_IQR = IQR(RMSE_beta),
                Sigma_M = median(RMSE_SigmaG), Sigma_IQR = IQR(RMSE_SigmaG)
            )

            rm(result); gc()
        }
    }
}

if (length(results_list) > 0) {
    results_summary <- do.call(rbind, results_list)
    results_summary <- round(results_summary, 3)
    cat("\n=== Table S1 (available settings) ===\n\n")
    print(results_summary)
    write.csv(as.data.frame(results_summary), "../output/TableS1_fast.csv")
    message("\nSaved: ../output/TableS1_fast.csv")
} else {
    message("No result files found in ../output/. Run Simulation_Main.R first.")
}
