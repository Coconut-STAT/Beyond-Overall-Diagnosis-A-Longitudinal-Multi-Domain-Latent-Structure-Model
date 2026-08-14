# TableS1_MainSimulation.R
# Produces Table S1 and Figure S13 from the main simulation results
# Input:  ../output/{s}+{CR}+{N}.Rdata (from Simulation_Main.R)
# Output: ../output/TableS1.csv and ../output/Figure_S13.pdf

library(clue)
library(mclust)

base_path <- "../output"
DS <- 100
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
results_summary <- matrix(NA, nrow = 16, ncol = 8)
colnames(results_summary) <- c(
    "ARI_M", "ARI_IQR", "Theta_M", "Theta_IQR",
    "Beta_M", "Beta_IQR", "Sigma_M", "Sigma_IQR"
)
row_names <- c()
row_idx <- 1

for (s in scenarios) {
    if (s %in% c(1, 4)) {
        folder <- ""
        G <- 3
        J <- 15
        B <- 2
        beta_real <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:5], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    } else if (s == 2) {
        folder <- ""
        G <- 5
        J <- 25
        B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:5], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    } else if (s == 3) {
        folder <- ""
        G <- 5
        J <- 35
        B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
        t_list <- lapply(all_theta_vecs[1:7], function(v) matrix(v, nrow = 3, ncol = 2))
        Theta_real <- do.call(cbind, rep(t_list, times = G))
    }
    s_real <- rep(1:G, each = J / G)

    for (cen in censoring) {
        for (n in samples) {
            file_path <- file.path(base_path, folder, paste0(s, "+", cen, "+", n, ".Rdata"))
            if (!file.exists(file_path)) {
                message("Skipping: ", file_path)
                next
            }
            load(file_path)
            row_names <- c(row_names, paste0("S", s, "_", cen, "_N", n))

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

            results_summary[row_idx, ] <- c(
                median(S_ARI), IQR(S_ARI),
                median(RMSE_Theta), IQR(RMSE_Theta),
                median(RMSE_beta), IQR(RMSE_beta),
                median(RMSE_SigmaG), IQR(RMSE_SigmaG)
            )

            row_idx <- row_idx + 1
            rm(result)
            gc()
        }
    }
}

rownames(results_summary) <- row_names
results_summary <- round(results_summary, 3)
print(results_summary)
write.csv(as.data.frame(results_summary), "../output/TableS1.csv")
message("Saved: ../output/TableS1.csv")

calc_R_metrics <- function(file_path, G = 5, DS = 100) {

  load(file_path) 
  R_accuracy <- matrix(0, G, DS)
  R_abs <- matrix(0, G, DS)

  for (ds in 1:DS) {
    R_real <- result[[ds]]$R_real
    R_est <- result[[ds]]$mode_R
    inner_product <- t(R_real) %*% R_est
    permutation <- apply(inner_product, 1, which.max)
    for (g in 1:G) {
      R_accuracy[g, ds] <- mean(R_real[, g] == R_est[, permutation[g]])
      R_abs[g, ds] <- mean(abs(R_real[, g] - R_est[, permutation[g]]))
    }
  }

  return(c(
    Accuracy = mean(apply(R_accuracy, 1, mean)),
    MAE = mean(apply(R_abs, 1, mean))
  ))
}

paths <- c(
  S2 = "2+40%+200.Rdata",
  S3 = "3+40%+200.Rdata"
)
R_results_list <- lapply(file.path(base_path, paths), calc_R_metrics)
R_summary_table <- do.call(rbind, R_results_list)
R_summary_table <- as.data.frame(R_summary_table)
R_summary_table$Setting <- c("(25, 5, 1)", "(35, 5, 1)")
R_summary_table <- R_summary_table[, c("Setting", "Accuracy", "MAE")]
R_summary_table$Accuracy = round(R_summary_table$Accuracy, 2)
R_summary_table$MAE = round(R_summary_table$MAE, 3)
print(R_summary_table)

library(ggplot2)
library(gridExtra)

plot_list <- list()
plot_idx_200 <- c()
plot_idx_500 <- c()

for (s in scenarios) {
    if (s %in% c(1, 4)) {
        G <- 3
        J <- 15
        B <- 2
        beta_real <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
    } else if (s == 2) {
        G <- 5
        J <- 25
        B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
    } else if (s == 3) {
        G <- 5
        J <- 35
        B <- 2
        beta_real <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
    }

    for (n in samples) {
        if (s == 1) {
            setting_id <- paste0("S1_N", n)
            plot_title <- substitute("(J, G, " ~ lambda["0g"] ~ "(t)) = (15, 3, 1), N=" ~ N, list(N = n))
        } else if (s == 2) {
            setting_id <- paste0("S2_N", n)
            plot_title <- substitute("(J, G, " ~ lambda["0g"] ~ "(t)) = (25, 5, 1), N=" ~ N, list(N = n))
        } else if (s == 3) {
            setting_id <- paste0("S3_N", n)
            plot_title <- substitute("(J, G, " ~ lambda["0g"] ~ "(t)) = (35, 5, 1), N=" ~ N, list(N = n))
        } else if (s == 4) {
            setting_id <- paste0("S4_N", n)
            plot_title <- substitute("(J, G, " ~ lambda["0g"] ~ "(t)) = (15, 3, 2t), N=" ~ N, list(N = n))
        }
        rmse_data <- data.frame()

        for (cen in censoring) {
            file_path <- file.path(base_path, paste0(s, "+", cen, "+", n, ".Rdata"))
            if (!file.exists(file_path)) next
            load(file_path)

            for (ds in 1:DS) {
                res_ds <- result[[ds]]
                inner_prod <- t(beta_real) %*% res_ds$mean_beta
                cost_matrix <- max(inner_prod) - inner_prod
                perm <- as.vector(solve_LSAP(cost_matrix))
                beta_aligned <- res_ds$mean_beta[, perm]
                rmse_beta_val <- sqrt(mean((beta_aligned - beta_real)^2))

                rmse_data <- rbind(rmse_data, data.frame(
                    RMSE = rmse_beta_val,
                    Censoring = cen,
                    Setting = setting_id,
                    stringsAsFactors = FALSE
                ))
            }
            rm(result)
        }

        if (nrow(rmse_data) > 0) {
            p <- ggplot(rmse_data, aes(x = factor(Setting, levels = unique(Setting)), y = RMSE, fill = Censoring)) +
                geom_boxplot(alpha = 0.7, outlier.size = 0.8, width = 0.75) +
                scale_fill_manual(values = c("20%" = "#3498db", "40%" = "#e74c3c")) +
                scale_x_discrete(expand = expansion(add = 0.4)) +
                scale_y_continuous(limits = c(0, 0.4), expand = expansion(mult = c(0.01, 0.01))) +
                labs(x = NULL, y = expression(bold("RMSE of ") ~ bolditalic(beta)), fill = NULL) +
                ggtitle(plot_title) +
                theme_minimal() +
                theme(
                    panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
                    axis.text = element_text(size = 11, face = "bold"),
                    axis.text.x = element_blank(),
                    axis.title.y = element_text(size = 12, face = "bold"),
                    axis.ticks.x = element_blank(),
                    plot.title = element_text(size = 11, hjust = 0.5, face = "bold",
                                              margin = margin(b = 2)),
                    legend.position = "bottom",
                    legend.text = element_text(size = 11, face = "bold"),
                    legend.key.size = unit(0.5, "cm"),
                    legend.margin = margin(t = 1, b = 0),
                    legend.box.margin = margin(t = -4, b = 2),
                    plot.margin = margin(t = 4, r = 8, b = 2, l = 8)
                )

            if (n == 200) {
                plot_idx_200 <- c(plot_idx_200, length(plot_list) + 1)
            } else {
                plot_idx_500 <- c(plot_idx_500, length(plot_list) + 1)
            }
            plot_list <- append(plot_list, list(p))
        }
    }
}

plot_list_ordered <- list()
for (i in 1:4) {
    if (length(plot_idx_200) >= i) plot_list_ordered <- append(plot_list_ordered, plot_list[plot_idx_200[i]])
    if (length(plot_idx_500) >= i) plot_list_ordered <- append(plot_list_ordered, plot_list[plot_idx_500[i]])
}

if (length(plot_list_ordered) > 0) {
    pdf("../output/Figure_S13.pdf", width = 12, height = 16)
    grid.arrange(grobs = plot_list_ordered, ncol = 2, nrow = 4)
    dev.off()
    message("Saved: ../output/Figure_S13.pdf")

}
