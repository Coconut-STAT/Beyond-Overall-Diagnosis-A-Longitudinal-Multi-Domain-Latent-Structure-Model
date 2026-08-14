# TableS3_MimicPPMI.R
# Produces Table S3: Mimic-PPMI simulation with misspecified G (3-7)
# Input: ../output/MimicPPMI_G_{G}.Rdata (from Simulation_MimicPPMI.R)
#        Fitting_PPMI.RData

rm(list = ls())
library(clue)
library(mclust)

load("Fitting_PPMI.RData")

s_real <- Fitting_PPMI$mode_S
J <- length(s_real)
G_true <- length(unique(s_real))
beta_real <- Fitting_PPMI$mean_beta 
sigma_real <- Fitting_PPMI$mean_SigmaG
theta_real <- Fitting_PPMI$mean_Theta

G_values <- c(3, 4, 5, 6, 7)
rdata_files <- paste0("../output/MimicPPMI_G_", G_values, ".Rdata")

results_list <- vector("list", length(G_values))
names(results_list) <- paste0("G", G_values)

for (i in seq_along(G_values)) {
  if (!file.exists(rdata_files[i])) {
    message("File not found: ", rdata_files[i])
    next
  }
  env <- new.env()
  load(rdata_files[i], envir = env)
  results_list[[i]] <- env$result
  cat(sprintf("Loading %s  (DS = %d)\n", rdata_files[i], length(env$result)))
}

align_labels_flexible <- function(s_real, s_est) {
  G_true <- length(unique(s_real))
  G_est <- length(unique(s_est))
  confusion_matrix <- table(Real = s_real, Est = s_est)

  if (G_est == G_true) {
    match_result <- as.vector(solve_LSAP(max(confusion_matrix) - confusion_matrix))
    mapping <- numeric(G_est)
    mapping[match_result] <- 1:G_true
    names(mapping) <- 1:G_est
  } else if (G_est < G_true) {
    mapping <- numeric(G_est)
    names(mapping) <- 1:G_est
    for (g_est in 1:G_est) {
      mapping[g_est] <- which.max(confusion_matrix[, g_est])
    }
  } else {
    mapping <- numeric(G_est)
    names(mapping) <- 1:G_est
    for (g_est in 1:G_est) {
      mapping[g_est] <- which.max(confusion_matrix[, g_est])
    }
  }

  s_aligned <- numeric(length(s_est))
  for (i in 1:length(s_est)) {
    s_aligned[i] <- mapping[s_est[i]]
  }

  return(list(
    s_aligned        = s_aligned,
    mapping          = mapping,
    confusion_matrix = confusion_matrix
  ))
}

compute_beta_rmse <- function(beta_real, beta_est, mapping, s_real, s_est) {
  G_true <- ncol(beta_real)
  all_diffs <- c()
  for (g_true in 1:G_true) {
    est_groups <- which(mapping == g_true)
    if (length(est_groups) == 0) next
    beta_est_avg <- if (length(est_groups) > 1) {
      rowMeans(beta_est[, est_groups, drop = FALSE])
    } else {
      beta_est[, est_groups]
    }
    all_diffs <- c(all_diffs, beta_est_avg - beta_real[, g_true])
  }
  sqrt(mean(all_diffs^2))
}

results_summary <- matrix(NA, nrow = length(G_values), ncol = 8)
colnames(results_summary) <- c(
  "ARI_M", "ARI_IQR", "Theta_M", "Theta_IQR",
  "Beta_M", "Beta_IQR", "Sigma_M", "Sigma_IQR"
)
rownames(results_summary) <- paste0("G=", G_values)

for (i in seq_along(G_values)) {
  result <- results_list[[i]]
  if (is.null(result)) next

  DS_cur <- length(result)

  all_ari <- numeric(DS_cur)
  all_theta_rmse <- numeric(DS_cur)
  all_beta_rmse <- numeric(DS_cur)
  all_sigma_rmse <- numeric(DS_cur)

  for (ds in 1:DS_cur) {
    res_ds <- result[[ds]]

    s_est_vec <- as.vector(res_ds$mode_S)

    all_ari[ds] <- adjustedRandIndex(s_real, s_est_vec)

    align_res <- align_labels_flexible(s_real, s_est_vec)
    mapping <- align_res$mapping

    all_theta_rmse[ds] <- sqrt(mean((res_ds$mean_Theta - theta_real)^2))

    all_beta_rmse[ds] <- compute_beta_rmse(beta_real, res_ds$mean_beta, mapping, s_real, s_est_vec)

    all_sigma_rmse[ds] <- abs(sigma_real - res_ds$mean_SigmaG)
  }

  results_summary[i, ] <- c(
    median(all_ari),        IQR(all_ari),
    median(all_theta_rmse), IQR(all_theta_rmse),
    median(all_beta_rmse),  IQR(all_beta_rmse),
    median(all_sigma_rmse), IQR(all_sigma_rmse)
  )
}

cat(" \n")
cat("Misspecification of the number of domains:\n")
print(round(results_summary, 3))
write.csv(as.data.frame(round(results_summary, 3)), "../output/TableS3.csv")
message("Saved: ../output/TableS3.csv")

cat(" \n")
cat("G = 5 domain misclassification:")

result_G5 <- results_list[["G5"]]

if (!is.null(result_G5)) {
  DS5 <- length(result_G5)

  all_ari_G5 <- numeric(DS5)
  incorrect_cases <- c()

  for (ds in 1:DS5) {
    all_ari_G5[ds] <- adjustedRandIndex(s_real, as.vector(result_G5[[ds]]$mode_S))
    if (all_ari_G5[ds] < 1) incorrect_cases <- c(incorrect_cases, ds)
  }

  cat(sprintf(
    "\nAccuracy (ARI = 1): %.4f  (%d / %d)\n",
    mean(all_ari_G5 == 1), sum(all_ari_G5 == 1), DS5
  ))
  cat(sprintf("# of Misspecification: %d / %d\n", length(incorrect_cases), DS5))

  if (length(incorrect_cases) > 0) {
    for (ds in incorrect_cases) {
      res_ds <- result_G5[[ds]]
      s_est_vec <- as.vector(res_ds$mode_S)

      align_res <- align_labels_flexible(s_real, s_est_vec)
      mapping <- align_res$mapping

      ari_val <- all_ari_G5[ds]
      theta_val <- sqrt(mean((res_ds$mean_Theta - theta_real)^2))
      beta_val <- compute_beta_rmse(beta_real, res_ds$mean_beta, mapping, s_real, s_est_vec)
      sigma_val <- abs(sigma_real - res_ds$mean_SigmaG)

      cat(sprintf("\n  Case ds = %d\n", ds))
      cat(sprintf("    ARI        = %.2f\n", round(ari_val, 2)))
      cat(sprintf("    RMSE Theta = %.3f\n", round(theta_val, 3)))
      cat(sprintf("    RMSE Beta  = %.3f\n", round(beta_val, 3)))
      cat(sprintf("    RMSE Sigma = %.3f\n", round(sigma_val, 3)))

      cat("  Sign of Beta (-: correct, x: incorrect):\n")
      for (g_true in 1:G_true) {
        est_groups <- which(mapping == g_true)
        if (length(est_groups) == 0) next
        beta_avg <- if (length(est_groups) > 1) {
          rowMeans(res_ds$mean_beta[, est_groups, drop = FALSE])
        } else {
          res_ds$mean_beta[, est_groups]
        }
        sign_match <- ifelse(sign(beta_avg) == sign(beta_real[, g_true]), "-", "x")
        cat(sprintf("      Group %d: [%s]\n", g_true, paste(sign_match, collapse = ", ")))
      }
    }
  } else {
    cat("All cases are correctly classified (ARI = 1)\n")
  }
}
