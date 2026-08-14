# TableS2_ModelSelection.R
# Produces Table S2: Model selection via AIC/BIC/WAIC
# Input: ../output/ModelSelection_{setting}_N{N}.Rdata (from Simulation_ModelSelection.R)

rm(list = ls())

library(dplyr)
library(tidyr)

setting_info <- list(
  list(label = "K=15, G=3", G_true = 3, file_N100 = "../output/ModelSelection_1_N100.Rdata", file_N200 = "../output/ModelSelection_1_N200.Rdata"),
  list(label = "K=25, G=5", G_true = 5, file_N100 = "../output/ModelSelection_2_N100.Rdata", file_N200 = "../output/ModelSelection_2_N200.Rdata"),
  list(label = "K=35, G=5", G_true = 5, file_N100 = "../output/ModelSelection_3_N100.Rdata", file_N200 = "../output/ModelSelection_3_N200.Rdata")
)

count_correct <- function(result, G_true) {
  if (G_true == 3) {
    G_candidates <- 2:6
    correct_idx <- 2
  } else if (G_true == 5) {
    G_candidates <- 3:7
    correct_idx <- 3
  } else {
    stop("only for G_true=3 or G_true=5")
  }

  n_candidates <- length(G_candidates)
  DS <- length(result)

  AIC_matrix  <- matrix(NA, n_candidates, DS)
  BIC_matrix  <- matrix(NA, n_candidates, DS)
  WAIC_matrix <- matrix(NA, n_candidates, DS)

  for (ds in 1:DS) {
    for (g_idx in 1:n_candidates) {
      g_name <- paste0("G", G_candidates[g_idx])
      res_g <- result[[ds]][[g_name]]
      if (!is.null(res_g)) {
        AIC_matrix[g_idx, ds]  <- ifelse(is.null(res_g$AIC),  NA, res_g$AIC)
        BIC_matrix[g_idx, ds]  <- ifelse(is.null(res_g$BIC),  NA, res_g$BIC)
        WAIC_matrix[g_idx, ds] <- ifelse(is.null(res_g$WAIC), NA, res_g$WAIC)
      }
    }
  }

  sel <- function(mat) apply(mat, 2, function(col) if (all(is.na(col))) NA else which.min(col))
  AIC_sel  <- sel(AIC_matrix)
  BIC_sel  <- sel(BIC_matrix)
  WAIC_sel <- sel(WAIC_matrix)

  c(
    AIC  = sum(AIC_sel  == correct_idx, na.rm = TRUE),
    BIC  = sum(BIC_sel  == correct_idx, na.rm = TRUE),
    WAIC = sum(WAIC_sel == correct_idx, na.rm = TRUE)
  )
}

rows <- vector("list", length(setting_info))

for (i in seq_along(setting_info)) {
  s <- setting_info[[i]]

  if (file.exists(s$file_N100)) {
    load(s$file_N100)
    cnt100 <- count_correct(result, s$G_true)
  } else {
    warning(sprintf("file not found: %s", s$file_N100))
    cnt100 <- c(AIC = NA, BIC = NA, WAIC = NA)
  }

  if (file.exists(s$file_N200)) {
    load(s$file_N200)
    cnt200 <- count_correct(result, s$G_true)
  } else {
    warning(sprintf("file not found: %s", s$file_N200))
    cnt200 <- c(AIC = NA, BIC = NA, WAIC = NA)
  }

  rows[[i]] <- data.frame(
    Setting        = s$label,
    N100_AIC       = cnt100["AIC"],
    N100_BIC       = cnt100["BIC"],
    N100_WAIC      = cnt100["WAIC"],
    N200_AIC       = cnt200["AIC"],
    N200_BIC       = cnt200["BIC"],
    N200_WAIC      = cnt200["WAIC"],
    row.names      = NULL
  )
}

table_s2 <- do.call(rbind, rows)
colnames(table_s2) <- c("Setting",
                         "AIC (100)", "BIC (100)", "WAIC (100)",
                         "AIC (200)", "BIC (200)", "WAIC (200)")
print(table_s2, row.names = FALSE)
write.csv(table_s2, "../output/TableS2.csv", row.names = FALSE)
message("Saved: ../output/TableS2.csv")
