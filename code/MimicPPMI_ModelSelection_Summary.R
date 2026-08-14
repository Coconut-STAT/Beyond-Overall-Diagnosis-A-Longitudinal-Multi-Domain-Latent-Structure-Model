# MimicPPMI_ModelSelection_Summary.R
# Produces the Mimic-PPMI model-selection summary reported in the text
# (this is not a numbered table in the paper)
# Input: ../output/MimicPPMI_ModelSelection.Rdata (from Simulation_MimicPPMI_ModelSelection.R)

rm(list = ls())

library(dplyr)
library(tidyr)

G_true       <- 5
G_candidates <- 3:7
correct_idx  <- which(G_candidates == G_true)

data_file <- "../output/MimicPPMI_ModelSelection.Rdata"

if (!file.exists(data_file)) {
  stop(sprintf("File not found: %s", data_file))
}

load(data_file)

DS             <- length(result_all)
n_candidates   <- length(G_candidates)

AIC_matrix  <- matrix(NA, n_candidates, DS)
BIC_matrix  <- matrix(NA, n_candidates, DS)
WAIC_matrix <- matrix(NA, n_candidates, DS)

for (ds in 1:DS) {
  for (g_idx in 1:n_candidates) {
    g_name <- paste0("G_", G_candidates[g_idx])
    res_g  <- result_all[[ds]][[g_name]]
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

cnt <- c(
  AIC  = sum(AIC_sel  == correct_idx, na.rm = TRUE),
  BIC  = sum(BIC_sel  == correct_idx, na.rm = TRUE),
  WAIC = sum(WAIC_sel == correct_idx, na.rm = TRUE)
)

summary_table <- data.frame(
  Setting      = sprintf("G_true=%d, N=849, DS=%d", G_true, DS),
  AIC_correct  = cnt["AIC"],
  BIC_correct  = cnt["BIC"],
  WAIC_correct = cnt["WAIC"],
  row.names    = NULL
)

colnames(summary_table) <- c("Setting", "AIC", "BIC", "WAIC")

cat("\n===== Mimic-PPMI Model Selection Summary (correct G selections out of", DS, ") =====\n\n")
print(summary_table, row.names = FALSE)

cat("\n--- Detailed selection frequency table ---\n")
freq_df <- data.frame(
  G_candidate = G_candidates,
  AIC_freq    = tabulate(AIC_sel,  nbins = n_candidates),
  BIC_freq    = tabulate(BIC_sel,  nbins = n_candidates),
  WAIC_freq   = tabulate(WAIC_sel, nbins = n_candidates)
)
freq_df$G_candidate <- paste0("G=", freq_df$G_candidate)
freq_df$G_candidate[correct_idx] <- paste0(freq_df$G_candidate[correct_idx], " (true)")
colnames(freq_df) <- c("G", "AIC", "BIC", "WAIC")
print(freq_df, row.names = FALSE)
write.csv(freq_df, "../output/MimicPPMI_ModelSelection_Summary.csv", row.names = FALSE)
message("Saved: ../output/MimicPPMI_ModelSelection_Summary.csv")
