# TableS5_EstimationStrategies.R
# Produces Table S5: Clustering accuracy comparison (Additional vs Suggestion)
# Input: ../output/Additional_{s}+{CR}+{N}+2.Rdata (from Simulation_Addition.R)
#        ../output/Suggestion_{s}+{CR}+{N}+2.Rdata (from Simulation_Suggestion.R)

library(clue)
library(mclust)

base_path   <- "../output"
DS          <- 100
G           <- 3
J           <- 15

s_real <- rep(1:G, each = J / G)

settings_t2 <- list(
    list(s = 1, cen = "20%", n = 200),
    list(s = 1, cen = "20%", n = 500),
    list(s = 1, cen = "40%", n = 200),
    list(s = 1, cen = "40%", n = 500),
    list(s = 4, cen = "20%", n = 200),
    list(s = 4, cen = "20%", n = 500),
    list(s = 4, cen = "40%", n = 200),
    list(s = 4, cen = "40%", n = 500)
)

compute_acc <- function(result) {
    S_ARI <- sapply(1:DS, function(d) {
        if (is.null(result[[d]]$mode_S)) return(0)
        adjustedRandIndex(result[[d]]$mode_S, s_real)
    })
    mean(S_ARI == 1)
}

n_set  <- length(settings_t2)
rnames <- character(n_set)

A_acc <- numeric(n_set)
S_acc <- numeric(n_set)
A_ok  <- S_ok  <- logical(n_set)

for (idx in seq_along(settings_t2)) {
    cfg <- settings_t2[[idx]]
    rnames[idx] <- sprintf("S%d, CR=%s, N=%d", cfg$s, cfg$cen, cfg$n)

    file_A <- file.path(base_path,
        paste0("Additional_", cfg$s, "+", cfg$cen, "+", cfg$n, "+2.Rdata"))
    file_S <- file.path(base_path,
        paste0("Suggestion_", cfg$s, "+", cfg$cen, "+", cfg$n, "+2.Rdata"))

    if (file.exists(file_A)) {
        load(file_A)
        A_acc[idx] <- compute_acc(result)
        A_ok[idx]  <- TRUE
    }

    if (file.exists(file_S)) {
        load(file_S)
        S_acc[idx] <- compute_acc(result)
        S_ok[idx]  <- TRUE
    }
}

f1 <- function(x) formatC(round(x * 100, 1), format = "f", digits = 1)
up <- function(d) ifelse(d > 0.0005, "+", ifelse(d < -0.0005, "-", "="))

tbl <- data.frame(
    Setting         = rnames,
    Add_Acc         = ifelse(A_ok, f1(A_acc), "NA"),
    Sug_Acc         = ifelse(S_ok, f1(S_acc), "NA"),
    Delta           = ifelse(A_ok & S_ok,
                             paste0(f1(S_acc - A_acc), "%  ", up(S_acc - A_acc)), "NA"),
    stringsAsFactors = FALSE
)
colnames(tbl) <- c("Setting", "Baseline Accuracy (%)", "Suggestion Accuracy (%)", "Delta Accuracy (%)")

cat("\n====== Table S5: Clustering Accuracy Comparison ======\n")
cat("Accuracy = Proportion of 100 replications with perfect grouping (ARI == 1)\n")
cat("Delta = Suggestion - Baseline  (+/- = improved / worsened)\n\n")
print(tbl, row.names = FALSE)
write.csv(tbl, "../output/TableS5.csv", row.names = FALSE)
message("Saved: ../output/TableS5.csv")

valid <- A_ok & S_ok
if (any(valid)) {
    cat(sprintf(
        "\nMean across %d settings:\n  Baseline:   %.1f%%\n  Suggestion: %.1f%%\n  Mean Delta: %.1f%%\n",
        sum(valid),
        mean(A_acc[valid]) * 100, mean(S_acc[valid]) * 100,
        mean(S_acc[valid] - A_acc[valid]) * 100
    ))
    cat(sprintf("  Settings with improved accuracy: %d / %d\n",
                sum((S_acc - A_acc)[valid] > 0.005), sum(valid)))
}
