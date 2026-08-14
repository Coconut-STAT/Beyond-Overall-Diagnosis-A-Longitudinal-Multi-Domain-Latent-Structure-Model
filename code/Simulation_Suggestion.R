# Simulation_Suggestion.R
# Clustering-prior initialization with 5 runs per replication (best log-likelihood)
# Output: ../output/Suggestion_{s}+{CR}+{N}+2.Rdata
# Usage: Rscript Simulation_Suggestion.R <1-8>

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
library(vcd)

n_core_set <- as.integer(Sys.getenv("N_CORES", unset = "102"))

seed_set = 31
DS = 100
set_para_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
n_runs_per_rep = 5

get_para_by_num = function(num) {
  all_settings = list(
    list(1, "20%", 200),
    list(1, "20%", 500),
    list(1, "40%", 200),
    list(1, "40%", 500),
    list(4, "20%", 200),
    list(4, "20%", 500),
    list(4, "40%", 200),
    list(4, "40%", 500)
  )
  return(all_settings[[num]])
}

para = get_para_by_num(set_para_num)
set_para1 = para[[1]]
set_para2 = para[[2]]
set_para3 = para[[3]]

cat(sprintf("Running Setting #%d: setting_id=%d, CR=%s, N=%d\n", 
            set_para_num, set_para1, set_para2, set_para3))

dqset.seed(seed_set) 
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(min(n_core_set, n_cores)) 
registerDoParallel(cl)
registerDoRNG(seed_set)

IR <- 4000
BI <- 2000

get_setting = function(J = 15, G = 3, lambda_label = "1", CR_label = "20%", N = 200) {
  B = 2
  if (G == 5) {
    beta <- matrix(c(1, 1, 1.5, 0.5, 1, -1, -1, -1, -1, 1), B, G, byrow = FALSE)
  } else if (G == 3) {
    beta <- matrix(c(-1, -1, 1, -1, 1, 1), B, G, byrow = FALSE)
  } else {
    stop("Invalid G value; only 3 or 5 supported.")
  }
  
  if (J / G == 5) {
    theta_vecs <- list(
      c(0.7, 0.2, 0.1, 0.5, 0.3, 0.2), 
      c(0.3, 0.6, 0.1, 0.2, 0.5, 0.3),
      c(0.6, 0.2, 0.2, 0.4, 0.3, 0.3),
      c(0.9, 0.05, 0.05, 0.7, 0.15, 0.15),
      c(0.5, 0.4, 0.1, 0.4, 0.3, 0.3)
    )
  } else {
    stop("Invalid J/G ratio; only J/G=5 supported for setting 1 and 4.")
  }
  
  if (CR_label == "20%" && lambda_label == "1") {
    baseline_nu <- 1
    expRate <- 3
  } else if (CR_label == "40%" && lambda_label == "1") {
    baseline_nu <- 1
    expRate <- 9.2
  } else if (CR_label == "60%" && lambda_label == "1") {
    baseline_nu <- 1
    expRate <- 25
  } else if (CR_label == "20%" && lambda_label == "2t") {
    baseline_nu <- 2
    expRate <- 5
  } else if (CR_label == "40%" && lambda_label == "2t") {
    baseline_nu <- 2
    expRate <- 9
  } else {
    stop("Invalid CR_label or lambda_label combination.")
  }
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
get_setting_num = function(id, CR_label = "20%", N = 200) {
  if (id == 1) {
    get_setting(J = 15, G = 3, lambda_label = "1", CR_label = CR_label, N = N)
  } else if (id == 4) {
    get_setting(J = 15, G = 3, lambda_label = "2t", CR_label = CR_label, N = N)
  } else {
    stop("Invalid id; only setting 1 and 4 are supported.")
  }
}

cfg = get_setting_num(set_para1, set_para2, set_para3)

result <- foreach(ds = 1:DS, .packages = c("Rcpp", "RcppArmadillo", "MASS", "RcppParallel", "MCMCpack", "dqrng", "vcd")) %dorng% {
  Rcpp::sourceCpp("f_MCMC.cpp")
  RcppParallel::setThreadOptions(numThreads = 1)
  
  iter_seed <- seed_set * 1000 + ds
  set.seed(iter_seed)
  dqrng::dqset.seed(iter_seed)

  {
    N = cfg$N
    G = cfg$G
    J = cfg$J
    beta = cfg$beta
    theta_vecs = cfg$theta_vecs
    baseline_nu = cfg$baseline_nu
    expRate = cfg$expRate

    B = 2
    X = matrix(0, B, N * G)
    for (g in 1:G) {
      if (g == 1) {
        for (i in 1:N) {
          X[, i + N * (g - 1)] = rnorm(B, 0, 1)
        }
      }
      else {
        for (i in 1:N) {
          X[, i + N * (g - 1)] = X[, i + N * (1 - 1)]
        }
      }
    } 
    SigmaG = 0.5
    Ui = rnorm(N, 0, sqrt(SigmaG))
    Unif_gener = matrix(runif(N * G, 0, 1), nrow = N, ncol = G)
    Tig = matrix(0, N, G)
    lambda0 = 1 
    for (i in 1:N) {
      for (g in 1:G) {
        Tig[i, g] = ((-log(Unif_gener[i, g]) * exp(-t(beta[, g]) %*% X[, i + N * (g - 1)] - Ui[i] ))^(1/baseline_nu)) / lambda0
      }
    }
    Tig = round(Tig,3)
  
    lower_len = 7
    upper_len = 12
    ri = sample(lower_len:upper_len, N, replace = TRUE)
    TT = max(ri)
    tig = matrix(0, TT, N)
    for (i in 1:N) {
      tig[1,i] = rexp(1, rate = expRate)
      for (r in 2:ri[i]) {
        tig[r,i] = rexp(1, rate = expRate) + tig[r - 1,i]
      }
    }
    tig = round(tig, 3)
    
    InterL = matrix(0, N, G)
    InterR = matrix(0, N, G)
    status = matrix(0, N, G)
    R = matrix(0, N, G)
    for (i in 1:N) {
      for (g in 1:G) {
        if (findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE) == 0) {
          InterL[i,g] = 0
          InterR[i,g] = tig[1,i]
          status[i,g] = 0
          R[i,g] = 1
        }
        if (findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE) == ri[i]) {
          InterL[i,g] = tig[ri[i],i]
          InterR[i,g] = NA
          status[i,g] = 2
          R[i,g] = ri[i] + 1
        }
        if (findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE) > 0 & findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE) < ri[i]) {
          Interval_index = findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE)
          InterL[i,g] = tig[Interval_index,i]
          InterR[i,g] = tig[Interval_index + 1,i]
          status[i,g] = 1
          R[i,g] = Interval_index + 1
        }
      }
    }
    CR = rep(0, G)
    for (g in 1:G) {
      CR[g] = sum(status[,g] == 2) / (N)
    }

    Z = matrix(0, G, max(ri) * N)
    for (i in 1:N) {
      for (g in 1:G) {
        Z[g, ((1 + TT * (i - 1)):(R[i, g] - 1 + TT * (i - 1)))] = 1
        if (R[i, g] < ri[i] + 1) {
          Z[g, ((R[i, g] + TT * (i - 1)):(ri[i] + TT * (i - 1)))] = 2
        }
      }
    }
    
    L = matrix(0, J, G)
    for (g in 1:G) {
      L[((g - 1) * (J/G) + 1):(g * (J/G)), g] = 1
    }
    dj = rep(3, J) 
    Theta = matrix(0, max(dj), 2 * J)
    Theta_list = lapply(seq_along(theta_vecs), function(k) {matrix(theta_vecs[[k]], nrow = dj[k], ncol = 2, byrow = FALSE)})
    Theta = do.call(cbind, rep(Theta_list, times = G))

    Y = matrix(0, J, TT * N)
    LL = which(L == 1, arr.ind = TRUE)
    for (i in 1:N) {
      for (j in 1:J) {
        g = unname(LL[which(LL[, 1] == j), 2])
        for (t in 1:ri[i]) {
          Y[j, t + TT * (i - 1)] = sample(1:dj[j],1,TRUE,Theta[, Z[g, t + TT * (i - 1)] + 2 * (j - 1)])
        }
      }
    }
    
    s = numeric(J)
    for (j in 1:J) {
      for (g in 1:G) {
        if (L[j, g] == 1)
          s[j] = g
      }
    }
    
    Yijtc = matrix(0, max(ri) * J, max(dj) * N)
    for (i in 1:N) {
      for (j in 1:J) {
        for (t in 1:ri[i]) {
          for (c in 1:dj[j]) {
            if (Y[j, t + TT * (i - 1)] == c) {
              Yijtc[t + max(ri) * (j - 1), c + max(dj) * (i - 1)] = 1
            } else{
              Yijtc[t + max(ri) * (j - 1), c + max(dj) * (i - 1)] = 0
            }
          }
        }
      }
    }
    
    beta_real = beta
    R_real = R
    L_real = L
    Theta_real = Theta
    status_real = status
    InterL_real = InterL
    InterR_real = InterR
    Ui_real = Ui
    SigmaG_real = SigmaG
    s_real = s
    X_real = X
  }

  {
    Y_transposed <- t(Y)
    cor_matrix <- matrix(1, J, J)
    for (i in 1:(J - 1)) {
      for (j in (i + 1):J) {
        tab <- table(Y_transposed[, i], Y_transposed[, j])
        v_val <- tryCatch(assocstats(tab)$cramer, error = function(e) 0)
        cor_matrix[i, j] <- cor_matrix[j, i] <- v_val
      }
    }
    dist_matrix  <- as.dist(1 - cor_matrix)
    hc           <- hclust(dist_matrix, method = "complete")
    prior_groups <- cutree(hc, k = G)

    strength     <- 0.9
    prob_matrix  <- matrix((1 - strength) / (G - 1), nrow = J, ncol = G)
    for (j in 1:J) prob_matrix[j, prior_groups[j]] <- strength
  }

  all_loglik <- numeric(n_runs_per_rep)
  first_result <- NULL
  best_result <- NULL
  best_loglik <- -Inf
  
  for (run_id in 1:n_runs_per_rep) {
    run_seed <- iter_seed * 100 + run_id
    set.seed(run_seed)
    dqrng::dqset.seed(run_seed)
    
    {
      repeat {
        s <- numeric(J)
        for (j in 1:J) s[j] <- sample(1:G, 1, prob = prob_matrix[j, ])
        if (length(unique(s)) == G && all(table(s) >= 2)) break
      }
      
      for (j in 1:J) {
        for (k in 1:2) {
          repeat {
            sample = (rdirichlet(1, rep(1, dj[j])))
            if (any(sample < 0.1))
              next
            Theta[1:dj[j], k + 2 * (j - 1)] = sample
            if (k == 2)
              if (Theta[1, 1 + 2 * (j - 1)] < Theta[1, 2 + 2 * (j - 1)])
                next
            break
          } 
        }
      }
      
      order = 3
      sig0 = 2
      coef_range_beta = 10
      coef_range_ui = 10
      a_eta = 0.1
      b_eta = 0.1
      LengthVKnots = 20
      knots <- seq(min(tig), max(tig), length = LengthVKnots)
      Lknots <- LengthVKnots - 2 + order

      eta = rgamma(G, a_eta, rate = b_eta)
      beta = matrix(0, B, G)
      gammaL = matrix(rgamma(Lknots * G, 1, 1), nrow = Lknots, ncol = G)
      
      sigma_mean = 1
      sigma_variance = 0.5^2
      sigma_alpha = sigma_mean^2 / sigma_variance + 2   
      sigma_beta = sigma_mean * (sigma_alpha - 1)
      sigma_beta = sigma_alpha
      SigmaG = rinvgamma(1, sigma_alpha, sigma_beta)
      Ui = rnorm(N, 0, sqrt(SigmaG))
    } 

    res <- f_mcmc(IR, BI, N, J, G, B, order, Lknots, ri, dj, Y, s, knots, Theta, Yijtc, beta, gammaL, eta, X, tig, Ui, sig0, coef_range_beta, coef_range_ui, a_eta, b_eta, SigmaG, sigma_alpha, sigma_beta, ds) 
    
    current_loglik <- res$LogLikelihood
    all_loglik[run_id] <- current_loglik
    
    res$Ui_e <- NULL        
    res$gammaL_e <- NULL        
    res$PRigt_e <- NULL 
    res$beta_e <- NULL        
    res$Theta_e <- NULL        
    res$SigmaG_e <- NULL 
    res$s_e <- NULL        
    res$R_real <- R_real
    
    if (current_loglik > best_loglik) {
      best_loglik <- current_loglik
      best_result <- res
      best_result$best_run_id <- run_id
    }
  }

  best_result$all_loglik    <- all_loglik
  best_result$prior_groups  <- prior_groups

  return(best_result)
}

stopCluster(cl)

out_file <- file.path("../output", paste0("Suggestion_", set_para1, "+", set_para2, "+", set_para3, "+2.Rdata"))
save(result, file = out_file)

cat("\n=== Summary ===\n")

cat(sprintf(
    "Done. Saved %s\n",
    out_file
))
