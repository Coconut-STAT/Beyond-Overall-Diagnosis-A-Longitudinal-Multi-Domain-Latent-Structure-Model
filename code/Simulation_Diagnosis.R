# Simulation_Diagnosis.R
# MCMC diagnostics: baseline function estimation + convergence diagnostics
# Part 1: 100 replications (Setting 1, CR=20%, N=200) -> baseline survival/hazard plots
# Part 2: 3 chains (Setting 1, CR=20%, N=500) -> traceplots + Gelman-Rubin
# Output: ../output/diag_part1+20%+200.Rdata, ../output/diag_part2+20%+500.Rdata,
#         ../output/Figure_S1.pdf, ../output/Figure_S2.pdf
# Usage: Rscript Simulation_Diagnosis.R

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
library(mclust)
library(clue)
library(coda)
library(ggplot2)
library(patchwork)

seedss = 2026 
theme_diag <- function(base_size = 9) {
  theme_minimal(base_size = base_size) %+replace% theme(
    axis.text         = element_text(size = base_size,     color = "black", face = "bold"),
    axis.title        = element_text(size = base_size + 1, face = "bold"),
    axis.ticks        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.15, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position   = "top",
    legend.text       = element_text(size = base_size, face = "bold"),
    legend.margin     = margin(t = 0, b = -4),
    plot.margin       = margin(4, 6, 3, 4)
  )
}

CHAIN_COLS <- c("Chain 1" = "#237B9F", "Chain 2" = "#E07B39", "Chain 3" = "#5BB55A")

LengthVKnots <- 20; order <- 3
sig0 <- 10; coef_range_beta <- 10; coef_range_ui <- 10
a_eta <- 0.1; b_eta <- 0.1

sigma_mean <- 1; sigma_variance <- 0.5^2
sigma_alpha <- sigma_mean^2 / sigma_variance + 2
sigma_beta  <- sigma_mean * (sigma_alpha - 1)

build_data <- function(cfg, iter_seed) {
  set.seed(iter_seed); dqrng::dqset.seed(iter_seed)
  N <- cfg$N; G <- cfg$G; J <- cfg$J
  beta <- cfg$beta; theta_vecs <- cfg$theta_vecs
  baseline_nu <- cfg$baseline_nu; expRate <- cfg$expRate
  B <- 2

  X <- matrix(0, B, N * G)
  for (i in 1:N) X[, i] <- rnorm(B, 0, 1)
  for (g in 2:G) for (i in 1:N) X[, i + N*(g-1)] <- X[, i]

  SigmaG <- 0.5
  Ui     <- rnorm(N, 0, sqrt(SigmaG))
  Unif   <- matrix(runif(N * G), N, G)
  lambda0 <- 1
  Tig <- matrix(0, N, G)
  for (i in 1:N)
    for (g in 1:G)
      Tig[i, g] <- ((-log(Unif[i,g]) *
                       exp(-t(beta[,g]) %*% X[, i+N*(g-1)] - Ui[i]))^(1/baseline_nu)) / lambda0
  Tig <- round(Tig, 4)

  ri  <- sample(7:12, N, replace = TRUE)
  TT  <- max(ri)
  tig <- matrix(0, TT, N)
  for (i in 1:N) {
    tig[1, i] <- rexp(1, rate = expRate)
    for (r in 2:ri[i]) tig[r, i] <- rexp(1, rate = expRate) + tig[r-1, i]
  }
  tig <- round(tig, 4)

  InterL <- matrix(0, N, G); InterR <- matrix(0, N, G)
  status <- matrix(0, N, G); R <- matrix(0, N, G)
  for (i in 1:N) for (g in 1:G) {
    idx <- findInterval(Tig[i,g], vec = tig[1:ri[i],i], rightmost.closed = TRUE)
    if (idx == 0) {
      InterL[i,g] <- 0; InterR[i,g] <- tig[1,i]; status[i,g] <- 0; R[i,g] <- 1
    } else if (idx == ri[i]) {
      InterL[i,g] <- tig[ri[i],i]; InterR[i,g] <- NA; status[i,g] <- 2; R[i,g] <- ri[i]+1
    } else {
      InterL[i,g] <- tig[idx,i]; InterR[i,g] <- tig[idx+1,i]; status[i,g] <- 1; R[i,g] <- idx+1
    }
  }

  Z <- matrix(0, G, TT*N)
  for (i in 1:N) for (g in 1:G) {
    Z[g, (1+TT*(i-1)):(R[i,g]-1+TT*(i-1))] <- 1
    if (R[i,g] < ri[i]+1) Z[g, (R[i,g]+TT*(i-1)):(ri[i]+TT*(i-1))] <- 2
  }

  L <- matrix(0, J, G)
  for (g in 1:G) L[((g-1)*(J/G)+1):(g*(J/G)), g] <- 1
  dj        <- rep(3, J)
  Theta_list <- lapply(seq_along(theta_vecs), function(k)
    matrix(theta_vecs[[k]], nrow = 3, ncol = 2, byrow = FALSE))
  Theta <- do.call(cbind, rep(Theta_list, times = G))

  Y  <- matrix(0, J, TT*N)
  LL <- which(L == 1, arr.ind = TRUE)
  for (i in 1:N) for (j in 1:J) {
    gj <- LL[LL[,1] == j, 2]
    for (t in 1:ri[i])
      Y[j, t+TT*(i-1)] <- sample(1:3, 1, TRUE, Theta[, Z[gj, t+TT*(i-1)] + 2*(j-1)])
  }

  s_real <- numeric(J)
  for (j in 1:J) for (g in 1:G) if (L[j,g]==1) s_real[j] <- g

  Yijtc <- matrix(0, max(ri)*J, max(dj)*N)
  for (i in 1:N) for (j in 1:J) for (t in 1:ri[i]) for (c in 1:dj[j])
    if (Y[j, t+TT*(i-1)] == c) Yijtc[t+max(ri)*(j-1), c+max(dj)*(i-1)] <- 1

  list(N=N, G=G, J=J, B=B, dj=dj, TT=TT, ri=ri, tig=tig, Y=Y, Yijtc=Yijtc,
       Theta=Theta, L=L, X=X, s_real=s_real, R_real=R, Tig=Tig,
       InterR=InterR,
       SigmaG_real=SigmaG, beta_real=beta, Ui=Ui, lambda0=lambda0)
}

rand_init <- function(d, G, order, LengthVKnots, a_eta, b_eta,
                      sigma_alpha, sigma_beta) {
  repeat {
    XI <- rdirichlet(1, rep(1, G)); s <- sample(1:G, d$J, TRUE, XI)
    if (length(unique(s)) == G && all(table(s) >= 2)) break
  }
  Theta <- d$Theta
  for (j in 1:d$J) for (k in 1:2) repeat {
    samp <- as.vector(rdirichlet(1, rep(1, d$dj[j])))
    if (any(samp < 0.1)) next
    Theta[1:d$dj[j], k+2*(j-1)] <- samp
    if (k == 2 && Theta[1, 1+2*(j-1)] < Theta[1, 2+2*(j-1)]) next
    break
  }
  Lknots <- LengthVKnots - 2 + order
  eta    <- rgamma(G, a_eta, rate = b_eta)
  beta   <- matrix(runif(d$B * G, -1, 1), d$B, G)
  gammaL <- matrix(rgamma(Lknots * G, 1, 1), Lknots, G)
  SigmaG <- rinvgamma(1, sigma_alpha, sigma_beta)
  Ui     <- rnorm(d$N, 0, sqrt(SigmaG))
  list(s=s, Theta=Theta, eta=eta, beta=beta, gammaL=gammaL, SigmaG=SigmaG, Ui=Ui)
}

extract_beta <- function(beta_e, b, g, G) {
  nc <- ncol(beta_e)
  beta_e[b, seq(g, nc, by = G)]
}

extract_theta <- function(Theta_e, cat, j, state, J) {
  nc <- ncol(Theta_e)
  col_within <- 2*(j-1) + state
  Theta_e[cat, seq(col_within, nc, by = 2*J)]
}

get_label_perm <- function(mode_s, s_real, G) {
  conf <- table(factor(mode_s, levels = 1:G), factor(s_real, levels = 1:G))
  perm <- as.integer(clue::solve_LSAP(conf, maximum = TRUE))
  order(perm)
}

seed_p1 <- seedss
DS_p1   <- 100
IR_p1   <- 10000
BI_p1   <- 5000

cfg_p1 <- list(
  N = 200, G = 3, J = 15,
  beta = matrix(c(-1, -1, 1, -1, 1, 1), 2, 3, byrow = FALSE),
  theta_vecs = list(
    c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
    c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
    c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
    c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
    c(0.8, 0.1,  0.1, 0.2, 0.6, 0.2)
  ),
  baseline_nu = 1, expRate = 3.0
)

seed_p2 <- seedss
IR_p2   <- 10000
BI_p2   <- 5000

cfg_p2 <- list(
  N = 500, G = 3, J = 15,
  beta = matrix(c(-1, -1, 1, -1, 1, 1), 2, 3, byrow = FALSE),
  theta_vecs = list(
    c(0.7, 0.2, 0.1, 0.1, 0.8, 0.1),
    c(0.3, 0.6, 0.1, 0.1, 0.2, 0.7),
    c(0.6, 0.2, 0.2, 0.05, 0.05, 0.9),
    c(0.9, 0.05, 0.05, 0.6, 0.3, 0.1),
    c(0.8, 0.1,  0.1, 0.2, 0.6, 0.2)
  ),
  baseline_nu = 1, expRate = 3.0
)

p1_file <- "../output/diag_part1+20%+200.Rdata"
p2_file <- "../output/diag_part2+20%+500.Rdata"

d_ref_p1 <- build_data(cfg_p1, iter_seed = seed_p1 * 1000)
knots_p1  <- seq(0, max(d_ref_p1$tig), length = LengthVKnots)
Lknots_p1 <- LengthVKnots - 2 + order

d_p2         <- build_data(cfg_p2, iter_seed = seed_p2 * 1000)
knots_p2     <- seq(min(d_p2$tig), max(d_p2$tig), length = LengthVKnots)
Lknots_p2    <- LengthVKnots - 2 + order
chain_seeds  <- c(seed_p2, seed_p2 + 500, seed_p2 + 1000)

chain_inits <- lapply(1:3, function(ch) {
  set.seed(chain_seeds[ch]); dqrng::dqset.seed(chain_seeds[ch])
  rand_init(d_p2, d_p2$G, order, LengthVKnots, a_eta, b_eta, sigma_alpha, sigma_beta)
})

  n_jobs    <- DS_p1 + 3
  n_core_max <- as.integer(Sys.getenv("N_CORES", unset = "105"))
  n_cores   <- min(n_jobs + 2, n_core_max, parallel::detectCores() - 1)

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(seed_p1)

  all_res <- foreach(
    job = 1:n_jobs,
    .packages = c("Rcpp","RcppArmadillo","MASS","RcppParallel","MCMCpack","dqrng")
  ) %dorng% {
    Rcpp::sourceCpp("f_MCMC.cpp")
    RcppParallel::setThreadOptions(numThreads = 1)

    if (job <= DS_p1) {
      ds        <- job
      iter_seed <- seed_p1 * 1000 + ds
      d         <- build_data(cfg_p1, iter_seed)
      set.seed(iter_seed + 1); dqrng::dqset.seed(iter_seed + 1)
      init      <- rand_init(d, d$G, order, LengthVKnots, a_eta, b_eta,
                             sigma_alpha, sigma_beta)
      set.seed(iter_seed + 2); dqrng::dqset.seed(iter_seed + 2)
      res <- f_mcmc(IR_p1, BI_p1, d$N, d$J, d$G, d$B, order, Lknots_p1,
                    d$ri, d$dj, d$Y, init$s, knots_p1, init$Theta, d$Yijtc,
                    init$beta, init$gammaL, init$eta, d$X, d$tig, init$Ui,
                    sig0, coef_range_beta, coef_range_ui, a_eta, b_eta,
                    init$SigmaG, sigma_alpha, sigma_beta, ds)
      list(part = 1L, mean_gammaL = res$mean_gammaL,
           mode_S = res$mode_S, s_real = d$s_real)

    } else {
      ch <- job - DS_p1
      set.seed(chain_seeds[ch]); dqrng::dqset.seed(chain_seeds[ch])
      init_ch <- chain_inits[[ch]]
      res_ch <- f_mcmc(IR_p2, BI_p2, d_p2$N, d_p2$J, d_p2$G, d_p2$B,
                       order, Lknots_p2,
                       d_p2$ri, d_p2$dj, d_p2$Y, init_ch$s, knots_p2,
                       init_ch$Theta, d_p2$Yijtc,
                       init_ch$beta, init_ch$gammaL, init_ch$eta,
                       d_p2$X, d_p2$tig, init_ch$Ui,
                       sig0, coef_range_beta, coef_range_ui, a_eta, b_eta,
                       init_ch$SigmaG, sigma_alpha, sigma_beta, ch)
      res_ch$Ui_e     <- NULL
      res_ch$gammaL_e <- NULL
      res_ch$PRigt_e  <- NULL
      list(part = 2L, chain = ch, res = res_ch)
    }
  }
  stopCluster(cl)

  result_p1 <- all_res[1:DS_p1]
  save(result_p1, file = p1_file)

  chains_p2 <- lapply(all_res[(DS_p1+1):n_jobs], function(x) x$res)
  save(chains_p2, d_p2, file = p2_file)

G_p1   <- cfg_p1$G
Lgrids <- 1000

valid_mask  <- sapply(result_p1, function(r) !is.null(r$mean_gammaL))
gammaL_list <- lapply(result_p1[valid_mask], function(r) {
  inv_perm <- get_label_perm(r$mode_S, r$s_real, G_p1)
  r$mean_gammaL[, inv_perm]
})
gammaL_arr  <- array(unlist(gammaL_list),
                     dim = c(dim(gammaL_list[[1]]), length(gammaL_list)))
gammaL_avg  <- apply(gammaL_arr, c(1, 2), median)
t_max_p <- 6
grids   <- seq(0, t_max_p, length.out = Lgrids)

cl_tmp <- makeCluster(1); registerDoParallel(cl_tmp)
spl <- foreach(dummy = 1:1,
               .packages = c("Rcpp","RcppArmadillo")) %dopar% {
  Rcpp::sourceCpp("f_MCMC.cpp")
  list(Isp = f_Ispline(order, grids, knots_p1))
}
stopCluster(cl_tmp)
Isp <- spl[[1]]$Isp

surv0 <- Haz0 <- matrix(0, G_p1, Lgrids)
for (g in 1:G_p1) {
  Haz0[g,]  <- gammaL_avg[,g] %*% Isp
  surv0[g,] <- exp(-Haz0[g,])
}
true_surv <- exp(-grids); true_cumh <- grids

df_p1 <- do.call(rbind, lapply(1:G_p1, function(g) {
  rbind(
    data.frame(Time=grids, Est=surv0[g,], True=true_surv,
               Func="Baseline Survival",          Domain=paste0("Domain ",g)),
    data.frame(Time=grids, Est=Haz0[g,],  True=true_cumh,
               Func="Cumulative Baseline Hazard", Domain=paste0("Domain ",g))
  )
}))
df_p1$Func   <- factor(df_p1$Func,
  levels = c("Baseline Survival","Cumulative Baseline Hazard"))
df_p1$Domain <- factor(df_p1$Domain, levels = paste0("Domain ", 1:G_p1))

all_funcs    <- levels(df_p1$Func)
n_funcs      <- length(all_funcs)
plot_p1_list <- lapply(seq_along(all_funcs), function(fi) {
  fn       <- all_funcs[fi]
  sub      <- df_p1[df_p1$Func == fn, ]
  show_leg <- if (fi == 1) "top" else "none"
  x_lab    <- if (fi == n_funcs) "Time" else NULL

  sub$Domain <- factor(sub$Domain, levels = levels(df_p1$Domain))

  ggplot(sub, aes(x = Time)) +
    geom_line(aes(y = Est,  color = "Average Estimate", linetype = "Average Estimate"),
              linewidth = 1.0) +
    geom_line(aes(y = True, color = "True Function",    linetype = "True Function"),
              linewidth = 1.0) +
    scale_color_manual(
      name   = "",
      values = c("Average Estimate" = "orange", "True Function" = "#2CA02C")
    ) +
    scale_linetype_manual(
      name   = "",
      values = c("Average Estimate" = "solid", "True Function" = "dotted")
    ) +
    coord_cartesian(xlim = c(0, 6)) +
    labs(x = x_lab, y = fn) +
    facet_wrap(~ Domain, nrow = 1) +
    theme_diag(9) +
    theme(
      strip.text        = element_text(face = "bold", size = 9),
      strip.background  = element_rect(fill = "grey95", color = "black", linewidth = 0.7),
      legend.key.width  = unit(1.5, "cm"),
      legend.position   = show_leg,
      plot.margin       = margin(3, 6, 2, 4)
    )
})

p1_combined <- wrap_plots(plot_p1_list, ncol = 1)

G_p2 <- cfg_p2$G; J_p2 <- cfg_p2$J
n_stored   <- length(chains_p2[[1]]$SigmaG_e)
post_start <- BI_p2 + 1

get_label_perm <- get_label_perm

traces <- lapply(chains_p2, function(res) {
  inv_perm <- get_label_perm(res$mode_S, d_p2$s_real, G_p2)
  list(
    sigma2  = res$SigmaG_e,
    beta11  = extract_beta(res$beta_e,  b=1, g=inv_perm[1], G=G_p2),
    beta32  = extract_beta(res$beta_e,  b=2, g=inv_perm[3], G=G_p2),
    theta11 = extract_theta(res$Theta_e, cat=1, j=1, state=1, J=J_p2),
    theta21 = extract_theta(res$Theta_e, cat=1, j=2, state=1, J=J_p2)
  )
})

true_vals <- list(sigma2=0.5, beta11=-1, beta32=1, theta11=0.7, theta21=0.3)
param_labels <- list(
  sigma2  = expression(sigma^2),
  beta11  = expression(beta[11]),
  beta32  = expression(beta[32]),
  theta11 = expression(theta["1,1,1"]),
  theta21 = expression(theta["2,1,1"])
)

fixed_ylim <- function(pn) {
  if (pn == "sigma2")             c(0, 2)
  else if (pn %in% c("beta11", "beta32")) c(-2, 2)
  else                            c(0, 1)
}

make_pair <- function(pn, lab_full, lab_post, show_legend = FALSE, show_xlab = FALSE) {
  n_iter <- length(traces[[1]][[pn]])
  df_all <- do.call(rbind, lapply(1:3, function(ch)
    data.frame(iter  = 1:n_iter, value = traces[[ch]][[pn]],
               Chain = factor(paste("Chain", ch), levels = paste("Chain", 1:3)))))
  yl     <- fixed_ylim(pn)
  tv     <- true_vals[[pn]]
  yl_lab <- param_labels[[pn]]
  shade  <- data.frame(xmin=1, xmax=post_start, ymin=-Inf, ymax=Inf)
  leg    <- if (show_legend) "top" else "none"
  xl_full <- if (show_xlab) "Iteration" else NULL
  xl_post <- if (show_xlab) "Post Burn-in Iteration" else NULL

  p_full <- ggplot(df_all, aes(x=iter, y=value, color=Chain)) +
    geom_rect(data=shade, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill="grey92", alpha=0.5, inherit.aes=FALSE) +
    geom_line(linewidth=0.35, alpha=0.85) +
    geom_hline(yintercept=tv, linetype="dashed", color="black", linewidth=0.7) +
    scale_color_manual(name="", values=CHAIN_COLS) +
    coord_cartesian(xlim=c(1, n_iter), ylim=yl) +
    labs(x=xl_full, y=yl_lab, title=lab_full) +
    theme_diag(8) +
    theme(plot.title      = element_text(face="bold", size=9, hjust=0),
          axis.title.y    = element_text(size=9),
          axis.text       = element_text(size=8),
          legend.position = leg,
          legend.text     = element_text(size=8))

  df_post        <- df_all[df_all$iter > post_start, ]
  df_post$iter_p <- df_post$iter - post_start

  p_post <- ggplot(df_post, aes(x=iter_p, y=value, color=Chain)) +
    geom_line(linewidth=0.35, alpha=0.85) +
    geom_hline(yintercept=tv, linetype="dashed", color="black", linewidth=0.7) +
    scale_color_manual(name="", values=CHAIN_COLS) +
    coord_cartesian(xlim=c(0, n_iter-post_start), ylim=yl) +
    labs(x=xl_post, y=yl_lab, title=lab_post) +
    theme_diag(8) +
    theme(plot.title      = element_text(face="bold", size=9, hjust=0),
          axis.title.y    = element_text(size=9),
          axis.text       = element_text(size=8),
          legend.position = leg,
          legend.text     = element_text(size=8))

  list(full=p_full, post=p_post)
}

panel_pairs <- list(
  list("sigma2",  "(a)", "(b)", TRUE,  FALSE),
  list("beta11",  "(c)", "(d)", FALSE, FALSE),
  list("beta32",  "(e)", "(f)", FALSE, FALSE),
  list("theta11", "(g)", "(h)", FALSE, FALSE),
  list("theta21", "(i)", "(j)", FALSE, TRUE)
)

all_trace_plots <- list()
for (pp in panel_pairs) {
  out <- make_pair(pp[[1]], pp[[2]], pp[[3]], show_legend = pp[[4]], show_xlab = pp[[5]])
  all_trace_plots <- c(all_trace_plots, list(out$full, out$post))
}

p2_combined <- wrap_plots(all_trace_plots, ncol=2, byrow=TRUE)

pdf("../output/Figure_S1.pdf", width = 10, height = 6)
print(p1_combined)
dev.off()

pdf("../output/Figure_S2.pdf", width = 10, height = 12)
print(p2_combined)
dev.off()

cat("Gelman-Rubin Diagnostics (post burn-in only)\n")
cat(sprintf("  %-10s  %s\n", "Parameter", "Rhat"))
cat(strrep("-", 30), "\n")

param_names_gr <- c("sigma2","beta11","beta32","theta11","theta21")
for (pn in param_names_gr) {
  post_chains <- lapply(1:3, function(ch)
    coda::mcmc(traces[[ch]][[pn]][post_start:n_stored]))
  ml  <- coda::mcmc.list(post_chains)
  gr  <- coda::gelman.diag(ml, multivariate = FALSE)
  cat(sprintf("  %-10s  %.4f\n", pn, gr$psrf[1,1]))
}

