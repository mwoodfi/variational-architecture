############################################################
## 03_run_W3_recovery.R
## W3: Parameter recovery / identification (Catastrophic Recovery)
## - Two-block design with correct likelihood: beliefs reset per block
## - Block 2 uses asymmetric identity priors to identify lambda
## - Fits with fixed beta (params$beta_choice)
##
## AUTO-EXPORT:
## - Table: tables/Table1_W3_recovery_summary.csv
## - Figure 14: figures/Figure 14_W3_recovery_boxplots.png
## - Appendix figs: figures/Figure 15-... (see below)
## Raw outputs:
## - data/W3_recovery_twoblock.RData
## - tables/W3_res_small.csv
## - tables/W3_res_medium.csv
## - tables/W3_res_large.csv
############################################################


rm(list = ls())

stopifnot(file.exists("Part-A_EFE Simulations.R"))

# RNG standardisation
RNGkind(kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")

# ----------------------------
# Dependencies (NO side-effect sourcing)
# ----------------------------
source("Part-A_EFE Simulations.R")

# Hard dependency checks (fail fast, explicit)
stopifnot(
  exists("params_catastrophic_recovery"),
  exists("cr_efe_policy"),
  exists("cr_choice_prob"),
  exists("cr_epistemic_value"),
  exists("cr_instr_cost")
)

# ----------------------------
# Config (define AFTER dependencies)
# ----------------------------
RUN_HEAVY <- TRUE  # set FALSE to skip recomputing and only render from existing CSVs

N_SMALL  <- 30
N_MEDIUM <- 100
N_LARGE  <- 300

TRUE_OMEGA <- list(gamma_p = 1, gamma_o = 6, lambda = 1)

SEED0_SMALL  <- 1000
SEED0_MEDIUM <- 2000
SEED0_LARGE  <- 3000

# ----------------------------
# Output directories
# ----------------------------
OUT_ROOT   <- "outputs"
OUT_DATA   <- file.path(OUT_ROOT, "data")
OUT_FIG    <- file.path(OUT_ROOT, "figures")
OUT_TABLES <- file.path(OUT_ROOT, "tables")

dir.create(OUT_ROOT,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DATA,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TABLES, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Local helpers
# (Inlining removes rm(list=ls()) side-effects from source() calls.)
# ============================================================

# ---- W1 simulator (verbatim behaviour: outputs t,s,p_B,choice_B)
simulate_cr_constant_omega <- function(params,
                                       omega = list(gamma_p = 1, gamma_o = 6, lambda = 1),
                                       seed = 1L) {
  set.seed(seed)
  
  s <- numeric(params$T_total)
  s[1] <- params$s0
  
  beliefs <- list(
    harm_aB = params$harm_prior_a, harm_bB = params$harm_prior_b,
    id_a    = params$id_prior_a,   id_b    = params$id_prior_b
  )
  
  out <- vector("list", params$T_total)
  
  for (t in seq_len(params$T_total)) {
    
    GA <- cr_efe_policy(s[t], "A", beliefs, omega, params)
    GB <- cr_efe_policy(s[t], "B", beliefs, omega, params)
    
    pB <- cr_choice_prob(GA, GB, params$beta_choice)
    choiceB <- rbinom(1, 1, pB)
    policy_chosen <- if (choiceB == 1) "B" else "A"
    
    p_harm_true <- if (policy_chosen == "A") params$p_harm_A_true else params$p_harm_B_true
    harm_obs <- rbinom(1, 1, p_harm_true)
    
    if (policy_chosen == "B") {
      beliefs$harm_aB <- beliefs$harm_aB + harm_obs
      beliefs$harm_bB <- beliefs$harm_bB + (1 - harm_obs)
    }
    
    if (policy_chosen == "A") {
      beliefs$id_a <- beliefs$id_a + params$id_update_strength
    } else {
      beliefs$id_b <- beliefs$id_b + params$id_update_strength * params$id_recovery_ratio
    }
    
    if (t < params$T_total) {
      drift <- if (policy_chosen == "A") params$s_drift_A else -params$s_drift_B
      s[t + 1] <- pmin(pmax(s[t] + drift + rnorm(1, 0, params$s_noise), 0), 1)
    }
    
    out[[t]] <- data.frame(
      t = t,
      s = s[t],
      p_B = pB,
      choice_B = choiceB,
      harm_obs = harm_obs
    )
  }
  
  do.call(rbind, out)
}

# ---- W2 likelihood helpers (minimal set used by W3)
.clip01 <- function(p, eps = 1e-9) pmin(pmax(p, eps), 1 - eps)
.BIG <- 1e12

.as_finite_scalar <- function(x, fallback = .BIG) {
  if (length(x) != 1L) return(fallback)
  if (!is.finite(x)) return(fallback)
  x
}

.assert_params_cr <- function(params) {
  need <- c(
    "harm_prior_a","harm_prior_b","id_prior_a","id_prior_b",
    "id_update_strength","id_recovery_ratio",
    "p_harm_A_true","p_harm_B_true",
    "s_drift_A","s_drift_B","s_noise",
    "k_A0","k_B0","k_state",
    "beta_choice"
  )
  missing <- need[!vapply(need, function(nm) !is.null(params[[nm]]), logical(1))]
  if (length(missing)) {
    stop(paste0("params is missing required fields: ", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

loglik_cr <- function(theta, df, params, model = "EFE_full") {
  
  .assert_params_cr(params)
  
  gamma_p     <- theta[1]
  gamma_o     <- theta[2]
  lambda      <- theta[3]
  beta_choice <- theta[4]
  
  if (!all(is.finite(c(gamma_p, gamma_o, lambda, beta_choice)))) return(-Inf)
  if (beta_choice <= 0) return(-Inf)
  
  params$beta_choice <- beta_choice
  
  omega <- list(gamma_p = gamma_p, gamma_o = gamma_o, lambda = lambda)
  
  T_total <- nrow(df)
  if (T_total <= 0) return(-Inf)
  
  s <- df$s
  if (length(s) != T_total) return(-Inf)
  
  beliefs <- list(
    harm_aB = params$harm_prior_a, harm_bB = params$harm_prior_b,
    id_a    = params$id_prior_a,   id_b    = params$id_prior_b
  )
  
  ll <- 0
  
  for (t in seq_len(T_total)) {
    
    GA <- cr_efe_policy(s[t], "A", beliefs, omega, params)
    GB <- cr_efe_policy(s[t], "B", beliefs, omega, params)
    
    if (model == "no_identity") {
      omega_tmp <- omega; omega_tmp$lambda <- 0
      GA <- cr_efe_policy(s[t], "A", beliefs, omega_tmp, params)
      GB <- cr_efe_policy(s[t], "B", beliefs, omega_tmp, params)
    }
    
    if (model == "no_epistemic") {
      omega_tmp <- omega; omega_tmp$gamma_o <- 0
      GA <- cr_efe_policy(s[t], "A", beliefs, omega_tmp, params)
      GB <- cr_efe_policy(s[t], "B", beliefs, omega_tmp, params)
    }
    
    if (model == "RI_cost") {
      kappa <- abs(gamma_o)
      
      IG_B <- cr_epistemic_value(beliefs$harm_aB, beliefs$harm_bB)
      IG_A <- 0
      
      instr_A <- gamma_p * cr_instr_cost(s[t], "A", params)
      instr_B <- gamma_p * cr_instr_cost(s[t], "B", params)
      
      GA <- instr_A + kappa * IG_A
      GB <- instr_B + kappa * IG_B
    }
    
    pB <- cr_choice_prob(GA, GB, params$beta_choice)
    pB <- .clip01(pB)
    
    choiceB <- df$choice_B[t]
    if (!is.finite(choiceB)) return(-Inf)
    
    ll_step <- dbinom(choiceB, 1, pB, log = TRUE)
    if (length(ll_step) != 1L || !is.finite(ll_step)) return(-Inf)
    ll <- ll + ll_step
    
    policy_chosen <- if (choiceB == 1) "B" else "A"
    harm_obs <- df$harm_obs[t]
    if (!is.finite(harm_obs)) return(-Inf)
    
    if (policy_chosen == "B") {
      beliefs$harm_aB <- beliefs$harm_aB + harm_obs
      beliefs$harm_bB <- beliefs$harm_bB + (1 - harm_obs)
    }
    
    if (policy_chosen == "A") {
      beliefs$id_a <- beliefs$id_a + params$id_update_strength
    } else {
      beliefs$id_b <- beliefs$id_b + params$id_update_strength * params$id_recovery_ratio
    }
  }
  
  ll
}

fit_cr_fixbeta <- function(df, params, model = "EFE_full",
                           init  = c(1, 1, 1),
                           lower = c(0.001, 0.001, 0.0001),
                           upper = c(50, 50, 50)) {
  
  .assert_params_cr(params)
  
  nll <- function(th) {
    theta4 <- c(th[1], th[2], th[3], params$beta_choice)
    val <- -loglik_cr(theta4, df, params, model = model)
    .as_finite_scalar(val, fallback = .BIG)
  }
  
  v0 <- nll(init)
  if (!is.finite(v0) || v0 >= .BIG) {
    init <- (lower + upper) / 2
    v1 <- nll(init)
    if (!is.finite(v1) || v1 >= .BIG) {
      stop("fit_cr_fixbeta: cannot find a finite starting point. Check df and params consistency.")
    }
  }
  
  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)
  
  list(
    par   = c(opt$par, beta_choice = params$beta_choice),
    value = opt$value,
    conv  = opt$convergence,
    msg   = opt$message
  )
}

# -------------------------------------------------------------
# Two-block generator: return blocks separately (do not rbind for fitting)
# -------------------------------------------------------------
simulate_two_blocks_list <- function(params, omega, seed) {
  
  # helper: enforce required column names for likelihood
  harmonise_df_for_ll <- function(df, params, seed_tag) {
    
    # choice column
    if (!"choice_B" %in% names(df)) {
      if ("choiceB" %in% names(df)) df$choice_B <- df$choiceB
      if ("choice"  %in% names(df)) df$choice_B <- df$choice
    }
    
    # harm observation column (try aliases first)
    if (!"harm_obs" %in% names(df)) {
      if ("harm"     %in% names(df)) df$harm_obs <- df$harm
      if ("harmObs"  %in% names(df)) df$harm_obs <- df$harmObs
      if ("harm_B"   %in% names(df)) df$harm_obs <- df$harm_B
      if ("obs"      %in% names(df)) df$harm_obs <- df$obs
      if ("y"        %in% names(df)) df$harm_obs <- df$y
      if ("outcome"  %in% names(df)) df$harm_obs <- df$outcome
    }
    
    if (!"harm_obs" %in% names(df)) {
      stop(paste0(
        "simulate_* returned no harm_obs column (and no recognised alias). ",
        "Available columns: ", paste(names(df), collapse = ", "), "\n",
        "Head(df):\n", paste(capture.output(utils::head(df, 3)), collapse = "\n")
      ))
    }
    
    if (!"choice_B" %in% names(df)) {
      stop(paste0(
        "simulate_* returned no choice_B column (and no recognised alias). ",
        "Available columns: ", paste(names(df), collapse = ", "), "\n",
        "Head(df):\n", paste(capture.output(utils::head(df, 3)), collapse = "\n")
      ))
    }
    
    if (!"harm_obs" %in% names(df)) {
      stop(paste0(
        "simulate_* returned no harm_obs column (and no recognised alias). ",
        "Available columns: ", paste(names(df), collapse = ", "), "\n",
        "Head(df):\n", paste(capture.output(utils::head(df, 3)), collapse = "\n")
      ))
    }
    
    df
  }
  
  # Block 1: baseline identity prior
  df1 <- simulate_cr_constant_omega(params, omega = omega, seed = seed)
  df1 <- harmonise_df_for_ll(df1, params = params, seed_tag = seed + 11001)
  
  # Block 2: asymmetric identity priors to load lambda more strongly
  params2 <- params
  params2$id_prior_a <- 12
  params2$id_prior_b <- 2
  
  df2 <- simulate_cr_constant_omega(params2, omega = omega, seed = seed + 999)
  df2 <- harmonise_df_for_ll(df2, params = params2, seed_tag = seed + 22002)
  
  list(df1 = df1, df2 = df2, params2 = params2)
}

# -------------------------------------------------------------
# Fixed-beta MLE over TWO blocks (beliefs reset, correct priors per block)
# theta3 = c(gamma_p, gamma_o, lambda); beta fixed to params$beta_choice
# -------------------------------------------------------------
fit_cr_fixbeta_twoblock <- function(blocks, params, model = "EFE_full",
                                    init = c(1, 1, 1),
                                    lower = c(0.001, 0.001, 0.0001),
                                    upper = c(50, 50, 50)) {
  
  df1 <- blocks$df1
  df2 <- blocks$df2
  params2 <- blocks$params2
  
  nll <- function(th) {
    theta4 <- c(th[1], th[2], th[3], params$beta_choice)
    
    ll1 <- loglik_cr(theta4, df1, params, model = model)
    theta4_2 <- c(th[1], th[2], th[3], params2$beta_choice)
    ll2 <- loglik_cr(theta4_2, df2, params2, model = model)
    
    -(ll1 + ll2)
  }

  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)
  
  list(par = opt$par, value = opt$value, conv = opt$convergence, msg = opt$message)
}

# -------------------------------------------------------------
# Recovery runner
# -------------------------------------------------------------
run_recovery <- function(N = 10,
                         true_omega = TRUE_OMEGA,
                         params = params_catastrophic_recovery,
                         seed0 = 1000,
                         model_fit = "EFE_full",
                         init = c(1, 1, 1)) {
  
  est <- matrix(NA_real_, nrow = N, ncol = 3)
  colnames(est) <- c("gamma_p", "gamma_o", "lambda")
  conv <- integer(N)
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (i in seq_len(N)) {
    blocks <- simulate_two_blocks_list(params, omega = true_omega, seed = seed0 + i)
    
    fit <- fit_cr_fixbeta_twoblock(blocks, params, model = model_fit, init = init)
    
    est[i, ] <- fit$par
    conv[i] <- fit$conv
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  out <- data.frame(est)
  out$convergence <- conv
  out
}

# ----------------------------
# Summary stats (Table 1)
# ----------------------------
recovery_stats <- function(res, truth = c(gamma_p=1, gamma_o=6, lambda=1)) {
  X <- as.matrix(res[, c("gamma_p","gamma_o","lambda")])
  bias <- colMeans(X) - truth
  rmse <- sqrt(colMeans((X - matrix(truth, nrow=nrow(X), ncol=3, byrow=TRUE))^2))
  data.frame(
    param = names(truth),
    mean  = colMeans(X),
    bias  = as.numeric(bias),
    rmse  = as.numeric(rmse),
    conv_fail_rate = mean(res$convergence != 0)
  )
}

# ----------------------------
# Plotting helpers (base R)
# ----------------------------
plot_w3_grouped_boxplots <- function(res_small, res_medium, res_large,
                                     file_png) {
  
  if (!exists("set_plot_style") || !exists("draw_pub_grid")) {
    stop("Manuscript style functions not found. Ensure you source the EFE simulation file first.")
  }
  
  blue <- rgb(44, 100, 156, maxColorValue = 255)
  
  res_small$N  <- N_SMALL
  res_medium$N <- N_MEDIUM
  res_large$N  <- N_LARGE
  all <- rbind(res_small, res_medium, res_large)
  
  Ns <- c(N_SMALL, N_MEDIUM, N_LARGE)
  x_pos <- seq_along(Ns)
  
  draw_panel <- function(y, truth, ylab_expr, main_expr, ylim = NULL) {
    
    if (is.null(ylim)) {
      rng <- range(y, finite = TRUE)
      pad <- 0.08 * diff(rng)
      ylim <- c(rng[1] - pad, rng[2] + pad)
    }
    
    plot(NA, xlim = c(0.5, 3.5), ylim = ylim, axes = FALSE,
         xlab = "N", ylab = ylab_expr, main = main_expr, type = "n")
    
    x_ticks <- x_pos
    y_ticks_major <- pretty(ylim, n = 5)
    y_ticks_minor <- pretty(ylim, n = 9)
    
    draw_pub_grid(x_ticks = x_ticks,
                  y_ticks_major = y_ticks_major,
                  y_ticks_minor = y_ticks_minor)
    
    b <- boxplot(y ~ factor(all$N, levels = Ns), plot = FALSE)
    
    bxp(b,
        at = x_pos,
        add = TRUE,
        boxlwd = 1.6,
        whisklwd = 1.3,
        staplelwd = 1.3,
        outpch = 1,
        outcex = 0.75,
        boxfill = "grey88",
        border = "grey20",
        whiskcol = "grey20",
        staplecol = "grey20",
        outcol = "grey20",
        medlwd = 3.0,
        medcol = blue)
    
    abline(h = truth, lty = 2, lwd = 2.2, col = "grey40")
    
    axis(1, at = x_pos, labels = Ns, lwd = 0, lwd.ticks = 1.35)
    axis(2, at = y_ticks_major, labels = y_ticks_major, las = 1, lwd = 0, lwd.ticks = 1.35)
    
    box(lwd = 1.15, col = "grey20")
  }
  
  png(file_png, width = 900, height = 600)
  op_outer <- set_plot_style()
  par(font = 1, font.axis = 1, font.lab = 1, font.main = 1)
  par(xaxs = "i", yaxs = "i")
  on.exit(par(op_outer), add = TRUE)
  
  par(mfrow = c(1,3))
  
  draw_panel(all$gamma_p, 1,
             ylab_expr = expression(gamma[p]),
             main_expr = expression(gamma[p]))
  
  draw_panel(all$gamma_o, 6,
             ylab_expr = expression(gamma[o]),
             main_expr = expression(gamma[o]))
  
  draw_panel(all$lambda, 1,
             ylab_expr = expression(lambda),
             main_expr = expression(lambda))
  
  dev.off()
}

plot_w3_scatter_confounding <- function(res, Nlabel, file_png) {
  png(file_png, width = 900, height = 600)
  op_outer <- set_plot_style()
  par(font = 1, font.axis = 1, font.lab = 1, font.main = 1)
  par(xaxs = "i", yaxs = "i")
  on.exit(par(op_outer), add = TRUE)
  
  par(mfrow = c(1,2))
  
  # Panel 1: gamma_o vs lambda
  x1 <- res$gamma_o
  y1 <- res$lambda
  x_ticks1 <- pretty(x1, n = 5)
  y_ticks1_major <- pretty(y1, n = 5)
  y_ticks1_minor <- pretty(y1, n = 9)
  
  plot(x1, y1,
       type = "n",
       axes = FALSE,
       xlab = expression(gamma[o]),
       ylab = expression(lambda),
       main = bquote("N=" * .(Nlabel) * ": " * gamma[o] * " vs " * lambda))
  
  draw_pub_grid(x_ticks1, y_ticks1_major, y_ticks1_minor)
  points(x1, y1, pch = 16, cex = 0.7, col = rgb(44, 100, 156, maxColorValue = 255))
  abline(h = 1, lty = 2, lwd = 2.2, col = "grey40")
  abline(v = 6, lty = 2, lwd = 2.2, col = "grey40")
  axis(1, at = x_ticks1, labels = x_ticks1, lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks1_major, labels = y_ticks1_major, las = 1, lwd = 0, lwd.ticks = 1.35)
  box(lwd = 1.15, col = "grey20")
  
  # Panel 2: gamma_p vs lambda
  x2 <- res$gamma_p
  y2 <- res$lambda
  x_ticks2 <- pretty(x2, n = 5)
  y_ticks2_major <- pretty(y2, n = 5)
  y_ticks2_minor <- pretty(y2, n = 9)
  
  plot(x2, y2,
       type = "n",
       axes = FALSE,
       xlab = expression(gamma[p]),
       ylab = expression(lambda),
       main = bquote("N=" * .(Nlabel) * ": " * gamma[p] * " vs " * lambda))
  
  draw_pub_grid(x_ticks2, y_ticks2_major, y_ticks2_minor)
  points(x2, y2, pch = 16, cex = 0.7, col = rgb(44, 100, 156, maxColorValue = 255))
  abline(h = 1, lty = 2, lwd = 2.2, col = "grey40")
  abline(v = 1, lty = 2, lwd = 2.2, col = "grey40")
  axis(1, at = x_ticks2, labels = x_ticks2, lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks2_major, labels = y_ticks2_major, las = 1, lwd = 0, lwd.ticks = 1.35)
  box(lwd = 1.15, col = "grey20")
  
  dev.off()
}

# ----------------------------
# Run (or load) W3 results
# ----------------------------
if (RUN_HEAVY) {
  
  cat("\nRunning W3 (small)...\n")
  res_small  <- run_recovery(N = N_SMALL,  true_omega = TRUE_OMEGA, seed0 = SEED0_SMALL)
  
  cat("\nRunning W3 (medium)...\n")
  res_medium <- run_recovery(N = N_MEDIUM, true_omega = TRUE_OMEGA, seed0 = SEED0_MEDIUM)
  
  cat("\nRunning W3 (large)...\n")
  res_large  <- run_recovery(N = N_LARGE,  true_omega = TRUE_OMEGA, seed0 = SEED0_LARGE)
  
  save(res_small, res_medium, res_large, file = file.path(OUT_DATA, "W3_recovery_twoblock.RData"))
  
  write.csv(res_small,  file.path(OUT_TABLES, "W3_res_small.csv"),  row.names = FALSE)
  write.csv(res_medium, file.path(OUT_TABLES, "W3_res_medium.csv"), row.names = FALSE)
  write.csv(res_large,  file.path(OUT_TABLES, "W3_res_large.csv"),  row.names = FALSE)
  
} else {
  
  res_small  <- read.csv(file.path(OUT_TABLES, "W3_res_small.csv"))
  res_medium <- read.csv(file.path(OUT_TABLES, "W3_res_medium.csv"))
  res_large  <- read.csv(file.path(OUT_TABLES, "W3_res_large.csv"))
}

# ----------------------------
# Table (W3 summary) + export
# ----------------------------
tab_small  <- recovery_stats(res_small);  tab_small$N  <- N_SMALL
tab_medium <- recovery_stats(res_medium); tab_medium$N <- N_MEDIUM
tab_large  <- recovery_stats(res_large);  tab_large$N  <- N_LARGE

tab_w3 <- rbind(tab_small, tab_medium, tab_large)
tab_w3 <- tab_w3[, c("N","param","mean","bias","rmse","conv_fail_rate")]

print(tab_w3)
write.csv(tab_w3, file.path(OUT_TABLES, "W3_recovery_summary.csv"), row.names = FALSE)

# ----------------------------
# Figure 14 (main text): grouped boxplots across N (PNG ONLY)
# ----------------------------
plot_w3_grouped_boxplots(
  res_small, res_medium, res_large,
  file_png = file.path(OUT_FIG, "Figure 14_W3_Recovery Boxplots.png")
)

# ----------------------------
# Appendix diagnostics (PNG ONLY)
# ----------------------------
plot_w3_scatter_confounding(
  res_large, N_LARGE,
  file_png = file.path(OUT_FIG, "Figure 15_W3 Confounding Scatter vN300.png")
)

cat("\nW3 outputs written to outputs/, tables/, figures/main/, figures/appendix/\n")