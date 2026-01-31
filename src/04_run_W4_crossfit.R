############################################################
## 04_run_W4_crossfit.R
## W4: Cross-fitting / model discrimination (Catastrophic Recovery)
## - Two-/three-block design to discriminate models out-of-sample
## - Fit each candidate model by MLE (beta fixed) on TRAIN set (Block 1)
## - Evaluate held-out predictive log-likelihood on TEST set (Blocks 2+3)
##
## AUTO-EXPORT (src/outputs/{data,figures,tables}):
## - data/<run_id>.rds
## - tables/<run_id>_confusion.csv
## - tables/<run_id>_subject_diagnostics.csv
## - tables/<run_id>_na_summary.csv
## - tables/<run_id>_referee_summary_no_identity.csv
## - tables/<run_id>_meta.csv
## - figures/Figure 16_W4_Crossfit_Bar.png  (only if make_figure = TRUE)
## - figures/Figure 17_W4_Crossfit_Confusion.png  (only if make_figure = TRUE)
## - (optional replication outputs, only if do_replication = TRUE):
## - data/<run_id>_rep_seed<rep_seed>.rds
## - tables/<run_id>_rep_seed<rep_seed>_confusion_diff.csv
## - tables/<run_id>_rep_seed<rep_seed>_replication_summary.csv
############################################################

rm(list = ls())

# RNG standardisation (referee-safe determinism; match W1–W3)
RNGkind(kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")

# ----------------------------
# Dependencies (NO side-effect sourcing)
# ----------------------------
stopifnot(file.exists("Part-A_EFE Simulations.R"))
source("Part-A_EFE Simulations.R")

# Hard dependency checks (fail fast, explicit)
stopifnot(
  exists("params_catastrophic_recovery"),
  exists("cr_efe_policy"),
  exists("cr_choice_prob"),
  exists("cr_epistemic_value"),
  exists("cr_instr_cost"),
  exists("set_plot_style"),
  exists("draw_pub_grid")
)

# ----------------------------
# Paths
# ----------------------------
OUT_ROOT     <- "outputs"
OUT_DATA     <- file.path(OUT_ROOT, "data")
OUT_FIG      <- file.path(OUT_ROOT, "figures")
OUT_TABLES   <- file.path(OUT_ROOT, "tables")

dir.create(OUT_ROOT,     showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DATA,     showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG,      showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TABLES,   showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Local likelihood helpers (verbatim W2/W3 minimal set)
# ============================================================

# Safe clipping + finite objective guards (optim() robustness)
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

# Self-contained likelihood
# theta = c(gamma_p, gamma_o, lambda, beta_choice)
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

stopifnot(is.function(loglik_cr))

# -------------------------------------------------------------
# Simulator that matches the likelihood exactly (per model)
# -------------------------------------------------------------
simulate_cr_constant_omega_model <- function(params,
                                             omega = list(gamma_p = 1, gamma_o = 6, lambda = 1),
                                             model = "EFE_full",
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
    
    if (model == "EFE_full") {
      GA <- cr_efe_policy(s[t], "A", beliefs, omega, params)
      GB <- cr_efe_policy(s[t], "B", beliefs, omega, params)
    }
    
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
      kappa <- abs(omega$gamma_o)
      IG_B <- cr_epistemic_value(beliefs$harm_aB, beliefs$harm_bB)
      instr_A <- omega$gamma_p * cr_instr_cost(s[t], "A", params)
      instr_B <- omega$gamma_p * cr_instr_cost(s[t], "B", params)
      GA <- instr_A
      GB <- instr_B + kappa * IG_B
    }
    
    pB <- cr_choice_prob(GA, GB, params$beta_choice)
    pB <- pmin(pmax(pB, 1e-9), 1 - 1e-9)
    
    choiceB <- rbinom(1, 1, pB)
    policy_chosen <- if (choiceB == 1) "B" else "A"
    
    p_harm_true <- if (policy_chosen == "A") params$p_harm_A_true else params$p_harm_B_true
    harm_obs <- rbinom(1, 1, p_harm_true)
    
    # Belief updates (verbatim)
    if (policy_chosen == "B") {
      beliefs$harm_aB <- beliefs$harm_aB + harm_obs
      beliefs$harm_bB <- beliefs$harm_bB + (1 - harm_obs)
    }
    
    if (policy_chosen == "A") {
      beliefs$id_a <- beliefs$id_a + params$id_update_strength
    } else {
      beliefs$id_b <- beliefs$id_b + params$id_update_strength * params$id_recovery_ratio
    }
    
    # State evolution (verbatim)
    if (t < params$T_total) {
      drift <- if (policy_chosen == "A") params$s_drift_A else -params$s_drift_B
      s[t + 1] <- pmin(pmax(s[t] + drift + rnorm(1, 0, params$s_noise), 0), 1)
    }
    
    out[[t]] <- data.frame(
      t = t,
      s = s[t],
      choice_B = choiceB,
      harm_obs = harm_obs
    )
  }
  
  do.call(rbind, out)
}

# -------------------------------------------------------------
# Three-block generator (IDENTIFICATION-SAFE TEST DESIGN):
# Block 1: Baseline (Training)
# Block 2: Lambda identification block (ONLY identity/rigidity breaks symmetry)
# Block 3: Epistemic identification block (identity neutralised; instrumental parity)
#
# Principle (referee-safe): omega (agent phenotype) is invariant across blocks;
# only task parameters (params2/params3) change.
# -------------------------------------------------------------
simulate_three_blocks_list_model <- function(params, omega, seed, dgp_model) {
  
  # --- Block 1 (Baseline) ---
  df1 <- simulate_cr_constant_omega_model(
    params, omega = omega, model = dgp_model, seed = seed
  )
  
  # --- Block 2 (Lambda identification / NO epistemic attenuation; instrumental parity) ---
  params2 <- params
  
  # 1) Identity salience via asymmetric identity priors (directional pull toward A)
  total_id <- params$id_prior_a + params$id_prior_b
  params2$id_prior_a <- total_id * 0.90
  params2$id_prior_b <- total_id * 0.10
  
  # 2) Instrumental parity: remove preference differences so lambda is the only driver
  params2$k_B0    <- params2$k_A0
  params2$k_state <- 0
  
  # 3) Epistemic attenuation REMOVED: keep harm priors unchanged (no +1000 inflation)
  prior_inflation <- 0
  params2$harm_prior_a <- params$harm_prior_a
  params2$harm_prior_b <- params$harm_prior_b
  
  df2 <- simulate_cr_constant_omega_model(
    params2, omega = omega, model = dgp_model, seed = seed + 9999
  )
  
  # --- Block 3 (Epistemic identification / Identity neutralised) ---
  params3 <- params
  
  # 1) Instrumental parity (remove instrumental asymmetry)
  params3$k_B0 <- params3$k_A0
  params3$k_state <- 0
  
  # 2) Neutralise identity task-side:
  #    (a) stop identity learning updates
  params3$id_update_strength <- 0
  #    (b) remove directional identity bias via symmetric priors
  total_id3 <- params$id_prior_a + params$id_prior_b
  params3$id_prior_a <- total_id3 / 2
  params3$id_prior_b <- total_id3 / 2
  
  df3 <- simulate_cr_constant_omega_model(
    params3, omega = omega, model = dgp_model, seed = seed + 19999
  )
  
  # Minimal numeric diagnostics for metadata (no side effects; returned to caller):
  IG_base <- cr_epistemic_value(params$harm_prior_a,  params$harm_prior_b)
  IG_b2   <- cr_epistemic_value(params2$harm_prior_a, params2$harm_prior_b)
  
  list(
    df1 = df1, df2 = df2, df3 = df3,
    params2 = params2, params3 = params3,
    diag = list(
      prior_inflation = prior_inflation,
      IG_base = IG_base,
      IG_block2 = IG_b2
    )
  )
}

# -------------------------------------------------------------
# Fit on one block (fixed beta). Model-specific free parameters.
# Returns free parameter vector for that model.
# -------------------------------------------------------------
fit_block_fixedbeta <- function(df, params, model, init, lower, upper) {
  
  # Finite objective wrapper (optim() robustness)
  nll <- function(th_free) {
    
    if (model == "EFE_full") {
      theta4 <- c(th_free[1], th_free[2], th_free[3], params$beta_choice)
    }
    
    if (model == "no_identity") {
      theta4 <- c(th_free[1], th_free[2], 0, params$beta_choice)
    }
    
    if (model == "no_epistemic") {
      theta4 <- c(th_free[1], 0, th_free[2], params$beta_choice)
    }
    
    if (model == "RI_cost") {
      theta4 <- c(th_free[1], abs(th_free[2]), 0, params$beta_choice)
    }
    
    val <- -loglik_cr(theta4, df, params, model = model)
    .as_finite_scalar(val, fallback = .BIG)
  }
  
  # Ensure finite start (platform-stable)
  v0 <- nll(init)
  if (!is.finite(v0) || v0 >= .BIG) {
    init <- (lower + upper) / 2
    v1 <- nll(init)
    if (!is.finite(v1) || v1 >= .BIG) {
      return(list(par = init, value = .BIG, conv = 999L))
    }
  }
  
  opt <- optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper)
  list(par = opt$par, value = opt$value, conv = opt$convergence)
}

predict_ll_block <- function(df, params, model, fit_free_par) {
  
  if (model == "EFE_full")     theta4 <- c(fit_free_par[1], fit_free_par[2], fit_free_par[3], params$beta_choice)
  if (model == "no_identity")  theta4 <- c(fit_free_par[1], fit_free_par[2], 0, params$beta_choice)
  if (model == "no_epistemic") theta4 <- c(fit_free_par[1], 0, fit_free_par[2], params$beta_choice)
  if (model == "RI_cost")      theta4 <- c(fit_free_par[1], abs(fit_free_par[2]), 0, params$beta_choice)
  
  loglik_cr(theta4, df, params, model = model)
}

# -------------------------------------------------------------
# One subject cross-fit:
# train on Block 1; test on Block 2 + Block 3 (sum predictive LL)
# -------------------------------------------------------------
crossfit_one_subject <- function(params, omega, seed,
                                 dgp_model,
                                 candidate_models = c("EFE_full","no_identity","no_epistemic","RI_cost")) {
  
  blocks <- simulate_three_blocks_list_model(params, omega, seed, dgp_model)
  
  df1 <- blocks$df1
  df2 <- blocks$df2
  df3 <- blocks$df3
  params2 <- blocks$params2
  params3 <- blocks$params3
  
  delta_ll_min <- 2.0
  
  # raw out-of-sample predictive log score (no penalty)
  raw_scores <- setNames(rep(NA_real_, length(candidate_models)), candidate_models)
  
  # penalised score (decision rule)
  pen_scores <- setNames(rep(NA_real_, length(candidate_models)), candidate_models)
  
  # penalty components
  k_free <- c(EFE_full = 3, no_identity = 2, no_epistemic = 2, RI_cost = 2)
  
  n_test <- nrow(df2) + nrow(df3)
  
  for (m in candidate_models) {
    
    # IMPROVED INITIALIZATION (Interior Starts)
    if (m == "EFE_full") {
      init  <- c(1, 6, 1)        # Start lambda at 1, not 0
      lower <- c(0.001, 0.001, 0)
      upper <- c(50, 50, 50)
    }
    if (m == "no_identity") {
      init  <- c(1, 6)
      lower <- c(0.001, 0.001)
      upper <- c(50, 50)
    }
    if (m == "no_epistemic") {
      init  <- c(1, 1)          # Start lambda at 1, not 0
      lower <- c(0.001, 0)
      upper <- c(50, 50)
    }
    if (m == "RI_cost") {
      init  <- c(1, 6)
      lower <- c(0.001, 0)
      upper <- c(50, 50)
    }
    
    fit <- fit_block_fixedbeta(df1, params, m, init, lower, upper)
    
    if (!is.null(fit$conv) && fit$conv != 0) {
      raw_scores[m] <- -Inf
      pen_scores[m] <- -Inf
    } else {
      ll2 <- predict_ll_block(df2, params2, m, fit$par)
      ll3 <- predict_ll_block(df3, params3, m, fit$par)
      
      raw <- ll2 + ll3
      pen <- 0.5 * k_free[m] * log(n_test)
      
      raw_scores[m] <- raw
      pen_scores[m] <- raw - pen
    }
  }
  
  winner_raw <- names(which.max(raw_scores))
  
  ord <- order(pen_scores, decreasing = TRUE)
  best   <- ord[1]
  second <- ord[2]
  
  if (pen_scores[best] - pen_scores[second] < delta_ll_min) {
    winner_pen <- NA_character_
  } else {
    winner_pen <- names(pen_scores)[best]
  }
  
  delta_raw <- raw_scores["EFE_full"] - raw_scores["no_identity"]
  delta_pen <- pen_scores["EFE_full"] - pen_scores["no_identity"]
  
  list(
    winner = winner_pen,
    winner_raw = winner_raw,
    winner_pen = winner_pen,
    scores_raw = raw_scores,
    scores_pen = pen_scores,
    delta_raw = delta_raw,
    delta_pen = delta_pen,
    n_test = n_test,
    dgp = dgp_model,
    seed = seed
  )
}

# -------------------------------------------------------------
# Run W4 across subjects and DGPs -> confusion matrix
# + diagnostics dataframe for referee-safe margins and stability checks
# -------------------------------------------------------------
run_W4 <- function(N = 200,
                   params = params_catastrophic_recovery,
                   omega_truth = list(gamma_p = 1, gamma_o = 6, lambda = 1),
                   dgps = c("EFE_full","no_identity","no_epistemic","RI_cost"),
                   candidate_models = c("EFE_full","no_identity","no_epistemic","RI_cost"),
                   seed0 = 5000) {
  
  results <- list()
  diag_rows <- vector("list", length(dgps) * N)
  
  pb <- txtProgressBar(min = 0, max = length(dgps) * N, style = 3)
  k <- 0
  r <- 0
  
  for (dgp in dgps) {
    winners <- character(N)
    
    for (i in seq_len(N)) {
      
      this_seed <- seed0 + i + 100000 * match(dgp, dgps)
      
      out <- crossfit_one_subject(
        params, omega_truth, this_seed,
        dgp_model = dgp,
        candidate_models = candidate_models
      )
      
      winners[i] <- out$winner
      
      r <- r + 1
      diag_rows[[r]] <- data.frame(
        dgp = out$dgp,
        seed = out$seed,
        n_test = out$n_test,
        
        winner_raw = out$winner_raw,
        winner_pen = out$winner_pen,
        
        raw_EFE_full = out$scores_raw["EFE_full"],
        raw_no_identity = out$scores_raw["no_identity"],
        raw_no_epistemic = out$scores_raw["no_epistemic"],
        raw_RI_cost = out$scores_raw["RI_cost"],
        
        pen_EFE_full = out$scores_pen["EFE_full"],
        pen_no_identity = out$scores_pen["no_identity"],
        pen_no_epistemic = out$scores_pen["no_epistemic"],
        pen_RI_cost = out$scores_pen["RI_cost"],
        
        delta_raw_EFE_minus_noid = out$delta_raw,
        delta_pen_EFE_minus_noid = out$delta_pen,
        
        stringsAsFactors = FALSE
      )
      
      k <- k + 1
      setTxtProgressBar(pb, k)
    }
    
    results[[dgp]] <- winners
  }
  
  close(pb)
  
  conf <- matrix(0L, nrow = length(dgps), ncol = length(candidate_models),
                 dimnames = list(true = dgps, selected = candidate_models))
  
  for (dgp in dgps) {
    tab <- table(factor(results[[dgp]][!is.na(results[[dgp]])],
                        levels = candidate_models))
    conf[dgp, ] <- as.integer(tab)
  }
  
  diagnostics <- do.call(rbind, diag_rows)
  
  list(confusion = conf, winners = results, diagnostics = diagnostics)
}

# -------------------------------------------------------------
# Plot confusion matrix (base R, manuscript-aligned)
# -------------------------------------------------------------

MODEL_ORDER <- c("no_epistemic", "no_identity", "RI_cost", "EFE_full")

plot_confusion <- function(M, main = "Confusion matrix") {
  
  M <- M[MODEL_ORDER, MODEL_ORDER, drop = FALSE]
  
  blue <- rgb(44, 100, 156, maxColorValue = 255)
  
  row_tot <- rowSums(M)
  P <- sweep(M, 1, row_tot, "/")
  
  x_pos <- seq_len(ncol(M))
  y_pos <- seq_len(nrow(M))
  
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.18, 0.82))
  
  # TOP PANEL: title + legend
  par(mar = c(0, 0, 1.2, 0))
  plot.new()
  title(main, line = 0.2, cex.main = 1)
  
  pal_nz <- grDevices::colorRampPalette(
    c(rgb(220, 233, 246, maxColorValue = 255), blue)
  )(5)
  
  legend(
    "center",
    legend = c("0%", "0–20%", "20–40%", "40–60%", "60–80%", "80–100%"),
    fill   = c("white", pal_nz),
    title  = "Row-normalised share",
    horiz  = TRUE,
    bty    = "o",
    bg     = rgb(1, 1, 1, 0.72),
    box.col = "grey35",
    box.lwd = 1,
    cex = 1,
    x.intersp = 1.2,
    y.intersp = 1.05,
    seg.len = 2.6
  )
  
  # BOTTOM PANEL: matrix
  par(
    mar = c(4.6, 8.2, 0.5, 1.5),
    xaxs = "i",
    yaxs = "i"
  )
  
  plot(NA,
       xlim = c(0.5, ncol(M) + 0.5),
       ylim = c(0.5, nrow(M) + 0.5),
       type = "n",
       axes = FALSE,
       xlab = "Selected model",
       ylab = "")
  
  mtext("True DGP", side = 2, line = 6.0, cex = 1.18, font = 1)
  
  nz <- (M > 0)
  brks_nz <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0000001)
  
  z <- matrix(0L, nrow = nrow(M), ncol = ncol(M))
  if (any(nz)) {
    z[nz] <- as.integer(cut(P[nz], breaks = brks_nz, include.lowest = TRUE))
  }
  
  fill_cols <- matrix("white", nrow = nrow(M), ncol = ncol(M))
  if (any(nz)) fill_cols[nz] <- pal_nz[z[nz]]
  
  grid_col <- "grey90"
  
  for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
           col = fill_cols[i, j],
           border = grid_col, lwd = 1.1)
      
      txt <- sprintf("%d\n(%.0f%%)", M[i, j], 100 * P[i, j])
      text_col <- if (nz[i, j] && z[i, j] >= 4) "white" else "black"
      text(j, i, txt, cex = 0.88, col = text_col)
    }
  }
  
  axis(1, at = x_pos, labels = colnames(M),
       lwd = 0, lwd.ticks = 1.35)
  
  axis(2, at = y_pos, labels = rownames(M),
       las = 1, lwd = 0, lwd.ticks = 1.35,
       line = -0.6)
  
  box(lwd = 1.15, col = "grey20")
}

# -------------------------------------------------------------
# Manuscript-required exports (fixed filenames)
# -------------------------------------------------------------

# Simple selection-rate bar plot (penalised winner, NA excluded)
plot_w4_bar <- function(conf, file_png, main = "W4: Out-of-sample model selection rates") {
  
  png(file_png, width = 900, height = 600)
  op <- set_plot_style()
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  conf <- conf[MODEL_ORDER, MODEL_ORDER, drop = FALSE]
  
  # row-normalised selection shares
  row_tot <- rowSums(conf)
  P <- sweep(conf, 1, row_tot, "/")
  
  # plot as grouped bars by true DGP
  y <- t(P)
  
  # (a) House Standard Blue Palette (gradient for 4 models)
  blue_dark  <- rgb(44, 100, 156, maxColorValue = 255)
  blue_light <- rgb(220, 233, 246, maxColorValue = 255)
  bar_cols   <- grDevices::colorRampPalette(c(blue_light, blue_dark))(nrow(y))
  
  # (b) Layout: Top Panel (Legend), Bottom Panel (Plot) - matches Fig 17
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.18, 0.82))
  
  # --- Top Panel: Title & Legend ---
  par(mar = c(0, 0, 1.2, 0))
  plot.new()
  title(main, line = 0.2, cex.main = 1)
  
  legend("center",
         legend = rownames(y),
         fill   = bar_cols,
         horiz  = TRUE,
         bty    = "o",
         bg     = rgb(1, 1, 1, 0.72),
         box.col = "grey35",
         box.lwd = 1,
         cex    = 1,
         x.intersp = 1.2,
         y.intersp = 1.05)
  
  # --- Bottom Panel: Bar Plot ---
  par(mar = c(4.6, 5.0, 0.5, 1.5), xaxs = "i", yaxs = "i")
  
  # Determine bar geometry without plotting to calculate limits
  bp_locs <- barplot(y, beside = TRUE, plot = FALSE)
  
  # (e) Breathing room on x-axis
  x_min <- min(bp_locs)
  x_max <- max(bp_locs)
  xlim_val <- c(x_min - 2, x_max + 2)
  
  # Setup coordinate system
  plot(NA, xlim = xlim_val, ylim = c(0, 1), 
       axes = FALSE, xlab = "", ylab = "")
  
  # (d) Horizontal Grid Only (no vertical grid)
  y_ticks <- seq(0, 1, by = 0.2)
  abline(h = y_ticks, col = "grey90", lwd = 1)
  
  # Draw Bars
  barplot(y, beside = TRUE, add = TRUE, 
          col = bar_cols, 
          border = "grey20",
          axes = FALSE)
  
  # X-axis labels
  x_ticks <- colMeans(bp_locs)
  axis(1, at = x_ticks, labels = colnames(conf), 
       lwd = 0, lwd.ticks = 1.35)
  mtext("True DGP", side = 1, line = 3.0, cex = 1)
  
  # (c) Y-axis as Percentages
  axis(2, at = y_ticks, labels = paste0(y_ticks * 100, "%"),
       las = 1, lwd = 0, lwd.ticks = 1.35)
  mtext("Share selected (row-normalised; NA excluded)", side = 2, line = 3.5, cex = 1)
  
  box(lwd = 1.15, col = "grey20")
}

# -------------------------------------------------------------
# One-shot runner that produces a complete output bundle.
# -------------------------------------------------------------
run_W4_and_write_outputs <- function(N = 200,
                                     params = params_catastrophic_recovery,
                                     omega_truth = list(gamma_p = 1, gamma_o = 6, lambda = 1),
                                     dgps = c("EFE_full","no_identity","no_epistemic","RI_cost"),
                                     candidate_models = c("EFE_full","no_identity","no_epistemic","RI_cost"),
                                     seed0 = 5000,
                                     make_figures = TRUE) {
  
  # Ensure recursive dirs even if out_root changes in future edits
  dir.create(OUT_ROOT,   showWarnings = FALSE, recursive = TRUE)
  dir.create(OUT_DATA,   showWarnings = FALSE, recursive = TRUE)
  dir.create(OUT_FIG,    showWarnings = FALSE, recursive = TRUE)
  dir.create(OUT_TABLES, showWarnings = FALSE, recursive = TRUE)

  w4 <- run_W4(N = N, params = params, omega_truth = omega_truth,
               dgps = dgps, candidate_models = candidate_models, seed0 = seed0)
  
  # -------------------------------------------------------------
  # Forensic run_id bundle (hostile-referee transparency)
  # -------------------------------------------------------------
  run_id <- paste0("W4_N", N, "_seed0_", seed0, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  
  saveRDS(w4, file = file.path(OUT_DATA, paste0(run_id, ".rds")))
  
  write.csv(w4$confusion,
            file.path(OUT_TABLES, paste0(run_id, "_confusion.csv")),
            row.names = TRUE)
  
  write.csv(w4$diagnostics,
            file.path(OUT_TABLES, paste0(run_id, "_subject_diagnostics.csv")),
            row.names = FALSE)
  
  na_summary <- aggregate(is.na(w4$diagnostics$winner_pen),
                          by = list(dgp = w4$diagnostics$dgp),
                          FUN = function(z) c(n_na = sum(z), N = length(z), share = mean(z)))
  na_summary <- data.frame(
    dgp = na_summary$dgp,
    n_na = na_summary$x[, "n_na"],
    N = na_summary$x[, "N"],
    share_na = na_summary$x[, "share"],
    row.names = NULL
  )
  write.csv(na_summary,
            file.path(OUT_TABLES, paste0(run_id, "_na_summary.csv")),
            row.names = FALSE)
  
  boot_ci_mean <- function(x, B = 5000L, conf = 0.95) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 2) return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
    m <- mean(x)
    idx <- replicate(B, sample.int(n, size = n, replace = TRUE))
    boots <- colMeans(matrix(x[idx], nrow = n))
    alpha <- (1 - conf) / 2
    qs <- as.numeric(quantile(boots, probs = c(alpha, 1 - alpha), names = FALSE))
    c(mean = m, lo = qs[1], hi = qs[2])
  }
  
  d_noid <- subset(w4$diagnostics, dgp == "no_identity")
  ci_raw <- boot_ci_mean(d_noid$delta_raw_EFE_minus_noid)
  ci_pen <- boot_ci_mean(d_noid$delta_pen_EFE_minus_noid)
  
  summary_ref <- data.frame(
    dgp = "no_identity",
    N = nrow(d_noid),
    n_test_unique = paste(sort(unique(d_noid$n_test)), collapse = ";"),
    mean_delta_raw_EFE_minus_noid = ci_raw["mean"],
    ci95_lo_delta_raw = ci_raw["lo"],
    ci95_hi_delta_raw = ci_raw["hi"],
    mean_delta_pen_EFE_minus_noid = ci_pen["mean"],
    ci95_lo_delta_pen = ci_pen["lo"],
    ci95_hi_delta_pen = ci_pen["hi"],
    stringsAsFactors = FALSE
  )
  write.csv(summary_ref,
            file.path(OUT_TABLES, paste0(run_id, "_referee_summary_no_identity.csv")),
            row.names = FALSE)
  
  # -------------------------------------------------------------
  # META export block
  # -------------------------------------------------------------
  # Pull a representative diagnostic from the *first* subject of the run,
  # deterministically keyed by seed0 and the first DGP label.
  # This is a referee-safe numeric anchor for “IG is negligible in Block 2”.
  diag_seed <- seed0 + 1 + 100000 * match(dgps[1], dgps)
  diag_blocks <- simulate_three_blocks_list_model(params, omega_truth, diag_seed, dgp_model = dgps[1])
  
  meta <- data.frame(
    run_id = run_id,
    N = N,
    seed0 = seed0,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    dgps = paste(dgps, collapse = ";"),
    candidates = paste(candidate_models, collapse = ";"),
    beta_choice = params$beta_choice,
    delta_ll_min = 2.0,
    block_seed_offsets = "B2:+9999;B3:+19999",
    
    # Block 2 (Identity salience + epistemic attenuation)
    block2_identity_salience =
      "total_id <- params$id_prior_a + params$id_prior_b; params2$id_prior_a <- total_id*0.90; params2$id_prior_b <- total_id*0.10; params2$k_B0 <- params2$k_A0; params2$k_state <- 0",
    block2_epistemic_attenuation =
      "None (no harm-prior inflation; Block 2 is a lambda-identification block)",
    
    # Block 3 (Instrumental parity + identity neutralisation)
    block3_instrumental_parity =
      "params3$k_B0 <- params3$k_A0; params3$k_state <- 0",
    block3_identity_neutralisation =
      "params3$id_update_strength <- 0; total_id <- params$id_prior_a + params$id_prior_b; params3$id_prior_a <- total_id/2; params3$id_prior_b <- total_id/2",
    
    # Numeric diagnostic (single-run anchor; makes “IG ~ 0” checkable)
    IG_base_priors =
      as.numeric(diag_blocks$diag$IG_base),
    IG_block2_inflated_priors =
      as.numeric(diag_blocks$diag$IG_block2),
    IG_block2_to_base_ratio =
      as.numeric(ifelse(diag_blocks$diag$IG_base > 0,
                        diag_blocks$diag$IG_block2 / diag_blocks$diag$IG_base,
                        NA_real_)),
    
    stringsAsFactors = FALSE
  )
  
  write.csv(
    meta,
    file.path(OUT_TABLES, paste0(run_id, "_meta.csv")),
    row.names = FALSE
  )
  
  # Export primary bundle (fixed names; manuscript-aligned)
  save(w4, file = file.path(OUT_DATA, "W4_crossfit_results.RData"))
  
  # ll matrix: mean held-out predictive LL by true DGP (rows) and candidate model (cols)
  d <- w4$diagnostics
  ll_mat <- matrix(NA_real_, nrow = length(dgps), ncol = length(candidate_models),
                   dimnames = list(true = dgps, model = candidate_models))
  
  for (g in dgps) {
    dg <- d[d$dgp == g, , drop = FALSE]
    ll_mat[g, "EFE_full"]     <- mean(dg$raw_EFE_full[is.finite(dg$raw_EFE_full)])
    ll_mat[g, "no_identity"]  <- mean(dg$raw_no_identity[is.finite(dg$raw_no_identity)])
    ll_mat[g, "no_epistemic"] <- mean(dg$raw_no_epistemic[is.finite(dg$raw_no_epistemic)])
    ll_mat[g, "RI_cost"]      <- mean(dg$raw_RI_cost[is.finite(dg$raw_RI_cost)])
  }
  
  write.csv(ll_mat, file.path(OUT_TABLES, "W4_ll_matrix.csv"), row.names = TRUE)
  
  # summary: confusion counts + row-normalised shares + NA rate
  conf <- w4$confusion
  row_tot <- rowSums(conf)
  share <- sweep(conf, 1, row_tot, "/")
  
  na_rate <- tapply(is.na(d$winner_pen), d$dgp, mean)
  na_rate <- na_rate[dgps]
  
  w4_summary <- data.frame(
    dgp = rep(dgps, each = length(candidate_models)),
    selected = rep(candidate_models, times = length(dgps)),
    count = as.integer(t(conf)),
    share_row_norm = as.numeric(t(share)),
    stringsAsFactors = FALSE
  )
  
  # attach NA rate per row of dgp (repeated)
  w4_summary$share_na_winner_pen <- rep(as.numeric(na_rate), each = length(candidate_models))
  
  write.csv(w4_summary, file.path(OUT_TABLES, "W4_crossfit_summary.csv"), row.names = FALSE)
  
  # Figures
  if (isTRUE(make_figures)) {
    
    # Figure 17 (appendix): confusion matrix (counts + shares)
    local({
      png(file.path(OUT_FIG, "Figure 17_W4_Crossfit_Confusion.png"),
          width = 900, height = 600)
      op <- set_plot_style()
      on.exit({ par(op); dev.off() }, add = TRUE)
      
      par(font = 1, font.axis = 1, font.lab = 1, font.main = 1)
      par(xaxs = "i", yaxs = "i")
      
      plot_confusion(conf, main = "W4: Cross-fit confusion matrix (penalised; NA excluded)")
    })
    
    # Figure 16 (main): grouped bar plot of row-normalised selection shares
    plot_w4_bar(conf,
                file_png = file.path(OUT_FIG, "Figure 16_W4_Crossfit_Bar.png"),
                main = "W4: Cross-fit selection shares (penalised; row-normalised)")
  }
  
  invisible(w4)
}

# -------------------------------------------------------------
# OPTIONAL: Run once when sourcing (comment out for “defs-only”)
# -------------------------------------------------------------
# Example:
# w4 <- run_W4_and_write_outputs(N = 200, seed0 = 5000, make_figures = TRUE)