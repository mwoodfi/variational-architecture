############################################################
## 02_fit_cr_params.R
## W2: Likelihood + estimation helpers  +  single-subject runner (combined)
## - Provides: loglik_cr(), fit_cr(), fit_cr_fixbeta()
## - Also runs a minimal single-subject sanity check and writes outputs
##
## Outputs (src/outputs/{data,figures,tables}):
##  - data/W2_single_subject_data.rds
##  - tables/W2_single_subject_fit.csv
##  - data/W2_single_subject_fit.RData
##  - figures/Figure 13_W2_Single Subject Fit.png
############################################################

rm(list = ls())

# ----------------------------
# RNG standardisation
# ----------------------------
RNGkind(kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection')

# ----------------------------
# Dependencies
# ----------------------------

stopifnot(file.exists("Part-A_EFE Simulations.R"))
source("Part-A_EFE Simulations.R")

# ----------------------------
# Paths (relative to src/)
# ----------------------------
OUT_ROOT   <- "outputs"
OUT_DATA   <- file.path(OUT_ROOT, "data")
OUT_FIG    <- file.path(OUT_ROOT, "figures")
OUT_TABLES <- file.path(OUT_ROOT, "tables")

dir.create(OUT_ROOT,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DATA,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TABLES, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------

# Safe clipping to avoid log(0)
.clip01 <- function(p, eps = 1e-9) pmin(pmax(p, eps), 1 - eps)

# Large finite penalty (must be finite for L-BFGS-B)
.BIG <- 1e12

# Helper: return finite scalar or .BIG
.as_finite_scalar <- function(x, fallback = .BIG) {
  if (length(x) != 1L) return(fallback)
  if (!is.finite(x)) return(fallback)
  x
}

# Helper: basic parameter sanity (avoids NULL propagation -> numeric(0))
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

# Compute log-likelihood for theta on one subject dataset df
# theta = c(gamma_p, gamma_o, lambda, beta_choice)
loglik_cr <- function(theta, df, params, model = "EFE_full") {
  
  .assert_params_cr(params)
  
  gamma_p     <- theta[1]
  gamma_o     <- theta[2]
  lambda      <- theta[3]
  beta_choice <- theta[4]
  
  # Guard against pathological inputs
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

fit_cr <- function(df, params, model = "EFE_full",
                   init  = c(1, 1, 1, params$beta_choice),
                   inits = NULL,
                   lower = c(0.001, 0.001, 0.0001, 0.5),
                   upper = c(50, 50, 50, 30)) {
  
  .assert_params_cr(params)
  
  # Build candidate start points
  if (is.null(inits)) {
    inits <- rbind(
      init,
      c(0.5, 5, 2, 8),
      c(2, 10, 0.5, 4),
      c(1, 6, 1, params$beta_choice)
    )
  }
  inits <- as.matrix(inits)
  
  # Finite objective wrapper for L-BFGS-B
  nll <- function(th) {
    val <- -loglik_cr(th, df, params, model = model)
    .as_finite_scalar(val, fallback = .BIG)
  }
  
  best_val <- .BIG
  best <- NULL
  
  for (k in seq_len(nrow(inits))) {
    init_k <- as.numeric(inits[k, ])
    
    # If fn(init) is not finite, replace by safe midpoint of bounds
    v0 <- nll(init_k)
    if (!is.finite(v0) || v0 >= .BIG) {
      init_k <- (lower + upper) / 2
      v1 <- nll(init_k)
      if (!is.finite(v1) || v1 >= .BIG) next
    }
    
    opt <- optim(init_k, nll, method = "L-BFGS-B", lower = lower, upper = upper)
    
    if (is.finite(opt$value) && opt$value < best_val) {
      best_val <- opt$value
      best <- list(par = opt$par, value = opt$value, conv = opt$convergence, msg = opt$message)
    }
  }
  
  if (is.null(best)) {
    stop("fit_cr: all starting values produced non-finite objective. Check df and params consistency.")
  }
  
  best
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
  
  # Ensure finite start for L-BFGS-B
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

############################################################
## W2 runner (single-subject sanity check) -- appended here
############################################################

# helper: enforce required columns for likelihood and plotting
harmonise_df_for_ll <- function(df, params, seed_tag) {
  
  # choice column
  if (!"choice_B" %in% names(df)) {
    if ("choiceB" %in% names(df)) df$choice_B <- df$choiceB
    if ("choice"  %in% names(df)) df$choice_B <- df$choice
  }
  
  # harm observation column (aliases)
  if (!"harm_obs" %in% names(df)) {
    if ("harm"     %in% names(df)) df$harm_obs <- df$harm
    if ("harmObs"  %in% names(df)) df$harm_obs <- df$harmObs
    if ("harm_B"   %in% names(df)) df$harm_obs <- df$harm_B
    if ("obs"      %in% names(df)) df$harm_obs <- df$obs
    if ("y"        %in% names(df)) df$harm_obs <- df$y
    if ("outcome"  %in% names(df)) df$harm_obs <- df$outcome
  }
  
  # If still missing, stop: fabricating harm_obs breaks simulator–estimator coherence.
  if (!"harm_obs" %in% names(df)) {
    stop("harmonise_df_for_ll: harm_obs missing. W2 requires simulator-exported harm_obs for estimator–simulator coherence.")
  }
  
  if (!("choice_B" %in% names(df))) stop("harmonise_df_for_ll: choice_B missing after harmonisation.")
  if (!("harm_obs" %in% names(df))) stop("harmonise_df_for_ll: harm_obs missing after harmonisation.")
  
  df
}

# ----------------------------
# Simulation (single subject for estimation)
# - Must export harm_obs from the same run; otherwise likelihood belief-updating is incoherent.
# ----------------------------
simulate_cr_single_subject_w2 <- function(params,
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
      harm_obs = harm_obs  # **CHANGED**
    )
  }
  
  do.call(rbind, out)
}

# ----------------------------
# Simulation (single subject for estimation)
# ----------------------------

# 1) Simulate one subject (must include realised harm_obs)
df <- simulate_cr_single_subject_w2(
  params = params_catastrophic_recovery,
  omega  = list(gamma_p = 1, gamma_o = 6, lambda = 1),
  seed   = 2026
)

df <- harmonise_df_for_ll(df, params = params_catastrophic_recovery, seed_tag = 2026 + 22002)

saveRDS(df, file.path(OUT_DATA, "W2_single_subject_data.rds"))

# 2) Fit (fixed beta)
fit <- fit_cr_fixbeta(df, params_catastrophic_recovery, model = "EFE_full",
                      init = c(1, 6, 1))

tab <- data.frame(
  param    = c("gamma_p", "gamma_o", "lambda", "beta_choice"),
  estimate = as.numeric(fit$par),
  row.names = NULL
)

write.csv(tab, file.path(OUT_TABLES, "W2_single_subject_fit.csv"), row.names = FALSE)
save(fit, tab, file = file.path(OUT_DATA, "W2_single_subject_fit.RData"))

# 3) compute fitted p_B(t) implied by recovered parameters (not the simulated p_B)
# ----------------------------
compute_pB_path_given_choices <- function(df, params, omega_hat) {
  .assert_params_cr(params)
  
  T_total <- nrow(df)
  s <- df$s
  
  beliefs <- list(
    harm_aB = params$harm_prior_a, harm_bB = params$harm_prior_b,
    id_a    = params$id_prior_a,   id_b    = params$id_prior_b
  )
  
  pB_fit <- numeric(T_total)
  
  for (t in seq_len(T_total)) {
    GA <- cr_efe_policy(s[t], "A", beliefs, omega_hat, params)
    GB <- cr_efe_policy(s[t], "B", beliefs, omega_hat, params)
    
    pB_fit[t] <- cr_choice_prob(GA, GB, params$beta_choice)
    
    choiceB <- df$choice_B[t]
    harm_obs <- df$harm_obs[t]
    policy_chosen <- if (choiceB == 1) "B" else "A"
    
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
  
  pB_fit
}

omega_hat <- list(
  gamma_p = fit$par[1],
  gamma_o = fit$par[2],
  lambda  = fit$par[3]
)

df$p_B_fit <- compute_pB_path_given_choices(df, params_catastrophic_recovery, omega_hat)


# ----------------------------
# Plot (choice probability + realised choice)
# ----------------------------

png(file.path(OUT_FIG, "Figure 13_W2_Single Subject Fit.png"),
    width = 900, height = 600)

op <- set_plot_style()
par(font = 1, font.axis = 1, font.lab = 1, font.main = 1)
par(xaxs = "i", yaxs = "i")
on.exit(par(op), add = TRUE)

blue <- rgb(44, 100, 156, maxColorValue = 255)

x_ticks <- pretty(df$t, 5)
y_ticks_major <- seq(0, 1, by = 0.2)
y_ticks_minor <- seq(0, 1, by = 0.1)

plot(df$t, df$p_B_fit,
     type = "n",
     xlim = range(df$t),
     ylim = c(0, 1),
     axes = FALSE,
     xlab = "Time t",
     ylab = expression(paste("Selection Probability ", Pr(pi[B]))),
     main = "W2: Single-subject fit (fixed beta)")

draw_pub_grid(x_ticks = x_ticks,
              y_ticks_major = y_ticks_major,
              y_ticks_minor = y_ticks_minor)

lines(df$t, df$p_B_fit, lwd = 3.4, col = blue)
points(df$t, df$choice_B, pch = 16, cex = 0.7, col = blue)

axis(1, at = x_ticks, labels = x_ticks, lwd = 0, lwd.ticks = 1.35)
axis(2, at = y_ticks_major, labels = format(y_ticks_major, nsmall = 1),
     las = 1, lwd = 0, lwd.ticks = 1.35)

box(lwd = 1.15, col = "grey20")

dev.off()

cat("\nW2 outputs written to outputs/{data,figures,tables}/\n")