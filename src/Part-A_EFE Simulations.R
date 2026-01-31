############################################################
## efe_simulations.R
## Reference implementation for Part A: "Beyond Utility"
##
## Deterministic replication:
## - For full reproduction of all figures and data as used in
##   the paper, call:
##
##       run_all_simulations(seed = 2026)
##
## - The non-interactive entry point (Rscript efe_simulations.R)
##   automatically uses seed = 2026 and writes outputs/...
##
## - Experiments 4.2 and 4.3 are seed-isolated; the same seed reproduces
##   each experiment independently of call order.
##
## Structure:
## 0. Utility functions (Maths, Style & Divergences)
## 1. Generic EFE implementation (The Canonical Decomposition)
## 2. Experiment 4.1: Identity Hysteresis (Theorem 4 / App B.4)
## 3. Experiment 4.2: Emergent Probability Weighting (Prop 2 / App A.8)
## 4. Experiment 4.3: Catastrophic Recovery (App B.3 / App B.5)
## 5. Master Execution Function
## 6. Command Line Entry Point
############################################################

############################################################
## 0. Utility functions ------------------------------------
## Mathematical primitives used throughout the proofs.
############################################################

## Numerically stable softmax
## ALIGNMENT: Section 3.1
## Implements the Boltzmann-Luce choice rule for smooth attractors.
softmax <- function(x) {
  z <- x - max(x)
  ex <- exp(z)
  ex / sum(ex)
}

## Standard Sigmoid function
## Used for smooth transitions in dynamics (Exp 4.3)
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## Shannon entropy H(p) in nats
entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0L) return(0)
  -sum(p * log(p))
}

## General KL divergence D_KL(p || q)
kl_divergence <- function(p, q) {
  if (length(p) != length(q)) {
    stop("kl_divergence: p and q must have same length.")
  }
  idx <- (p > 0) & (q > 0)
  if (!any(idx)) return(0)
  sum(p[idx] * log(p[idx] / q[idx]))
}

## KL divergence for binary distributions with epsilon-guard
## Used for bounded entropy calculations in Appendix A.8
kl_binary <- function(p, q) {
  eps <- 1e-12
  p <- min(max(p, eps), 1 - eps)
  q <- min(max(q, eps), 1 - eps)
  p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
}

## KL divergence between univariate Gaussians
## ALIGNMENT: Appendix A.6 (Theorem 4)
## Used to calculate Identity Rigidity costs locally.
kl_gauss_1d <- function(mu_q, sig_q, mu_p, sig_p) {
  if (sig_q <= 0 || sig_p <= 0) {
    stop("kl_gauss_1d: sigmas must be > 0.")
  }
  log(sig_p / sig_q) +
    (sig_q^2 + (mu_q - mu_p)^2) / (2 * sig_p^2) -
    0.5
}

## Bernoulli entropy with epsilon-guard (nats)
entropy_binary_guarded <- function(p, eps = 1e-12) {
  p <- min(max(p, eps), 1 - eps)
  -p * log(p) - (1 - p) * log(1 - p)
}

############################################################
## Global plotting style (publication quality)
## ALIGNMENT: Figures 9, 10, 11
## Defines the visual standard for the manuscript.
############################################################

## Font fallback helper to ensure portability across systems
get_font_family <- function() {
  # Uses generic 'sans' for cross-platform portability.
  # Manuscript figures were rendered with Calibri on the author's system.
  "sans" 
}

set_plot_style <- function() {
  op <- par(no.readonly = TRUE)
  par(
    family  = get_font_family(),
    bg      = "white",
    fg      = "black",
    col     = "black",
    
    ## Default line aesthetics (curves can override)
    lwd     = 2.2,
    lend    = "round",
    ljoin   = "round",
    
    ## Typography / sizing
    cex     = 1.08,
    cex.axis= 1.05,
    cex.lab = 1.18,
    cex.main= 1.20,
    
    ## Margins and axis label placement
    mar     = c(4.6, 5.3, 1.5, 1.5),
    mgp     = c(2.9, 0.9, 0),
    tcl     = -0.30
  )
  op
}

## Draw a publication-style grid
draw_pub_grid <- function(x_ticks, y_ticks_major, y_ticks_minor = NULL,
                          col_major = "grey90", col_minor = "grey94",
                          lwd_major = 1.1, lwd_minor = 1.0) {
  
  ## Minor grid first (very faint)
  if (!is.null(y_ticks_minor)) {
    abline(h = y_ticks_minor, col = col_minor, lwd = lwd_minor, lty = 1)
  }
  
  ## Major horizontals
  abline(h = y_ticks_major, col = col_major, lwd = lwd_major, lty = 1)
  
  ## Major verticals exactly at x tick marks
  abline(v = x_ticks, col = col_minor, lwd = lwd_minor, lty = 1)
}

############################################################
## 1. Generic EFE components -------------------------------
##
## Implements the "Canonical Decomposition" (Section 2.2).
## All experiments are specific limiting regimes of this
## generic generative model structure.
############################################################

## 1.1 Environment specification ---------------------------

## Base outcome support: o ??? {0, 1}
outcomes <- c(0, 1)

## ALIGNMENT: Appendix A.2 / Theorem 1
## Constructs the Preference Prior via exponential tilting of utility:
## P(o) = exp(U(o)) / Z  (Eq A.2.2)
make_preference_prior <- function(pref_strength = 2.0) {
  u <- c(0, pref_strength)
  p <- softmax(u)
  names(p) <- outcomes
  p
}

make_state_prior <- function(bias = 0.5) {
  if (bias < 0 || bias > 1) {
    stop("make_state_prior: bias must be in [0,1].")
  }
  p <- c(1 - bias, bias)
  names(p) <- c("s0", "s1")
  p
}

make_obs_model <- function(p_succ_policy1 = c(0.2, 0.8),
                           p_succ_policy2 = c(0.8, 0.2)) {
  if (length(p_succ_policy1) != 2L || length(p_succ_policy2) != 2L) {
    stop("make_obs_model: p_succ_policy1 and p_succ_policy2 must be length 2.")
  }
  P_succ <- cbind(p_succ_policy1, p_succ_policy2)
  rownames(P_succ) <- c("s0", "s1")
  colnames(P_succ) <- c("pi1", "pi2")
  P_succ
}

## 1.2 Predictive outcome distribution Q(o | ??) ------------

predict_Q_o_given_pi <- function(pi_idx, P_s, P_succ) {
  if (!pi_idx %in% c(1L, 2L)) stop("predict_Q_o_given_pi: pi_idx must be 1 or 2.")
  if (length(P_s) != 2L) stop("predict_Q_o_given_pi: P_s must be length 2.")
  p1 <- sum(P_s * P_succ[, pi_idx])
  c(`0` = 1 - p1, `1` = p1)
}

## 1.3 Expected Information Gain IG(??) ---------------------
## ALIGNMENT: Appendix A.4 / Theorem 2
## Implements Eq (A.4.2): E_o [D_KL(Q(s|o,??) || Q(s|??))]
## Also equivalent to Value of Information (VoI).

compute_IG <- function(pi_idx, P_s, P_succ) {
  prior_s <- P_s
  
  Q_o <- predict_Q_o_given_pi(pi_idx, P_s, P_succ)
  
  IG <- 0
  for (o_val in outcomes) {
    p_o <- Q_o[as.character(o_val)]
    if (p_o <= 0) next
    
    if (o_val == 1) {
      lik <- P_succ[, pi_idx]
    } else {
      lik <- 1 - P_succ[, pi_idx]
    }
    
    post_unnorm <- prior_s * lik
    post_s <- post_unnorm / sum(post_unnorm)
    
    IG <- IG + p_o * kl_divergence(post_s, prior_s)
  }
  IG
}

## 1.4 Complexity term C(??) ??? generic placeholder ----------
## ALIGNMENT: Appendix A.5 (Rational Inattention) vs A.6 (Identity)
## Acts as a switch between 'Information Cost' (RI) and 'Identity Cost'.

compute_complexity_generic <- function(pi_idx,
                                       P_s,
                                       P_succ,
                                       complexity_mode = c("info", "none")) {
  mode <- match.arg(complexity_mode)
  if (mode == "none") {
    return(0)
  }
  if (mode == "info") {
    return(compute_IG(pi_idx, P_s, P_succ))
  }
  stop("compute_complexity_generic: Unknown complexity mode")
}

## 1.5 Expected Free Energy G(??; ??) ------------------------
## ALIGNMENT: Section 2.2 (The Canonical Decomposition)
## G(??) = Instr - ??_o * IG + ?? * C
##
## This function is the structural basis for Theorem 6 (Structural Completeness).

compute_EFE <- function(pi_idx,
                        P_o_pref,
                        P_s,
                        P_succ,
                        gamma_p,
                        gamma_o,
                        lambda,
                        complexity_mode = c("info", "none"),
                        identity_cost = 0) {
  complexity_mode <- match.arg(complexity_mode)
  
  Q_o <- predict_Q_o_given_pi(pi_idx, P_s, P_succ)
  
  E_neglogP <- sum(Q_o * (-log(P_o_pref)))
  H_Q <- entropy(Q_o)
  IG <- compute_IG(pi_idx, P_s, P_succ)
  C_pi <- compute_complexity_generic(pi_idx, P_s, P_succ, complexity_mode)
  
  G_core <- gamma_p * E_neglogP - gamma_p * H_Q - gamma_o * IG + lambda * C_pi
  G_core + identity_cost
}

## Special case for Identity Regime (Experiment 4.1)
## ALIGNMENT: Theorem 4 (Identity Economics)
compute_EFE_identity_regime <- function(
    pi_idx,
    P_o_pref,
    P_s,
    P_succ,
    gamma_p,
    identity_cost = 0
) {
  Q_o <- predict_Q_o_given_pi(pi_idx, P_s, P_succ)
  E_neglogP <- sum(Q_o * (-log(P_o_pref)))
  gamma_p * E_neglogP + identity_cost
}

identity_cost_for_policy <- function(policy_label, theta, lambda_id, sigma_q, sigma_p) {
  target <- if (policy_label == "pi1") -1.0 else 1.0
  lambda_id * kl_gauss_1d(
    mu_q  = theta,  sig_q = sigma_q,
    mu_p  = target, sig_p = sigma_p
  )
}

#############################################################
## 2. Experiment 4.1 ??? Identity hysteresis
##
## Implementation of Theorem 4 (Identity Economics Representation)
## and Appendix B.4 (Hysteresis and Path Dependence).
##
## Mathematical Alignment:
## - Identity Cost:   Appendix A.6 (Eq A.6.2 - KL Divergence)
## - Objective:       Section 3.6 (Instr. + Identity, gamma_o=0)
## - Dynamics:        Appendix B.4 (Slow prior adaptation)
############################################################

simulate_identity_hysteresis <- function(
    phi_seq_forward  = seq(0.0, 1.5, length.out = 100),
    phi_seq_backward = seq(1.5, 0.0, length.out = 100),
    gamma_p   = 8.0,
    lambda_id = 6.0,
    id_adapt_rate = 0.02,
    beta_plot = 10.0,
    sigma_q = 0.40,
    sigma_p = 1.60
) {
  P_o_pref <- make_preference_prior(pref_strength = 2.0)
  P_s      <- make_state_prior(bias = 0.5)
  
  ## Incentive calibration:
  ## Policy 1: "status quo" stable at 0.75.
  ## Policy 2: "challenger" improves with phi; saturates at 0.99.  
  make_P_succ_phi <- function(phi) {
    p1_val <- 0.75
    p2_val <- pmin(pmax(0.35 + 0.55 * phi, 0.01), 0.99)
    make_obs_model(c(p1_val, p1_val), c(p2_val, p2_val))
  }
  
  ## Identity state: -1 = ??1 basin, +1 = ??2 basin
  ## Represents the 'Deep Prior' P(theta) in Appendix A.6
  theta <- -1.0
  
  n_total <- length(phi_seq_forward) + length(phi_seq_backward)
  
  phi_path   <- numeric(n_total)
  theta_path <- numeric(n_total)
  p_piB_path <- numeric(n_total)
  direction  <- character(n_total)
  
  idx <- 1L
  
  ## One simulation step for a given phi
  run_step <- function(phi, current_theta) {
    P_succ <- make_P_succ_phi(phi)
    
    ## Policy-contingent posterior means for identity (simple contraction dynamics)
    ## Approximates the posterior update Q(theta|pi) described in Appendix A.6
    kappa <- 0.35
    mu_1 <- current_theta + kappa * (-1.0 - current_theta)
    mu_2 <- current_theta + kappa * ( 1.0 - current_theta)
    
    ## Identity rigidity costs: lambda_id * KL(Q || P)
    ## STRICT ALIGNMENT: Appendix A.6, Eq (A.6.2)
    ## Implements D_KL(Q(theta|pi) || P(theta)) using Gaussian proxy
    ## where 'current_theta' acts as the slowly adapting prior P(theta).
    ID1 <- lambda_id * kl_gauss_1d(mu_q = mu_1, sig_q = sigma_q,
                                   mu_p = current_theta, sig_p = sigma_p)
    ID2 <- lambda_id * kl_gauss_1d(mu_q = mu_2, sig_q = sigma_q,
                                   mu_p = current_theta, sig_p = sigma_p)
    
    ## Instrumental term + identity cost only
    ## STRICT ALIGNMENT: Section 3.6 / Theorem 4
    ## Enforces Condition ID4 (Epistemic Neutralisation, gamma_o = 0).
    ## G reduces to: -gamma_p * E[U] + lambda * D_KL
    G1 <- compute_EFE_identity_regime(
      pi_idx = 1,
      P_o_pref = P_o_pref,
      P_s = P_s,
      P_succ = P_succ,
      gamma_p = gamma_p,
      identity_cost = ID1
    )
    G2 <- compute_EFE_identity_regime(
      pi_idx = 2,
      P_o_pref = P_o_pref,
      P_s = P_s,
      P_succ = P_succ,
      gamma_p = gamma_p,
      identity_cost = ID2
    )
    
    ## Choice probabilities
    ## ALIGNMENT: Boltzmann???Luce/Gibbs choice mapping; decision precision ??
    probs  <- softmax(-beta_plot * c(G1, G2))
    prob_B <- probs[2]
    
    ## Identity update: slow drift towards the expected posterior mean
    ## ALIGNMENT: Appendix B.4 (Hysteresis and Path Dependence)
    ## The 'id_adapt_rate' creates the timescale separation required for hysteresis.
    mu_post_expect <- probs[1] * mu_1 + probs[2] * mu_2
    new_theta <- current_theta + id_adapt_rate * (mu_post_expect - current_theta)
    
    list(p_B = prob_B, theta = new_theta)
  }
  
  ## Forward path (increasing phi)
  for (phi in phi_seq_forward) {
    res <- run_step(phi, theta)
    theta <- res$theta
    
    phi_path[idx]   <- phi
    theta_path[idx] <- theta
    p_piB_path[idx] <- res$p_B
    direction[idx]  <- "Increasing (entry path)"
    idx <- idx + 1
  }
  
  ## Backward path (decreasing phi)
  for (phi in phi_seq_backward) {
    res <- run_step(phi, theta)
    theta <- res$theta
    
    phi_path[idx]   <- phi
    theta_path[idx] <- theta
    p_piB_path[idx] <- res$p_B
    direction[idx]  <- "Decreasing (exit path)"
    idx <- idx + 1
  }
  
  data.frame(
    phi       = phi_path,
    theta     = theta_path,
    prob_piB  = p_piB_path,
    direction = factor(direction,
                       levels = c("Increasing (entry path)", "Decreasing (exit path)"))
  )
}

## Polished Figure 9 style plot-------
plot_identity_hysteresis <- function(df,
                                     blue = rgb(44, 100, 156, maxColorValue = 255)) {
  
  op <- set_plot_style()
  on.exit(par(op), add = TRUE)
  
  ## Deterministic tick locations (publication feel)
  x_ticks <- c(0.0, 0.5, 1.0, 1.5)
  y_ticks_major <- seq(0, 1, by = 0.2)
  y_ticks_minor <- seq(0, 1, by = 0.1)
  
  ## Split data
  entry <- df[df$direction == "Increasing (entry path)", ]
  exit  <- df[df$direction == "Decreasing (exit path)", ]
  
  ## Canvas (hard cap at 1.0)
  plot(df$phi, df$prob_piB,
       type = "n",
       xlim = c(0, 1.5),
       ylim = c(0, 1.0),
       axes = FALSE,
       xlab = expression(paste("Incentive parameter ", phi)),
       ylab = expression(paste("Identity adherence (Pr(", pi[B], "))"))
  )
  
  ## Grid and frame
  draw_pub_grid(x_ticks, y_ticks_major, y_ticks_minor)
  box(lwd = 1.15, col = "grey20")
  
  ## Axes (black, clean)
  axis(1, at = x_ticks, labels = format(x_ticks, nsmall = 1),
       lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks_major, labels = format(y_ticks_major, nsmall = 1),
       las = 1, lwd = 0, lwd.ticks = 1.35)
  
  ## Curves
  lines(entry$phi, entry$prob_piB, col = blue, lwd = 3.4, lty = 1)
  lines(exit$phi,  exit$prob_piB,  col = blue, lwd = 3.0, lty = 2)
  
  ## Legend: let legend() manage layout (precise alignment), but style it
  par(xpd = NA)
  legend("topleft",
         inset = 0.02,
         legend = c(expression(paste("Increasing ", phi, " (entry path)")),
                    expression(paste("Decreasing ", phi, " (exit path)"))),
         col = blue,
         lty = c(1, 2),
         lwd = c(2.6, 2.4),
         bty = "o",
         bg = rgb(1, 1, 1, 0.72),
         box.col = "grey35",
         box.lwd = 1.0,
         cex = 1.00,
         x.intersp = 1.2,
         y.intersp = 1.05,
         seg.len = 2.6,
         text.col = "black")
  par(xpd = FALSE)
}

run_identity_hysteresis <- function() {
  df <- simulate_identity_hysteresis()
  plot_identity_hysteresis(df)
  invisible(df)
}

############################################################
## 3. Experiment 4.2 ??? Emergent probability weighting
##
## Implementation of Proposition 2 (Prospect Theory Asymmetries)
## and Appendix A.8 (Certainty Collapse & Mutual Information).
##
## Mathematical Alignment:
## - Epistemic Value: Appendix A.8 (Eq A.8.13 - Entropy Gap)
## - Objective:       Section 3.8 (Net epistemic coefficient)
## - Dynamics:        Prop 2(iii) (Transient weighting via learning)
##
## Determinism Note:
## - If seed is provided, full determinism (Calibration + MC) is guaranteed.
## - If seed is NULL, execution continues the current RNG stream.
############################################################

## Bernoulli entropy with epsilon-guard (nats): H_B(x)
H_B_guarded <- function(x, eps = 1e-12) {
  x <- min(max(x, eps), 1 - eps)
  -x * log(x) - (1 - x) * log(1 - x)
}

## E_{p ~ Beta(a,b)}[H_B(p)] in closed form using digamma
## STRICT ALIGNMENT: Appendix A.8, Eq (A.8.14)
E_H_Beta <- function(a, b) {
  if (a <= 0 || b <= 0) stop("E_H_Beta: a and b must be > 0.")
  ab <- a + b
  
  ## E[p ln p]
  E_p_log_p <- (a / ab) * (base::digamma(a + 1) - base::digamma(ab + 1))
  
  ## E[(1-p) ln(1-p)]
  E_1mp_log_1mp <- (b / ab) * (base::digamma(b + 1) - base::digamma(ab + 1))
  
  ## E[H_B(p)] = -E[p ln p] - E[(1-p) ln(1-p)]
  -(E_p_log_p + E_1mp_log_1mp)
}

## Mutual information IG_t(R) for Beta???Bernoulli
## STRICT ALIGNMENT: Appendix A.8, Eq (A.8.13)
## Calculates the exact "Jensen Gap" of the entropy function.
IG_BetaBernoulli_MI <- function(a, b, eps = 1e-12) {
  p_hat <- a / (a + b)
  ig <- H_B_guarded(p_hat, eps = eps) - E_H_Beta(a, b)
  max(0, ig)
}

## Internal: EU-clean instrumental term
## ALIGNMENT: Proposition 2 (Instrumental Baseline)
.compute_G_inst_R_EU <- function(a, b, m_safe, R_win, gamma_p) {
  p_hat <- a / (a + b)
  EU_R <- p_hat * R_win
  -gamma_p * EU_R
}

.compute_G_inst_S_EU <- function(m_safe, gamma_p) {
  -gamma_p * m_safe
}

## Internal: compute risky EFE given Beta(a,b)
## ALIGNMENT: Section 3.8 (Net Epistemic Coefficient)
.compute_G_R_pw <- function(a, b,
                            m_safe, R_win,
                            gamma_p, gamma_o, lambda,
                            IG_cache, eps = 1e-12) {
  
  G_inst_R <- .compute_G_inst_R_EU(
    a = a, b = b,
    m_safe = m_safe,
    R_win  = R_win,
    gamma_p = gamma_p
  )
  
  ## Cache key robust to non-integer a,b
  key <- paste0(format(a, digits = 12), "_", format(b, digits = 12))
  
  IG_R <- IG_cache[[key]]
  if (is.null(IG_R)) {
    IG_R <- IG_BetaBernoulli_MI(a, b, eps = eps)
    IG_cache[[key]] <- IG_R
  }
  
  ## The Canonical Reduction: G = Inst - gamma*IG + lambda*C
  ## Since C = IG here (Appendix A.8), term becomes -(gamma_o - lambda)*IG
  G_inst_R - (gamma_o - lambda) * IG_R
}

## Internal: Calibration offset Delta
## Ensures valid comparison at p_ref (Figure 10 reporting level)
.compute_calibration_delta_pw_choice_MI <- function(p_ref,
                                                    prior_a, prior_b,
                                                    m_safe, R_win,
                                                    gamma_p, gamma_o, lambda,
                                                    beta_choice,
                                                    IG_cache) {
  eps <- 1e-12
  p_ref <- min(max(p_ref, eps), 1 - eps)
  
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  
  ## Ensure caller RNG state is restored even if something errors inside calibration
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  
  ## Fixed local calibration seed (Appendix D.2)
  set.seed(2026)
  ## --- END PATCH ---
  
  G_S_const <- .compute_G_inst_S_EU(m_safe = m_safe, gamma_p = gamma_p)
  
  ## Calibration runs use a fixed local seed and restore the caller RNG state.
  n_cal_runs   <- 300L
  n_cal_trials <- 300L
  
  freq_given_delta <- function(Delta) {
    freq_runs <- numeric(n_cal_runs)
    
    for (r in seq_len(n_cal_runs)) {
      a <- prior_a
      b <- prior_b
      n_choose_R <- 0L
      
      for (t in seq_len(n_cal_trials)) {
        G_R <- .compute_G_R_pw(
          a = a, b = b,
          m_safe = m_safe,
          R_win  = R_win,
          gamma_p = gamma_p,
          gamma_o = gamma_o,
          lambda  = lambda,
          IG_cache = IG_cache
        ) + Delta
        
        pr <- softmax(-beta_choice * c(G_S_const, G_R))
        choose_R <- (runif(1) < pr[2])
        
        if (choose_R) {
          n_choose_R <- n_choose_R + 1L
          win <- (runif(1) < p_ref)
          if (win) a <- a + 1 else b <- b + 1
        }
      }
      
      freq_runs[r] <- n_choose_R / n_cal_trials
    }
    
    mean(freq_runs)
  }
  
  f_root <- function(Delta) freq_given_delta(Delta) - p_ref
  
  lo <- -50
  hi <-  50
  f_lo <- f_root(lo)
  f_hi <- f_root(hi)
  
  if (f_lo * f_hi > 0) {
    return(0)
  }
  
  uniroot(f_root, lower = lo, upper = hi, tol = 1e-3)$root
}

simulate_probability_weighting <- function(
    p_grid       = seq(0.01, 0.99, length.out = 120),
    m_safe       = 0.5,
    R_win        = 1.0,
    gamma_p      = 3.0,
    gamma_o      = 1.15,
    lambda       = 1.05,
    beta_choice = 6.0,
    
    ## Monte Carlo controls
    n_trials  = 800,
    n_mc_runs = 1000,
    
    ## Belief initialisation (Beta prior)
    prior_a = 2.0,
    prior_b = 2.0,
    
    ## Diagnostics toggle
    return_belief_diagnostics = TRUE,
    
    ## Progress reporting
    show_progress  = TRUE,
    progress_inner = FALSE,
    
    ## Explicit seeding for determinism
    seed = NULL
) {
  # Determinism: if seed is provided, set it here.
  # Otherwise, this function continues the RNG stream from where it was called.
  if (!is.null(seed)) set.seed(seed)
  
  if (R_win <= 0) stop("R_win must be > 0.")
  if (m_safe < 0) stop("m_safe must be >= 0.")
  if (any(p_grid <= 0 | p_grid >= 1)) stop("p_grid must be strictly within (0,1).")
  if (n_trials <= 0) stop("n_trials must be > 0.")
  if (n_mc_runs <= 0) stop("n_mc_runs must be > 0.")
  if (prior_a <= 0 || prior_b <= 0) stop("prior_a and prior_b must be > 0.")
  if (beta_choice <= 0) stop("beta_choice must be > 0.")
  
  p_ref <- m_safe / R_win
  if (!(p_ref > 0 && p_ref < 1)) {
    stop("EU threshold p_ref = m_safe / R_win must lie strictly in (0,1).")
  }
  
  mean_freq_R <- numeric(length(p_grid))
  sd_freq_R   <- numeric(length(p_grid))
  
  mean_p_hat_end <- rep(NA_real_, length(p_grid))
  sd_p_hat_end   <- rep(NA_real_, length(p_grid))
  
  IG_cache <- new.env(parent = emptyenv())
  
  ## Safe EFE constant (EU-clean; IG_S = 0)
  G_S_const <- .compute_G_inst_S_EU(m_safe = m_safe, gamma_p = gamma_p)
  
  ## Calibration offset
  ## Calibration uses a fixed local seed and restores the caller RNG state.
  Delta <- .compute_calibration_delta_pw_choice_MI(
    p_ref = p_ref,
    prior_a = prior_a, prior_b = prior_b,
    m_safe = m_safe, R_win = R_win,
    gamma_p = gamma_p,
    gamma_o = gamma_o,
    lambda  = lambda,
    beta_choice = beta_choice,
    IG_cache = IG_cache
  )
  
  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = length(p_grid), style = 3)
    on.exit({
      if (!is.null(pb)) close(pb)
    }, add = TRUE)
  }
  
  for (i in seq_along(p_grid)) {
    p_true <- p_grid[i]
    
    pb_in <- NULL
    if (isTRUE(show_progress) && isTRUE(progress_inner)) {
      pb_in <- utils::txtProgressBar(min = 0, max = n_mc_runs, style = 3)
      on.exit({
        if (!is.null(pb_in)) close(pb_in)
      }, add = TRUE)
    }
    
    freq_runs       <- numeric(n_mc_runs)
    p_hat_end_runs <- numeric(n_mc_runs)
    
    for (r in seq_len(n_mc_runs)) {
      a <- prior_a
      b <- prior_b
      n_choose_R <- 0L
      
      ## Inner Loop: Sequential Learning
      ## ALIGNMENT: Proposition 2(iii) - Certainty Collapse
      ## As 't' increases, 'a' and 'b' grow, causing IG -> 0
      ## and behavior to converge to EU baseline.
      for (t in seq_len(n_trials)) {
        G_R <- .compute_G_R_pw(
          a = a, b = b,
          m_safe = m_safe,
          R_win  = R_win,
          gamma_p = gamma_p,
          gamma_o = gamma_o,
          lambda  = lambda,
          IG_cache = IG_cache
        ) + Delta
        
        pr <- softmax(-beta_choice * c(G_S_const, G_R))
        choose_R <- (runif(1) < pr[2])
        
        if (choose_R) {
          n_choose_R <- n_choose_R + 1L
          win <- (runif(1) < p_true)
          if (win) a <- a + 1 else b <- b + 1
        }
      }
      
      freq_runs[r]       <- n_choose_R / n_trials
      p_hat_end_runs[r] <- a / (a + b)
      
      if (!is.null(pb_in)) utils::setTxtProgressBar(pb_in, r)
    }
    
    if (!is.null(pb_in)) close(pb_in)
    
    mean_freq_R[i] <- mean(freq_runs)
    sd_freq_R[i]   <- stats::sd(freq_runs)
    
    if (return_belief_diagnostics) {
      mean_p_hat_end[i] <- mean(p_hat_end_runs)
      sd_p_hat_end[i]   <- stats::sd(p_hat_end_runs)
    }
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
  }
  
  out <- data.frame(
    p           = p_grid,
    freq_R      = mean_freq_R,
    sd_R        = sd_freq_R,
    p_thresh_EU = p_ref
  )
  
  if (return_belief_diagnostics) {
    out$p_hat_end_mean <- mean_p_hat_end
    out$p_hat_end_sd   <- sd_p_hat_end
  }
  
  out
}

plot_probability_weighting <- function(df,
                                       show_p_hat_end = FALSE,
                                       blue = rgb(44, 100, 156, maxColorValue = 255)) {
  if (!all(c("p", "freq_R", "p_thresh_EU") %in% names(df))) {
    stop("plot_probability_weighting: df must contain columns p, freq_R, p_thresh_EU.")
  }
  if (show_p_hat_end && !"p_hat_end_mean" %in% names(df)) {
    stop("plot_probability_weighting: p_hat_end_mean not found in df. Re-run simulate_probability_weighting(return_belief_diagnostics = TRUE).")
  }
  
  op <- set_plot_style()
  on.exit(par(op), add = TRUE)
  
  ## Deterministic tick locations (publication feel)
  x_ticks <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
  y_ticks_major <- seq(0, 1, by = 0.2)
  y_ticks_minor <- seq(0, 1, by = 0.1)
  
  ## Canvas
  plot(df$p, df$freq_R,
       type = "n",
       xlim = c(0, 1),
       ylim = c(0, 1),
       axes = FALSE,
       xlab = "Objective probability p",
       ylab = "Frequency of choosing risky option",
       main = "")
  
  ## Grid and frame
  draw_pub_grid(x_ticks, y_ticks_major, y_ticks_minor)
  box(lwd = 1.15, col = "grey20")
  
  ## Axes (black, clean)
  axis(1, at = x_ticks, labels = format(x_ticks, nsmall = 1),
       lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks_major, labels = format(y_ticks_major, nsmall = 1),
       las = 1, lwd = 0, lwd.ticks = 1.35)
  
  ## Curves and reference lines (annotations unchanged)
  lines(df$p, df$freq_R, col = blue, lwd = 3.4, lty = 1)
  lines(df$p, df$p, lty = 2, lwd = 2.2)
  abline(v = df$p_thresh_EU[1L], lty = 3, lwd = 2.2)
  
  if (show_p_hat_end) {
    par(new = TRUE)
    plot(df$p, df$p_hat_end_mean,
         type = "l", lty = 4,
         axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4, at = seq(0, 1, by = 0.2))
    mtext("Mean end-belief p_hat_end", side = 4, line = 2.5)
  }
  
  ## Legend styling and placement matched to Experiment 4.1
  par(xpd = NA)
  legend("topleft",
         inset = 0.02,
         legend = c(
           "EFE simulated",
           "Linear benchmark (y = p)",
           "EU threshold",
           if (show_p_hat_end) "Mean end-belief p_hat_end" else NULL
         ),
         col = c(blue, "black", "black", if (show_p_hat_end) "black" else NULL),
         lty = c(1, 2, 3, if (show_p_hat_end) 4 else NULL),
         lwd = c(2.6, 2.4, 2.4, if (show_p_hat_end) 2.2 else NULL),
         bty = "o",
         bg = rgb(1, 1, 1, 0.72),
         box.col = "grey35",
         box.lwd = 1.0,
         cex = 1.00,
         x.intersp = 1.2,
         y.intersp = 1.05,
         seg.len = 2.6,
         text.col = "black")
  par(xpd = FALSE)
}

run_probability_weighting <- function(show_p_hat_end = FALSE,
                                      show_progress = TRUE,
                                      progress_inner = FALSE,
                                      seed = 2026) {
  
  # By default, use a fixed seed to ensure the entire run (Calibration + MC)
  # is fully deterministic and independent of prior RNG states.
  df <- simulate_probability_weighting(
    return_belief_diagnostics = show_p_hat_end,
    show_progress  = show_progress,
    progress_inner = progress_inner,
    seed = seed 
  )
  plot_probability_weighting(df, show_p_hat_end = show_p_hat_end)
  invisible(df)
}

############################################################
## 4. Experiment 4.3 ??? Catastrophic Recovery
##
## Implements the "Phase Transition" mechanism described in 
## Section 3 and formally characterized in Appendix B.3 (Bifurcations).
##
## Mathematical Alignment:
## - Epistemic Value: Eq (A.8.13) "Entropy Gap" (exact digamma solution)
## - Identity Cost:   Appendix A.6 (KL-divergence rigidity)
## - Dynamics:        Appendix B.5 (Slow hyperparameter adaptation)
##
## IMPLEMENTATION PATCH (behavioural activation of identity term):
## - Introduces params$id_cost_scale and applies it to the identity KL term.
## - Keeps all other logic unchanged.
## - Adds optional output columns for term magnitudes (instr_B, epi_B, comp_B, lambda_dyn).
############################################################

# ----------------------------
# Helper maths (Namespaced to Exp 4.3)
# ----------------------------

# Exact KL divergence for Beta distributions
# Used for Identity Rigidity Term (Appendix A.6)
cr_kl_beta <- function(a, b, c, d) {
  lbeta(c, d) - lbeta(a, b) +
    (a - c) * (digamma(a) - digamma(a + b)) +
    (b - d) * (digamma(b) - digamma(a + b))
}

# Mutual information for Beta???Bernoulli in "entropy gap" form
# EXACT implementation of Appendix A.8, Eq (A.8.13) and (A.8.14)
cr_ig_beta_bernoulli <- function(a, b) {
  phat <- a / (a + b)
  
  # Bernoulli entropy H_B(p)
  H_B <- function(x) {
    x <- pmin(pmax(x, 1e-12), 1 - 1e-12)
    -x * log(x) - (1 - x) * log(1 - x)
  }
  
  # Expected Entropy E[H_B(p)] using digamma (Eq A.8.14)
  Ep_logp   <- (a / (a + b)) * (digamma(a + 1) - digamma(a + b + 1))
  E1m_log1m <- (b / (a + b)) * (digamma(b + 1) - digamma(a + b + 1))
  EH_B      <- -Ep_logp - E1m_log1m
  
  # IG = Entropy - Expected Entropy (The Jensen Gap)
  H_B(phat) - EH_B
}

# Softmax/Sigmoid Choice Rule (Appendix A.8.18)
cr_choice_prob <- function(GA, GB, beta) {
  sigmoid(beta * (GA - GB))
}

# ----------------------------
# Model components
# ----------------------------

# Instrumental Cost (The Utility Interface - Appendix A.2)
cr_instr_cost <- function(s, policy, params) {
  if (policy == "A") {
    # Policy A (Addiction): Fixed low cost (short-term utility)
    params$k_A0
  } else {
    # Policy B (Health): Cost depends on state s
    params$k_B0 + params$k_state * s
  }
}

# Identity Rigidity (Appendix A.6)
# PATCH: scale identity KL so lambda has behavioural leverage in the choice kernel.
cr_identity_cost <- function(id_a, id_b, params) {
  params$id_cost_scale * cr_kl_beta(id_a, id_b, params$id_prior_a, params$id_prior_b)
}

cr_epistemic_value <- function(a, b) {
  cr_ig_beta_bernoulli(a, b)
}

# Hyperparameter Dynamics (Appendix B.5 - Slow Adaptation)
cr_omega_path <- function(t, params) {
  # Exogenous rise in Insight (Gamma_o)
  gamma_o_t <- params$gamma_o_min +
    (params$gamma_o_max - params$gamma_o_min) *
    sigmoid((t - params$t_event) / params$gamma_o_ramp)
  
  list(
    gamma_p = params$gamma_p,
    gamma_o = gamma_o_t,
    lambda  = params$lambda_base
  )
}

# Expected Free Energy Evaluation
# Implements the Canonical Decomposition: G = Instr - gamma*IG + lambda*C
cr_efe_policy <- function(s, policy, beliefs, omega, params) {
  
  # 1. Instrumental Term (Preference)
  instr <- omega$gamma_p * cr_instr_cost(s, policy, params)
  
  # 2. Epistemic Term (Information Gain)
  IG <- if (policy == "A") 0 else cr_epistemic_value(beliefs$harm_aB, beliefs$harm_bB)
  epi <- -omega$gamma_o * IG
  
  # 3. Complexity Term (Identity Rigidity)
  comp <- if (policy == "B") {
    omega$lambda * cr_identity_cost(beliefs$id_a, beliefs$id_b, params)
  } else {
    0
  }
  
  instr + epi + comp
}

# ----------------------------
# Core simulator (internal)
# ----------------------------

.simulate_catastrophic_recovery_core <- function(params) {
  set.seed(params$seed)
  
  s <- numeric(params$T_total)
  s[1] <- params$s0
  
  beliefs <- list(
    harm_aB = params$harm_prior_a, harm_bB = params$harm_prior_b,
    id_a    = params$id_prior_a,   id_b    = params$id_prior_b
  )
  
  out <- vector("list", params$T_total)
  lambda_dyn <- params$lambda0
  
  for (t in seq_len(params$T_total)) {
    om <- cr_omega_path(t, params)
    
    # PATCH (bookkeeping only): keep static lambda separate from dynamic barrier
    om$lambda_static <- om$lambda
    om$lambda <- om$lambda_static + lambda_dyn
    
    # Calculate EFE for both policies
    GA <- cr_efe_policy(s[t], "A", beliefs, om, params)
    GB <- cr_efe_policy(s[t], "B", beliefs, om, params)
    
    # Selection (Phase Transition Mechanism - Appendix B.3)
    pB <- cr_choice_prob(GA, GB, params$beta_choice)
    choiceB <- rbinom(1, 1, pB)  # 1 => choose B, 0 => choose A
    policy_chosen <- if (choiceB == 1) "B" else "A"
    
    # Environment Feedback
    p_harm_true <- if (policy_chosen == "A") params$p_harm_A_true else params$p_harm_B_true
    harm_obs <- rbinom(1, 1, p_harm_true)
    
    # Belief Updating (Bayesian Learning)
    if (policy_chosen == "B") {
      beliefs$harm_aB <- beliefs$harm_aB + harm_obs
      beliefs$harm_bB <- beliefs$harm_bB + (1 - harm_obs)
    }
    
    # Identity Accumulation (Hysteresis effect - Appendix B.4)
    if (policy_chosen == "A") {
      beliefs$id_a <- beliefs$id_a + params$id_update_strength
    } else {
      beliefs$id_b <- beliefs$id_b + params$id_update_strength * params$id_recovery_ratio
    }
    
    # Dynamic Rigidity (Lambda)
    lambda_dyn <- lambda_dyn +
      (policy_chosen == "A") * params$lambda_gain_A -
      (policy_chosen == "B") * params$lambda_relax_B
    lambda_dyn <- lambda_dyn -
      params$lambda_event_drop * exp(-0.5 * ((t - params$t_event) / params$lambda_event_width)^2)
    lambda_dyn <- max(lambda_dyn, 0)
    
    # State Evolution (Drift-Diffusion)
    if (t < params$T_total) {
      drift <- if (policy_chosen == "A") params$s_drift_A else -params$s_drift_B
      s[t + 1] <- pmin(pmax(s[t] + drift + rnorm(1, 0, params$s_noise), 0), 1)
    }
    
    # Optional term magnitudes for referee diagnostics (no behavioural effect)
    IG_B <- cr_epistemic_value(beliefs$harm_aB, beliefs$harm_bB)
    instr_B <- om$gamma_p * cr_instr_cost(s[t], "B", params)
    epi_B   <- -om$gamma_o * IG_B
    comp_B  <- om$lambda * cr_identity_cost(beliefs$id_a, beliefs$id_b, params)
    
    out[[t]] <- data.frame(
      t = t,
      s = s[t],
      gamma_p = om$gamma_p,
      gamma_o = om$gamma_o,
      lambda  = om$lambda,
      lambda_static = om$lambda_static,
      lambda_dyn = lambda_dyn,
      G_A = GA,
      G_B = GB,
      DeltaG = GA - GB,
      p_B = pB,
      choice_B = choiceB,
      harm_obs = harm_obs,
      harm_mean_B = beliefs$harm_aB / (beliefs$harm_aB + beliefs$harm_bB),
      id_a = beliefs$id_a,
      id_b = beliefs$id_b,
      id_KL = cr_identity_cost(beliefs$id_a, beliefs$id_b, params),
      instr_B = instr_B,
      epi_B = epi_B,
      comp_B = comp_B
    )
  }
  
  do.call(rbind, out)
}

# ----------------------------
# Default parameters
# ----------------------------

params_catastrophic_recovery <- list(
  seed = 2026,
  T_total = 220L,
  s0 = 0.95,
  s_pref = 0.0,
  beta_choice = 8.0,
  
  # Initial Rigidity (Barrier Height)
  lambda0 = 0.05,
  
  # Insight Dynamics (The Trigger)
  gamma_o_min = 0.1,
  gamma_o_max = 18.0,
  
  # Instrumental Landscape (The Trap)
  k_A0 = -1.0,
  k_B0 = -1.5,
  k_state = 1.1,
  
  gamma_p = 1.0,
  t_event = 80L,
  gamma_o_ramp = 15,
  
  # Environment Statistics
  p_harm_A_true = 0.28,
  p_harm_B_true = 0.08,
  
  # Priors
  harm_prior_a = 1,
  harm_prior_b = 1,
  id_prior_a = 6,
  id_prior_b = 6,
  
  # Dynamics constants
  id_update_strength = 0.01,
  id_recovery_ratio = 0.5,
  
  # Identity-cost normalisation (PATCH)
  # Increases the behavioural leverage of the identity rigidity term without redefining lambda.
  id_cost_scale = 25,
  
  lambda_base = 0.1,
  lambda_gain_A = 0.01,
  lambda_relax_B = 0.05,
  lambda_event_drop = 0.1,
  lambda_event_width = 15,
  
  s_drift_A = 0.01,
  s_drift_B = 0.05,
  s_noise = 0.01
)

# Public entry point
simulate_catastrophic_recovery <- function(params = params_catastrophic_recovery) {
  .simulate_catastrophic_recovery_core(params)
}

# ----------------------------
# Plotting (Fully self-contained, Modular, and Publication Quality)
# ----------------------------

# Individual Panel Plotters (Factored out as requested for cleaner code)
plot_cr_panel_a <- function(df, blue) {
  # --- Panel 11.A: Gamma_o (Insight) ---
  y_ticks_A <- seq(0, 20, by = 5)
  ylim_A <- c(0, max(df$gamma_o) * 1.05)
  t_ticks <- seq(0, 250, by = 50)
  
  plot(df$t, df$gamma_o, type = "n",
       xlim = range(df$t), ylim = ylim_A,
       axes = FALSE,
       xlab = "Time t",
       ylab = expression(paste("Epistemic Precision ", gamma[o])),
       main = expression(paste("11.A: Exogenous Epistemic Precision (Insight) ", gamma[o](t))))
  
  draw_pub_grid(t_ticks, y_ticks_A)
  lines(df$t, df$gamma_o, lwd = 3.0, col = blue)
  axis(1, at = t_ticks, labels = t_ticks, lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks_A, labels = y_ticks_A, las = 1, lwd = 0, lwd.ticks = 1.35)
  box(lwd = 1.15, col = "grey20")
}

plot_cr_panel_b <- function(df, blue) {
  # --- Panel 11.B: Latent Mechanism ---
  y_ticks_prob <- seq(0, 1, by = 0.2)
  t_ticks <- seq(0, 250, by = 50)
  
  plot(df$t, df$p_B, type = "n", ylim = c(0, 1),
       axes = FALSE,
       xlab = "Time t",
       ylab = expression(paste("Selection Probability ", P(pi[B]))),
       main = expression(paste("11.B: Latent Mechanism (Phase Transition) ", Pr(pi[B] == "Health")(t))))
  
  draw_pub_grid(t_ticks, y_ticks_prob)
  abline(h = 0.5, lty = 2, col = "grey40", lwd = 1.5)
  lines(df$t, df$p_B, lwd = 3.0, col = blue)
  axis(1, at = t_ticks, labels = t_ticks, lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks_prob, labels = format(y_ticks_prob, nsmall = 1), las = 1, lwd = 0, lwd.ticks = 1.35)
  box(lwd = 1.15, col = "grey20")
}

plot_cr_panel_c <- function(df, blue) {
  # --- Panel 11.C: Observable Behaviour ---
  y_ticks_prob <- seq(0, 1, by = 0.2)
  t_ticks <- seq(0, 250, by = 50)
  
  plot(df$t, df$p_B, type = "n", ylim = c(0, 1),
       axes = FALSE,
       xlab = "Time (Learning Episodes) t",
       ylab = "Probability / realised choice",
       main = expression(paste("11.C: Observable Behavioural Recovery - Choices on ", Pr(pi[B])(t))))
  
  draw_pub_grid(t_ticks, y_ticks_prob)
  abline(h = 0.5, lty = 2, col = "grey40", lwd = 1.5)
  lines(df$t, df$p_B, lwd = 2.8, col = "grey60")
  points(df$t, df$choice_B, pch = 16, cex = 0.7, col = blue)
  axis(1, at = t_ticks, labels = t_ticks, lwd = 0, lwd.ticks = 1.35)
  axis(2, at = y_ticks_prob, labels = format(y_ticks_prob, nsmall = 1), las = 1, lwd = 0, lwd.ticks = 1.35)
  box(lwd = 1.15, col = "grey20")
}

# Composite Plotter (Orchestrates the panels into Figure 11)
plot_catastrophic_recovery <- function(df,
                                       blue = rgb(44, 100, 156, maxColorValue = 255)) {
  
  op <- set_plot_style()
  par(mfrow = c(3, 1), mar = c(4.0, 5.2, 2.2, 1.2))
  on.exit(par(op), add = TRUE)
  
  plot_cr_panel_a(df, blue)
  plot_cr_panel_b(df, blue)
  plot_cr_panel_c(df, blue)
}

############################################################
## 5. Master Execution Function ---------------------------
############################################################

run_all_simulations <- function(output_dir = "outputs", seed = 2026) {
  # Standardize RNG for exact cross-platform replication
  # Requires R >= 3.6 (for RNGkind(sample.kind="Rejection"))
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  figs_dir <- file.path(output_dir, "figures")
  if (!dir.exists(figs_dir)) dir.create(figs_dir)
  data_dir <- file.path(output_dir, "data")
  if (!dir.exists(data_dir)) dir.create(data_dir)
  
  message("Starting Experiment 4.1: Identity Hysteresis...")
  # Exp 4.1 is naturally deterministic (no random draws)
  df_hys <- simulate_identity_hysteresis()
  saveRDS(df_hys, file.path(data_dir, "identity_hysteresis.rds"))
  png(file.path(figs_dir, "Figure 9_Identity Hysteresis.png"), width = 900, height = 600)
  plot_identity_hysteresis(df_hys)
  dev.off()
  
  message("Starting Experiment 4.2: Probability Weighting...")
  # Explicitly passing seed ensures MC runs are seed-isolated modules
  df_pw <- simulate_probability_weighting(seed = seed)
  saveRDS(df_pw, file.path(data_dir, "probability_weighting.rds"))
  png(file.path(figs_dir, "Figure 10_Probability Weighting.png"), width = 900, height = 600)
  plot_probability_weighting(df_pw)
  dev.off()
  
  message("Starting Experiment 4.3: Catastrophic Recovery...")
  # Explicitly passing seed ensures seed-isolated determinism
  p_cr <- params_catastrophic_recovery
  if (!is.null(seed)) p_cr$seed <- seed
  df_cr <- simulate_catastrophic_recovery(params = p_cr)
  saveRDS(df_cr, file.path(data_dir, "catastrophic_recovery.rds"))
  
  # --- Exp 4.3: save three separate panels (preferred) ---
  blue <- rgb(44, 100, 156, maxColorValue = 255)
  
  png(file.path(figs_dir, "Figure 11A_Catastrophic Recovery.png"), width = 900, height = 600)
  op <- set_plot_style(); 
  plot_cr_panel_a(df_cr, blue)
  par(op)
  dev.off()
  
  png(file.path(figs_dir, "Figure 11B_Catastrophic Recovery.png"), width = 900, height = 600)
  op <- set_plot_style(); 
  plot_cr_panel_b(df_cr, blue)
  par(op)
  dev.off()
  
  png(file.path(figs_dir, "Figure 11C_Catastrophic Recovery.png"), width = 900, height = 600)
  op <- set_plot_style(); 
  plot_cr_panel_c(df_cr, blue)
  par(op)
  dev.off()
  
  # Optional: also save the manuscript-style composite
  png(file.path(figs_dir, "Figure 11_Catastrophic Recovery vComposite.png"), width = 900, height = 1013)
  plot_catastrophic_recovery(df_cr, blue = blue)
  dev.off()
  
  invisible(list(
    hysteresis = df_hys,
    prob_weighting = df_pw,
    catastrophic_recovery = df_cr
  ))
}

############################################################
## 6. Execution Entry Point -------------------------------
############################################################

if (!interactive()) {
  message("Running all simulations with fixed seed (2026)...")
  run_all_simulations(seed = 2026)
  message("Done. Results saved to ./outputs/")
}