############################################################
## 01_run_W1_constant_omega.R
## W1: Constant-Omega catastrophic recovery (illustrative)
############################################################

rm(list = ls())

stopifnot(file.exists("Part-A_EFE Simulations.R"))

# RNG standardisation explicitly (referee-safe determinism).
RNGkind(kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")

# ----------------------------
# Paths (relative to src/)
# ----------------------------
OUT_ROOT <- "outputs"
OUT_DATA <- file.path(OUT_ROOT, "data")
OUT_FIG  <- file.path(OUT_ROOT, "figures")
OUT_TAB  <- file.path(OUT_ROOT, "tables")

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB,  showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Dependencies
# ----------------------------
source("Part-A_EFE Simulations.R")

stopifnot(
  exists("params_catastrophic_recovery"),
  exists("cr_efe_policy"),
  exists("cr_choice_prob"),
  exists("set_plot_style"),
  exists("draw_pub_grid")
)

# ----------------------------
# Simulation
# ----------------------------
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
      choice_B = choiceB
    )
  }
  
  do.call(rbind, out)
}

# ----------------------------
# Run W1
# ----------------------------
df_w1 <- simulate_cr_constant_omega(
  params = params_catastrophic_recovery,
  omega  = list(gamma_p = 1, gamma_o = 6, lambda = 1),
  seed   = 2026
)

saveRDS(df_w1, file.path(OUT_DATA, "W1_constant_omega.rds"))

# ----------------------------
# Export referee-facing summaries (W1)
# ----------------------------

# 1) Numeric summary matching Appendix E.1 exactly (6 rows x 4 variables)
w1_summary <- data.frame(
  Parameter = c("Min", "1st Quartile", "Median", "Mean", "3rd Quartile", "Max"),
  t        = c(min(df_w1$t),
               as.numeric(quantile(df_w1$t, 0.25, names = FALSE)),
               median(df_w1$t),
               mean(df_w1$t),
               as.numeric(quantile(df_w1$t, 0.75, names = FALSE)),
               max(df_w1$t)),
  s        = c(min(df_w1$s),
               as.numeric(quantile(df_w1$s, 0.25, names = FALSE)),
               median(df_w1$s),
               mean(df_w1$s),
               as.numeric(quantile(df_w1$s, 0.75, names = FALSE)),
               max(df_w1$s)),
  p_B      = c(min(df_w1$p_B),
               as.numeric(quantile(df_w1$p_B, 0.25, names = FALSE)),
               median(df_w1$p_B),
               mean(df_w1$p_B),
               as.numeric(quantile(df_w1$p_B, 0.75, names = FALSE)),
               max(df_w1$p_B)),
  choice_B = c(min(df_w1$choice_B),
               as.numeric(quantile(df_w1$choice_B, 0.25, names = FALSE)),
               median(df_w1$choice_B),
               mean(df_w1$choice_B),
               as.numeric(quantile(df_w1$choice_B, 0.75, names = FALSE)),
               max(df_w1$choice_B)),
  row.names = NULL
)

write.csv(
  w1_summary,
  file = file.path(OUT_TAB, "W1_constant_omega_summary.csv"),
  row.names = FALSE
)

# 2) Structure/metadata table (referee-facing)
w1_structure <- data.frame(
  variable      = names(df_w1),
  class         = vapply(df_w1, function(x) paste(class(x), collapse = "/"), character(1)),
  n_obs         = rep.int(nrow(df_w1), length(df_w1)),
  n_nonmissing  = vapply(df_w1, function(x) sum(!is.na(x)), integer(1)),
  n_missing     = vapply(df_w1, function(x) sum(is.na(x)), integer(1)),
  stringsAsFactors = FALSE
)

# Add min/max where meaningful (numeric/integer/logical)
is_num_like <- vapply(df_w1, function(x) is.numeric(x) || is.integer(x) || is.logical(x), logical(1))
w1_structure$min <- NA_real_
w1_structure$max <- NA_real_
w1_structure$min[is_num_like] <- vapply(df_w1[is_num_like], function(x) suppressWarnings(min(as.numeric(x), na.rm = TRUE)), numeric(1))
w1_structure$max[is_num_like] <- vapply(df_w1[is_num_like], function(x) suppressWarnings(max(as.numeric(x), na.rm = TRUE)), numeric(1))

write.csv(
  w1_structure,
  file = file.path(OUT_TAB, "W1_constant_omega_structure.csv"),
  row.names = FALSE
)

# ----------------------------
# Plot
# ----------------------------
png(file.path(OUT_FIG, "Figure 12_W1_Constant Omega.png"),
    width = 900, height = 600)

op <- set_plot_style()
par(font = 1, font.axis = 1, font.lab = 1, font.main = 1)
par(xaxs = "i", yaxs = "i")
on.exit(par(op), add = TRUE)

blue <- rgb(44, 100, 156, maxColorValue = 255)

x_ticks <- pretty(df_w1$t, 5)
y_ticks_major <- seq(0, 1, by = 0.2)
y_ticks_minor <- seq(0, 1, by = 0.1)

plot(df_w1$t, df_w1$p_B,
     type = "n",
     axes = FALSE,
     xlab = "Time t",
     ylab = expression(paste("Selection Probability ", Pr(pi[B]))),
     main = "W1: Constant-Ω recovery dynamics")

draw_pub_grid(x_ticks = x_ticks,
              y_ticks_major = y_ticks_major,
              y_ticks_minor = y_ticks_minor)

lines(df_w1$t, df_w1$p_B, lwd = 3.4, col = blue)
points(df_w1$t, df_w1$choice_B, pch = 16, cex = 0.7, col = blue)

axis(1, at = x_ticks, labels = x_ticks, lwd = 0, lwd.ticks = 1.35)
axis(2, at = y_ticks_major, labels = format(y_ticks_major, nsmall = 1),
     las = 1, lwd = 0, lwd.ticks = 1.35)

box(lwd = 1.15, col = "grey20")

dev.off()