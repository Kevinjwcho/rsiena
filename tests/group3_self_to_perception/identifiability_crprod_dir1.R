############################################################
## simulation_crprod_identifiability.R
##
## Purpose: Empirical identifiability test for the two
##   cross-network directions.
##
##   Direction 1 (D1)  — perceptual accuracy, lives in
##     Y[shared]'s objective (perceived slice objective):
##       s_j(Y[i]) = Σ_k x_{ijk} * x_{jjk}
##     includeEffects(eff, crprod,
##                    interaction1 = "Y[self]",
##                    name         = "Y[shared]")
##
##   Direction 2 (D2)  — perception shaping reality, lives in
##     Y[self]'s objective (aggregate over all K perceivers):
##       s_j(Y[self]) = Σ_m Σ_k x_{jjk} * x_{mjk}
##     includeEffects(eff, crprod,
##                    interaction1 = "Y[shared]")   # name=Y[self] auto
##
##   DGP used here (asymmetric):
##     β_D1_true  = 0.6      (same as simulation_crprod_recovery.R)
##     β_D2_true  = 0        ← IDENTIFIABILITY PROBE
##
##   Estimation includes BOTH D1 and D2.  Success criteria:
##     (i)  β̂_D1 recovers ~0.6  (unbiased)
##     (ii) β̂_D2 centered on ~0 (no spurious activation)
##     (iii) cor(β̂_D1, β̂_D2) across reps ≈ small
##
##   If (ii) fails (β̂_D2 systematically nonzero) the two
##   directions are not separately identifiable under the
##   current effect-set.
############################################################

library(devtools)
load_all()
library(ggplot2)

## =========================================================
## 0. True parameter values
## =========================================================

n               <- 20
lambdaSelf      <- 4.75
lambdaPerceived <- 4.66

beta_true <- list(
  dens_perceived   = -1.68,
  reciprocity      =  1.98,
  self_outdeg      = -0.98,
  self_recip       =  1.50,
  ## Direction 1 — in DGP, true = 0.6
  crprod_D1        =  0.6,
  ## Direction 2 aggregate — NOT in DGP, true = 0 (identifiability probe)
  crprod_D2_agg    =  0.0
)

## =========================================================
## 1. Load real Wave 1 as fixed starting state
## =========================================================

data_path <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData/Advice"

make_css_wave <- function(df,
                          perceiver_col = "Participant_ID",
                          sender_col    = "Subject_ID",
                          receiver_cols = NULL,
                          forbid_sender_to_self = TRUE) {
  if (is.null(receiver_cols))
    receiver_cols <- grep("^[0-9]+$", names(df), value = TRUE)
  actor_ids <- receiver_cols
  n_local   <- length(actor_ids)
  i_idx <- match(as.character(df[[perceiver_col]]), actor_ids)
  j_idx <- match(as.character(df[[sender_col]]),    actor_ids)
  Xk <- as.matrix(df[, receiver_cols])
  storage.mode(Xk) <- "integer"
  css <- array(0L, dim = c(n_local, n_local, n_local))
  for (r in seq_len(nrow(df))) css[i_idx[r], j_idx[r], ] <- Xk[r, ]
  if (forbid_sender_to_self)
    for (ii in seq_len(n_local)) for (jj in seq_len(n_local)) css[ii, jj, jj] <- 0L
  css
}

f_list   <- list.files(data_path, pattern = "Advice_Week [0-9]+\\.csv", full.names = TRUE)
df_waves <- lapply(f_list, read.csv, check.names = FALSE)
css_w1   <- make_css_wave(df_waves[[1]])
for (i in 1:n) for (j in 1:n) css_w1[i, j, j] <- 0L

## =========================================================
## 2. Delta functions — ASYMMETRIC DGP (D1 only)
##
##   Y[shared] toggle x[i,j,k]  :  includes β_D1 * x[j,j,k]
##   Y[self]   toggle x[i,i,k]  :  NO D2 term  (β_D2 = 0 hard-coded)
## =========================================================

## Y[self]: x[i, i, k] toggle.  No Direction 2 in this DGP.
delta_f_self_iik <- function(x, i, k, beta) {
  sgn <- (1 - 2 * x[i, i, k])
  sgn * (beta$self_outdeg +
         beta$self_recip * x[k, k, i])
  ## NOTE: no β_D2_agg * Σ_m x[m, i, k] term here (β_D2 = 0).
}

## Y[shared] (perceived slice): x[i, j, k] toggle.  Direction 1 active.
delta_f_perceived_ijk <- function(x, i, j, k, beta) {
  sgn <- (1 - 2 * x[i, j, k])
  recip_val <- if (k == i) 0 else x[i, k, j]
  sgn * (beta$dens_perceived +
         beta$reciprocity * recip_val +
         beta$crprod_D1   * x[j, j, k])
}

softmax_choice <- function(delta) {
  m <- max(c(0, delta))
  w <- c(exp(0 - m), exp(delta - m))
  p <- w / sum(w)
  sample.int(length(p), 1, prob = p) - 1L
}

saom_step_self <- function(x, i, beta) {
  k_set  <- setdiff(1:n, i)
  delta  <- sapply(k_set, function(k) delta_f_self_iik(x, i, k, beta))
  choice <- softmax_choice(delta)
  if (choice == 0L) return(x)
  k <- k_set[choice]
  x[i, i, k] <- 1L - x[i, i, k]
  x
}

saom_step_perceived <- function(x, i, j, beta) {
  k_set  <- setdiff(1:n, j)
  delta  <- sapply(k_set,
                   function(k) delta_f_perceived_ijk(x, i, j, k, beta))
  choice <- softmax_choice(delta)
  if (choice == 0L) return(x)
  k <- k_set[choice]
  x[i, j, k] <- 1L - x[i, j, k]
  x
}

simulate_interval <- function(x0, t_end = 1) {
  x          <- x0
  t          <- 0
  total_rate <- n * lambdaSelf + n * (n - 1) * lambdaPerceived
  p_self     <- (n * lambdaSelf) / total_rate
  while (t < t_end) {
    t <- t + rexp(1, rate = total_rate)
    if (t >= t_end) break
    if (runif(1) < p_self) {
      i <- sample.int(n, 1)
      x <- saom_step_self(x, i, beta_true)
    } else {
      i <- sample.int(n, 1)
      j <- sample(setdiff(1:n, i), 1)
      x <- saom_step_perceived(x, i, j, beta_true)
    }
  }
  x
}

## =========================================================
## 3. RSiena algorithm for estimation
## =========================================================

alg_est <- sienaAlgorithmCreate(
  projname = "crprod_identifiability",
  n3       = 500,
  nsub     = 3,
  seed     = NULL
)

## =========================================================
## 4. Helper: build RSiena objects with BOTH directions
## =========================================================

build_model <- function(css_arr) {
  Y   <- sienaDependent(css_arr, type = "threeway",
                        shareParameters = TRUE, sharedCov = TRUE)
  dat <- sienaDataCreate(Y)
  eff <- getEffects(dat)

  ## Perceived-slice baseline (matches DGP)
  eff <- includeEffects(eff, recip, name = "Y[shared]")

  ## Y[self] baseline (matches DGP)
  eff <- includeEffects(eff, recip, name = "Y[self]")

  ## --- Direction 1 — perceived slice ← Y[self] (IN DGP, true = 0.6) ---
  eff <- includeEffects(eff, crprod,
                        interaction1 = "Y[self]",
                        name         = "Y[shared]")

  ## --- Direction 2 aggregate — Y[self] ← Σ_m Y[m]  (NOT IN DGP, true = 0) ---
  ## name auto-defaults to "Y[self]" via sienaeffects.r auto-default.
  eff <- includeEffects(eff, crprod,
                        interaction1 = "Y[shared]")

  list(dat = dat, eff = eff)
}

## =========================================================
## 5. Quick sanity check
## =========================================================

cat("=== Sanity check: building model with D1 + D2 ===\n")
css_test        <- array(0L, dim = c(n, n, n, 2))
css_test[,,,1]  <- css_w1
css_test[,,,2]  <- css_w1
mdl_test <- tryCatch(build_model(css_test), error = function(e) {
  cat("ERROR building model:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(mdl_test)) {
  cat("\nIncluded effects:\n")
  inc <- mdl_test$eff[mdl_test$eff$include, ]
  print.data.frame(inc[, c("name", "effectName", "shortName",
                            "interaction1", "type",
                            if ("sharedDup" %in% names(inc)) "sharedDup"
                            else NULL)])

  cat("\nDirection 1 row(s) (name ends in [shared], i1=Y[self]):\n")
  d1 <- inc[grepl("\\[shared\\]$", inc$name) &
            inc$shortName == "crprod" &
            inc$interaction1 == "Y[self]" &
            (if (is.null(inc$sharedDup)) TRUE else !inc$sharedDup), ]
  print.data.frame(d1[, c("name", "effectName", "shortName", "interaction1")])

  cat("\nDirection 2 row(s) (name ends in [self], i1=Y[shared]):\n")
  d2 <- inc[grepl("\\[self\\]$", inc$name) &
            inc$shortName == "crprod" &
            inc$interaction1 == "Y[shared]", ]
  print.data.frame(d2[, c("name", "effectName", "shortName", "interaction1")])
}

## =========================================================
## 6. Parameter recovery loop
## =========================================================

N_SIM <- 10
# N_SIM <- 1

res_list      <- vector("list", N_SIM)
cov_list      <- vector("list", N_SIM)  # estimator covariance per rep
converged     <- logical(N_SIM)

cat(sprintf("\nStarting parameter recovery: %d replications\n", N_SIM))
cat(sprintf("DGP:  β_D1 = %.2f (active),  β_D2 = %.2f (null)\n\n",
            beta_true$crprod_D1, beta_true$crprod_D2_agg))

set.seed(2027)

for (sim in seq_len(N_SIM)) {

  cat(sprintf("[%d/%d] Generating data... ", sim, N_SIM))

  css_sim           <- simulate_interval(css_w1, t_end = 1)
  css_array_sim     <- array(0L, dim = c(n, n, n, 2))
  css_array_sim[,,,1] <- css_w1
  css_array_sim[,,,2] <- css_sim

  mdl <- tryCatch(build_model(css_array_sim), error = function(e) {
    cat("SKIP (model build error:", conditionMessage(e), ")\n"); NULL
  })
  if (is.null(mdl)) next

  cat("estimating... ")
  ans_rep <- tryCatch(
    siena07(alg_est,
            data       = mdl$dat,
            effects    = mdl$eff,
            batch      = TRUE,
            silent     = TRUE,
            thetaBound = 500),
    error = function(e) { cat("\n  ERR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(ans_rep)) { cat("SKIP (siena07 error)\n"); next }

  conv_ok <- !is.null(ans_rep$tconv.max) && ans_rep$tconv.max < 0.35
  converged[sim] <- conv_ok

  ef_df <- tryCatch(as.data.frame(ans_rep$requestedEffects),
                    error = function(e) NULL)
  if (is.null(ef_df)) { cat("SKIP (effects error)\n"); next }

  dup_col <- if ("sharedDup" %in% names(ef_df)) ef_df$sharedDup
             else rep(FALSE, nrow(ef_df))

  out <- cbind(
    data.frame(rep = sim, converged = conv_ok),
    ef_df[, c("name", "effectName", "shortName", "interaction1")],
    theta     = ans_rep$theta,
    sharedDup = dup_col
  )
  res_list[[sim]] <- out

  ## Capture D1 ↔ D2 estimator covariance entry from covtheta
  ##   D1: name ends [shared], shortName crprod, i1=Y[self],   !sharedDup
  ##   D2: name ends [self],   shortName crprod, i1=Y[shared], !sharedDup
  ##       (canonical aggregate row before C++ expansion renames i1→Y[1])
  i_d1 <- which(grepl("\\[shared\\]$", out$name) &
                out$shortName == "crprod" &
                out$interaction1 == "Y[self]" &
                !out$sharedDup)
  ## After initializeFRAN expansion, D2 canonical aggregate row has i1=Y[1]
  ## AND effectName contains "(agg)". Use the (agg) marker to disambiguate.
  i_d2 <- which(grepl("\\[self\\]$", out$name) &
                out$shortName == "crprod" &
                grepl("\\(agg\\)", out$effectName) &
                !out$sharedDup)

  if (length(i_d1) == 1 && length(i_d2) == 1 &&
      !is.null(ans_rep$covtheta)) {
    v11 <- ans_rep$covtheta[i_d1, i_d1]
    v22 <- ans_rep$covtheta[i_d2, i_d2]
    v12 <- ans_rep$covtheta[i_d1, i_d2]
    cov_list[[sim]] <- data.frame(
      rep      = sim,
      theta_d1 = ans_rep$theta[i_d1],
      theta_d2 = ans_rep$theta[i_d2],
      se_d1    = sqrt(v11),
      se_d2    = sqrt(v22),
      cov_d12  = v12,
      cor_d12  = v12 / sqrt(v11 * v22)
    )
  }

  cat(sprintf("done  (tconv.max = %.3f%s)\n",
              ans_rep$tconv.max,
              if (conv_ok) "" else "  [NOT CONVERGED]"))

  if (sim == 1) {
    cat("\n=== EFFECT TABLE (rep 1) ===\n")
    cat(sprintf("  %-20s %-38s %-12s %-12s %-6s %10s\n",
                "name", "effectName", "shortName",
                "interaction1", "dup", "theta"))
    cat(strrep("-", 112), "\n")
    for (rr in seq_len(nrow(out))) {
      cat(sprintf("  %-20s %-38s %-12s %-12s %-6s %10.4f\n",
                  as.character(out$name[rr]),
                  substr(as.character(out$effectName[rr]), 1, 38),
                  as.character(out$shortName[rr]),
                  as.character(out$interaction1[rr]),
                  as.character(out$sharedDup[rr]),
                  out$theta[rr]))
    }
    cat("=== END ===\n\n")
  }
}

res_df <- do.call(rbind, res_list[!sapply(res_list, is.null)])
res_canon <- res_df[!res_df$sharedDup, ]

cov_df <- do.call(rbind, cov_list[!sapply(cov_list, is.null)])

N_reps <- length(unique(res_canon$rep))
cat(sprintf("\nReplications used: %d / %d  (converged: %d)\n\n",
            N_reps, N_SIM, sum(converged, na.rm=TRUE)))

## =========================================================
## 7. Summary table: true vs. estimated
## =========================================================

get_true <- function(nm, ef, sn, i1) {
  is_self <- grepl("\\[self\\]$", nm)
  if (grepl("^basic rate",  ef)) return(if (is_self) lambdaSelf else lambdaPerceived)
  if (grepl("^outdegree",   ef)) return(if (is_self) beta_true$self_outdeg else beta_true$dens_perceived)
  if (grepl("^reciprocity", ef)) return(if (is_self) beta_true$self_recip else beta_true$reciprocity)
  if (sn == "crprod") {
    ## D1: Y[shared] with i1=Y[self]
    if (!is_self && i1 == "Y[self]") return(beta_true$crprod_D1)
    ## D2 aggregate canonical: Y[self] + (agg) in effectName
    if (is_self && grepl("\\(agg\\)", ef)) return(beta_true$crprod_D2_agg)
  }
  NA_real_
}

param_keys <- unique(res_canon[, c("name", "effectName", "shortName", "interaction1")])

cat(sprintf("  %-55s  %8s  %8s  %8s  %8s\n", "Parameter", "True", "Mean", "SD", "Bias"))
cat(strrep("-", 91), "\n")
for (ii in seq_len(nrow(param_keys))) {
  nm  <- as.character(param_keys$name[ii])
  ef  <- as.character(param_keys$effectName[ii])
  sn  <- as.character(param_keys$shortName[ii])
  i1  <- as.character(param_keys$interaction1[ii])
  tv  <- get_true(nm, ef, sn, i1)
  vals <- res_canon$theta[res_canon$name == nm &
                          res_canon$effectName == ef &
                          res_canon$interaction1 == i1]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0 || !is.finite(tv)) next
  lbl <- paste0(nm, " / ", sn, " (", i1, ")")
  cat(sprintf("  %-55s  %8.3f  %8.3f  %8.3f  %+8.3f\n",
              lbl, tv, mean(vals), sd(vals), mean(vals) - tv))
}

## =========================================================
## 7b. Identifiability diagnostics
## =========================================================

cat("\n============================================================\n")
cat("IDENTIFIABILITY DIAGNOSTICS  (D1 vs D2 aggregate)\n")
cat("============================================================\n")

if (!is.null(cov_df) && nrow(cov_df) > 0) {
  cat(sprintf("\nEstimator correlation cor(β̂_D1, β̂_D2)  (median across %d reps)\n",
              nrow(cov_df)))
  cat(sprintf("  median: %+.3f   mean: %+.3f   range: [%+.3f, %+.3f]\n\n",
              median(cov_df$cor_d12), mean(cov_df$cor_d12),
              min(cov_df$cor_d12),    max(cov_df$cor_d12)))

  ## Cross-rep correlation between point estimates θ̂_D1 and θ̂_D2
  cross_cor <- cor(cov_df$theta_d1, cov_df$theta_d2)
  cat(sprintf("Cross-rep correlation  cor_{reps}(θ̂_D1, θ̂_D2) = %+.3f\n",
              cross_cor))
  cat("  (strong negative → weakly identified: D1 and D2 trade off.)\n\n")

  ## Summary stats per direction
  cat(sprintf("D1:  mean θ̂ = %+.3f   SD = %.3f   target = %.2f\n",
              mean(cov_df$theta_d1), sd(cov_df$theta_d1),
              beta_true$crprod_D1))
  cat(sprintf("D2:  mean θ̂ = %+.3f   SD = %.3f   target = %.2f\n",
              mean(cov_df$theta_d2), sd(cov_df$theta_d2),
              beta_true$crprod_D2_agg))

  ## t-test: is D2 centered at 0?
  t_d2 <- t.test(cov_df$theta_d2, mu = 0)
  cat(sprintf("\nH0: β_D2 = 0   t = %.2f   p = %.4f\n",
              t_d2$statistic, t_d2$p.value))
  if (t_d2$p.value < 0.01) {
    cat("  ** D2 SIGNIFICANTLY NONZERO — spurious activation; weak identification.\n")
  } else {
    cat("  ✓ D2 not distinguishable from 0 — identification OK under this DGP.\n")
  }
} else {
  cat("  (cov_df empty — no rep produced both D1 and D2 rows.)\n")
}

## =========================================================
## 8. Plotting helpers
## =========================================================

is_self    <- grepl("\\[self\\]$", res_canon$name)
is_rate    <- grepl("^basic rate",  res_canon$effectName)
is_outdeg  <- grepl("^outdegree",   res_canon$effectName)
is_recip   <- grepl("^reciprocity", res_canon$effectName)
is_d1      <- !is_self & res_canon$shortName == "crprod" &
              res_canon$interaction1 == "Y[self]"
is_d2_agg  <- is_self  & res_canon$shortName == "crprod" &
              grepl("\\(agg\\)", res_canon$effectName)

theme_recovery <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position = "bottom"
    )
}

make_histogram <- function(df_sub, true_val, title, subtitle = NULL,
                           fill_col = "#4A90D9", vline_col = "darkblue",
                           hw = 1) {
  if (is.null(subtitle))
    subtitle <- sprintf("N = %d reps  (true = %.3f)", nrow(df_sub), true_val)

  ggplot(df_sub, aes(x = theta)) +
    geom_histogram(bins = 30, fill = fill_col, color = "white") +
    geom_vline(xintercept = true_val,
               color = "red", linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = mean(df_sub$theta, na.rm = TRUE),
               color = vline_col, linewidth = 0.9) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("true = %.3f\nmean = %.3f   bias = %+.3f",
                             true_val,
                             mean(df_sub$theta, na.rm = TRUE),
                             mean(df_sub$theta, na.rm = TRUE) - true_val),
             hjust = 1.05, vjust = 1.3, size = 3.5, color = "grey30") +
    coord_cartesian(xlim = c(true_val - hw, true_val + hw)) +
    labs(title = title, subtitle = subtitle,
         x = "Estimated theta", y = "Frequency") +
    theme_recovery()
}

safe_histogram <- function(df_sub, ...) {
  if (nrow(df_sub) == 0) { cat("  (skipped: 0 rows)\n"); return(NULL) }
  make_histogram(df_sub, ...)
}

## =========================================================
## 9. Build plots
## =========================================================

## Plot 1 - Rate (self + perceived)
rate_self_df <- res_canon[is_self  & is_rate & is.finite(res_canon$theta), ]
rate_perc_df <- res_canon[!is_self & is_rate & is.finite(res_canon$theta), ]
if (nrow(rate_self_df) > 0) rate_self_df$rate_type <- "Self"
if (nrow(rate_perc_df) > 0) rate_perc_df$rate_type <- "Perceived"
rate_df <- rbind(rate_self_df, rate_perc_df)
rate_true_lines <- data.frame(
  rate_type = c("Self", "Perceived"),
  true_val  = c(lambdaSelf, lambdaPerceived)
)
p_rate <- if (nrow(rate_df) > 0) {
  ggplot(rate_df, aes(x = theta, fill = rate_type)) +
    geom_histogram(bins = 30, color = "white", alpha = 0.75, position = "identity") +
    geom_vline(data = rate_true_lines,
               aes(xintercept = true_val, color = rate_type),
               linetype = "dashed", linewidth = 1.2) +
    scale_fill_manual(values  = c("Self" = "#90CAF9", "Perceived" = "#FFCC80")) +
    scale_color_manual(values = c("Self" = "#1565C0", "Perceived" = "#E65100")) +
    labs(title = "Rate Parameters",
         subtitle = sprintf("Self: true=%.2f  |  Perceived: true=%.2f",
                            lambdaSelf, lambdaPerceived),
         x = "Estimated rate", y = "Frequency",
         fill = "Slice type", color = "True value") +
    theme_recovery()
} else NULL

## Plots 2-5: baseline
cat("Building plot: Outdegree Shared\n")
p_od_shared <- safe_histogram(
  res_canon[!is_self & is_outdeg & is.finite(res_canon$theta), ],
  beta_true$dens_perceived, "Outdegree  -  Shared")

cat("Building plot: Reciprocity Shared\n")
p_recip_shared <- safe_histogram(
  res_canon[!is_self & is_recip & is.finite(res_canon$theta), ],
  beta_true$reciprocity, "Reciprocity  -  Shared (Perceived)")

cat("Building plot: Outdegree Self\n")
p_od_self <- safe_histogram(
  res_canon[is_self & is_outdeg & is.finite(res_canon$theta), ],
  beta_true$self_outdeg, "Outdegree  -  Self",
  fill_col = "#5BAD6F", vline_col = "darkgreen")

cat("Building plot: Reciprocity Self\n")
p_recip_self <- safe_histogram(
  res_canon[is_self & is_recip & is.finite(res_canon$theta), ],
  beta_true$self_recip, "Reciprocity  -  Self",
  fill_col = "#5BAD6F", vline_col = "darkgreen")

## Plot 6 - Direction 1 recovery  (true = 0.6, active in DGP)
cat("Building plot: Direction 1  (active, true=0.6)\n")
p_d1 <- safe_histogram(
  res_canon[is_d1 & is.finite(res_canon$theta), ],
  beta_true$crprod_D1,
  "crprod  -  Direction 1  (Y[shared] ← Y[self])",
  subtitle = sprintf(
    paste0("N = %d reps  (DGP true = %.2f)\n",
           "Perceptual accuracy: Σ_k x_{ijk} * x_{jjk}"),
    sum(is_d1 & is.finite(res_canon$theta)),
    beta_true$crprod_D1),
  fill_col = "#E57373", vline_col = "#B71C1C")

## Plot 7 - Direction 2 aggregate recovery  (true = 0, nuisance)
cat("Building plot: Direction 2 aggregate  (null, true=0)\n")
p_d2 <- safe_histogram(
  res_canon[is_d2_agg & is.finite(res_canon$theta), ],
  beta_true$crprod_D2_agg,
  "crprod  -  Direction 2 AGG  (Y[self] ← Σ_m Y[m])",
  subtitle = sprintf(
    paste0("N = %d reps  (DGP true = %.2f; null probe)\n",
           "If mean ≈ 0: D1 vs D2 are identifiable under this effect-set."),
    sum(is_d2_agg & is.finite(res_canon$theta)),
    beta_true$crprod_D2_agg),
  fill_col = "#BA68C8", vline_col = "#4A148C")

## Plot 8 - Joint scatter of (θ̂_D1, θ̂_D2) across reps
p_joint <- NULL
if (!is.null(cov_df) && nrow(cov_df) >= 2) {
  cat("Building plot: joint (D1, D2) scatter\n")
  p_joint <- ggplot(cov_df, aes(x = theta_d1, y = theta_d2)) +
    geom_vline(xintercept = beta_true$crprod_D1,
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = beta_true$crprod_D2_agg,
               color = "red", linetype = "dashed") +
    geom_point(size = 2.3, alpha = 0.75, color = "#1F2B7A") +
    stat_ellipse(level = 0.95, linetype = "dotted",
                 color = "#1F2B7A") +
    labs(title = "Joint estimator:  (θ̂_D1, θ̂_D2) across replications",
         subtitle = sprintf(
           "DGP: β_D1 = %.2f, β_D2 = 0.  Cross-rep cor(θ̂_D1, θ̂_D2) = %+.3f",
           beta_true$crprod_D1,
           cor(cov_df$theta_d1, cov_df$theta_d2)),
         x = "θ̂_D1  (Y[shared] ← Y[self])",
         y = "θ̂_D2  (Y[self] ← Σ_m Y[m])") +
    theme_recovery()
}

## Plot 9 - Within-rep estimator correlation
p_within_cor <- NULL
if (!is.null(cov_df) && nrow(cov_df) >= 2) {
  cat("Building plot: within-rep estimator correlation\n")
  p_within_cor <- ggplot(cov_df, aes(x = cor_d12)) +
    geom_histogram(bins = 20, fill = "#78909C", color = "white") +
    geom_vline(xintercept = 0,
               color = "red", linetype = "dashed") +
    geom_vline(xintercept = mean(cov_df$cor_d12),
               color = "darkblue") +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(title = "Within-rep correlation  cor(β̂_D1, β̂_D2)",
         subtitle = sprintf(
           paste0("N = %d reps  |  median = %+.3f\n",
                  "Values near ±1 → D1 and D2 are linearly dependent",
                  " in the local likelihood."),
           nrow(cov_df), median(cov_df$cor_d12)),
         x = "cor(β̂_D1, β̂_D2)  from ans$covtheta",
         y = "Frequency") +
    theme_recovery()
}

## =========================================================
## 10. Save plots
## =========================================================

save_dir <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData"
out_dir  <- file.path(save_dir, "hist_crprod_identifiability")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list(
  "01_rate"               = p_rate,
  "02_outdeg_shared"      = p_od_shared,
  "03_recip_shared"       = p_recip_shared,
  "04_outdeg_self"        = p_od_self,
  "05_recip_self"         = p_recip_self,
  "06_d1_active"          = p_d1,
  "07_d2_null"            = p_d2,
  "08_joint_d1_d2"        = p_joint,
  "09_within_rep_corr"    = p_within_cor
)

plots <- plots[!sapply(plots, is.null)]

for (nm in names(plots)) {
  path_png <- file.path(out_dir, paste0(nm, ".png"))
  ggsave(path_png, plots[[nm]], width = 7, height = 5, dpi = 200)
  cat("Saved:", path_png, "\n")
}

save(res_df, res_canon, cov_df, beta_true, lambdaSelf, lambdaPerceived,
     file = file.path(save_dir, "sim_crprod_identifiability.RData"))

for (p in plots) print(p)

cat("\nDone. Results saved to:\n  ", out_dir, "\n")
