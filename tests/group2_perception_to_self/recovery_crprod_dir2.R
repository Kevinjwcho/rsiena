############################################################
## simulation_crprod_direction2_recovery.R
##
## Purpose: Parameter recovery simulation for the Direction 2
##   cross-network effect — perceived slice influence on
##   the self-reported network Y[self].
##
##   Direction 1 (existing, see simulation_crprod_recovery.R):
##     Lives in Y[shared]'s (perceived slice) objective.
##       s_j(Y[i]) = sum_k Y[i]_{jk} * Y[self]_{jk}
##                 = sum_k x_{ijk} * x_{jjk}
##     Interpretation: perceptual accuracy.
##
##   Direction 2 (this file):
##     Lives in Y[self]'s objective. Two variants tested:
##
##     (C) aggregate  (interaction1 = "Y[shared]"):
##       s_j(Y[self]) = sum_{m=1..K} sum_k Y[self]_{jk} * Y[m]_{jk}
##                    = sum_m sum_k x_{jjk} * x_{mjk}
##       All K perceivers' views drive Y[self]; one β shared by K
##       expanded rows (initializeFRAN.r does the expansion).
##       Interpretation: perception shaping reality, aggregate.
##
##     (A) single-slice (interaction1 = "Y[1]"):
##       s_j(Y[self]) = sum_k Y[self]_{jk} * Y[1]_{jk}
##                    = sum_k x_{jjk} * x_{1jk}
##       Only perceiver 1's view drives Y[self]. One row, one β.
##       (Included as a sanity check alongside the aggregate.)
##
##   For each replication:
##     1. Generate synthetic CSS Wave 2 starting from real Wave 1,
##        using a DGP that includes the Direction 2 aggregate.
##     2. Estimate with siena07() including both (A) and (C).
##   After N_SIM replications, plot histograms of estimated
##   parameters overlaid with true values.
##
##   Key DGP derivation for Direction 2 aggregate:
##     When actor i toggles self-report x[i, i, k], the change
##     in the Y[self] aggregate statistic summed over all m is
##
##       Δ = sgn * sum_{m=1..n} x[m, i, k]
##
##     i.e. "how many perceivers currently believe i -> k".
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
  dens_perceived    = -1.68,
  reciprocity       =  1.98,
  self_outdeg       = -0.98,
  self_recip        =  1.50,
  ## Direction 2 aggregate crprod (Y[self] <- Σ_m Y[m]):
  ##   influence of the *aggregate* perceiver view on self-reports.
  crprod_self_perc  =  0.3,
  ## Direction 2 single-slice crprod (Y[self] <- Y[1] only):
  ##   NOT part of the DGP.  We include it as an estimated nuisance
  ##   to verify that when the DGP is aggregate-only, the single-slice
  ##   estimator is not systematically biased (should recover ~0).
  crprod_self_perc1 =  0.0
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
## 2. Delta functions & SAOM steps
## =========================================================

## ---- Y[self] objective: x[i, i, k] toggle ----
## Y[self]'s Direction 2 aggregate stat:
##   s_i(Y[self]) = sum_{m=1..n} sum_{k'} x[i, i, k'] * x[m, i, k']
## Derivative w.r.t. x[i, i, k]:
##   d/dx[i,i,k] = sum_m x[m, i, k]
## i.e. the cross-perceiver count of "does m see i -> k".
delta_f_self_iik <- function(x, i, k, beta) {
  sgn <- (1 - 2 * x[i, i, k])
  agg_perc <- sum(x[, i, k])   # Σ_m x[m, i, k]
  sgn * (beta$self_outdeg +
         beta$self_recip       * x[k, k, i] +
         beta$crprod_self_perc * agg_perc)
}

## ---- Y[shared] objective: x[i, j, k] toggle ----
## Purely perceived dynamics (reciprocity only, no Direction 1 crprod
## in this DGP so we isolate Direction 2 recovery).
## Structural zero guard for k == i (consistent with Direction 1 sim).
delta_f_perceived_ijk <- function(x, i, j, k, beta) {
  sgn <- (1 - 2 * x[i, j, k])
  recip_val <- if (k == i) 0 else x[i, k, j]
  sgn * (beta$dens_perceived +
         beta$reciprocity * recip_val)
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
  projname = "crprod_direction2_recovery",
  n3       = 500,
  nsub     = 3,
  seed     = NULL
)

## =========================================================
## 4. Helper: build RSiena objects for one css_array
## =========================================================

build_model <- function(css_arr) {
  Y   <- sienaDependent(css_arr, type = "threeway",
                        shareParameters = TRUE, sharedCov = TRUE)
  dat <- sienaDataCreate(Y)
  eff <- getEffects(dat)

  ## Perceived-slice baseline (matches the DGP)
  eff <- includeEffects(eff, recip, name = "Y[shared]")

  ## Y[self] baseline (matches the DGP)
  eff <- includeEffects(eff, recip, name = "Y[self]")

  ## Direction 2 — AGGREGATE (Y[self] <- Σ_m Y[m])
  ## name auto-defaults to Y[self] from interaction1="Y[shared]".
  eff <- includeEffects(eff, crprod,
                        interaction1 = "Y[shared]")

  ## Direction 2 — SINGLE-SLICE (Y[self] <- Y[1] only)
  ## Added as nuisance; true value = 0 under aggregate-only DGP.
  ## name auto-defaults to Y[self] from interaction1="Y[1]".
  eff <- includeEffects(eff, crprod,
                        interaction1 = "Y[1]", name = "Y[self]")

  list(dat = dat, eff = eff)
}

## =========================================================
## 5. Quick sanity check
## =========================================================

cat("=== Sanity check: building model with Direction 2 crprod ===\n")
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

  cat("\nDirection 2 rows (name ends in [self], shortName starts with crprod):\n")
  d2 <- mdl_test$eff[grepl("\\[self\\]$", mdl_test$eff$name) &
                     grepl("^crprod", mdl_test$eff$shortName) &
                     mdl_test$eff$include, ]
  if (nrow(d2) > 0) {
    print.data.frame(d2[, c("name", "effectName", "shortName",
                             "interaction1",
                             if ("sharedDup" %in% names(d2)) "sharedDup"
                             else NULL)])
  } else {
    cat("  *** No Direction 2 crprod effects found! ***\n")
  }
}

## =========================================================
## 6. Parameter recovery loop
## =========================================================

N_SIM <- 30
# N_SIM <- 1

res_list  <- vector("list", N_SIM)
converged <- logical(N_SIM)

cat(sprintf("\nStarting parameter recovery: %d replications\n\n", N_SIM))

set.seed(2026)

for (sim in seq_len(N_SIM)) {

  cat(sprintf("[%d/%d] Generating data... ", sim, N_SIM))

  ## Step 1: generate synthetic Wave 2 under Direction 2 DGP
  css_sim           <- simulate_interval(css_w1, t_end = 1)
  css_array_sim     <- array(0L, dim = c(n, n, n, 2))
  css_array_sim[,,,1] <- css_w1
  css_array_sim[,,,2] <- css_sim

  ## Step 2: build RSiena model (same spec as DGP + single-slice nuisance)
  mdl <- tryCatch(build_model(css_array_sim), error = function(e) {
    cat("SKIP (model build error:", conditionMessage(e), ")\n"); NULL
  })
  if (is.null(mdl)) next

  ## Step 3: estimate
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

  ## Step 4: store as long-format rows
  ef_df <- tryCatch(as.data.frame(ans_rep$requestedEffects),
                    error = function(e) NULL)
  if (is.null(ef_df)) { cat("SKIP (effects error)\n"); next }

  ## The requestedEffects table after initializeFRAN has K expanded rows
  ## for the aggregate (sharedDup=TRUE on K-1 of them). We keep the full
  ## table here; downstream we filter to non-dup rows for plotting.
  dup_col <- if ("sharedDup" %in% names(ef_df)) ef_df$sharedDup
             else rep(FALSE, nrow(ef_df))

  out <- cbind(
    data.frame(rep = sim, converged = conv_ok),
    ef_df[, c("name", "effectName", "shortName", "interaction1")],
    theta     = ans_rep$theta,
    sharedDup = dup_col
  )
  res_list[[sim]] <- out

  cat(sprintf("done  (tconv.max = %.3f%s)\n",
              ans_rep$tconv.max,
              if (conv_ok) "" else "  [NOT CONVERGED]"))

  ## Print first replication's effect table for verification
  if (sim == 1) {
    cat("\n=== EFFECT TABLE (rep 1) ===\n")
    cat(sprintf("  %-20s %-35s %-15s %-12s %-6s %10s\n",
                "name", "effectName", "shortName",
                "interaction1", "dup", "theta"))
    cat(strrep("-", 108), "\n")
    for (rr in seq_len(nrow(out))) {
      cat(sprintf("  %-20s %-35s %-15s %-12s %-6s %10.4f\n",
                  as.character(out$name[rr]),
                  substr(as.character(out$effectName[rr]), 1, 35),
                  as.character(out$shortName[rr]),
                  as.character(out$interaction1[rr]),
                  as.character(out$sharedDup[rr]),
                  out$theta[rr]))
    }
    cat("=== END ===\n\n")
  }
}

res_df      <- do.call(rbind, res_list[!sapply(res_list, is.null)])
res_df_conv <- res_df
# res_df_conv <- res_df[res_df$converged, ]   # uncomment to filter converged only

N_reps <- length(unique(res_df_conv$rep))
cat(sprintf("\nReplications used for plots: %d / %d\n\n", N_reps, N_SIM))

## Diagnostic: show what names and effects are in the result
cat("=== Diagnostic: unique (name, effectName, shortName, i1, dup) in results ===\n")
if (nrow(res_df_conv) > 0) {
  diag_df <- unique(res_df_conv[, c("name", "effectName", "shortName",
                                     "interaction1", "sharedDup")])
  print.data.frame(diag_df)
} else {
  cat("  (no results)\n")
}

## =========================================================
## 7. Summary table: true vs. estimated
##
## Row-classification logic:
##   - "rate"       → basic rate effect
##   - "outdegree"  → density
##   - "reciprocity"→ reciprocity
##   - "crprod, agg, Y[self]"         → Direction 2 aggregate (canonical row only)
##   - "crprod, single Y[1], Y[self]" → Direction 2 single-slice
## =========================================================

## Keep only canonical rows (drop sharedDup=TRUE duplicates) for summary/plots.
res_canon <- res_df_conv[!res_df_conv$sharedDup, ]

get_true <- function(nm, ef, sn, i1) {
  is_self <- grepl("\\[self\\]$", nm)
  if (grepl("^basic rate",  ef)) return(if (is_self) lambdaSelf else lambdaPerceived)
  if (grepl("^outdegree",   ef)) return(if (is_self) beta_true$self_outdeg else beta_true$dens_perceived)
  if (grepl("^reciprocity", ef)) return(if (is_self) beta_true$self_recip else beta_true$reciprocity)
  ## Direction 2: name=Y[self], shortName=crprod, i1 distinguishes variant
  if (sn == "crprod" && is_self) {
    ## Aggregate canonical: effectName contains "(agg)" OR (post-expansion)
    ## i1 == "Y[1]" with sharedDup==FALSE comes from aggregate canonical
    ## (but single-slice also uses i1="Y[1]"). We distinguish by effectName:
    if (grepl("\\(agg\\)", ef)) return(beta_true$crprod_self_perc)
    ## single-slice row: effectName has no "(agg)" tag
    return(beta_true$crprod_self_perc1)
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
## 8. Plotting helpers
## =========================================================

is_self    <- grepl("\\[self\\]$", res_canon$name)
is_rate    <- grepl("^basic rate",  res_canon$effectName)
is_outdeg  <- grepl("^outdegree",   res_canon$effectName)
is_recip   <- grepl("^reciprocity", res_canon$effectName)

## Direction 2 aggregate (canonical row): Y[self] + shortName=crprod + effectName has (agg)
is_crprod_d2_agg <- is_self &
                    res_canon$shortName == "crprod" &
                    grepl("\\(agg\\)", res_canon$effectName)

## Direction 2 single-slice: Y[self] + shortName=crprod + effectName has NO (agg)
is_crprod_d2_one <- is_self &
                    res_canon$shortName == "crprod" &
                    !grepl("\\(agg\\)", res_canon$effectName)

theme_recovery <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 9),
      legend.position = "bottom"
    )
}

mean_bias_label <- function(vals, true_val) {
  sprintf("mean = %.3f  |  bias = %+.3f",
          mean(vals, na.rm = TRUE),
          mean(vals, na.rm = TRUE) - true_val)
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
             label = sprintf("true = %.3f\n%s", true_val,
                             mean_bias_label(df_sub$theta, true_val)),
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
## 9. Build all plots
## =========================================================

## Plot 1 - Rate (self + perceived overlaid)
rate_self_df <- res_canon[is_self  & is_rate & is.finite(res_canon$theta), ]
rate_perc_df <- res_canon[!is_self & is_rate & is.finite(res_canon$theta), ]
if (nrow(rate_self_df) > 0) rate_self_df$rate_type <- "Self"
if (nrow(rate_perc_df) > 0) rate_perc_df$rate_type <- "Perceived"
rate_df <- rbind(rate_self_df, rate_perc_df)

rate_true_lines <- data.frame(
  rate_type = c("Self", "Perceived"),
  true_val  = c(lambdaSelf, lambdaPerceived)
)

if (nrow(rate_df) > 0) {
  p_rate <- ggplot(rate_df, aes(x = theta, fill = rate_type)) +
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
} else {
  cat("WARNING: no rate parameters found in results.\n")
  p_rate <- NULL
}

## Plot 2 - Outdegree, Shared
cat("Building plot: Outdegree Shared\n")
od_shared_df <- res_canon[!is_self & is_outdeg & is.finite(res_canon$theta), ]
p_od_shared  <- safe_histogram(od_shared_df, beta_true$dens_perceived,
                               "Outdegree  -  Shared")

## Plot 3 - Reciprocity, Shared (perceived)
cat("Building plot: Reciprocity Shared\n")
recip_df <- res_canon[!is_self & is_recip & is.finite(res_canon$theta), ]
p_recip  <- safe_histogram(recip_df, beta_true$reciprocity,
                           "Reciprocity  -  Shared (Perceived)")

## Plot 4 - Outdegree, Self
cat("Building plot: Outdegree Self\n")
od_self_df <- res_canon[is_self & is_outdeg & is.finite(res_canon$theta), ]
p_od_self  <- safe_histogram(od_self_df, beta_true$self_outdeg,
                             "Outdegree  -  Self",
                             fill_col = "#5BAD6F", vline_col = "darkgreen")

## Plot 5 - Reciprocity, Self
cat("Building plot: Reciprocity Self\n")
recip_self_df <- res_canon[is_self & is_recip & is.finite(res_canon$theta), ]
p_recip_self  <- safe_histogram(recip_self_df, beta_true$self_recip,
                                "Reciprocity  -  Self",
                                fill_col = "#5BAD6F", vline_col = "darkgreen")

## Plot 6 - Direction 2 AGGREGATE crprod  *** THE KEY TEST ***
cat("Building plot: crprod Direction 2 aggregate\n")
crprod_agg_df <- res_canon[is_crprod_d2_agg & is.finite(res_canon$theta), ]
p_crprod_agg  <- safe_histogram(crprod_agg_df, beta_true$crprod_self_perc,
                                "crprod  -  Direction 2 AGGREGATE  (Y[self] <- Σ_m Y[m])",
                                subtitle = sprintf(
                                  paste0("N = %d reps  (true = %.2f)\n",
                                         "Σ_m Σ_k x_{jjk} * x_{mjk}: does aggregate ",
                                         "perceiver view shape actor j's self-report?"),
                                  nrow(crprod_agg_df),
                                  beta_true$crprod_self_perc),
                                fill_col  = "#E57373",
                                vline_col = "#B71C1C")

## Plot 7 - Direction 2 SINGLE-SLICE crprod  (nuisance; true = 0)
cat("Building plot: crprod Direction 2 single-slice\n")
crprod_one_df <- res_canon[is_crprod_d2_one & is.finite(res_canon$theta), ]
p_crprod_one  <- safe_histogram(crprod_one_df, beta_true$crprod_self_perc1,
                                "crprod  -  Direction 2 SINGLE-SLICE  (Y[self] <- Y[1])",
                                subtitle = sprintf(
                                  paste0("N = %d reps  (true = %.2f; aggregate-only DGP)\n",
                                         "Nuisance term: should not be systematically biased."),
                                  nrow(crprod_one_df),
                                  beta_true$crprod_self_perc1),
                                fill_col  = "#BA68C8",
                                vline_col = "#4A148C")

## =========================================================
## 10. Save plots
## =========================================================

save_dir <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData"
out_dir  <- file.path(save_dir, "hist_crprod_direction2_recovery")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list(
  "01_rate"             = p_rate,
  "02_outdeg_shared"    = p_od_shared,
  "03_recip_shared"     = p_recip,
  "04_outdeg_self"      = p_od_self,
  "05_recip_self"       = p_recip_self,
  "06_crprod_d2_agg"    = p_crprod_agg,
  "07_crprod_d2_single" = p_crprod_one
)

plots <- plots[!sapply(plots, is.null)]

for (nm in names(plots)) {
  path_png <- file.path(out_dir, paste0(nm, ".png"))
  ggsave(path_png, plots[[nm]], width = 7, height = 5, dpi = 200)
  cat("Saved:", path_png, "\n")
}

## Save raw results for post-hoc inspection
save(res_df, res_canon, beta_true, lambdaSelf, lambdaPerceived,
     file = file.path(save_dir, "sim_crprod_direction2_recovery.RData"))

## Print to screen
for (p in plots) print(p)

cat("\nDone. Results saved to:\n  ", out_dir, "\n")
