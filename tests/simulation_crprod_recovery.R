############################################################
## simulation_crprod_recovery.R
##
## Purpose: Parameter recovery simulation for cross-network
##   effects between perceived slices and Y[self], using
##   the existing RSiena crprod infrastructure.
##
##   The key statistic is:
##     crprod (Y[shared] with Y[self]):
##       s_j(x) = sum_k Y[i]_{jk} * Y[self]_{jk}
##              = sum_k x_{ijk} * x_{jjk}
##   This measures perceptual accuracy: how much does
##   perceiver i's view of j's ties match j's actual ties.
##
##   For each replication:
##     1. Generate synthetic CSS data using known true
##        parameters (including crprod for perc-self cross)
##     2. Estimate the model with siena07()
##   After N_SIM replications, plot histograms of
##   estimated parameters overlaid with true values.
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
  ## Cross-network effect: perceptual accuracy
  ## crprod Y[shared] with Y[self]: sum_k x_{ijk} * x_{jjk}
  crprod_perc_self =  0.6
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

## Self-reported tie change: x[i, i, k]
## Reciprocity in Y[self]: does k self-report a tie to i?  x[k, k, i]
delta_f_self_iik <- function(x, i, k, beta) {
  sgn <- (1 - 2 * x[i, i, k])
  sgn * (beta$self_outdeg +
           beta$self_recip * x[k, k, i])
}

## Perceived tie change: x[i, j, k]
## Includes crprod effect: when toggling x[i,j,k],
## the change in sum_k' x[i,j,k'] * x[j,j,k'] is:
##   delta = sgn * x[j, j, k]
## (because only the k-th term changes, and sgn = +1/-1)
delta_f_perceived_ijk <- function(x, i, j, k, beta) {
  sgn <- (1 - 2 * x[i, j, k])
  ## When k == i, the reciprocity term x[i, k, j] = x[i, i, j] lives
  ## in row i of Y[i], which is structural zero in RSiena.  The C++
  ## engine always sees 0 for that cell, so the DGP must agree;
  ## otherwise the data are generated under a different model than what
  ## RSiena estimates, causing systematic bias.
  recip_val <- if (k == i) 0 else x[i, k, j]
  sgn * (beta$dens_perceived +
           beta$reciprocity      * recip_val +
           beta$crprod_perc_self * x[j, j, k])
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
  projname = "crprod_recovery",
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

  ## Include crprod: Y[shared] with Y[self]
  ## This is the cross-network effect we're testing
  eff <- includeEffects(eff, crprod,
                        interaction1 = "Y[self]",
                        name = "Y[shared]")

  ## Self-reported reciprocity (in DGP)
  eff <- includeEffects(eff, recip, name = "Y[self]")

  list(dat = dat, eff = eff)
}

## =========================================================
## 5. Quick sanity check: does crprod appear?
## =========================================================

cat("=== Sanity check: building model with crprod ===\n")
css_test        <- array(0L, dim = c(n, n, n, 2))
css_test[,,,1]  <- css_w1
css_test[,,,2]  <- css_w1  # same wave, just for structure check
mdl_test <- tryCatch(build_model(css_test), error = function(e) {
  cat("ERROR building model:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(mdl_test)) {
  cat("\nEffects table (included only):\n")
  inc <- mdl_test$eff[mdl_test$eff$include, ]
  print.data.frame(inc[, c("name", "effectName", "shortName", "interaction1", "type")])

  cat("\nAll crprod-family effects available:\n")
  crp <- mdl_test$eff[grep("crprod", mdl_test$eff$shortName), ]
  if (nrow(crp) > 0) {
    print.data.frame(crp[, c("name", "effectName", "shortName", "interaction1", "include", "sharedDup")])
  } else {
    cat("  *** No crprod effects found! Check threeWayNet() changes. ***\n")
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

set.seed(2025)

for (sim in seq_len(N_SIM)) {

  cat(sprintf("[%d/%d] Generating data... ", sim, N_SIM))

  ## Step 1: generate synthetic Wave 2
  css_sim           <- simulate_interval(css_w1, t_end = 1)
  css_array_sim     <- array(0L, dim = c(n, n, n, 2))
  css_array_sim[,,,1] <- css_w1
  css_array_sim[,,,2] <- css_sim

  ## Step 2: build RSiena model
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
    error = function(e) NULL
  )
  if (is.null(ans_rep)) { cat("SKIP (siena07 error)\n"); next }

  conv_ok <- !is.null(ans_rep$tconv.max) && ans_rep$tconv.max < 0.35
  converged[sim] <- conv_ok

  ## Step 4: store as long-format rows
  ef_df <- tryCatch(as.data.frame(ans_rep$requestedEffects),
                    error = function(e) NULL)
  if (is.null(ef_df)) { cat("SKIP (effects error)\n"); next }

  out <- cbind(
    data.frame(rep = sim, converged = conv_ok),
    ef_df[, c("name", "effectName", "shortName", "interaction1")],
    theta = ans_rep$theta
  )
  res_list[[sim]] <- out

  cat(sprintf("done  (tconv.max = %.3f%s)\n",
              ans_rep$tconv.max,
              if (conv_ok) "" else "  [NOT CONVERGED]"))

  ## Print first replication's effect table for verification
  if (sim == 1) {
    cat("\n=== EFFECT TABLE (rep 1) ===\n")
    cat(sprintf("  %-20s %-30s %-12s %-15s %10s\n",
                "name", "effectName", "shortName", "interaction1", "theta"))
    cat(strrep("-", 92), "\n")
    for (rr in seq_len(nrow(out))) {
      cat(sprintf("  %-20s %-30s %-12s %-15s %10.4f\n",
                  as.character(out$name[rr]),
                  as.character(out$effectName[rr]),
                  as.character(out$shortName[rr]),
                  as.character(out$interaction1[rr]),
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
cat("=== Diagnostic: unique (name, effectName, shortName) in results ===\n")
if (nrow(res_df_conv) > 0) {
  diag_df <- unique(res_df_conv[, c("name", "effectName", "shortName", "interaction1")])
  print.data.frame(diag_df)
} else {
  cat("  (no results)\n")
}

## =========================================================
## 7. Summary table: true vs. estimated
## =========================================================

get_true <- function(nm, ef, sn, i1) {
  is_self <- grepl("\\[self\\]", nm)
  if (grepl("^basic rate",  ef)) return(if (is_self) lambdaSelf else lambdaPerceived)
  if (grepl("^outdegree",   ef)) return(if (is_self) beta_true$self_outdeg else beta_true$dens_perceived)
  if (grepl("^reciprocity", ef)) return(if (is_self) beta_true$self_recip else beta_true$reciprocity)
  if (sn == "crprod" && grepl("self", i1)) return(beta_true$crprod_perc_self)
  NA_real_
}

param_keys <- unique(res_df_conv[, c("name", "effectName", "shortName", "interaction1")])

cat(sprintf("  %-50s  %8s  %8s  %8s  %8s\n", "Parameter", "True", "Mean", "SD", "Bias"))
cat(strrep("-", 86), "\n")
for (ii in seq_len(nrow(param_keys))) {
  nm  <- as.character(param_keys$name[ii])
  ef  <- as.character(param_keys$effectName[ii])
  sn  <- as.character(param_keys$shortName[ii])
  i1  <- as.character(param_keys$interaction1[ii])
  tv  <- get_true(nm, ef, sn, i1)
  vals <- res_df_conv$theta[res_df_conv$name == nm & res_df_conv$effectName == ef]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0 || !is.finite(tv)) next
  cat(sprintf("  %-50s  %8.3f  %8.3f  %8.3f  %+8.3f\n",
              paste0(nm, " / ", sn), tv,
              mean(vals), sd(vals), mean(vals) - tv))
}

## =========================================================
## 8. Plotting helpers
## =========================================================

is_self    <- grepl("\\[self\\]",  res_df_conv$name)
is_rate    <- grepl("^basic rate",  res_df_conv$effectName)
is_outdeg  <- grepl("^outdegree",   res_df_conv$effectName)
is_recip   <- grepl("^reciprocity", res_df_conv$effectName)
is_crprod  <- (res_df_conv$shortName == "crprod") & grepl("self", res_df_conv$interaction1)

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

## =========================================================
## 9. Build all plots
## =========================================================

## Plot 1 - Rate (self + perceived overlaid)
rate_self_df <- res_df_conv[is_self  & is_rate & is.finite(res_df_conv$theta), ]
rate_perc_df <- res_df_conv[!is_self & is_rate & is.finite(res_df_conv$theta), ]
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
  cat("WARNING: no rate parameters found in results. Check name patterns.\n")
  p_rate <- NULL
}

## Helper: safe histogram (returns NULL if no data)
safe_histogram <- function(df_sub, ...) {
  if (nrow(df_sub) == 0) { cat("  (skipped: 0 rows)\n"); return(NULL) }
  make_histogram(df_sub, ...)
}

## Plot 2 - Outdegree, Shared
cat("Building plot: Outdegree Shared\n")
od_shared_df <- res_df_conv[!is_self & is_outdeg & is.finite(res_df_conv$theta), ]
p_od_shared  <- safe_histogram(od_shared_df, beta_true$dens_perceived,
                               "Outdegree  -  Shared")

## Plot 3 - Reciprocity, Shared (perceived)
cat("Building plot: Reciprocity Shared\n")
recip_df <- res_df_conv[!is_self & is_recip & is.finite(res_df_conv$theta), ]
p_recip  <- safe_histogram(recip_df, beta_true$reciprocity,
                           "Reciprocity  -  Shared (Perceived)")

## Plot 4 - Outdegree, Self
cat("Building plot: Outdegree Self\n")
od_self_df <- res_df_conv[is_self & is_outdeg & is.finite(res_df_conv$theta), ]
p_od_self  <- safe_histogram(od_self_df, beta_true$self_outdeg,
                             "Outdegree  -  Self",
                             fill_col = "#5BAD6F", vline_col = "darkgreen")

## Plot 5 - Reciprocity, Self
cat("Building plot: Reciprocity Self\n")
recip_self_df <- res_df_conv[is_self & is_recip & is.finite(res_df_conv$theta), ]
p_recip_self  <- safe_histogram(recip_self_df, beta_true$self_recip,
                                "Reciprocity  -  Self",
                                fill_col = "#5BAD6F", vline_col = "darkgreen")

## Plot 6 - crprod (Y[shared] with Y[self])  *** THE KEY TEST ***
cat("Building plot: crprod\n")
crprod_df <- res_df_conv[is_crprod & is.finite(res_df_conv$theta), ]
p_crprod  <- safe_histogram(crprod_df, beta_true$crprod_perc_self,
                            "crprod  -  Perceptual Accuracy (Y[shared] x Y[self])",
                            subtitle = sprintf(
                              "N = %d reps  (true = %.2f)\nsum_k x_{ijk} * x_{jjk}: does i's perception of j match j's reality?",
                              nrow(crprod_df), beta_true$crprod_perc_self),
                            fill_col  = "#E57373",
                            vline_col = "#B71C1C")

## =========================================================
## 10. Save plots
## =========================================================

save_dir <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData"
out_dir  <- file.path(save_dir, "hist_crprod_recovery")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list(
  "01_rate"           = p_rate,
  "02_outdeg_shared"  = p_od_shared,
  "03_recip_shared"   = p_recip,
  "04_outdeg_self"    = p_od_self,
  "05_recip_self"     = p_recip_self,
  "06_crprod"         = p_crprod
)

## Remove NULL entries (skipped plots)
plots <- plots[!sapply(plots, is.null)]

for (nm in names(plots)) {
  path_png <- file.path(out_dir, paste0(nm, ".png"))
  ggsave(path_png, plots[[nm]], width = 7, height = 5, dpi = 200)
  cat("Saved:", path_png, "\n")
}

## Print to screen
for (p in plots) print(p)
