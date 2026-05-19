###############################################################################
##
##  recovery_percRecip.R
##
##  Recovery study for the new `percRecip` effect (perceived reciprocity)
##  using the fork's three-way SAOM mode.
##
##  Pilot setup:
##    - n = 10 actors  (so K = 10 perceiver slices, plus Y[self])
##    - 2 waves
##    - 50 replications
##    - shareParameters = TRUE  (one shared beta per structural effect across
##                               perceiver slices Y[1]..Y[K])
##
##  Effects in the model:
##    * Shared across perceiver slices (single beta each):
##        density, reciprocity, percRecip
##    * Self-reported network Y[self]:
##        density, reciprocity
##    * No cross-network (crprod) effects in the pilot.
##    * Rate parameters: per-network defaults; treated as nuisance.
##
##  PREREQ: the dev RSiena fork at /Users/jinwoocho/Desktop/rsiena must
##          already include the percRecip C++ effect.  This script loads
##          it via devtools::load_all() so no install step is needed.
##
###############################################################################

## ---- 0.  Load dev RSiena (not the CRAN install) -------------------------

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

cat("Loaded RSiena from:", find.package("RSiena"), "\n")

set.seed(20260430)

## ---- 1.  Hyper-parameters ----------------------------------------------

## Minimal pilot setup, smallest possible:
##   * n=10 actors (K=10 perceiver slices + Y[self])
##   * Only structural effects kept: density on Y[shared] and Y[self],
##     plus the percRecip effect being tested.
##   * recip (Y[shared] and Y[self]) DROPPED to cut moments
##     from 6 to 3.
##   * Rates retained (one per network) but fixed at truth in
##     estimation so they don't contribute free parameters.
##   * In estimation everything except percRecip is fixed at truth;
##     percRecip starts from default initialValue=0 and must reach 0.8.
## n_actors bumped to 20 to give percRecip enough signal: at n=10 the
## expected statistic Σ_i Σ_j Y[self]_{i,j} Y[i]_{j,i} ≈ 0.9 — too small
## a difference between β=±ε for MoM to grab.  At n=20 the expectation
## ≈ 3.8 (4× larger), comparable to typical SAOM moment magnitudes.
n_actors    <- 20          # K = n in this three-way mode
n_reps      <- 1
density_p1  <- 0.10        # initial-wave Bernoulli tie probability

## Truth (parameters we will try to recover).  Only percRecip is free
## in estimation; the rest are pinned at truth via fix=TRUE.
truth <- c(
  rate_perc      =  6.0,   # rate of Y[1]..Y[K] (per slice, per period)
  rate_self      =  6.0,   # rate of Y[self]
  density_perc   = -2.0,   # shared across Y[1]..Y[K]
  percRecip      =  0.8,   # shared (the parameter of interest)
  density_self   = -2.0    # Y[self]
)

## ---- 2.  Build wave-1 four-way array -----------------------------------

make_wave1_4way <- function(n, p) {
  arr <- array(rbinom(n^3, 1, p), c(n, n, n))   # K=n perceivers
  ## j -> j self-loops not allowed: zero out diagonal positions in each slice
  for (k in seq_len(n)) {
    diag(arr[k, , ]) <- 0L
  }
  arr
}

wave1 <- make_wave1_4way(n_actors, density_p1)

## Build the 4D (n, n, n, T=2) array used by sienaDependent(type="threeway").
## Wave 2 is a placeholder (= wave 1) — we overwrite it after simulation.
build_4d <- function(arr_t1, arr_t2) {
  T_ <- 2
  out <- array(0L, dim = c(dim(arr_t1), T_))
  out[ , , , 1] <- arr_t1
  out[ , , , 2] <- arr_t2
  out
}

base_arr <- build_4d(wave1, wave1)

## ---- 3.  Build base sienaData with shared parameters --------------------

## sharedCov = TRUE is CRITICAL: it triggers the covariate-like share
## marking in effects.r:1029-1032 for cross-network rows (interaction1
## != "").  Without it, cross-network effects fall through to a
## fallback share branch (effects.r:1037-1043) that breaks the
## simulator's evaluation function wiring — β values planted into
## theta never reach the choice probability computation.  Confirmed by
## tests/test_crprod_effects.R, which uses sharedCov=TRUE and recovers
## crprod-family parameters successfully.
Y         <- sienaDependent(base_arr, type = "threeway",
                            shareParameters = TRUE, sharedCov = TRUE)
base_data <- sienaDataCreate(Y)

## ---- 4.  Build effects + truth (for simulation only) -------------------

sim_eff <- getEffects(base_data)

## Prune the auto-added cross-network (crprod-family) effects so we keep
## only the structural rows we explicitly want.  Then add percRecip back
## via includeEffects below with its share=Y[shared] tag.
.cross_idx <- which(sim_eff$include & sim_eff$interaction1 != "")
sim_eff[.cross_idx, "include"] <- FALSE
rm(.cross_idx)

## With shareParameters = TRUE, the canonical row is renamed to "Y[shared]"
## (see R/effects.r line 1007).  Y[2]..Y[K] still exist but are flagged
## sharedDup = TRUE.  So we add and set effects on "Y[shared]" / "Y[self]".

## Minimal effect set:
##   * density (auto-included by default) on Y[shared] and Y[self]
##   * percRecip as a Direction-1 CROSS-NETWORK effect on each perceived
##     slice Y[k]'s objective, referencing Y[self] via interaction1.
##     With shareParameters=TRUE the K cross-network rows are share-marked
##     (effects.r ~1034-1042) so the canonical row's name is renamed to
##     "Y[shared]" and all K share one β.  We address the canonical via
##     name = "Y[shared]", interaction1 = "Y[self]".
##     (Y[self] cannot host the effect: it's structurally locked from
##     evolution in the three-way fork — sim_diag confirmed 0 transitions
##     on Y[self] regardless of rate.  An effect on Y[self]'s objective
##     therefore contributes to zero choice probabilities.)
##   * recip DROPPED on both Y[shared] and Y[self] to keep the moment
##     set as small as possible.
sim_eff <- includeEffects(sim_eff, percRecip,
                          name = "Y[shared]", interaction1 = "Y[self]")

## Set true initial values (no fix=TRUE needed for simOnly runs).
sim_eff <- setEffect(sim_eff, density,   name = "Y[shared]",
                     initialValue = truth["density_perc"])
sim_eff <- setEffect(sim_eff, percRecip,
                     name = "Y[shared]", interaction1 = "Y[self]",
                     initialValue = truth["percRecip"])

sim_eff <- setEffect(sim_eff, density,   name = "Y[self]",
                     initialValue = truth["density_self"])

## Rate parameters.  shareParameters does NOT share-mark rate rows in
## this fork, so we set the slice rate per Y[k] in a loop.  The self-net
## rate is its own row.
for (k in seq_len(n_actors)) {
  sim_eff <- setEffect(sim_eff, Rate, type = "rate",
                       name = paste0("Y[", k, "]"),
                       period = 1,
                       initialValue = truth["rate_perc"])
}
sim_eff <- setEffect(sim_eff, Rate, type = "rate", name = "Y[self]",
                     period = 1,
                     initialValue = truth["rate_self"])


## ---- 5.  Simulate n_reps wave-2 realizations under truth ---------------

## n_sim_draws is what we ask the simulator for; we always request at
## least 2 because n3=1 trips an upstream RSiena bug in phase3.r:80
## (nits = seq(1, 1, 1) has length 1, so nits[2] is NA and the iteration
## check "if (nit == nits[2])" blows up).  We then truncate to n_reps.
n_sim_draws <- max(n_reps, 2)
sim_alg <- sienaAlgorithmCreate(
  projname = "sim_percRecip",
  cond     = FALSE,
  simOnly  = TRUE,
  nsub     = 0,
  n3       = n_sim_draws,
  seed     = 20260430
)

cat("\n=== Drawing", n_sim_draws, "wave-2 realizations under truth ",
    "(harness will use first ", n_reps, ") ===\n", sep = "")
sim_out <- siena07(
  sim_alg,
  data       = base_data,
  effects    = sim_eff,
  returnDeps = TRUE,
  batch      = TRUE,
  verbose    = FALSE
)
stopifnot(length(sim_out$sims) >= n_reps)
## Keep only the first n_reps draws (discard the extra one(s) we asked
## for solely to dodge the n3=1 phase3 bug).
sim_out$sims <- sim_out$sims[seq_len(n_reps)]
stopifnot(length(sim_out$sims) == n_reps)

## sim_out$sims structure for threeway after FRAN expansion:
##   sims[[rep]][[group]][[depvar_name]][[period]] = edge list
## where depvar_name iterates over "Y[1]"..."Y[K]" and "Y[self]".
## We re-stack the K perceiver slices into a (K, n, n) array per rep
## so the next siena07 call receives a fresh 4D input.

edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) {
    m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  }
  m
}

extract_simulated_wave2 <- function(sims_one_rep, n) {
  inner <- sims_one_rep[[1]]                       # group 1
  ## Use positional indexing instead of name lookup.
  ## The threeway simulator returns sims[[rep]][[group]] as a 31-element
  ## list (K=30 slices + Y[self]) but the names are broken: the first
  ## slot is named "Y" and the rest are NA, so inner[["Y[1]"]] returns
  ## NULL and edgeList_to_adj(NULL, n) silently yields a zero matrix.
  ## Position 1..K is Y[1]..Y[K]; position K+1 is Y[self].
  K        <- n
  perc_arr <- array(0L, c(K, n, n))
  for (k in seq_len(K)) {
    perc_arr[k, , ] <- edgeList_to_adj(inner[[k]][[1]], n)
  }
  self_mat <- edgeList_to_adj(inner[[K + 1]][[1]], n)
  list(perc = perc_arr, self = self_mat)
}

simulated_wave2 <- lapply(sim_out$sims, extract_simulated_wave2, n = n_actors)

## Sanity: mean perceived-slice density and self density at wave 2.
mean_perc_density <- mean(sapply(simulated_wave2, function(x) mean(x$perc)))
mean_self_density <- mean(sapply(simulated_wave2, function(x) mean(x$self)))
cat(sprintf("\nMean simulated wave-2 density   slices: %.3f   self: %.3f\n",
            mean_perc_density, mean_self_density))

## ---- 6.  Helper: build a fresh sienaData from one simulated rep --------

build_rep_data <- function(sim_rep) {
  ## Reassemble the 4D (n, n, n, T=2) array.
  arr_t2 <- sim_rep$perc
  ## Splice the simulated self-network onto the diagonal of each slice
  ## (this matches the threeway convention: array[k, k, j] = self k -> j).
  for (k in seq_len(n_actors)) {
    arr_t2[k, k, ] <- sim_rep$self[k, ]
  }
  full <- build_4d(wave1, arr_t2)
  ## IMPORTANT: keep the dependent-variable name as "Y" so that all the
  ## slice names (Y[shared], Y[self], Y[1]..Y[K]) match what the
  ## estimation loop's includeEffects/setEffect calls reference.  If we
  ## bind to a different symbol (e.g. Yrep), the slice names become
  ## Yrep[shared] etc. and the includeEffects calls silently no-op.
  Y <- sienaDependent(full, type = "threeway",
                      shareParameters = TRUE, sharedCov = TRUE)
  sienaDataCreate(Y)
}

## ---- 7.  Estimation loop ------------------------------------------------

## Estimation algorithm.  findiff was set TRUE earlier to dodge an
## auto-fix problem in the broken intra-slice percRecip; the new
## cross-network effect works through the standard score-based LR
## derivative just like crprodRecip, so we switch back to findiff=FALSE
## (default).  LR is faster and gives smoother derivatives.
est_alg <- sienaAlgorithmCreate(
  projname = "est_percRecip",
  cond     = FALSE,
  nsub     = 4,
  n3       = 1000,
  seed     = 20260430
)

## Minimal moment set: density (Y[shared] + Y[self]) and percRecip.
par_names <- c("density_perc", "percRecip", "density_self")

## Map par_names slot -> matching criteria in the included-effects table.
## percRecip is now a cross-network effect on Y[self]'s objective; after the
## Direction-2 Y[shared] expansion it appears as K rows (one per slice)
## sharing one β.  Only the canonical row (sharedDup=FALSE, interaction1
## ends in "[1]") carries the free parameter — the others are tied.  We
## locate the canonical via that extra interaction1 suffix.
## After Direction-1 share expansion the canonical percRecip row has
## name = "Y[shared]" (the user-facing name is restored from "Y[1]" by
## terminateFRAN.r:48-51) and interaction1 = "Y[self]".  Match on those.
param_locator <- list(
  density_perc = list(name = "Y[shared]", short = "density",   interaction1 = ""),
  percRecip    = list(name = "Y[shared]", short = "percRecip", interaction1 = "Y[self]"),
  density_self = list(name = "Y[self]",   short = "density",   interaction1 = "")
)

estimates <- matrix(NA_real_, n_reps, length(par_names),
                    dimnames = list(NULL, par_names))
ses       <- matrix(NA_real_, n_reps, length(par_names),
                    dimnames = list(NULL, par_names))
converged <- logical(n_reps)

for (rep in seq_len(n_reps)) {

  rep_data <- build_rep_data(simulated_wave2[[rep]])

  est_eff <- getEffects(rep_data)
  ## Prune auto-added cross-network rows (same as for sim_eff).
  .cx <- which(est_eff$include & est_eff$interaction1 != "")
  if (length(.cx) > 0) est_eff[.cx, "include"] <- FALSE
  rm(.cx)

  est_eff <- includeEffects(est_eff, percRecip,
                            name = "Y[shared]", interaction1 = "Y[self]")

  ## --- MINIMAL PILOT: fix everything except percRecip ---------------
  ## Free param: percRecip Y[shared] only (default initialValue=0,
  ## must move to 0.8).
  ## Fixed at truth: density Y[shared], density Y[self], all rates.
  ## Dropped entirely: recip on Y[shared] and Y[self].
  est_eff <- setEffect(est_eff, density,   name = "Y[shared]",
                       initialValue = truth["density_perc"], fix = TRUE)
  est_eff <- setEffect(est_eff, density,   name = "Y[self]",
                       initialValue = truth["density_self"], fix = TRUE)
  for (k in seq_len(n_actors)) {
    est_eff <- setEffect(est_eff, Rate, type = "rate",
                         name = paste0("Y[", k, "]"), period = 1,
                         initialValue = truth["rate_perc"], fix = TRUE)
  }
  est_eff <- setEffect(est_eff, Rate, type = "rate", name = "Y[self]",
                       period = 1, initialValue = truth["rate_self"],
                       fix = TRUE)
  ## --- end minimal pilot block --------------------------------------

  est_out <- tryCatch(
    siena07(est_alg, data = rep_data, effects = est_eff,
            batch = TRUE, verbose = FALSE),
    error = function(e) {
      message("rep ", rep, " failed: ", e$message); NULL
    }
  )
  if (is.null(est_out)) next

  ## Use est_out$requestedEffects: the post-expansion table that aligns
  ## with est_out$theta (length z$pp).  terminateFRAN.r:48-51 restores
  ## user-facing names ("Y[shared]" rather than "Y[1]") on this table.
  inc       <- est_out$requestedEffects
  est_theta <- est_out$theta
  est_se    <- sqrt(diag(est_out$covtheta))

  for (k in seq_along(par_names)) {
    spec <- param_locator[[ par_names[k] ]]
    cond <- inc$name         == spec$name        &
            inc$shortName    == spec$short       &
            inc$interaction1 == spec$interaction1
    idx <- which(cond)
    if (length(idx) == 1) {
      estimates[rep, k] <- est_theta[idx]
      ses[rep, k]       <- est_se[idx]
    }
  }
  converged[rep] <- isTRUE(est_out$tconv.max < 0.25)

  if (rep %% 5 == 0 || rep == 1) {
    cat(sprintf(
      "rep %3d  | tconv.max = %.3f | percRecip-hat = %+.3f (truth %.3f)\n",
      rep,
      est_out$tconv.max,
      estimates[rep, "percRecip"],
      truth["percRecip"]))
  }
}

## ---- 8.  Recovery summary ----------------------------------------------

ok <- which(converged & complete.cases(estimates))
cat(sprintf("\n=== Converged: %d / %d ===\n", length(ok), n_reps))

truth_mat <- matrix(truth[par_names], length(ok), length(par_names),
                    byrow = TRUE)

summary_tbl <- data.frame(
  parameter = par_names,
  truth     = unname(truth[par_names]),
  mean_hat  = colMeans(estimates[ok, , drop = FALSE]),
  bias      = colMeans(estimates[ok, , drop = FALSE]) -
              unname(truth[par_names]),
  rmse      = sqrt(colMeans((estimates[ok, , drop = FALSE] - truth_mat)^2)),
  mean_se   = colMeans(ses[ok, , drop = FALSE]),
  coverage  = colMeans(
    abs(estimates[ok, , drop = FALSE] - truth_mat) <
      qnorm(.975) * ses[ok, , drop = FALSE]
  )
)

cat("\n=== Recovery summary ===\n")
print(summary_tbl, row.names = FALSE, digits = 3)

## Outputs land next to this script (tests/), regardless of cwd.
.OUT_DIR <- "/Users/jinwoocho/Desktop/rsiena/tests/group4_perceiver_restricted/artifacts"
.summary_path <- file.path(.OUT_DIR, "recovery_percRecip_summary.csv")
.raw_path     <- file.path(.OUT_DIR, "recovery_percRecip_raw.rds")
write.csv(summary_tbl, file = .summary_path, row.names = FALSE)
saveRDS(list(estimates = estimates, ses = ses, truth = truth,
             converged = converged),
        file = .raw_path)

cat(sprintf("\nWrote: %s\n       %s\n", .summary_path, .raw_path))
