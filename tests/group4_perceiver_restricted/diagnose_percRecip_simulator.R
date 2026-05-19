###############################################################################
##
##  diagnose_percRecip_simulator.R
##
##  Direct test: does the simulator actually USE percRecip's contribution
##  when generating trajectories?  Runs sim_alg twice with the same data
##  and same RNG seed, varying ONLY β_percRecip (0 vs 5).  Computes the
##  observed percRecip statistic on each resulting wave-2 and compares.
##
##    diff > ~3   → simulator uses the effect ✓
##                  (then percRecip recovery failure is a signal/noise
##                  or MoM convergence issue, not a wiring issue)
##    diff ≈ 0    → simulator ignores the effect ✗
##                  (perceiver-restricted contribution isn't actually
##                  entering the choice probability — C++ bug)
##
##  Usage:
##    Rscript tests/diagnose_percRecip_simulator.R
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260430)

n_actors   <- 20
density_p1 <- 0.10

## ---- 1. Same base data construction as recovery_percRecip.R ------------

make_wave1_4way <- function(n, p) {
  arr <- array(rbinom(n^3, 1, p), c(n, n, n))
  for (k in seq_len(n)) diag(arr[k, , ]) <- 0L
  arr
}
build_4d <- function(t1, t2) {
  out <- array(0L, c(dim(t1), 2))
  out[ , , , 1] <- t1
  out[ , , , 2] <- t2
  out
}
wave1    <- make_wave1_4way(n_actors, density_p1)
base_arr <- build_4d(wave1, wave1)

Y         <- sienaDependent(base_arr, type = "threeway",
                            shareParameters = TRUE)
base_data <- sienaDataCreate(Y)

## ---- 2. Effects setup helper: pin everything, vary only β_percRecip ----

build_sim_eff <- function(beta_percRecip) {
  eff <- getEffects(base_data)
  ## Prune all auto-added cross-network rows (we'll add percRecip back).
  .cx <- which(eff$include & eff$interaction1 != "")
  if (length(.cx) > 0) eff[.cx, "include"] <- FALSE

  ## percRecip is now Direction 1: on each slice Y[k]'s objective with
  ## interaction1 = "Y[self]" (Y[self] doesn't evolve, so D2 won't work).
  eff <- includeEffects(eff, percRecip,
                        name = "Y[shared]", interaction1 = "Y[self]")

  ## Pin all structural and rate parameters at fixed values.
  eff <- setEffect(eff, density,   name = "Y[shared]",
                   initialValue = -2)
  eff <- setEffect(eff, density,   name = "Y[self]",
                   initialValue = -2)
  eff <- setEffect(eff, percRecip,
                   name = "Y[shared]", interaction1 = "Y[self]",
                   initialValue = beta_percRecip)
  ## CRITICAL: .markSharedDups in effects.r averages initialValues across
  ## the share group, then setEffect only updates the canonical row.  So
  ## the K-1 share-duplicate rows are stuck at 0 unless we explicitly
  ## cascade.  In sim_alg's simOnly mode the phase-2 theta propagation
  ## (theta[dup] <- theta[canon]) never runs, so without this cascade
  ## only slice Y[1] sees β=beta_percRecip; slices Y[2..K] see β=0.
  .dup_idx <- which(eff$shortName == "percRecip" & eff$sharedDup &
                    eff$interaction1 == "Y[self]")
  eff[.dup_idx, "initialValue"] <- beta_percRecip
  rm(.dup_idx)
  for (k in seq_len(n_actors)) {
    eff <- setEffect(eff, Rate, type = "rate",
                     name = paste0("Y[", k, "]"),
                     period = 1, initialValue = 6)
  }
  eff <- setEffect(eff, Rate, type = "rate", name = "Y[self]",
                   period = 1, initialValue = 6)
  eff
}

## Extra test: bump Y[self] rate to 100 (extreme) — if simulator respects
## the rate, wave-2 should change in MANY Y[self] cells.  If it stays at 0,
## the rate is being ignored entirely for Y[self].
build_extreme_eff <- function(beta_percRecip, self_rate) {
  eff <- build_sim_eff(beta_percRecip)
  eff <- setEffect(eff, Rate, type = "rate", name = "Y[self]",
                   period = 1, initialValue = self_rate)
  eff
}

## ---- 3. Run sim_alg twice ----------------------------------------------

sim_alg <- sienaAlgorithmCreate(projname = "diag_sim",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = 5,
                                seed = 20260430)

## Inspect what's actually in sim_eff: percRecip rows AND Y[self] rate row.
eff_inspect <- build_sim_eff(5)
.df <- as.data.frame(eff_inspect)
cat("\n--- percRecip rows in sim_eff (initialValue should be 5) ---\n")
print(.df[.df$include & .df$shortName == "percRecip",
          c("name", "shortName", "interaction1", "include", "initialValue")])
cat("\n--- All Y[self] rows in sim_eff (rate should be 6) ---\n")
print(.df[.df$name == "Y[self]",
          c("name", "shortName", "type", "period", "include", "initialValue", "fix")])

cat("\n=== sim 1: beta_percRecip = 0 ===\n")
sim_out_0 <- siena07(sim_alg, data = base_data, effects = build_sim_eff(0),
                     returnDeps = TRUE, batch = TRUE, verbose = FALSE)

cat("\n=== sim 2: beta_percRecip = 5 ===\n")
sim_out_5 <- siena07(sim_alg, data = base_data, effects = build_sim_eff(5),
                     returnDeps = TRUE, batch = TRUE, verbose = FALSE)

## Use rate=30 (under thetaBound=50).  If Y[self] evolves under rate=6 but
## not under rate=30, our rate setting is being applied normally.  If
## Y[self] differs in 0 cells under BOTH rate=6 AND rate=30, the rate is
## being silently ignored / overridden.
cat("\n=== sim 3: Y[self] rate = 30 (high rate test) ===\n")
sim_out_extreme <- siena07(sim_alg, data = base_data,
                           effects = build_extreme_eff(0, 30),
                           returnDeps = TRUE, batch = TRUE, verbose = FALSE)

## ---- 4. Compute observed percRecip stat on each draw -------------------

edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0)
    m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
extract_wave2 <- function(sims_one_rep, n) {
  inner   <- sims_one_rep[[1]]
  K       <- n
  perc    <- array(0L, c(K, n, n))
  for (k in seq_len(K)) perc[k, , ] <- edgeList_to_adj(inner[[k]][[1]], n)
  self    <- edgeList_to_adj(inner[[K + 1]][[1]], n)
  list(perc = perc, self = self)
}

## percRecip statistic: sum_i sum_j Y[self]_{i,j} * Y[i+1]_{j, i}
percRecip_stat <- function(perc_arr, self_mat) {
  s <- 0
  n <- nrow(self_mat)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (self_mat[i, j] && perc_arr[i, j, i]) s <- s + 1L
    }
  }
  s
}

stats_0 <- sapply(sim_out_0$sims, function(rep) {
  w2 <- extract_wave2(rep, n_actors)
  percRecip_stat(w2$perc, w2$self)
})
stats_5 <- sapply(sim_out_5$sims, function(rep) {
  w2 <- extract_wave2(rep, n_actors)
  percRecip_stat(w2$perc, w2$self)
})

## ---- 5. Verdict --------------------------------------------------------

## Also report whether wave-2 differs from wave-1 (if not, simulator
## didn't actually transition — rate parameters not registering).
wave1_self <- wave1     # use slice 1's self row? Actually wave1 is 3D (K,n,n)
## Compare first sim's wave-2 self matrix to wave-1's diagonal self matrix
## (Y[self][i,j] = wave1[i,i,j])
self_wave1 <- matrix(0L, n_actors, n_actors)
for (i in seq_len(n_actors)) self_wave1[i, ] <- wave1[i, i, ]
w2_first <- extract_wave2(sim_out_0$sims[[1]], n_actors)
n_diff_self <- sum(self_wave1 != w2_first$self)
n_diff_perc <- sum(wave1 != w2_first$perc)
cat(sprintf("\n--- Wave-1 vs Wave-2 (β=0, draw 1) diff cells ---\n"))
cat(sprintf("Y[self] differing cells: %d / %d\n",
            n_diff_self, n_actors*n_actors))
cat(sprintf("Y[perc] differing cells: %d / %d\n",
            n_diff_perc, n_actors*n_actors*n_actors))

w2_extreme <- extract_wave2(sim_out_extreme$sims[[1]], n_actors)
n_diff_self_ex <- sum(self_wave1 != w2_extreme$self)
n_diff_perc_ex <- sum(wave1 != w2_extreme$perc)
cat(sprintf("\n--- Wave-1 vs Wave-2 (Y[self] rate=30, draw 1) diff cells ---\n"))
cat(sprintf("Y[self] differing cells: %d / %d  (rate=30 → expect ~20-40 if rate respected)\n",
            n_diff_self_ex, n_actors*n_actors))
cat(sprintf("Y[perc] differing cells: %d / %d\n",
            n_diff_perc_ex, n_actors*n_actors*n_actors))

cat("\n=== Observed percRecip statistic on simulated wave-2 ===\n")
cat(sprintf("beta = 0  : mean=%.3f  sd=%.3f  values=%s\n",
            mean(stats_0), sd(stats_0), paste(stats_0, collapse=", ")))
cat(sprintf("beta = 5  : mean=%.3f  sd=%.3f  values=%s\n",
            mean(stats_5), sd(stats_5), paste(stats_5, collapse=", ")))
cat(sprintf("\nDifference (5 - 0): %.3f\n", mean(stats_5) - mean(stats_0)))

if (mean(stats_5) - mean(stats_0) > 2) {
  cat("\n✓ Simulator USES the effect (mean stat shifts with beta).\n")
  cat("  → recovery failure is signal/noise or MoM issue, NOT wiring.\n")
} else if (abs(mean(stats_5) - mean(stats_0)) < 0.5) {
  cat("\n✗ Simulator IGNORES the effect (mean stat unchanged at beta=0 vs 5).\n")
  cat("  → C++ wiring bug: contribution not actually entering choice probability.\n")
} else {
  cat("\n? Effect shifts mean stat slightly but not strongly.\n")
  cat("  Inconclusive — may need larger n or stronger beta to confirm.\n")
}
