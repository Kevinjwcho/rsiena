###############################################################################
##
##  diagnose_estimation_mode.R
##
##  Final isolation: does the simulator respect cross-network effects in
##  ESTIMATION mode (where the user's old test_crprod_effects.R succeeded),
##  or is the problem specific to sim_alg's simOnly mode?
##
##  Procedure:
##    1. Generate wave-2 ourselves in R (simple manual DGP — toggle some
##       cells in Y[i] based on Y[self] to encode a known crprodRecip
##       statistic, plus random density-driven toggles).
##    2. Run siena07 in ESTIMATION mode (cond=FALSE, nsub=4, n3=500).
##    3. See if the estimated crprodRecip-hat lands anywhere near truth.
##
##  Outcome A: crprodRecip-hat near 0.8 → estimation works; the bug is
##             specifically in sim_alg simOnly path for threeway
##             cross-network effects.  We replace our wave-2 generator.
##  Outcome B: crprodRecip-hat = 0 / NA → the whole pipeline is broken;
##             share/cross-network/threeway combination needs a deeper
##             fix.
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260430)

n_actors   <- 20
density_p1 <- 0.10

make_wave1_4way <- function(n, p) {
  arr <- array(rbinom(n^3, 1, p), c(n, n, n))
  for (k in seq_len(n)) diag(arr[k, , ]) <- 0L
  arr
}
wave1 <- make_wave1_4way(n_actors, density_p1)

## ---- Generate wave-2 manually with crprodRecip influence --------------
## For each slice i and cell (j, k), decide whether Y[i]_{j, k} flips
## based on: density baseline + crprodRecip boost when Y[self]_{k, j}=1.
##
## In wave 2 we want: cells (j, k) where Y[self]_{k, j}=1 are MORE likely
## to be 1 in Y[i].  This creates a positive crprodRecip statistic.

wave2 <- wave1   # start from wave-1
beta_density    <- -2
beta_crprodRecip <- 1.5   # truth we want to recover

for (i in seq_len(n_actors)) {
  ## Y[self]_{k, j} = wave1[k, k, j] (slice k's row k col j)
  self_kj <- function(k, j) wave1[k, k, j]
  for (j in seq_len(n_actors)) {
    if (j == i) next   # row i is structural-zero
    for (k in seq_len(n_actors)) {
      if (j == k) next
      ## Apply some toggles weighted by the crprodRecip contribution
      eta <- beta_density + beta_crprodRecip * self_kj(k, j)
      p_tie <- 1 / (1 + exp(-eta))
      wave2[i, j, k] <- as.integer(runif(1) < p_tie)
    }
  }
}
## Zero diagonals (j -> j self-loops)
for (k in seq_len(n_actors)) diag(wave2[k, , ]) <- 0L
## Make sure perceiver's own row is structural — leave as wave-1 (will be
## marked 10 by RSiena anyway)
for (i in seq_len(n_actors)) wave2[i, i, ] <- wave1[i, i, ]

cat(sprintf("Wave-1 vs Wave-2 perceived-cell diffs: %d / %d\n",
            sum(wave1 != wave2), length(wave1)))

base_arr <- array(0L, c(n_actors, n_actors, n_actors, 2))
base_arr[ , , , 1] <- wave1
base_arr[ , , , 2] <- wave2

Y <- sienaDependent(base_arr, type = "threeway",
                    shareParameters = TRUE, sharedCov = TRUE)
dat <- sienaDataCreate(Y)

## ---- Effects: crprodRecip on Y[shared] with interaction1 = Y[self] ----

eff <- getEffects(dat)
## DON'T prune — let the auto-cross-network rows stay so includeEffects
## can find them.  (This matches the user's old test_crprod_effects.R.)
eff <- includeEffects(eff, crprodRecip,
                      name = "Y[shared]", interaction1 = "Y[self]")

## ---- Estimation (NOT simOnly) -----------------------------------------

alg <- sienaAlgorithmCreate(projname = "diag_est",
                            cond = FALSE,
                            nsub = 4, n3 = 500,
                            seed = 20260430)

cat("\n=== Estimation: should recover crprodRecip near", beta_crprodRecip, "===\n")
ans <- tryCatch(
  siena07(alg, data = dat, effects = eff, batch = TRUE, verbose = FALSE,
          thetaBound = 500),
  error = function(e) { cat("siena07 error:", conditionMessage(e), "\n"); NULL })

if (!is.null(ans)) {
  reff <- as.data.frame(ans$requestedEffects)
  idx <- which(reff$shortName == "crprodRecip" & reff$name == "Y[shared]")
  if (length(idx) > 0) {
    cat(sprintf("\ncrprodRecip-hat: %s (truth %.3f)\n",
                paste(round(ans$theta[idx], 3), collapse=", "),
                beta_crprodRecip))
    cat(sprintf("tconv.max: %.3f\n", ans$tconv.max))

    if (any(abs(ans$theta[idx] - beta_crprodRecip) < 0.5)) {
      cat("\n✓ Estimation recovers truth — simulator IS using\n")
      cat("  cross-network effects in estimation mode.  The bug is\n")
      cat("  specifically in sim_alg simOnly path for threeway.\n")
    } else {
      cat("\n✗ Estimation also fails to recover.\n")
      cat("  Deeper issue: cross-network + threeway pipeline broken.\n")
    }
  } else {
    cat("\n✗ crprodRecip row not found in requestedEffects.\n")
  }
}
