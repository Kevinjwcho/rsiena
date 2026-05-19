###############################################################################
##
##  diagnose_baseline_crprodRecip.R
##
##  Sanity check: does the STANDARD, well-tested crprodRecip (no new C++,
##  no perceiver-restriction) recover its truth in our exact harness setup?
##
##  crprodRecip on Y[i] objective with interaction1 = "Y[self]" computes
##      s = Σ_i Σ_{j, k}  Y[i]_{j, k} · Y[self]_{k, j}
##  via the existing InTieFunction.  It's the UNRESTRICTED cross-network
##  reciprocity statistic — broader than the user's perceiver-restricted
##  target, but uses identical infrastructure: cross-network row layout,
##  shareParameters expansion, MoM share aggregation.
##
##  Two possible outcomes:
##
##    crprodRecip-hat ≈ 0.8  → recovery infrastructure OK
##                              → our PerceiverRestrictedInTieFunction
##                                has a wiring bug specifically (debug it)
##    crprodRecip-hat = NA / 0  → recovery infrastructure is broken
##                                regardless of our custom effect
##                                → fix harness / share / initialValue
##                                  propagation first
##
##  This narrows the problem space drastically.
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260430)

n_actors   <- 20
n_reps     <- 1
density_p1 <- 0.10

truth <- c(
  rate_perc      =  6.0,
  rate_self      =  6.0,
  density_perc   = -2.0,
  crprodRecip    =  0.8,    # the parameter we want to recover
  density_self   = -2.0
)

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
                            shareParameters = TRUE, sharedCov = TRUE)
base_data <- sienaDataCreate(Y)

## ---- sim_eff: crprodRecip on Y[i] with interaction1 = Y[self] ----------

sim_eff <- getEffects(base_data)
.cx <- which(sim_eff$include & sim_eff$interaction1 != "")
sim_eff[.cx, "include"] <- FALSE

## Standard, no-custom-code crprodRecip via the shared-canonical name.
sim_eff <- includeEffects(sim_eff, crprodRecip,
                          name = "Y[shared]", interaction1 = "Y[self]")

sim_eff <- setEffect(sim_eff, density,   name = "Y[shared]",
                     initialValue = truth["density_perc"])
sim_eff <- setEffect(sim_eff, crprodRecip,
                     name = "Y[shared]", interaction1 = "Y[self]",
                     initialValue = truth["crprodRecip"])
## Mirror the percRecip harness: cascade to all share-dup rows so every
## slice's effect sees β=0.8 in sim, not just slice Y[1].
.dup_idx <- which(sim_eff$shortName == "crprodRecip" & sim_eff$sharedDup &
                  sim_eff$interaction1 == "Y[self]")
sim_eff[.dup_idx, "initialValue"] <- truth["crprodRecip"]
rm(.dup_idx)

sim_eff <- setEffect(sim_eff, density,   name = "Y[self]",
                     initialValue = truth["density_self"])

for (k in seq_len(n_actors)) {
  sim_eff <- setEffect(sim_eff, Rate, type = "rate",
                       name = paste0("Y[", k, "]"), period = 1,
                       initialValue = truth["rate_perc"])
}
sim_eff <- setEffect(sim_eff, Rate, type = "rate", name = "Y[self]",
                     period = 1, initialValue = truth["rate_self"])

## ---- Simulate wave 2 ---------------------------------------------------

sim_alg <- sienaAlgorithmCreate(projname = "diag_crprodRecip",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = max(n_reps, 2),
                                seed = 20260430)

cat("=== sim with crprodRecip = 0.8 ===\n")
sim_out <- siena07(sim_alg, data = base_data, effects = sim_eff,
                   returnDeps = TRUE, batch = TRUE, verbose = FALSE)
sim_out$sims <- sim_out$sims[seq_len(n_reps)]

## Also run a quick β=0 vs β=5 sanity check to confirm the simulator
## actually feels crprodRecip.
build_at <- function(beta) {
  e <- getEffects(base_data)
  .cx <- which(e$include & e$interaction1 != "")
  e[.cx, "include"] <- FALSE
  e <- includeEffects(e, crprodRecip,
                      name = "Y[shared]", interaction1 = "Y[self]")
  e <- setEffect(e, density,     name = "Y[shared]", initialValue = -2)
  e <- setEffect(e, density,     name = "Y[self]",   initialValue = -2)
  e <- setEffect(e, crprodRecip, name = "Y[shared]", interaction1 = "Y[self]",
                 initialValue = beta)
  .di <- which(e$shortName == "crprodRecip" & e$sharedDup &
               e$interaction1 == "Y[self]")
  e[.di, "initialValue"] <- beta
  for (k in seq_len(n_actors)) {
    e <- setEffect(e, Rate, type = "rate",
                   name = paste0("Y[", k, "]"), period = 1, initialValue = 6)
  }
  e <- setEffect(e, Rate, type = "rate", name = "Y[self]",
                 period = 1, initialValue = 6)
  e
}
crprodRecip_stat <- function(perc_arr, self_mat) {
  s <- 0
  n <- nrow(self_mat)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      for (k in seq_len(n)) {
        if (perc_arr[i, j, k] && self_mat[k, j]) s <- s + 1L
      }
    }
  }
  s
}
edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
extract_wave2 <- function(sims_one_rep, n) {
  inner <- sims_one_rep[[1]]
  perc <- array(0L, c(n, n, n))
  for (k in seq_len(n)) perc[k, , ] <- edgeList_to_adj(inner[[k]][[1]], n)
  list(perc = perc, self = edgeList_to_adj(inner[[n + 1]][[1]], n))
}

cat("\n=== Sanity: β=0 vs β=5 for crprodRecip ===\n")
sim_at_0 <- siena07(sim_alg, data = base_data, effects = build_at(0),
                    returnDeps = TRUE, batch = TRUE, verbose = FALSE)
sim_at_5 <- siena07(sim_alg, data = base_data, effects = build_at(5),
                    returnDeps = TRUE, batch = TRUE, verbose = FALSE)
stats_0 <- sapply(sim_at_0$sims, function(r) {
  w <- extract_wave2(r, n_actors); crprodRecip_stat(w$perc, w$self) })
stats_5 <- sapply(sim_at_5$sims, function(r) {
  w <- extract_wave2(r, n_actors); crprodRecip_stat(w$perc, w$self) })
cat(sprintf("β=0 stats: %s\n", paste(stats_0, collapse=", ")))
cat(sprintf("β=5 stats: %s\n", paste(stats_5, collapse=", ")))
cat(sprintf("mean diff: %.3f\n", mean(stats_5) - mean(stats_0)))
if (abs(mean(stats_5) - mean(stats_0)) < 0.5) {
  cat("\n✗ EVEN crprodRecip is ignored by the simulator under this harness.\n")
  cat("  The bug is in the harness / share / setEffect chain,\n")
  cat("  NOT in our custom PerceiverRestrictedInTieFunction.\n")
} else {
  cat("\n✓ crprodRecip IS used by the simulator.\n")
  cat("  Our PerceiverRestrictedInTieFunction has a wiring issue specific\n")
  cat("  to it.  Compare side-by-side with InTieFunction's setup.\n")
}
