###############################################################################
##
##  diagnose_no_share_crprodRecip.R
##
##  Step A diagnostic — isolate the shareParameters machinery from the
##  cross-network mechanism.
##
##  Setup:
##    * threeway dependent variable with shareParameters = FALSE
##      (no Y[shared] canonical, no .markSharedDups averaging,
##       no share_groups in initializeFRAN.r)
##    * crprodRecip on a SINGLE perceived slice Y[1] with
##      interaction1 = "Y[self]"
##    * compare sim with β = 0 vs β = 5
##
##  Outcomes:
##    diff > 2   → cross-network mechanism works fine without share;
##                 the share machinery is what breaks Direction 1 sim
##    diff ≈ 0   → cross-network mechanism is broken regardless of
##                 share; the bug is even more fundamental
##                 (effect attachment to Y[i] NetworkVariable,
##                  Y[self] cache initialization, etc.)
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
build_4d <- function(t1, t2) {
  out <- array(0L, c(dim(t1), 2))
  out[ , , , 1] <- t1
  out[ , , , 2] <- t2
  out
}
wave1    <- make_wave1_4way(n_actors, density_p1)
base_arr <- build_4d(wave1, wave1)

## KEY DIFFERENCE: shareParameters = FALSE.
Y <- sienaDependent(base_arr, type = "threeway", shareParameters = FALSE)
base_data <- sienaDataCreate(Y)

cat("\n--- threeway depvar netnames after expansion ---\n")
oneModes <- names(base_data$depvars)
cat(oneModes, sep = "\n")

## ---- Build effects with crprodRecip ONLY on slice Y[1] ----------------

build_eff <- function(beta_crprodRecip) {
  eff <- getEffects(base_data)

  ## Prune all auto-added cross-network rows; we add back only the one
  ## we want.
  .cx <- which(eff$include & eff$interaction1 != "")
  eff[.cx, "include"] <- FALSE

  ## crprodRecip on slice Y[1] only (no share, no aggregate).
  eff <- includeEffects(eff, crprodRecip,
                        name = "Y[1]", interaction1 = "Y[self]")
  eff <- setEffect(eff, crprodRecip,
                   name = "Y[1]", interaction1 = "Y[self]",
                   initialValue = beta_crprodRecip)

  ## Set densities + rates so the network actually evolves.
  for (k in seq_len(n_actors)) {
    eff <- setEffect(eff, density, name = paste0("Y[", k, "]"),
                     initialValue = -2)
    eff <- setEffect(eff, Rate, type = "rate",
                     name = paste0("Y[", k, "]"), period = 1,
                     initialValue = 6)
  }
  eff <- setEffect(eff, density, name = "Y[self]", initialValue = -2)
  eff <- setEffect(eff, Rate, type = "rate", name = "Y[self]",
                   period = 1, initialValue = 6)
  eff
}

## ---- Verify the row is properly enabled --------------------------------

eff_inspect <- build_eff(5)
.df <- as.data.frame(eff_inspect)
cat("\n--- All Y[1] rows in eff (crprodRecip should be include=TRUE,",
    "initialValue=5) ---\n")
print(.df[.df$name == "Y[1]" & .df$include,
          c("shortName", "interaction1", "include", "initialValue", "fix")])

## ---- Simulate β=0 vs β=5 -----------------------------------------------

sim_alg <- sienaAlgorithmCreate(projname = "diag_no_share",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = 5,
                                seed = 20260430)

cat("\n=== sim β=0 ===\n")
sim_out_0 <- siena07(sim_alg, data = base_data, effects = build_eff(0),
                     returnDeps = TRUE, batch = TRUE, verbose = FALSE)

cat("\n=== sim β=5 ===\n")
sim_out_5 <- siena07(sim_alg, data = base_data, effects = build_eff(5),
                     returnDeps = TRUE, batch = TRUE, verbose = FALSE)

## ---- Compare stats on slice Y[1] only ----------------------------------

edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
extract_slice1_and_self <- function(sims_one_rep, n) {
  inner <- sims_one_rep[[1]]
  list(
    slice1 = edgeList_to_adj(inner[[1]][[1]], n),
    self   = edgeList_to_adj(inner[[n + 1]][[1]], n)
  )
}

## crprodRecip stat on a single slice Y[1] with interaction1 Y[self]:
##   s = Σ_{j, k} Y[1]_{j, k} · Y[self]_{k, j}
stat_one_slice <- function(slice_mat, self_mat) {
  s <- 0
  n <- nrow(self_mat)
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      if (slice_mat[j, k] && self_mat[k, j]) s <- s + 1L
    }
  }
  s
}

stats_0 <- sapply(sim_out_0$sims, function(rep) {
  w <- extract_slice1_and_self(rep, n_actors)
  stat_one_slice(w$slice1, w$self)
})
stats_5 <- sapply(sim_out_5$sims, function(rep) {
  w <- extract_slice1_and_self(rep, n_actors)
  stat_one_slice(w$slice1, w$self)
})

## Also report Y[1] differing cells to confirm slice 1 IS evolving.
w0 <- extract_slice1_and_self(sim_out_0$sims[[1]], n_actors)
w5 <- extract_slice1_and_self(sim_out_5$sims[[1]], n_actors)
slice1_wave1 <- wave1[1, , ]
cat(sprintf("\n--- Slice Y[1] wave-1 vs wave-2 diff cells ---\n"))
cat(sprintf("β=0 draw 1: %d / %d\n",
            sum(w0$slice1 != slice1_wave1), n_actors^2))
cat(sprintf("β=5 draw 1: %d / %d\n",
            sum(w5$slice1 != slice1_wave1), n_actors^2))

cat("\n=== crprodRecip stat on slice Y[1] ===\n")
cat(sprintf("β=0: %s   (mean=%.2f)\n",
            paste(stats_0, collapse=", "), mean(stats_0)))
cat(sprintf("β=5: %s   (mean=%.2f)\n",
            paste(stats_5, collapse=", "), mean(stats_5)))
cat(sprintf("\nmean diff (5 - 0): %.2f\n", mean(stats_5) - mean(stats_0)))

if (mean(stats_5) - mean(stats_0) > 2) {
  cat("\n✓ Cross-network mechanism works WITHOUT share machinery.\n")
  cat("  The bug is in share / setEffect cascade / aggregate setup.\n")
} else if (abs(mean(stats_5) - mean(stats_0)) < 0.5) {
  cat("\n✗ Cross-network mechanism is broken even WITHOUT share.\n")
  cat("  Deeper issue: Y[i] objective effects with interaction1=Y[self]\n")
  cat("  aren't reaching the simulator's choice probability computation.\n")
} else {
  cat("\n? Marginal difference; rerun with more draws to confirm.\n")
}
