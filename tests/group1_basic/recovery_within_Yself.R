###############################################################################
##
##  recovery_within_Yself.R
##
##  Group 1b: parameter recovery for standard structural effects on the
##  self-reported network Y[self].
##
##  STATUS: dead path expected.  Y[self] does not evolve in this fork
##  (confirmed in tests/group3_self_to_perception/diagnostics — rate=30
##  produces 0/400 cell changes).  Effects on Y[self]'s objective therefore
##  cannot influence trajectories, and MoM cannot recover their β.
##
##  This script exists to:
##    (a) document that the dead path is real (run and observe failure), and
##    (b) serve as a regression test once Y[self] evolution is enabled.
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260519)

n_actors    <- 20
n_reps      <- 3
density_p1  <- 0.10

truth <- c(
  rate_perc    = 6.0,
  rate_self    = 6.0,
  density_perc = -2.0,
  density_self = -2.0,
  recip_self   =  1.5   # the target parameter for recovery
)

make_wave1 <- function(n, p) {
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
wave1    <- make_wave1(n_actors, density_p1)
base_arr <- build_4d(wave1, wave1)
Y <- sienaDependent(base_arr, type = "threeway",
                    shareParameters = TRUE, sharedCov = TRUE)
base_data <- sienaDataCreate(Y)

## ---- sim_eff: truth values ---------------------------------------------

sim_eff <- getEffects(base_data)
.cx <- which(sim_eff$include & sim_eff$interaction1 != "")
sim_eff[.cx, "include"] <- FALSE

sim_eff <- includeEffects(sim_eff, recip, name = "Y[self]")

sim_eff <- setEffect(sim_eff, density, name = "Y[shared]",
                     initialValue = truth["density_perc"])
sim_eff <- setEffect(sim_eff, density, name = "Y[self]",
                     initialValue = truth["density_self"])
sim_eff <- setEffect(sim_eff, recip,   name = "Y[self]",
                     initialValue = truth["recip_self"])

for (k in seq_len(n_actors)) {
  sim_eff <- setEffect(sim_eff, Rate, type = "rate",
                       name = paste0("Y[", k, "]"), period = 1,
                       initialValue = truth["rate_perc"])
}
sim_eff <- setEffect(sim_eff, Rate, type = "rate", name = "Y[self]",
                     period = 1, initialValue = truth["rate_self"])

## ---- Simulate + sanity check Y[self] differs from wave-1 ---------------

sim_alg <- sienaAlgorithmCreate(projname = "sim_within_Yself",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = max(n_reps, 2),
                                seed = 20260519)
sim_out <- siena07(sim_alg, data = base_data, effects = sim_eff,
                   returnDeps = TRUE, batch = TRUE, verbose = FALSE)
sim_out$sims <- sim_out$sims[seq_len(n_reps)]

edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
extract_self <- function(sims_one_rep, n) {
  inner <- sims_one_rep[[1]]
  edgeList_to_adj(inner[[n + 1]][[1]], n)
}
self_wave1 <- matrix(0L, n_actors, n_actors)
for (i in seq_len(n_actors)) self_wave1[i, ] <- wave1[i, i, ]

cat("\n--- Sanity: does Y[self] evolve under truth recip_self = 1.5? ---\n")
for (r in seq_along(sim_out$sims)) {
  w2_self <- extract_self(sim_out$sims[[r]], n_actors)
  cat(sprintf("rep %d: Y[self] differing cells = %d / %d\n",
              r, sum(w2_self != self_wave1), n_actors^2))
}
cat("(If 0, Y[self] is structurally locked — recovery impossible.)\n")

cat("\n=== Recovery attempt would go here.\n")
cat("Skipping the estimation loop because Y[self] non-evolution means\n")
cat("recip_self stays at initial 0 regardless of MoM iterations.\n")
cat("Re-enable once Y[self] evolution is restored.\n")
