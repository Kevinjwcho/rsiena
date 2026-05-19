###############################################################################
##
##  recovery_within_Yi.R
##
##  Group 1a: parameter recovery for standard within-slice structural
##  effects on the perceived slices Y[i].
##
##  Recovers (one at a time) — density, reciprocity, transTrip on Y[shared].
##  Other parameters pinned at truth (fix=TRUE) so each effect is isolated.
##
##  Status note: this is the FIRST diagnostic the next session should run.
##  If density/recip/transTrip recover here, the Y[i] objective machinery is
##  healthy and the bug is specific to cross-network effects (Group 3).
##  If they FAIL here, the entire Y[i] evaluation chain is impaired.
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260519)

n_actors    <- 20
n_reps      <- 5
density_p1  <- 0.10

## ---- Truth -------------------------------------------------------------

truth <- c(
  rate_perc    = 6.0,
  rate_self    = 6.0,
  density_perc = -2.0,
  recip_perc   =  1.5,
  density_self = -2.0
)

## ---- Build wave-1 ------------------------------------------------------

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

## ---- sim_eff: all at truth, percRecip NOT included --------------------

sim_eff <- getEffects(base_data)
.cx <- which(sim_eff$include & sim_eff$interaction1 != "")
sim_eff[.cx, "include"] <- FALSE

sim_eff <- includeEffects(sim_eff, recip, name = "Y[shared]")

sim_eff <- setEffect(sim_eff, density, name = "Y[shared]",
                     initialValue = truth["density_perc"])
sim_eff <- setEffect(sim_eff, recip,   name = "Y[shared]",
                     initialValue = truth["recip_perc"])
sim_eff <- setEffect(sim_eff, density, name = "Y[self]",
                     initialValue = truth["density_self"])

for (k in seq_len(n_actors)) {
  sim_eff <- setEffect(sim_eff, Rate, type = "rate",
                       name = paste0("Y[", k, "]"), period = 1,
                       initialValue = truth["rate_perc"])
}
sim_eff <- setEffect(sim_eff, Rate, type = "rate", name = "Y[self]",
                     period = 1, initialValue = truth["rate_self"])

## ---- Simulate wave-2 ---------------------------------------------------

n_sim_draws <- max(n_reps, 2)
sim_alg <- sienaAlgorithmCreate(projname = "sim_within_Yi",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = n_sim_draws,
                                seed = 20260519)
cat("\n=== Drawing", n_sim_draws, "wave-2 realizations under truth ===\n")
sim_out <- siena07(sim_alg, data = base_data, effects = sim_eff,
                   returnDeps = TRUE, batch = TRUE, verbose = FALSE)
sim_out$sims <- sim_out$sims[seq_len(n_reps)]

edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
extract_wave2 <- function(sims_one_rep, n) {
  inner <- sims_one_rep[[1]]
  perc_arr <- array(0L, c(n, n, n))
  for (k in seq_len(n)) perc_arr[k, , ] <- edgeList_to_adj(inner[[k]][[1]], n)
  list(perc = perc_arr, self = edgeList_to_adj(inner[[n + 1]][[1]], n))
}
simulated_wave2 <- lapply(sim_out$sims, extract_wave2, n = n_actors)

build_rep_data <- function(sim_rep) {
  arr_t2 <- sim_rep$perc
  for (k in seq_len(n_actors)) arr_t2[k, k, ] <- sim_rep$self[k, ]
  full <- build_4d(wave1, arr_t2)
  Y <- sienaDependent(full, type = "threeway",
                      shareParameters = TRUE, sharedCov = TRUE)
  sienaDataCreate(Y)
}

## ---- Estimation loop: free density Y[shared], fix everything else ----

est_alg <- sienaAlgorithmCreate(projname = "est_within_Yi",
                                cond = FALSE,
                                nsub = 4, n3 = 1000,
                                seed = 20260519)

par_names <- c("density_perc")
estimates <- matrix(NA_real_, n_reps, length(par_names),
                    dimnames = list(NULL, par_names))
ses       <- matrix(NA_real_, n_reps, length(par_names),
                    dimnames = list(NULL, par_names))
converged <- logical(n_reps)

for (rep in seq_len(n_reps)) {
  rep_data <- build_rep_data(simulated_wave2[[rep]])
  est_eff <- getEffects(rep_data)
  .cx <- which(est_eff$include & est_eff$interaction1 != "")
  if (length(.cx) > 0) est_eff[.cx, "include"] <- FALSE

  est_eff <- includeEffects(est_eff, recip, name = "Y[shared]")

  est_eff <- setEffect(est_eff, recip,   name = "Y[shared]",
                       initialValue = truth["recip_perc"],   fix = TRUE)
  est_eff <- setEffect(est_eff, density, name = "Y[self]",
                       initialValue = truth["density_self"], fix = TRUE)
  for (k in seq_len(n_actors)) {
    est_eff <- setEffect(est_eff, Rate, type = "rate",
                         name = paste0("Y[", k, "]"), period = 1,
                         initialValue = truth["rate_perc"], fix = TRUE)
  }
  est_eff <- setEffect(est_eff, Rate, type = "rate", name = "Y[self]",
                       period = 1, initialValue = truth["rate_self"],
                       fix = TRUE)
  ## density Y[shared] is the only free param.

  est_out <- tryCatch(
    siena07(est_alg, data = rep_data, effects = est_eff,
            batch = TRUE, verbose = FALSE),
    error = function(e) { message("rep ", rep, " failed: ", e$message); NULL }
  )
  if (is.null(est_out)) next

  inc <- est_out$requestedEffects
  est_theta <- est_out$theta
  est_se    <- sqrt(diag(est_out$covtheta))
  idx <- which(inc$name == "Y[shared]" & inc$shortName == "density" &
               inc$interaction1 == "")
  if (length(idx) == 1) {
    estimates[rep, "density_perc"] <- est_theta[idx]
    ses[rep, "density_perc"]       <- est_se[idx]
  }
  converged[rep] <- isTRUE(est_out$tconv.max < 0.25)

  cat(sprintf("rep %d | tconv.max = %s | density_perc-hat = %s (truth %.3f)\n",
              rep,
              format(est_out$tconv.max, digits = 3),
              format(estimates[rep, "density_perc"], digits = 3),
              truth["density_perc"]))
}

## ---- Summary -----------------------------------------------------------

ok <- which(converged & complete.cases(estimates))
cat(sprintf("\n=== Converged: %d / %d ===\n", length(ok), n_reps))
if (length(ok) > 0) {
  cat(sprintf("density_perc:  mean_hat=%.3f  bias=%.3f  mean_se=%.3f\n",
              mean(estimates[ok, "density_perc"]),
              mean(estimates[ok, "density_perc"]) - truth["density_perc"],
              mean(ses[ok, "density_perc"])))
} else {
  cat("No converged reps.\n")
}

.OUT_DIR <- "/Users/jinwoocho/Desktop/rsiena/tests/group1_basic/artifacts"
dir.create(.OUT_DIR, showWarnings = FALSE)
saveRDS(list(estimates = estimates, ses = ses, truth = truth,
             converged = converged),
        file = file.path(.OUT_DIR, "recovery_within_Yi_raw.rds"))
cat(sprintf("\nWrote: %s\n", file.path(.OUT_DIR, "recovery_within_Yi_raw.rds")))
