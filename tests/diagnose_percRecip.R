###############################################################################
##
##  diagnose_percRecip.R
##
##  Single-rep diagnostic for the percRecip recovery harness.
##  Run this from /Users/jinwoocho/Desktop/rsiena to figure out *where*
##  the full recovery_percRecip.R is silently failing (last full run gave
##  50 identical estimates and 100% NA SEs, which means siena07 returned
##  initialValues without iterating).
##
##  This script:
##    1. Builds the same base sienaData as recovery_percRecip.R
##    2. Inspects the effects table BEFORE and AFTER includeEffects()
##       so we can confirm the percRecip row really is there with the
##       expected (name, shortName, include, sharedDup) values.
##    3. Runs siena07 ONCE with verbose=TRUE so the simulator's internal
##       progress is visible.
##    4. Dumps the full structure of the result, including theta vs.
##       initialValue, the dimnames of covtheta, tconv.max, and any
##       warnings/messages generated.
##
##  USAGE:  Rscript tests/diagnose_percRecip.R 2>&1 | tee diagnose.log
##
###############################################################################

suppressPackageStartupMessages({
  library(devtools)
  load_all("/Users/jinwoocho/Desktop/rsiena")
})

set.seed(20260430)

n_actors   <- 10
density_p1 <- 0.10

truth <- c(
  rate_perc    =  6.0,
  rate_self    =  6.0,
  density_perc = -2.0,
  recip_perc   =  1.5,
  percRecip    =  0.8,
  density_self = -2.0,
  recip_self   =  1.5
)

## ---- 1.  Build wave-1 array (same as full harness) ---------------------

make_wave1_4way <- function(n, p) {
  arr <- array(rbinom(n^3, 1, p), c(n, n, n))
  for (k in seq_len(n)) diag(arr[k, , ]) <- 0L
  arr
}

build_4d <- function(arr_t1, arr_t2) {
  out <- array(0L, dim = c(dim(arr_t1), 2))
  out[ , , , 1] <- arr_t1
  out[ , , , 2] <- arr_t2
  out
}

wave1    <- make_wave1_4way(n_actors, density_p1)
base_arr <- build_4d(wave1, wave1)

Y    <- sienaDependent(base_arr, type = "threeway", shareParameters = TRUE)
data <- sienaDataCreate(Y)

cat("\n=== Step 1: data structure ===\n")
print(names(data$depvars))
cat("type:",       attr(data$depvars[[1]], "type"), "\n")
cat("netdims:",    attr(data$depvars[[1]], "netdims"), "\n")
cat("share par:",  isTRUE(attr(data$depvars[[1]], "shareParameters")), "\n")

## ---- 2.  Effects table BEFORE includeEffects ---------------------------

eff <- getEffects(data)

cat("\n=== Step 2a: effects table BEFORE includeEffects ===\n")
cat("nrow =", nrow(eff), "  sum(include) =", sum(eff$include), "\n")
cat("\n-- all rows whose shortName is one of recip/density/percRecip/Rate --\n")
key <- c("density","recip","percRecip","Rate")
.cols <- intersect(c("name","shortName","interaction1","include","sharedDup",
                     "type","initialValue"),
                   colnames(eff))
## NOTE: cast to plain data.frame to bypass print.sienaEffects(), which
## requires a fixed column set (name/effectName/include/fix/test/...).
print(as.data.frame(eff)[eff$shortName %in% key, .cols])

## ---- Prune auto-added cross-network rows -------------------------------
.cx <- which(eff$include & eff$interaction1 != "")
cat("\nPruning", length(.cx), "auto-added cross-network rows\n")
if (length(.cx) > 0) eff[.cx, "include"] <- FALSE

## ---- 3.  includeEffects + setEffect (same as harness) -----------------

eff <- includeEffects(eff, recip,     name = "Y[shared]")
eff <- includeEffects(eff, percRecip, name = "Y[shared]")
eff <- includeEffects(eff, recip,     name = "Y[self]")

cat("\n=== Step 2b: AFTER includeEffects (verify rows are there) ===\n")
print(as.data.frame(eff)[eff$include & eff$shortName %in% key, .cols])

## Set true values via initialValue (sim only).
eff <- setEffect(eff, density,   name = "Y[shared]", initialValue = truth["density_perc"])
eff <- setEffect(eff, recip,     name = "Y[shared]", initialValue = truth["recip_perc"])
eff <- setEffect(eff, percRecip, name = "Y[shared]", initialValue = truth["percRecip"])
eff <- setEffect(eff, density,   name = "Y[self]",   initialValue = truth["density_self"])
eff <- setEffect(eff, recip,     name = "Y[self]",   initialValue = truth["recip_self"])
for (k in seq_len(n_actors)) {
  eff <- setEffect(eff, Rate, type = "rate",
                   name = paste0("Y[", k, "]"),
                   period = 1,
                   initialValue = truth["rate_perc"])
}
eff <- setEffect(eff, Rate, type = "rate", name = "Y[self]",
                 period = 1, initialValue = truth["rate_self"])

cat("\n=== Step 2c: included rows with initialValues ===\n")
print(as.data.frame(eff)[eff$include,
                         intersect(c("name","shortName","initialValue","type"),
                                   colnames(eff))])

## ---- 4.  Simulate ONE wave-2 realization -------------------------------

## NOTE: n3 must be >= 2 to avoid an RSiena bug in phase3.r:80
##   nits <- seq(1, endNit, int)
## With endNit=1, nits has length 1 and nits[2] is NA, blowing up the
## "if (nit == nits[2])" check inside doPhase1or3Iterations.  This is a
## bug in upstream RSiena, not in our percRecip code.  Use n3=5 here
## and keep the first replicate.
sim_alg <- sienaAlgorithmCreate(projname = "diag_sim",
                                cond = FALSE, simOnly = TRUE,
                                nsub = 0, n3 = 5,
                                seed = 20260430)

cat("\n=== Step 3: simulating five wave-2 draws (we use rep 1) ===\n")
sim_out <- siena07(sim_alg, data = data, effects = eff,
                   returnDeps = TRUE, batch = TRUE, verbose = FALSE)
cat("simulator OK, length(sims) =", length(sim_out$sims), "\n")
cat("\nstr(sim_out$sims[[1]], max.level = 3):\n")
str(sim_out$sims[[1]], max.level = 3)

## extract wave-2
edgeList_to_adj <- function(edges, n) {
  m <- matrix(0L, n, n)
  if (NROW(edges) > 0) m[ cbind(edges[, 1], edges[, 2]) ] <- as.integer(edges[, 3])
  m
}
inner    <- sim_out$sims[[1]][[1]]
perc_arr <- array(0L, c(n_actors, n_actors, n_actors))
for (k in seq_len(n_actors)) {
  perc_arr[k, , ] <- edgeList_to_adj(inner[[paste0("Y[", k, "]")]][[1]], n_actors)
}
self_mat <- edgeList_to_adj(inner[["Y[self]"]][[1]], n_actors)
cat(sprintf("\nwave-2 mean perceived density = %.3f, self density = %.3f\n",
            mean(perc_arr), mean(self_mat)))

## ---- 5.  Build the rep's sienaData -------------------------------------

arr_t2 <- perc_arr
for (k in seq_len(n_actors)) arr_t2[k, k, ] <- self_mat[k, ]
rep_arr <- build_4d(wave1, arr_t2)
Y       <- sienaDependent(rep_arr, type = "threeway", shareParameters = TRUE)
rep_dat <- sienaDataCreate(Y)

## ---- 6.  Estimation (FRESH effects, same includeEffects path) ----------

est_eff <- getEffects(rep_dat)
.cx <- which(est_eff$include & est_eff$interaction1 != "")
if (length(.cx) > 0) est_eff[.cx, "include"] <- FALSE
est_eff <- includeEffects(est_eff, recip,     name = "Y[shared]")
est_eff <- includeEffects(est_eff, percRecip, name = "Y[shared]")
est_eff <- includeEffects(est_eff, recip,     name = "Y[self]")

cat("\n=== Step 4: est_eff included rows (estimation input) ===\n")
print(as.data.frame(est_eff)[
  est_eff$include,
  intersect(c("name","shortName","initialValue","fix","include"),
            colnames(est_eff))])
cat("nrow(included):", sum(est_eff$include), "\n")

est_alg <- sienaAlgorithmCreate(projname = "diag_est",
                                cond = FALSE,
                                nsub = 4, n3 = 500,
                                seed = NULL)   # <-- fresh randomness

cat("\n=== Step 5: siena07 (verbose) ===\n")
est_out <- tryCatch(
  withCallingHandlers(
    siena07(est_alg, data = rep_dat, effects = est_eff,
            batch = TRUE, verbose = TRUE),
    warning = function(w) {
      cat("[WARN] ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat("[MSG] ", conditionMessage(m), sep = "")
      invokeRestart("muffleMessage")
    }
  ),
  error = function(e) {
    cat("[ERROR] siena07 threw: ", conditionMessage(e), "\n", sep = "")
    NULL
  }
)

cat("\n=== Step 6: result inspection ===\n")
if (is.null(est_out)) {
  cat("siena07 returned NULL — see error above\n")
} else {
  cat("class:", paste(class(est_out), collapse=","), "\n")
  cat("\ntheta:\n"); print(est_out$theta)
  cat("\nfixed:\n"); print(est_out$fixed)
  cat("\ntconv.max:", est_out$tconv.max, "\n")
  cat("\nnames(covtheta):\n"); print(dimnames(est_out$covtheta))
  cat("\ndiag(covtheta):\n"); print(diag(est_out$covtheta))
  cat("\nany NA in covtheta?", anyNA(est_out$covtheta), "\n")
  cat("\nrequestedEffects rows (name, shortName):\n")
  print(as.data.frame(est_out$requestedEffects)[, c("name","shortName")])
}

cat("\n=== diagnostic done ===\n")
