#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: threeWayDiagnostics.r
# *
# * Description: Pre-fit convergence diagnostics for three-way SAOM data.
# *   Reports wave-by-wave Jaccard and density change for each perceived
# *   slice Y[k] and for the self-reported network Y[self], and audits an
# *   effects object for combinations of three-way statistics that are known
# *   to be at risk of identifiability problems.
# *
# *   Three public entry points:
# *     threeWayJaccard()           - per (network, period) Jaccard / density
# *     threeWayEffectsAudit()      - flag risky statistic combinations
# *     threeWayConvergenceDiag()   - master wrapper with print method
# *
# *   Rule of thumb (Snijders, 2017, sec. 3.4): rate parameters are hard to
# *   estimate when J < 0.3 (too much net change between waves) and become
# *   uninformative when the number of changes is near zero. For three-way
# *   SAOM the rule applies independently to each perceiver slice and to the
# *   self-reported network, so a single low-Jaccard slice can prevent
# *   convergence of its slice-specific rate parameter even when the rest of
# *   the system is well-behaved.
# *****************************************************************************/


## =========================================================================
##  Public API
## =========================================================================


##@threeWayJaccard Convergence diagnostic
##
##  Compute wave-by-wave Jaccard index and density change for every
##  perceived slice Y[k] and the self-reported network Y[self] of a
##  three-way sienaDependent object. Flags slices / periods at risk of
##  rate non-convergence.
##
##  Arguments:
##    depvar           a sienaDependent threeway object, a sienaData object
##                     that contains one, or a raw 4D array with dimensions
##                     (slice, sender, receiver, time).
##    jaccard_low      threshold below which the slice / period is flagged
##                     "low_jaccard" (default 0.30, Snijders 2017).
##    jaccard_very_low threshold below which it is escalated to
##                     "very_low_jaccard" (default 0.20).
##    jaccard_high     threshold above which a slice / period with few
##                     changes is flagged "near_stationary" (default 0.95).
##    n_floor          maximum number of changes still considered
##                     "near_stationary" when jaccard > jaccard_high
##                     (default 5).
##
##  Returns a data.frame (class c("threeWayJaccardTable", "data.frame"))
##  with one row per (network, period) and columns:
##    network        "slice_k" (k = 1..N) or "self"
##    period         integer in 1..M-1
##    n_ties_t       number of 1-ties at wave m
##    n_ties_tp1     number of 1-ties at wave m+1
##    n_changed      symmetric difference (toggled entries)
##    n_stable_one   intersection (entries equal to 1 at both waves)
##    n_possible     non-structural-zero, non-missing entries considered
##    jaccard        n_stable_one / |union|
##    density_t      n_ties_t   / n_possible
##    density_tp1    n_ties_tp1 / n_possible
##    density_delta  density_tp1 - density_t
##    risk_flag      one of "OK", "low_jaccard", "very_low_jaccard",
##                   "near_stationary", "stationary", "no_data".
##
threeWayJaccard <- function(depvar,
                            jaccard_low = 0.30,
                            jaccard_very_low = 0.20,
                            jaccard_high = 0.95,
                            n_floor = 5L)
{
  arr <- .threeWayExtractArray(depvar)
  d   <- dim(arr)
  K   <- d[1]
  n   <- d[2]
  nT  <- d[4]
  if (d[3] != n)
    stop("array must be square in actor dimensions (got ",
         paste(d, collapse = "x"), ")")
  if (nT < 2L)
    stop("need at least 2 time points; got ", nT)

  out <- vector("list", (nT - 1L) * (K + 1L))
  idx <- 0L

  for (period in seq_len(nT - 1L)) {

    A <- arr[, , , period,      drop = FALSE]
    B <- arr[, , , period + 1L, drop = FALSE]
    A[A %in% c(10, 11)] <- NA
    B[B %in% c(10, 11)] <- NA

    ## --- Perceived slices ----------------------------------------------
    for (k in seq_len(K)) {
      Ak <- A[k, , , 1L]
      Bk <- B[k, , , 1L]
      diag(Ak) <- NA
      diag(Bk) <- NA
      Ak[k, ]  <- NA  # row k of slice k is structural zero
      Bk[k, ]  <- NA
      idx <- idx + 1L
      out[[idx]] <- .threeWayJaccardRow(
        Ak, Bk, network = paste0("slice_", k), period = period,
        jaccard_low, jaccard_very_low, jaccard_high, n_floor)
    }

    ## --- Self-reported network -----------------------------------------
    Aself <- matrix(NA_real_, n, n)
    Bself <- matrix(NA_real_, n, n)
    for (kk in seq_len(K)) {
      Aself[kk, ] <- A[kk, kk, , 1L]
      Bself[kk, ] <- B[kk, kk, , 1L]
    }
    diag(Aself) <- NA
    diag(Bself) <- NA
    idx <- idx + 1L
    out[[idx]] <- .threeWayJaccardRow(
      Aself, Bself, network = "self", period = period,
      jaccard_low, jaccard_very_low, jaccard_high, n_floor)
  }

  df <- do.call(rbind, out[seq_len(idx)])
  rownames(df) <- NULL
  class(df) <- c("threeWayJaccardTable", "data.frame")
  df
}


##@threeWayEffectsAudit Effects collinearity audit
##
##  Audit an effects object (as returned by getEffects()) for combinations of
##  three-way SAOM statistics that are known to be at risk of identifiability
##  problems. The rules encoded below cover collinearity within and across
##  the four three-way statistic groups (within-network, Y[self] referencing
##  perceived slices, Y[i] referencing Y[self], Y[self] with ego-bound slice
##  access).
##
##  Arguments:
##    effects   a sienaEffects object or a data.frame with the same column
##              names. Only rows with include == TRUE are inspected.
##
##  Returns a data.frame (class c("threeWayEffectsAudit", "data.frame")) with
##  columns effect_a, effect_b, network_a, network_b, risk_level, reason.
##
threeWayEffectsAudit <- function(effects)
{
  if (!is.data.frame(effects))
    stop("effects must be a sienaEffects data.frame (e.g. from getEffects())")
  required <- c("include", "shortName", "interaction1", "name")
  missing  <- setdiff(required, names(effects))
  if (length(missing) > 0L)
    stop("effects is missing required columns: ",
         paste(missing, collapse = ", "))

  inc <- effects[as.logical(effects$include), , drop = FALSE]
  if (nrow(inc) == 0L)
    return(.empty_audit())

  rules <- .threeWayKnownCollinearRules()
  hits  <- list()

  for (r in rules$pairs) {
    a_rows <- which(.match_effect(inc, r$a))
    b_rows <- which(.match_effect(inc, r$b))
    if (length(a_rows) == 0L || length(b_rows) == 0L) next

    for (ai in a_rows) for (bi in b_rows) {
      if (ai == bi) next
      same_net <- inc$name[ai] == inc$name[bi]
      if (isTRUE(r$same_network) && !same_net) next
      if (isTRUE(r$different_network) && same_net) next
      key <- paste(sort(c(ai, bi)), collapse = "_")
      if (!is.null(hits[[key]])) next  # de-dup pair
      hits[[key]] <- data.frame(
        effect_a   = inc$shortName[ai],
        effect_b   = inc$shortName[bi],
        network_a  = inc$name[ai],
        network_b  = inc$name[bi],
        risk_level = r$level,
        reason     = r$reason,
        stringsAsFactors = FALSE)
    }
  }

  if (length(hits) == 0L) return(.empty_audit())
  out <- do.call(rbind, hits)
  rownames(out) <- NULL
  ord <- order(factor(out$risk_level, levels = c("high", "medium", "low")),
               out$effect_a, out$effect_b)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("threeWayEffectsAudit", "data.frame")
  out
}


##@threeWayConvergenceDiag Master diagnostic wrapper
##
##  Run both the Jaccard / density diagnostic and the effects collinearity
##  audit and return an S3 object that prints a readable report.
##
##  Usage:
##    diag <- threeWayConvergenceDiag(myData)
##    diag <- threeWayConvergenceDiag(myData, myEffects)
##    print(diag)
##
##  The ... arguments are forwarded to threeWayJaccard().
##
threeWayConvergenceDiag <- function(data, effects = NULL, ...)
{
  jdf   <- threeWayJaccard(data, ...)
  audit <- if (!is.null(effects)) threeWayEffectsAudit(effects)
           else                   .empty_audit()
  out <- list(jaccard = jdf,
              audit   = audit,
              call    = sys.call())
  class(out) <- "threeWayConvergenceDiag"
  out
}


## =========================================================================
##  Print methods
## =========================================================================


print.threeWayJaccardTable <- function(x, max_rows = 40L, ...) {
  cat("threeWayJaccard table:",
      length(unique(x$network)), "networks,",
      length(unique(x$period)),  "periods\n\n")
  .print_table(as.data.frame(x), max_rows)
  invisible(x)
}


print.threeWayEffectsAudit <- function(x, max_rows = 40L, ...) {
  if (nrow(x) == 0L) {
    cat("threeWayEffectsAudit: no risky pairs identified.\n")
    return(invisible(x))
  }
  cat("threeWayEffectsAudit:", nrow(x), "risky pair(s) found.\n\n")
  .print_table(as.data.frame(x), max_rows)
  invisible(x)
}


print.threeWayConvergenceDiag <- function(x, max_rows = 30L, ...) {
  cat("Three-way SAOM convergence diagnostic\n")
  cat("=====================================\n\n")
  cat("Call: "); print(x$call); cat("\n")

  j <- x$jaccard
  slices  <- grep("^slice_", j$network, value = TRUE)
  n_slices  <- length(unique(slices))
  n_periods <- length(unique(j$period))

  cat(sprintf("Networks evaluated : %d perceiver slice(s) + Y[self]\n",
              n_slices))
  cat(sprintf("Periods evaluated  : %d (waves t -> t+1)\n\n", n_periods))

  cat("Rate-convergence risk summary (by risk_flag):\n")
  tab <- table(factor(j$risk_flag,
                      levels = c("OK", "near_stationary", "stationary",
                                 "low_jaccard", "very_low_jaccard",
                                 "no_data")))
  print(tab)
  cat("\n")

  risky <- j[j$risk_flag != "OK", , drop = FALSE]
  if (nrow(risky) > 0L) {
    cat(sprintf("Showing %d at-risk (network, period) combination(s):\n",
                nrow(risky)))
    .print_table(risky, max_rows)
  } else {
    cat("All (network, period) combinations passed Jaccard/density thresholds.\n")
  }
  cat("\n")

  if (nrow(x$audit) > 0L) {
    cat("Effects collinearity audit:\n")
    .print_table(as.data.frame(x$audit), max_rows)
  } else {
    cat("Effects collinearity audit: no risky pairs identified",
        "(or no effects supplied).\n")
  }

  invisible(x)
}


## =========================================================================
##  Internal helpers
## =========================================================================


## Extract a 4D threeway array from various input shapes.
.threeWayExtractArray <- function(depvar)
{
  ## Case 1: raw 4D array
  if (is.array(depvar) && length(dim(depvar)) == 4L) return(depvar)

  ## Case 2: a sienaDependent threeway object (carries 'type' attribute)
  ty <- attr(depvar, "type")
  if (!is.null(ty) && ty == "threeway") return(unclass(depvar))

  ## Case 3: a siena data object (list with $depvars)
  if (is.list(depvar) && !is.null(depvar$depvars)) {
    dv <- depvar$depvars
    typeMatches <- vapply(dv, function(d)
      identical(attr(d, "type"), "threeway"), logical(1L))
    twIdx <- which(typeMatches)
    if (length(twIdx) == 0L)
      stop("no three-way dependent variable found in siena data object")
    if (length(twIdx) > 1L)
      warning("multiple three-way dependents found; using the first (",
              names(dv)[twIdx[1L]], ")")
    return(unclass(dv[[twIdx[1L]]]))
  }

  stop("depvar must be a 4D array, a sienaDependent threeway object, ",
       "or a siena data object that contains one")
}


## Compute a single diagnostic row from two adjacency matrices A, B with
## the same NA mask (NA entries are excluded from possible / changes).
.threeWayJaccardRow <- function(A, B, network, period,
                                jaccard_low, jaccard_very_low,
                                jaccard_high, n_floor)
{
  ok <- !is.na(A) & !is.na(B)
  n_possible <- sum(ok)

  if (n_possible == 0L) {
    return(data.frame(network = network, period = period,
                      n_ties_t = NA_integer_, n_ties_tp1 = NA_integer_,
                      n_changed = NA_integer_, n_stable_one = NA_integer_,
                      n_possible = 0L,
                      jaccard = NA_real_,
                      density_t = NA_real_, density_tp1 = NA_real_,
                      density_delta = NA_real_,
                      risk_flag = "no_data",
                      stringsAsFactors = FALSE))
  }

  a <- A[ok]
  b <- B[ok]
  n_ties_t     <- sum(a == 1)
  n_ties_tp1   <- sum(b == 1)
  n_changed    <- sum(a != b)
  n_stable_one <- sum(a == 1 & b == 1)
  union_size   <- sum(a == 1 | b == 1)

  jaccard <- if (union_size > 0L) n_stable_one / union_size else NA_real_
  density_t   <- n_ties_t   / n_possible
  density_tp1 <- n_ties_tp1 / n_possible

  flag <- "OK"
  if (!is.na(jaccard)) {
    if (n_changed == 0L) {
      flag <- "stationary"
    } else if (jaccard < jaccard_very_low) {
      flag <- "very_low_jaccard"
    } else if (jaccard < jaccard_low) {
      flag <- "low_jaccard"
    } else if (jaccard > jaccard_high && n_changed < n_floor) {
      flag <- "near_stationary"
    }
  } else if (n_ties_t == 0L && n_ties_tp1 == 0L) {
    flag <- "stationary"
  }

  data.frame(network = network, period = period,
             n_ties_t = n_ties_t, n_ties_tp1 = n_ties_tp1,
             n_changed = n_changed, n_stable_one = n_stable_one,
             n_possible = n_possible,
             jaccard = jaccard,
             density_t = density_t, density_tp1 = density_tp1,
             density_delta = density_tp1 - density_t,
             risk_flag = flag,
             stringsAsFactors = FALSE)
}


## Empty audit shell, used when nothing is found.
.empty_audit <- function()
{
  out <- data.frame(effect_a = character(0), effect_b = character(0),
                    network_a = character(0), network_b = character(0),
                    risk_level = character(0), reason = character(0),
                    stringsAsFactors = FALSE)
  class(out) <- c("threeWayEffectsAudit", "data.frame")
  out
}


## Match an effect spec (a named list with at least shortName) against rows
## of an effects data.frame; returns a logical vector along rows.
.match_effect <- function(eff, spec)
{
  ok <- eff$shortName == spec$shortName
  if (!is.null(spec$interaction1))
    ok <- ok & (eff$interaction1 == spec$interaction1)
  if (!is.null(spec$name_pattern))
    ok <- ok & grepl(spec$name_pattern, eff$name)
  ok
}


## Knowledge base of risky three-way statistic combinations.
##
## Each entry of $pairs is a list with:
##   a, b           named lists specifying the effect (shortName, optionally
##                  interaction1 or name_pattern)
##   same_network    TRUE: only flag if both effects share the 'name' column
##                   (i.e., same dependent network)
##   different_network TRUE: only flag if the two effects belong to
##                   different dependent networks
##   level          "high" | "medium" | "low"
##   reason         human-readable explanation surfaced in the report
##
## The rules below intentionally err toward over-reporting: a flagged pair
## is a hypothesis to verify, not a final verdict. Users should consult the
## printed reason and decide whether the empirical correlation between
## change statistics warrants dropping one of the effects.
##
.threeWayKnownCollinearRules <- function()
{
  list(pairs = list(

    ## ---- Within Group 1: classic SAOM collinearities -------------------

    list(a = list(shortName = "transTrip"),
         b = list(shortName = "gwespFF"),
         same_network = TRUE, level = "high",
         reason = paste("transTrip and gwespFF capture overlapping",
                        "triadic-closure information; keep one.")),

    list(a = list(shortName = "density"),
         b = list(shortName = "outActSqrt"),
         same_network = TRUE, level = "medium",
         reason = paste("density and outActSqrt jointly confound mean",
                        "out-degree with degree-spread; convergence",
                        "often slow even though identification holds.")),

    list(a = list(shortName = "recip"),
         b = list(shortName = "crprodMutual"),
         same_network = FALSE, level = "medium",
         reason = paste("crprodMutual on a companion network reduces",
                        "to a constant times reciprocity when the two",
                        "networks rarely disagree on mutual ties.")),

    ## ---- Group 4 vs Group 1 (self) -------------------------------------

    list(a = list(shortName = "percRecip"),
         b = list(shortName = "recip", interaction1 = ""),
         different_network = FALSE, level = "high",
         reason = paste("Group 4 percRecip is a sub-count of self-network",
                        "reciprocity; the two become collinear when the",
                        "perceived reverse-tie rate is near 1.")),

    ## ---- Group 3 vs Group 1 (perceived density) -----------------------

    list(a = list(shortName = "crprodRecip"),
         b = list(shortName = "density"),
         different_network = FALSE, level = "medium",
         reason = paste("Group 3 crprodRecip approaches a constant times",
                        "perceived density when Y[self] reciprocity is",
                        "near-constant; identification weakens.")),

    list(a = list(shortName = "crprod"),
         b = list(shortName = "density"),
         different_network = FALSE, level = "medium",
         reason = paste("Group 2/3 crprod becomes ~ proportional to the",
                        "focal-network density when the companion",
                        "network's mean tie value is near 1 (or 0).")),

    ## ---- Group 4 vs Group 2 (overlapping cross-network signal) --------

    list(a = list(shortName = "percRecip"),
         b = list(shortName = "crprod"),
         different_network = FALSE, level = "medium",
         reason = paste("percRecip (Group 4) and crprod on Y[self]",
                        "(Group 2) both pull self-reported ties toward",
                        "perceived ties; if both reference a",
                        "reciprocity-like quantity they will compete."))
  ))
}


## Small helper to print a data.frame with an optional row cap.
.print_table <- function(df, max_rows)
{
  df <- as.data.frame(df)
  if (nrow(df) <= max_rows) {
    print(df, row.names = FALSE)
  } else {
    print(head(df, max_rows), row.names = FALSE)
    cat(sprintf("... (%d additional rows suppressed)\n",
                nrow(df) - as.integer(max_rows)))
  }
}
