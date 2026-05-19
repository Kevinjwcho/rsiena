## ============================================================
## test_crprod_direction2.R
## Direction 2 cross-network effects: Y[self] objective ← perceived slice.
##
## Verifies three includeEffects call patterns:
##   (A) Single-slice, explicit name:
##         includeEffects(eff, crprod, interaction1 = "Y[1]",
##                        name = "Y[self]")
##       → s_j = sum_k Y[self]_{jk} * Y[1]_{jk}
##         One statistic, one β, only perceiver 1's view used.
##
##   (B) Single-slice, auto-default name:
##         includeEffects(eff, crprod, interaction1 = "Y[1]")
##       → same as (A); name = "Y[self]" inferred from i1 pattern.
##
##   (C) Aggregate, auto-default name:
##         includeEffects(eff, crprod, interaction1 = "Y[shared]")
##       → s_j = sum_{i=1..K} sum_k Y[self]_{jk} * Y[i]_{jk}
##         One β, expanded by initializeFRAN.r into K rows with
##         sharedDup cascade so all K slice statistics share one param.
##
## Mirrors tests/test_crprod_effects.R data pipeline and assertion style.
## ============================================================

library(devtools)
load_all()

## ============================================================
## 0. Data path (same as test_crprod_effects.R)
## ============================================================

data_path <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData/Advice"

## ============================================================
## 1. Build CSS 3D array helper (copied from test_crprod_effects.R)
## ============================================================

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

## ============================================================
## 2. Load waves -> 4D array
## ============================================================

f_list   <- sort(list.files(data_path,
                             pattern = "Advice_Week [0-9]+\\.csv",
                             full.names = TRUE))
df_list  <- lapply(f_list, read.csv, check.names = FALSE)
css_list <- lapply(df_list, make_css_wave)

n  <- dim(css_list[[1]])[1]
TT <- length(css_list)

css_array <- array(0L, dim = c(n, n, n, TT))
for (tt in seq_len(TT)) css_array[,,,tt] <- css_list[[tt]]

## ============================================================
## 3. Create RSiena data
## ============================================================

Y   <- sienaDependent(css_array, type = "threeway",
                      shareParameters = TRUE, sharedCov = TRUE)
dat <- sienaDataCreate(Y)

## ============================================================
## 4. Specify effects — Direction 2, three variants
## ============================================================

eff <- getEffects(dat)

## Base model: density + reciprocity for both shared and self
eff <- includeEffects(eff, recip, name = "Y[shared]")
eff <- includeEffects(eff, recip, name = "Y[self]")

## --- Direction 1 (baseline, unchanged): perceived slice ← Y[self] ---
## Keep one Direction 1 effect so we can check the two directions coexist.
eff <- includeEffects(eff, crprod,
                      interaction1 = "Y[self]",
                      name         = "Y[shared]")

## --- Direction 2 (A): single-slice, explicit name = Y[self] ---
##     Should NOT trigger auto-default (name is explicit).
##     Expect: one row, name=Y[self], interaction1=Y[1], sharedDup=FALSE.
eff <- includeEffects(eff, crprodRecip,
                      interaction1 = "Y[1]",
                      name         = "Y[self]")

## --- Direction 2 (B): single-slice, auto-default name ---
##     missing(name) + interaction1="Y[2]" → auto name="Y[self]".
##     Expect: one row, name=Y[self], interaction1=Y[2], sharedDup=FALSE.
eff <- includeEffects(eff, crprodMutual,
                      interaction1 = "Y[2]")

## --- Direction 2 (C): aggregate, auto-default name ---
##     missing(name) + interaction1="Y[shared]" → auto name="Y[self]".
##     Matches the virtual aggregate row created in effects.r.
##     Before C++ dispatch, initializeFRAN.r will expand this into K rows
##     with interaction1 = Y[1]..Y[K], sharedDup cascade on rows 2..K.
eff <- includeEffects(eff, crprodInActIntn,
                      interaction1 = "Y[shared]")

## ============================================================
## 5. Print included effects for inspection
## ============================================================

cat("\n============================================================\n")
cat("INCLUDED EFFECTS (before estimation):\n")
cat("============================================================\n")
inc <- eff[eff$include, ]
print.data.frame(inc[, c("name", "effectName", "shortName",
                          "interaction1", "type", "sharedDup")])

cat("\n--- Direction 1 rows (name ends in [shared] or [k]) ---\n")
d1 <- inc[grepl("\\[(shared|[0-9]+)\\]$", inc$name) & !inc$sharedDup, ]
print.data.frame(d1[, c("name", "shortName", "interaction1")])

cat("\n--- Direction 2 rows (name ends in [self]) ---\n")
d2 <- inc[grepl("\\[self\\]$", inc$name) & !inc$sharedDup, ]
print.data.frame(d2[, c("name", "shortName", "interaction1")])

## ============================================================
## 6. Pre-estimation sanity checks
## ============================================================

cat("\n============================================================\n")
cat("PRE-ESTIMATION SANITY CHECKS:\n")
cat("============================================================\n")

## (A) Direction 2 single-slice explicit: crprodRecip with interaction1=Y[1]
sel_a <- inc$shortName == "crprodRecip" & inc$name == "Y[self]" &
         inc$interaction1 == "Y[1]"
cat(sprintf("(A) crprodRecip, name=Y[self], i1=Y[1]: %d row (expect 1)\n",
            sum(sel_a)))
stopifnot(sum(sel_a) == 1)

## (B) Direction 2 single-slice auto-default: crprodMutual with interaction1=Y[2]
sel_b <- inc$shortName == "crprodMutual" & inc$name == "Y[self]" &
         inc$interaction1 == "Y[2]"
cat(sprintf("(B) crprodMutual, name=Y[self] (auto), i1=Y[2]: %d row (expect 1)\n",
            sum(sel_b)))
stopifnot(sum(sel_b) == 1)

## (C) Direction 2 aggregate: crprodInActIntn with interaction1=Y[shared]
##     At this point (pre-initializeFRAN), the virtual row is still in the table.
sel_c <- inc$shortName == "crprodInActIntn" & inc$name == "Y[self]" &
         inc$interaction1 == "Y[shared]"
cat(sprintf("(C) crprodInActIntn, name=Y[self] (auto), i1=Y[shared]: %d row (expect 1)\n",
            sum(sel_c)))
stopifnot(sum(sel_c) == 1)

## Also: there must be NO crprod family row in Y[self] with interaction1=Y[self]
bad <- inc$name == "Y[self]" & grepl("^crprod", inc$shortName) &
       inc$interaction1 == "Y[self]"
cat(sprintf("(bad) crprod in Y[self] with i1=Y[self]: %d (expect 0)\n", sum(bad)))
stopifnot(sum(bad) == 0)

## ============================================================
## 7. Estimation (short run for testing)
## ============================================================

alg <- sienaAlgorithmCreate(
  projname = "test_crprod_direction2",
  n3       = 500,
  nsub     = 2,
  seed     = 23456
)

cat("\n============================================================\n")
cat("STARTING ESTIMATION...\n")
cat("============================================================\n")

ans <- siena07(alg, data = dat, effects = eff,
               batch = TRUE, silent = FALSE)

## ============================================================
## 8. Post-estimation verification
## ============================================================

re <- ans$requestedEffects
re_df <- data.frame(
  name       = re$name,
  effectName = re$effectName,
  shortName  = re$shortName,
  inter1     = re$interaction1,
  sharedDup  = if (is.null(re$sharedDup)) FALSE else re$sharedDup,
  theta      = ans$theta,
  se         = sqrt(diag(ans$covtheta)),
  stringsAsFactors = FALSE
)
re_df$tstat <- re_df$theta / re_df$se

cat("\n--- Post-initializeFRAN: requestedEffects for Y[self] objective ---\n")
self_rows <- re_df$name == "Y[self]"
print.data.frame(re_df[self_rows, c("name", "shortName", "inter1",
                                     "sharedDup", "theta", "se")])

## Infer K from requested effects (count unique Y[k] slice names)
K_inferred <- length(unique(grep("^Y\\[[0-9]+\\]$", re_df$name, value = TRUE)))
cat(sprintf("\nInferred K = %d perceived slices.\n", K_inferred))

## -- Assertion (A'): crprodRecip (single-slice explicit) stays a single row --
asrt_a <- re_df$shortName == "crprodRecip" & re_df$name == "Y[self]"
cat(sprintf("(A') crprodRecip Y[self] rows: %d (expect 1)\n", sum(asrt_a)))
stopifnot(sum(asrt_a) == 1)
stopifnot(!re_df$sharedDup[asrt_a])

## -- Assertion (B'): crprodMutual (single-slice auto) stays a single row --
asrt_b <- re_df$shortName == "crprodMutual" & re_df$name == "Y[self]"
cat(sprintf("(B') crprodMutual Y[self] rows: %d (expect 1)\n", sum(asrt_b)))
stopifnot(sum(asrt_b) == 1)
stopifnot(!re_df$sharedDup[asrt_b])

## -- Assertion (C'): crprodInActIntn aggregate expanded into K rows --
asrt_c <- re_df$shortName == "crprodInActIntn" & re_df$name == "Y[self]"
cat(sprintf("(C') crprodInActIntn Y[self] rows: %d (expect K = %d)\n",
            sum(asrt_c), K_inferred))
stopifnot(sum(asrt_c) == K_inferred)

## One canonical (sharedDup=FALSE), K-1 dups (sharedDup=TRUE)
n_canon <- sum(asrt_c & !re_df$sharedDup)
n_dup   <- sum(asrt_c &  re_df$sharedDup)
cat(sprintf("     canonical: %d (expect 1), sharedDup: %d (expect %d)\n",
            n_canon, n_dup, K_inferred - 1))
stopifnot(n_canon == 1)
stopifnot(n_dup == K_inferred - 1)

## All K rows reference distinct Y[1]..Y[K] in interaction1
i1_vals <- sort(re_df$inter1[asrt_c])
expected_i1 <- sort(paste0("Y[", 1:K_inferred, "]"))
cat(sprintf("     interaction1 values: %s\n", paste(i1_vals, collapse = ", ")))
stopifnot(identical(i1_vals, expected_i1))

## Duplicates must be fixed (z$fixed[dup] = TRUE via share_groups)
dup_idx <- which(asrt_c & re_df$sharedDup)
canon_theta <- re_df$theta[asrt_c & !re_df$sharedDup]
dup_thetas  <- re_df$theta[dup_idx]
cat(sprintf("     canonical θ = %.4f, dup θ = %s (should all equal canonical)\n",
            canon_theta, paste(round(dup_thetas, 4), collapse = ", ")))
stopifnot(all(abs(dup_thetas - canon_theta) < 1e-10))

## ============================================================
## 9. Final results table (non-dup rows only)
## ============================================================

cat("\n============================================================\n")
cat("RESULTS (non-dup rows only):\n")
cat("============================================================\n")
re_show <- re_df[!re_df$sharedDup, ]
cat(sprintf("\n  %-55s  %8s  %8s  %8s\n",
            "Effect", "theta", "SE", "t-stat"))
cat(strrep("-", 88), "\n")
for (i in seq_len(nrow(re_show))) {
  lbl <- paste0(re_show$name[i], " / ", re_show$shortName[i])
  if (nzchar(re_show$inter1[i])) lbl <- paste0(lbl, " (", re_show$inter1[i], ")")
  cat(sprintf("  %-55s  %8.3f  %8.3f  %8.3f\n",
              lbl, re_show$theta[i], re_show$se[i], re_show$tstat[i]))
}

cat(sprintf("\ntconv.max = %.4f  (< 0.25 is good)\n", ans$tconv.max))

## ============================================================
## 10. Save
## ============================================================

save_path <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData"
save(ans, file = file.path(save_path, "test_crprod_direction2.RData"))
cat("\nDone. Results saved.\n")
