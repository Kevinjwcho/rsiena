## ============================================================
## test_crprod_effects.R
## Check that crprod-family cross-network effects between
## perceived slices Y[shared] and Y[self] can be specified
## and estimated without errors.
##
## Tests 4 effects:
##   1. crprod        : s_j = sum_k Y[i]_{jk} * Y[self]_{jk}
##   2. crprodRecip   : s_j = sum_k Y[i]_{jk} * Y[self]_{kj}
##   3. crprodMutual  : s_j = sum_k Y[i]_{jk} * Y[i]_{kj} * Y[self]_{jk}
##   4. crprodInActIntn: s_j = sum_k Y[i]_{jk} * indeg(Y[self])_k^(1/p)
## ============================================================

library(devtools)
load_all()

## ============================================================
## 0. Data path
## ============================================================

data_path <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData/Advice"

## ============================================================
## 1. Build CSS 3D array helper
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

n <- dim(css_list[[1]])[1]
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
## 4. Specify effects — 4 crprod-family statistics
## ============================================================

eff <- getEffects(dat)

## Base model: density + reciprocity for both shared and self
eff <- includeEffects(eff, recip, name = "Y[shared]")
eff <- includeEffects(eff, recip, name = "Y[self]")

## --- Cross-network effects: Y[shared] <- Y[self] ---

## (1) crprod: perceptual accuracy
##     s_j = sum_k x_{ijk} * x_{jjk}
##     "Does perceiver i's view of j's ties match j's self-report?"
eff <- includeEffects(eff, crprod,
                      interaction1 = "Y[self]",
                      name = "Y[shared]")

## (2) crprodRecip: cross-product with reciprocal tie in Y[self]
##     s_j = sum_k x_{ijk} * x_{kkj}   (j->k perceived, k->j self-reported)
##     "Is perceiver i more likely to see j->k if k self-reports back to j?"
eff <- includeEffects(eff, crprodRecip,
                      interaction1 = "Y[self]",
                      name = "Y[shared]")

## (3) crprodMutual: cross-product gated by mutual tie
##     s_j = sum_k x_{ijk} * x_{ikj} * x_{jjk}
##     "Accuracy conditional on perceived reciprocation"
eff <- includeEffects(eff, crprodMutual,
                      interaction1 = "Y[self]",
                      name = "Y[shared]")

## (4) crprodInActIntn: cross-product weighted by Y[self] indegree
##     s_j = sum_k x_{ijk} * centr_indeg(Y[self])_k^(1/p)
##     "Perceiver i sees j tie to popular-in-Y[self] actors"
eff <- includeEffects(eff, crprodInActIntn,
                      interaction1 = "Y[self]",
                      name = "Y[shared]")

## Print the full included effects table
cat("\n============================================================\n")
cat("INCLUDED EFFECTS:\n")
cat("============================================================\n")
inc <- eff[eff$include, ]
print.data.frame(inc[, c("name", "effectName", "shortName",
                          "interaction1", "type", "sharedDup")])

## Verify: separate Y[self] vs Y[shared] effects
cat("\n--- Y[self] effects ---\n")
self_eff <- inc[grepl("\\[self\\]", inc$name) & !inc$sharedDup, ]
print.data.frame(self_eff[, c("name", "effectName", "shortName", "interaction1")])

cat("\n--- Y[shared] effects (canonical, non-dup) ---\n")
shared_eff <- inc[grepl("\\[shared\\]", inc$name) & !inc$sharedDup, ]
print.data.frame(shared_eff[, c("name", "effectName", "shortName", "interaction1")])

cat("\n--- Cross-network crprod effects (all, including dups) ---\n")
crp_all <- inc[grepl("crprod", inc$shortName), ]
cat(sprintf("  Total crprod rows: %d (canonical: %d, sharedDup: %d)\n",
            nrow(crp_all),
            sum(!crp_all$sharedDup),
            sum(crp_all$sharedDup)))

## ============================================================
## 5. Estimation (short run for testing)
## ============================================================

alg <- sienaAlgorithmCreate(
  projname = "test_crprod",
  n3       = 500,
  nsub     = 2,
  seed     = 12345
)

cat("\n============================================================\n")
cat("STARTING ESTIMATION...\n")
cat("============================================================\n")

ans <- siena07(alg, data = dat, effects = eff,
               batch = TRUE, silent = FALSE)

## ============================================================
## 6. Results
## ============================================================

cat("\n============================================================\n")
cat("RESULTS:\n")
cat("============================================================\n")

## Build results table from requestedEffects
re <- ans$requestedEffects
re_df <- data.frame(
  name       = re$name,
  effectName = re$effectName,
  shortName  = re$shortName,
  inter1     = re$interaction1,
  theta      = ans$theta,
  se         = sqrt(diag(ans$covtheta)),
  stringsAsFactors = FALSE
)
re_df$tstat <- re_df$theta / re_df$se

## Show only non-sharedDup effects
is_dup <- !is.na(re$sharedDup) & re$sharedDup
re_show <- re_df[!is_dup, ]

cat(sprintf("\n  %-45s  %8s  %8s  %8s\n",
            "Effect", "theta", "SE", "t-stat"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(re_show))) {
  lbl <- paste0(re_show$name[i], " / ", re_show$shortName[i])
  if (nzchar(re_show$inter1[i])) lbl <- paste0(lbl, " (", re_show$inter1[i], ")")
  cat(sprintf("  %-45s  %8.3f  %8.3f  %8.3f\n",
              lbl, re_show$theta[i], re_show$se[i], re_show$tstat[i]))
}

cat(sprintf("\ntconv.max = %.4f  (< 0.25 is good)\n", ans$tconv.max))

## ============================================================
## 7. Verification: check effectName alignment
## ============================================================

cat("\n--- Verification: Y[self] effectName should have [self] suffix ---\n")
self_rows <- grepl("\\[self\\]", re$name) & !is_dup
self_efnames <- re$effectName[self_rows]
bad_self <- self_efnames[!grepl("\\[self\\]", self_efnames)]
if (length(bad_self) > 0) {
  cat("  *** MISMATCH: these Y[self] effects have wrong effectName suffix:\n")
  cat(paste("    ", bad_self, collapse = "\n"), "\n")
} else {
  cat("  OK: all Y[self] effectNames correctly have [self] suffix.\n")
}

cat("\n--- Verification: crprod effects should only be in Y[shared] ---\n")
crp_in_self <- grepl("\\[self\\]", re$name) & grepl("crprod", re$shortName)
if (any(crp_in_self)) {
  cat("  *** WARNING: crprod found in Y[self] objective! ***\n")
  print.data.frame(re_df[crp_in_self, ])
} else {
  cat("  OK: no crprod effects in Y[self] objective.\n")
}

## ============================================================
## 8. Save
## ============================================================

save_path <- "/Users/jinwoocho/Library/CloudStorage/Dropbox/Nynke/LongitudinalKrackhardtData"
save(ans, file = file.path(save_path, "test_crprod_effects.RData"))
cat("\nDone. Results saved.\n")
