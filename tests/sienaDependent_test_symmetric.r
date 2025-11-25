library(devtools)
load_all()  # your modified RSiena-like package

set.seed(7)
X <- array(
  sample(c(0, 1, 10, 11, NA),
         8 * 8 * 8 * 3,
         replace = TRUE,
         prob = c(.66, .29, .02, .02, .01)),
  dim = c(8, 8, 8, 3)
)

Y <- sienaDependent(X, type = "threeway")
dat <- sienaDataCreate(Y)

## Make a fully symmetric 4D array: (slice, i, j, time)
## For each (slice, time), enforce X[,,k,t] == t(X[,,k,t])
make_all_symmetric <- function(X) {
  d <- dim(X)
  stopifnot(length(d) == 4)
  nslices <- d[1]; n <- d[2]; n2 <- d[3]; nt <- d[4]
  stopifnot(n == n2)  # square matrices
  
  Xsym <- X
  for (k in 1:nslices) {
    for (tt in 1:nt) {
      M <- Xsym[k, , , tt, drop = TRUE]
      ## Copy upper triangle to lower triangle (or vice versa)
      M[lower.tri(M)] <- t(M)[lower.tri(M)]
      Xsym[k, , , tt] <- M
    }
  }
  Xsym
}

## Starting from a fully symmetric array,
## break symmetry at a small number of dyads in selected slice/time.
## This gives a "mixed" case: most matrices look symmetric,
## but at least one is asymmetric.
make_mixed_asymmetric <- function(Xsym, slice = 2, time = 1, n_break = 3) {
  d <- dim(Xsym)
  nslices <- d[1]; n <- d[2]; nt <- d[4]
  stopifnot(slice >= 1 && slice <= nslices)
  stopifnot(time  >= 1 && time  <= nt)
  
  Xmix <- Xsym
  ## Choose some random dyads (i < j) and perturb only one direction
  idx <- combn(1:n, 2)
  n_pairs <- ncol(idx)
  n_break <- min(n_break, n_pairs)
  chosen <- sample(1:n_pairs, n_break)
  
  for (c in chosen) {
    i <- idx[1, c]
    j <- idx[2, c]
    ## Original values
    old_ij <- Xmix[slice, i, j, time]
    old_ji <- Xmix[slice, j, i, time]
    ## Ensure they become different (simple way: flip one side if not NA)
    if (!is.na(old_ij)) {
      Xmix[slice, i, j, time] <- old_ij
      Xmix[slice, j, i, time] <- if (is.na(old_ji) || old_ji == old_ij) {
        1L
      } else {
        old_ji
      }
    } else if (!is.na(old_ji)) {
      Xmix[slice, j, i, time] <- old_ji
      Xmix[slice, i, j, time] <- 0L
    } else {
      ## both NA -> force an asymmetric pair
      Xmix[slice, i, j, time] <- 1L
      Xmix[slice, j, i, time] <- 0L
    }
  }
  Xmix
}

## Force asymmetry in many dyads for every (slice, time),
## starting from a symmetric array.
## This gives a "fully" asymmetric case (practically no matrix is symmetric).
make_all_asymmetric <- function(Xsym, prob_break = 0.5) {
  d <- dim(Xsym)
  nslices <- d[1]; n <- d[2]; nt <- d[4]
  
  Xasym <- Xsym
  for (k in 1:nslices) {
    for (tt in 1:nt) {
      M <- Xasym[k, , , tt, drop = TRUE]
      ## For each dyad (i < j), with some probability, break symmetry
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          if (runif(1) < prob_break) {
            old_ij <- M[i, j]
            old_ji <- M[j, i]
            ## Change only one direction
            if (!is.na(old_ij)) {
              M[i, j] <- old_ij
              M[j, i] <- if (is.na(old_ji) || old_ji == old_ij) 1L else old_ji
            } else if (!is.na(old_ji)) {
              M[j, i] <- old_ji
              M[i, j] <- 0L
            } else {
              M[i, j] <- 1L
              M[j, i] <- 0L
            }
          }
        }
      }
      Xasym[k, , , tt] <- M
    }
  }
  Xasym
}

## Case 1: all slices and times are symmetric -----
X_sym <- make_all_symmetric(X)

Y_sym  <- sienaDependent(X_sym, type = "threeway")
dat_sym <- sienaDataCreate(Y_sym)

eff_sym <- getEffects(dat_sym)
print(eff_sym)

## Inspect: for each slice, it should be treated as symmetric
unique(eff_sym$shortName[eff_sym$type == "rate"])
head(as.data.frame(eff_sym)[grepl("Y\\[1\\]", eff_sym$name) & eff_sym$type == "rate", ])


## Case 2: start from symmetric and break a few dyads in one slice/time ----
X_sym2 <- make_all_symmetric(X)
X_mixed <- make_mixed_asymmetric(X_sym2, slice = 2, time = 1, n_break = 3)

Y_mixed  <- sienaDependent(X_mixed, type = "threeway")
dat_mixed <- sienaDataCreate(Y_mixed)

eff_mixed <- getEffects(dat_mixed)

## Check rate/objective types again
unique(eff_mixed$shortName[eff_mixed$type == "rate"])
head(as.data.frame(eff_mixed)[grepl("Y\\[1\\]", eff_mixed$name) & eff_mixed$type == "rate", ])
head(as.data.frame(eff_mixed)[grepl("Y\\[2\\]", eff_mixed$name) & eff_mixed$type == "rate", ])


## Case 3: force asymmetry widely for all slices and times -----
X_sym3   <- make_all_symmetric(X)
X_asym   <- make_all_asymmetric(X_sym3, prob_break = 0.7)

Y_asym   <- sienaDependent(X_asym, type = "threeway")
dat_asym <- sienaDataCreate(Y_asym)

eff_asym <- getEffects(dat_asym)

## Check again
unique(eff_asym$shortName[eff_asym$type == "rate"])
head(as.data.frame(eff_asym)[grepl("Y\\[1\\]", eff_asym$name) & eff_asym$type == "rate", ])
head(as.data.frame(eff_asym)[grepl("Y\\[2\\]", eff_asym$name) & eff_asym$type == "rate", ])
