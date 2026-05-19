set.seed(1)

n <- 16
density <- 0.02
css <- array(rbinom(n^3, size=1, p=density), c(n,n,n))
for (i in 1:n) for (j in 1:n) css[i, j, j] <- 0  # forbid j->j

# --- Rates ---
lambdaSelf      <- 0.5
lambdaPerceived <- 0.5

# ------------------------------------------------------------
# (A) Define slice/group membership for perceivers
#     Here: perceivers 1-3 in group 1, 4-6 in group 2
# ------------------------------------------------------------
perc_group <- rep(1:2, each = n/2)   # (1,1,1,2,2,2)
male <- rep(1:2, each = n/2)   # (1,1,1,2,2,2)

# ------------------------------------------------------------
# (B) Group-specific betas for perceived layer
#     Only a couple effects differ to create identifiability signal
# ------------------------------------------------------------
beta_perc <- list(
  list(
    dens_perceived = -1.2,
    reciprocity    =  0.8,   # group 1: strong reciprocity
    self_agree     =  0.0,
    perc_influence =  0.0,
    same_gender    =  0.0
  ),
  list(
    dens_perceived = -1.2,
    reciprocity    =  0.0,   # group 2: no reciprocity
    self_agree     =  0.0,
    perc_influence =  0.0,
    same_gender    =  0.0
  )
)

# ------------------------------------------------------------
# (C) Self-layer beta
# ------------------------------------------------------------
beta_self <- list(
  self_outdeg = -1.0
)

# --- delta for perceived tie x[i,j,k], using group-specific beta ---
delta_f_perceived_ijk <- function(y, i, j, k, beta_i, male) {
  y_cur <- y[i, j, k]
  diff  <- (1 - y_cur) - y_cur   # +1 if add, -1 if drop
  
  delta_dens   <- diff
  delta_recip  <- diff * y[i, k, j]
  delta_self   <- diff * y[j, j, k]
  delta_infl   <- diff * sum(y[i, i, ] * y[, j, k])
  delta_gender <- diff * as.numeric(male[j] == male[k])
  
  beta_i$dens_perceived * delta_dens +
    beta_i$reciprocity    * delta_recip +
    beta_i$self_agree     * delta_self +
    beta_i$perc_influence * delta_infl +
    beta_i$same_gender    * delta_gender
}

# --- delta for self tie x[i,i,k] ---
delta_f_self_iik <- function(y, i, k, beta_self) {
  y_cur <- y[i, i, k]
  diff  <- (1 - y_cur) - y_cur
  beta_self$self_outdeg * diff
}

# --- Bernoulli resample update (same as your toggle_rbinom_delta idea) ---
update_rbinom <- function(x, i, j, k, delta_f) {
  p1 <- exp(delta_f) / (1 + exp(delta_f))
  x[i, j, k] <- rbinom(1, 1, prob = p1)
  x
}

simulate_interval <- function(x0, t_end, lambdaSelf, lambdaPerceived,
                              beta_perc, beta_self, perc_group, male) {
  x <- x0
  t <- 0
  total_rate <- function(n) n*lambdaSelf + n*(n-1)*lambdaPerceived
  
  while (t < t_end) {
    dt <- rexp(1, rate = total_rate(n))
    t  <- t + dt
    if (t >= t_end) break
    
    if (runif(1) < (n*lambdaSelf) / total_rate(n)) {
      # --- self update ---
      i <- sample(1:n, 1)
      k <- sample(setdiff(1:n, i), 1)
      delta_f <- delta_f_self_iik(x, i, k, beta_self)
      x <- update_rbinom(x, i, i, k, delta_f)
      
    } else {
      # --- perceived update ---
      i <- sample(1:n, 1)
      j <- sample(setdiff(1:n, i), 1)
      k <- sample(setdiff(1:n, j), 1)
      
      g <- perc_group[i]
      beta_i <- beta_perc[[g]]
      
      delta_f <- delta_f_perceived_ijk(x, i, j, k, beta_i, male)
      x <- update_rbinom(x, i, j, k, delta_f)
    }
  }
  x
}

# --- Run waves ---
css0 <- css
css1 <- simulate_interval(css0, t_end = 1, lambdaSelf, lambdaPerceived,
                          beta_perc, beta_self, perc_group, male)
css2 <- simulate_interval(css1, t_end = 3, lambdaSelf, lambdaPerceived,
                          beta_perc, beta_self, perc_group, male)

# # # --- Plots (as in your code) ---
par(mfrow = c(1,3))
plot(network(css0[1,,]), vertex.col = perc_group, vertex.cex = 3, main = "T = 0")
plot(network(css1[1,,]), vertex.col = perc_group, vertex.cex = 3, main = "T = 1")
plot(network(css2[1,,]), vertex.col = perc_group, vertex.cex = 3, main = "T = 4")

## Run the rsiena with three-way option
css_array <- array(0L, dim = c(n, n, n, 3))
css_array[,,,1] <- css0
css_array[,,,2] <- css1
css_array[,,,3] <- css2


# install.packages("devtools")
library(devtools)
load_all()


Y <- sienaDependent(css_array, type="threeway")

attr(Y, "symmetric")                 # check length 2
attr(Y, "type")      # "threeway"
attr(Y, "netdims")   # 8 8 8 3
attr(Y, "nodeSet")   # "Actors"

dat = sienaDataCreate(Y) # undirected or directed density
dat$dycCovars
dat$vCovars
dat$observations
print(dat)


names(dat$depvars)
attr(dat$depvars[[1]], "type")             # type = threeway
attr(dat$depvars[[1]], "distance")         # numeric vector  (to avoid error)
attr(dat$depvars[[1]], "distanceSlices")   # list of distances by perceived slices : time point - 1 for each element in the list
attr(dat$depvars[[1]], "distanceSelf")     # distance of self slice : time point - 1
attr(dat$depvars[[1]], "symmetric")        # two length logical (oneMode)
attr(dat$depvars[[1]], "symmetricSlices")  # length K logical
attr(dat$depvars[[1]], "allUpOnly")        # FALSE/TRUE logical scalar (No NULL)
attr(dat$depvars[[1]], "allDownOnly")      # FALSE/TRUE logical scalar


eff <- getEffects(dat)
print(eff) # density  check the rate initial value
names(eff)          # what columns we have
head(as.data.frame(eff)) # view all effects
dim(as.data.frame(eff))  # 12 rows, 11 columns
unique(eff$name)



# eff <- includeEffects(eff, transTrip, name = "Y[1]")
# eff <- includeEffects(eff, transTrip, name = "Y[2]")
# eff <- includeEffects(eff, transTrip, name = "Y[self]")
# eff <- includeEffects(eff, transTrip, name = "Y")

siena_alg <- sienaAlgorithmCreate(projname = 'threeway_test')

ans <- siena07(siena_alg, data = dat, effects = eff, batch = TRUE, silent = FALSE)

ans$theta[4*(1:17)]
