# install.packages(c("network","sna"))   # once
library(network)
library(sna)

set.seed(1)

n <- 20
density <- 0.147
css <- array(rbinom(n^3, size=1, p=density), c(n,n,n))
for (i in 1:n) for (j in 1:n) css[i, j, j] <- 0  # diagonal elements zero

male <- c(rep(0,n/2), rep(1,n/2))

# --- Parameters ---
lambdaSelf      <- 4.75   # rate for self-reported micro-steps (i updates x[i,i,k])
lambdaPerceived <- 4.66   # rate for perceived micro-steps (i updates x[i,j,k], j!=i)
# Objective coefficients (β): tune these
beta <- list(
  # dens_perceived = -4.0,   # baseline density/activity for perceived j->k by i: Didn't work well because of sparsity?
  dens_perceived = -1.68,   # baseline density/activity for perceived j->k by i
  reciprocity    = 1.98,   # ∑_k x[i,j,k]*x[i,k,j]
  # self_agree     = 0.1,   # ∑_k x[i,j,k]*x[j,j,k]
  self_agree     = 0,   # ∑_k x[i,j,k]*x[j,j,k]
  # perc_influence = -0.5,   # ∑_k x[i,j,k] * (∑_h x[i,i,h]*x[h,j,k])
  perc_influence = 0,   # ∑_k x[i,j,k] * (∑_h x[i,i,h]*x[h,j,k])
  # same_gender    =  0.05,   # bonus if male[j]==male[k]
  same_gender    =  0,   # bonus if male[j]==male[k]
  self_outdeg    = -0.98    # baseline for self network outdegree ∑_k x[i,i,k]
)

# --- Helpers: objective pieces for a single perceiver i and sender j ---


# delta version: change in objective if y[i,j,k] is toggled
delta_f_perceived_ijk <- function(y, i, j, k, beta) {
  # current tie value
  y_cur <- y[i, j, k]
  y_new <- 1 - y_cur   # toggled value
  diff  <- y_new - y_cur   # = +1 if add, -1 if drop
  
  # (1) density/activity: perceived outdegree of j (as seen by i)
  # this changes by ±1 whenever we toggle y[i,j,k]
  delta_dens <- diff
  
  # (2) perceived reciprocity: x[i,j,k] * x[i,k,j]
  delta_recip <- diff * y[i, k, j]
  
  # (3) self-reported tie influence: x[i,j,k] * x[j,j,k]
  delta_self <- diff * y[j, j, k]
  
  # (4) perception-tie influence: x[i,j,k] * sum_h x[i,i,h]*x[h,j,k]
  delta_infl <- diff * sum(y[i, i, ] * y[, j, k])
  
  # (5) same-gender effect: x[i,j,k] * 1(male[j]==male[k])
  delta_gender <- diff * as.numeric(male[j] == male[k])
  
  # linear combination
  beta$dens_perceived * delta_dens +
    beta$reciprocity    * delta_recip +
    beta$self_agree     * delta_self  +
    beta$perc_influence * delta_infl  +
    beta$same_gender    * delta_gender
}

# objective for self network entry y[i,i,k]
# f_self_iik <- function(y, i, k, beta) {
#   # simple self-network outdegree penalty + (optional) placeholder for extras
#   s_out_self <- sum(y[i, i, ])
#   beta$self_outdeg * s_out_self
# }


# delta objective for self network entry y[i,i,k]
delta_f_self_iik <- function(y, i, j, k, beta, male) {
  y_cur <- y[i, i, k]
  y_new <- 1 - y_cur
  diff  <- y_new - y_cur
  beta$self_outdeg * diff
}

# --- NEW FUNCTION: Flip utility using rbinom ---
# This function calculates the probability of a tie existing and then uses
# rbinom() to draw a 0 or 1, updating the network state directly.
toggle_rbinom_delta <- function(x, i, j, k, beta, male, delta_fun, ...) {
  # delta_fun: Δf(x, i, j, k, ...) = f(y with tie=1) - f(y with tie=0)
  
  # Δf compute
  delta_f <- delta_fun(x, i, j, k, beta, ...)
  
  # logistic prob
  p1 <- exp(delta_f) / (1 + exp(delta_f))
  
  # binomial draw
  x[i, j, k] <- rbinom(1, 1, prob = p1)
  
  return(x)
}
# --- Continuous-time simulation with separate rates (Poisson thinning) ---
simulate_interval <- function(x0, t_end, lambdaSelf, lambdaPerceived, beta, male) {
  x <- x0
  t <- 0
  total_rate <- function(n) n*lambdaSelf + n*(n-1)*lambdaPerceived  # i picks self vs. perceived pools
  
  while (t < t_end) {
    dt <- rexp(1, rate = total_rate(n))
    t  <- t + dt
    if (t >= t_end) break
    
    # choose self vs. perceived event by their intensities
    if (runif(1) < (n*lambdaSelf) / total_rate(n)) {
      # --- self update ---
      i <- sample(1:n, 1)
      k <- sample(setdiff(1:n, i), 1)       # forbid i->i self-loop on receiver
      # j==i in the self layer by definition
      x <- toggle_rbinom_delta(  
        x, i, i, k, beta, male,
        delta_fun = delta_f_self_iik
      )
    } else {
      # --- perceived update ---
      i <- sample(1:n, 1)                  # perceiver
      j <- sample(setdiff(1:n, i), 1)      # sender != perceiver
      k <- sample(setdiff(1:n, c(j)), 1)   # receiver != sender; (can equal i, allowed)
      if (j != k) {
        x <- toggle_rbinom_delta( 
          x, i, j, k, beta, male,
          delta_fun = delta_f_perceived_ijk
        )
      }
    }
  }
  x
}


# --- Run your two intervals ---
css0 <- css
# css1 <- simulate_interval(css0, t_end = 1, lambdaSelf, lambdaPerceived, beta, male)
# css2 <- simulate_interval(css1, t_end = 1, lambdaSelf, lambdaPerceived, beta, male)
# css3 <- simulate_interval(css2, t_end = 1, lambdaSelf, lambdaPerceived, beta, male)
css_list <- vector("list", 6)

css_prev <- css0
for (tt in 1:6) {
  css_curr <- simulate_interval(
    css_prev,
    t_end = 1,
    lambdaSelf,
    lambdaPerceived,
    beta,
    male
  )
  css_list[[tt]] <- css_curr
  css_prev <- css_curr
}
# 
# # # --- Plots (as in your code) ---
# par(mfrow = c(1,3))
# plot(network(css0[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 0")
# plot(network(css1[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 1")
# plot(network(css2[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 2")

par(mfrow = c(2,3))
for (tt in 1:6) {
  plot(
    network(css_list[[tt]][1,,]),
    vertex.col = male,
    vertex.cex = 3,
    main = paste0("T = ", tt)
  )
}

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


# print(ans)
# summary(ans)
# ans$effects$effectName

