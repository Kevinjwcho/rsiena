library(network)

n <- 6

density <- 0.05
css <- array(rbinom(n^3, size=1, p=density), c(n,n,n))
# i : perceiver
# j : sender
# k : receiver

#  j to j ties should be 0 
for (i in 1:n)
  for (j in 1:n)
    css[i, j, j] <- 0

male <- c(rep(0,n/2), rep(1,n/2))

# people tend to perceive edges between others of the same gender
# (effect of self-reported edges on perceived edges)
# (perceived reciprocity may be higher than actual reciprocity)
  
# rate parameters
lambdaSelf <- 4
lambdaPerceived <- 4

css1 <- css
time <- 0
while (time < 1) {
  dt <- rexp(1, rate = n*2 + n*(n-1)*2)
  time <- time + dt
  
  selfChange <- rbinom(1, 1, n*2 / (n*2 + n*(n-1)*2))
  if (selfChange) {
    i <- sample(1:n, 1)
    # ...
  } else {
    i <- sample(1:n, 1)
    j <- sample(setdiff(1:n, i), 1)
    k <- sample(setdiff(1:n, c(i,j)), 1)
    if (male[j] == male[k]) {
      css1[i,j,k] <- 1
    }
  }
}

css2 <- css1
while (time < 3) {
  dt <- rexp(1, rate = n*2 + n*(n-1)*2)
  time <- time + dt
  
  selfChange <- rbinom(1, 1, n*2 / (n*2 + n*(n-1)*2))
  if (selfChange) {
    i <- sample(1:n, 1)
    # ...
  } else {
    i <- sample(1:n, 1)
    j <- sample(setdiff(1:n, i), 1)
    k <- sample(setdiff(1:n, c(i,j)), 1)
    if (male[j] == male[k]) {
      css2[i,j,k] <- 1
    }
  }
}

par(mfrow = c(1,3))
plot(network(css[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 0")
plot(network(css1[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 1")
plot(network(css2[1,,]), vertex.col = male, vertex.cex = 3, main = "T = 4")

# not in here: choice probability model

css_array <- array(0L, dim = c(n, n, n, 3))
css_array[,,,1] <- css
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



eff <- includeEffects(eff, transTrip, name = "Y[1]")
eff <- includeEffects(eff, transTrip, name = "Y[2]")
eff <- includeEffects(eff, transTrip, name = "Y[self]")
eff <- includeEffects(eff, transTrip, name = "Y")

siena_alg <- sienaAlgorithmCreate(projname = 'threeway_test')
ans <- siena07(siena_alg, data = dat, effects = eff, batch = TRUE, silent = FALSE)
# ans$
