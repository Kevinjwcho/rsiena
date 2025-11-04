# install.packages("devtools")
library(devtools)

load_all()

set.seed(7)
X <- array(sample(c(0,1,10,11,NA), 8*8*8*3, replace=TRUE, prob=c(.86,.09,.02,.02,.01)),
           dim=c(8,8,8,3))
Y <- sienaDependent(X, type = "threeway") 
attr(Y, "type")      # "threeway"
attr(Y, "netdims")   # 8 8 8 3
attr(Y, "nodeSet")   # "Actors"

dat = sienaDataCreate(Y) # undirected or directed density
dat$dycCovars
dat$vCovars
dat$observations
print(dat)
eff <- getEffects(dat)
print(eff) # density  check the rate initial value 
# Is the perceived network directed or undirected?
# Constraints on parameters (are some of them equal)
# getEffects(dat, homogenous = TRUE) / self-reported case (estimate are different with the perceived)
# Goodness of fit: outlier(include covariate effects)

# Individual characteristics
# Dyadic covariates: Distance btw where people live./ perceive more friendship ties for people who live close to you (Dyadic effect (perceiver and sender))
# Perceiver, sender, receiver effects (Right now I just included sender-receiver effects)

# Next includeEffects(): I think it would be similar to getEffects()
# sienaAlgorithmCreate() and siena07() would be the most intense work. 