# install.packages("devtools")
library(devtools)

load_all()

set.seed(7)
X <- array(sample(c(0,1,10,11,NA), 8*8*8*3, replace=TRUE, prob=c(.66,.29,.02,.02,.01)),
           dim=c(8,8,8,3))
# X <- array(sample(c(0,1,NA), 8*8*8*3, replace=TRUE, prob=c(.66,.29,.05)),
#            dim=c(8,8,8,3))
# # 
# X <- array(sample(c(0,1,10,11,NA), 8*8*8, replace=TRUE, prob=c(.66,.29,.02,.02,.01)),
#            dim=c(8,8,8))


# X = array(c(css0, css1, css2), dim = c(20, 20, 20, 3))
Y <- sienaDependent(X, type="threeway")
# Y <- sienaDependent(X, type="oneMode")
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


# Next includeEffects():

eff <- includeEffects(eff, transTrip, name = "Y[1]")
eff <- includeEffects(eff, transTrip, name = "Y[2]")
eff <- includeEffects(eff, transTrip, name = "Y[self]")
eff <- includeEffects(eff, transTrip, name = "Y")

siena_alg <- sienaAlgorithmCreate(projname = 'threeway_test')
ans <- siena07(siena_alg, data = dat, effects = eff, batch = TRUE, silent = TRUE)

