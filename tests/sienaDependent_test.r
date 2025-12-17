# install.packages("devtools")
library(devtools)

load_all()

set.seed(7)
X <- array(sample(c(0,1,10,11,NA), 8*8*8*3, replace=TRUE, prob=c(.66,.29,.02,.02,.01)),
           dim=c(8,8,8,3))
# 
X <- array(sample(c(0,1,10,11,NA), 8*8*8, replace=TRUE, prob=c(.66,.29,.02,.02,.01)),
           dim=c(8,8,8))


# X = array(c(css0, css1, css2), dim = c(20, 20, 20, 3))
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
names(eff)          # what columns we have
# head(as.data.frame(eff)) # view all effects


# Next includeEffects():

eff <- includeEffects(eff, transTrip, name = "Y[1]")
eff <- includeEffects(eff, transTrip, name = "Y[2]")

# eff <- includeEffects(eff, transTrip)# Include only in the first slice
# eff <- includeEffects(eff, transTrip, name = "Y")

slice_names <- unique(eff$name)

## we can define the group "parm" to control the overall parameters.
for (sn in slice_names) {
  eff <- includeEffects(eff,
                        density,
                        name = sn,
                        parm = 1) 
}
View(eff)


# sienaAlgorithmCreate() and siena07() would be the most intense work. 

siena_alg <- sienaAlgorithmCreate(projname = 'threeway_test')
