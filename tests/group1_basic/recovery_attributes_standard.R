###############################################################################
##
##  recovery_attributes_standard.R
##
##  Group 1c-1..3: recovery for STANDARD sender/receiver attribute effects
##  on Y[shared]:
##      egoX  : v_j * outdegree(j)
##      altX  : sum over k of v_k for j -> k ties
##      simX  : indicator[v_j = v_k] for j -> k ties (homophily)
##
##  These are standard SAOM effects that should work whenever the basic
##  within-Y[i] machinery works (Group 1a).
##
##  STUB — implementation deferred until Group 1a recovery is confirmed.
##  Pattern will follow recovery_within_Yi.R: pin everything except the
##  target attribute parameter at truth (fix=TRUE), simulate wave-2 under
##  truth, estimate the free parameter, compare with truth.
##
##  Required additions when implementing:
##    1. Define a constant or changing covariate via coCovar() / varCovar()
##    2. Attach it to the sienaData (sienaDataCreate(Y, covariate))
##    3. includeEffects(eff, egoX | altX | simX, name = "Y[shared]",
##                      interaction1 = "<covariate name>")
##    4. Set true beta and fix in sim; leave free in est
##
###############################################################################

stop("recovery_attributes_standard.R is a stub.\n",
     "Implement after recovery_within_Yi.R is confirmed working.\n",
     "See file header for the pattern.")
